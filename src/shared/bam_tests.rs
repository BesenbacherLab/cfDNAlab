use anyhow::{Context, Result, anyhow, bail};
use rust_htslib::bam::Read;
use std::{
    fs,
    io::{BufRead, BufReader, Write},
    net::{TcpListener, TcpStream},
    path::{Path, PathBuf},
    sync::{
        Arc, Mutex, MutexGuard, OnceLock,
        atomic::{AtomicBool, Ordering},
    },
    thread,
};
use tempfile::TempDir;

// These tests are kept in a separate file because they exercise private BAM-opening helpers from
// `bam.rs` and need more scaffolding than the usual small unit tests.
//
// Remote indexed BAM access is not just "open this URL". HTSlib expects the server to behave like a
// byte-addressable file store: it may send HEAD requests, ask for byte ranges, and discover the
// index by trying standard sidecar paths next to the BAM URL. The small server below implements only
// that subset so the test verifies cfDNAlab's URL path through rust-htslib without depending on an
// external public BAM file.
#[test]
fn indexed_reader_fetches_records_from_http_bam_url() -> Result<()> {
    let bam = crate::testing::single_contig_inward_pair_bam()?;
    let server = BamHttpServer::new(fs::read(bam.bam_path())?, fs::read(bam.bai_path())?)?;

    // When only a remote BAM URL is passed, HTSlib may cache the discovered index sidecar in the
    // current directory before opening the indexed reader. Without this guard, running the test from
    // the repository root can leave files such as `remote.bam.bai` behind. The guard also serializes
    // current-directory changes because the process current directory is global state.
    let _current_dir = TemporaryCurrentDir::new()?;

    let mut reader = open_indexed_bam_reader(Path::new(&server.bam_url()))?;
    reader.fetch(("chr1", 0, 200))?;
    let positions = reader
        .records()
        .map(|record| record.map(|record| record.pos()))
        .collect::<std::result::Result<Vec<_>, _>>()?;

    // The fixture contains one inward-facing pair spanning [20, 80):
    // forward read starts at 20, reverse read starts at 60.
    assert_eq!(positions, vec![20, 60]);
    Ok(())
}

#[test]
fn build_bam_bai_index_writes_explicit_bam_bai_sidecar() -> Result<()> {
    let bam = crate::testing::single_contig_inward_pair_bam()?;
    let temp_dir = TempDir::new()?;
    let output_bam = temp_dir.path().join("indexed-output.bam");
    fs::copy(bam.bam_path(), &output_bam).with_context(|| {
        format!(
            "copy fixture BAM {} to {}",
            bam.bam_path().display(),
            output_bam.display()
        )
    })?;

    let output_bai = build_bam_bai_index(&output_bam)?;

    // cfDNAlab writes BAM indexes next to generated BAMs using the standard `<name>.bam.bai`
    // sidecar name. That exact filename matters because users and HTSlib both discover indexes by
    // path convention when a command only receives the BAM path.
    assert_eq!(output_bai, temp_dir.path().join("indexed-output.bam.bai"));

    let mut reader = IndexedReader::from_path(&output_bam)?;
    reader.fetch(("chr1", 0, 200))?;
    let fetched_records = reader
        .records()
        .collect::<std::result::Result<Vec<_>, _>>()?;
    assert_eq!(
        fetched_records.len(),
        2,
        "the generated index should make the fixture's two records fetchable"
    );

    Ok(())
}

struct TemporaryCurrentDir {
    previous_dir: PathBuf,
    _temp_dir: TempDir,
    _guard: MutexGuard<'static, ()>,
}

impl TemporaryCurrentDir {
    fn new() -> Result<Self> {
        static LOCK: OnceLock<Mutex<()>> = OnceLock::new();
        let guard = LOCK
            .get_or_init(|| Mutex::new(()))
            .lock()
            .map_err(|err| anyhow!("temporary current directory lock poisoned: {err}"))?;
        let previous_dir = std::env::current_dir().context("get test current directory")?;
        let temp_dir = TempDir::new().context("create temporary current directory")?;
        std::env::set_current_dir(temp_dir.path()).with_context(|| {
            format!(
                "switch test current directory to {}",
                temp_dir.path().display()
            )
        })?;
        Ok(Self {
            previous_dir,
            _temp_dir: temp_dir,
            _guard: guard,
        })
    }
}

impl Drop for TemporaryCurrentDir {
    fn drop(&mut self) {
        // Drop cannot return `Result`. If restoring the current directory fails, there is no useful
        // recovery path inside the test helper; the next filesystem operation will surface the bad
        // process state.
        let _ = std::env::set_current_dir(&self.previous_dir);
    }
}

struct BamHttpServer {
    address: String,
    should_stop: Arc<AtomicBool>,
    handle: Option<thread::JoinHandle<()>>,
}

impl BamHttpServer {
    // The server owns the BAM and BAI bytes and serves them from fixed paths:
    //
    // - `/remote.bam` is the input URL passed to cfDNAlab's BAM opener.
    // - `/remote.bam.bai` and `/remote.bai` are the standard index sidecar names HTSlib probes when
    //   it has to infer the index URL from the BAM URL.
    //
    // Binding to port 0 asks the OS for a free local port, which avoids collisions between test runs.
    fn new(bam_bytes: Vec<u8>, bai_bytes: Vec<u8>) -> Result<Self> {
        let listener = TcpListener::bind("127.0.0.1:0").context("bind local HTTP server")?;
        let address = listener.local_addr()?.to_string();
        let should_stop = Arc::new(AtomicBool::new(false));
        let thread_should_stop = Arc::clone(&should_stop);

        let handle = thread::spawn(move || {
            for stream in listener.incoming() {
                if thread_should_stop.load(Ordering::SeqCst) {
                    break;
                }
                if let Ok(stream) = stream {
                    let _ = serve_bam_request(stream, &bam_bytes, &bai_bytes);
                }
            }
        });

        Ok(Self {
            address,
            should_stop,
            handle: Some(handle),
        })
    }

    fn bam_url(&self) -> String {
        format!("http://{}/remote.bam", self.address)
    }
}

impl Drop for BamHttpServer {
    fn drop(&mut self) {
        // `incoming()` blocks while waiting for the next connection. Store the stop flag, then open a
        // throwaway connection to wake the server thread so it can exit and be joined.
        self.should_stop.store(true, Ordering::SeqCst);
        let _ = TcpStream::connect(&self.address);
        if let Some(handle) = self.handle.take() {
            let _ = handle.join();
        }
    }
}

fn serve_bam_request(mut stream: TcpStream, bam_bytes: &[u8], bai_bytes: &[u8]) -> Result<()> {
    let request = HttpRequest::read(&stream)?;
    let body = match request.path.as_str() {
        "/remote.bam" => bam_bytes,
        "/remote.bam.bai" | "/remote.bai" => bai_bytes,
        _ => {
            write_http_response(&mut stream, "404 Not Found", None, &[], request.is_head)?;
            return Ok(());
        }
    };

    let (start, end) = requested_byte_range(&request.headers, body.len())?;
    let is_partial = start != 0 || end != body.len();
    let content_range =
        is_partial.then(|| format!("bytes {}-{}/{}", start, end - 1, body.len()));
    let status = if is_partial {
        "206 Partial Content"
    } else {
        "200 OK"
    };
    write_http_response(
        &mut stream,
        status,
        content_range.as_deref(),
        &body[start..end],
        request.is_head,
    )?;
    Ok(())
}

struct HttpRequest {
    path: String,
    headers: Vec<String>,
    is_head: bool,
}

impl HttpRequest {
    // This reads HTSlib's HTTP request, not the BAM file. The BAM and BAI are binary, but the HTTP
    // request framing that asks for them is text: a request line such as `GET /remote.bam HTTP/1.1`,
    // followed by headers such as `Range: bytes=123-456`, then a blank line. We only parse enough of
    // that text envelope to decide which fixture bytes to return and whether to omit the body for a
    // HEAD request.
    fn read(stream: &TcpStream) -> Result<Self> {
        let mut reader = BufReader::new(stream.try_clone()?);
        let mut request_line = String::new();
        reader.read_line(&mut request_line)?;
        let mut request_parts = request_line.split_whitespace();
        let method = request_parts.next().context("HTTP method")?;
        let path = request_parts.next().context("HTTP path")?.to_string();

        let mut headers = Vec::new();
        loop {
            let mut line = String::new();
            reader.read_line(&mut line)?;
            if line == "\r\n" || line.is_empty() {
                break;
            }
            headers.push(line);
        }

        Ok(Self {
            path,
            headers,
            is_head: method == "HEAD",
        })
    }
}

fn requested_byte_range(headers: &[String], body_len: usize) -> Result<(usize, usize)> {
    // HTSlib issues normal single-range requests such as `Range: bytes=123-456` while seeking in the
    // remote BAM/index. This parser intentionally rejects anything outside that subset instead of
    // pretending to be a general HTTP implementation.
    let Some(header) = headers
        .iter()
        .find(|line| line.to_ascii_lowercase().starts_with("range:"))
    else {
        return Ok((0, body_len));
    };

    let Some(raw_range) = header.split_once("bytes=").map(|(_, value)| value.trim()) else {
        bail!("unsupported HTTP range header: {}", header.trim());
    };
    let (raw_start, raw_end) = raw_range
        .split_once('-')
        .context("HTTP byte range separator")?;
    let start = raw_start.parse::<usize>().context("HTTP byte range start")?;
    let end = if raw_end.is_empty() {
        body_len
    } else {
        raw_end
            .parse::<usize>()
            .context("HTTP byte range end")?
            .saturating_add(1)
    };
    if start >= body_len || end > body_len || start >= end {
        bail!("HTTP byte range {start}-{end} outside body length {body_len}");
    }
    Ok((start, end))
}

fn write_http_response(
    stream: &mut TcpStream,
    status: &str,
    content_range: Option<&str>,
    body: &[u8],
    omit_body: bool,
) -> Result<()> {
    // `Accept-Ranges` advertises the behavior HTSlib needs for random access. `Content-Range` is
    // included only for partial responses; full-body and HEAD responses use the same helper.
    write!(
        stream,
        "HTTP/1.1 {status}\r\nContent-Length: {}\r\nAccept-Ranges: bytes\r\nConnection: close\r\n",
        body.len()
    )?;
    if let Some(content_range) = content_range {
        write!(stream, "Content-Range: {content_range}\r\n")?;
    }
    write!(stream, "\r\n")?;
    if !omit_body {
        stream.write_all(body)?;
    }
    stream.flush()?;
    Ok(())
}
