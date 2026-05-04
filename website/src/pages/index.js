import clsx from 'clsx';
import Layout from '@theme/Layout';
import Link from '@docusaurus/Link';
import useBaseUrl from '@docusaurus/useBaseUrl';

function HeroButtons() {
  return (
    <div className="heroActions">
      <Link className="button button--primary button--lg" to="/docs/get-started/installation">
        Get Started
      </Link>
      <Link className="button button--secondary button--lg" to="/docs/generated/cli/overview">
        Command Reference
      </Link>
      <Link className="button button--secondary button--lg" to="/docs/generated/release-notes">
        Release Notes
      </Link>
    </div>
  );
}

export default function Home() {
  const logoUrl = useBaseUrl('/img/cfdnalab_logo_750x500_150dpi.png');

  return (
    <Layout
      title="cfDNAlab"
      description="Fast and transparent cfDNA command-line analysis"
    >
      <main className={clsx('heroSection')}>
        <div className="heroInner">
          <img
            src={logoUrl}
            alt="cfDNAlab logo"
            className="heroLogo"
          />
          <div className="heroTaglineWrap">
            <p className="heroTagline">
              Extract <b>fragmentation features</b> from sequenced cell-free DNA
              with an ultra-fast, highly-flexible but easy-to-use command line tool.
            </p>
          </div>
          <HeroButtons />
        </div>
      </main>
    </Layout>
  );
}
