import clsx from 'clsx';
import Layout from '@theme/Layout';
import Link from '@docusaurus/Link';

function HeroButtons() {
  return (
    <div className="heroActions">
      <Link className="button button--primary button--lg" to="/docs/get-started/installation">
        Get Started
      </Link>
      <Link className="button button--secondary button--lg" to="/docs/generated/cli/index">
        Command Reference
      </Link>
      <Link className="button button--secondary button--lg" to="/docs/release-notes">
        Release Notes
      </Link>
    </div>
  );
}

export default function Home() {
  return (
    <Layout
      title="cfDNAlab"
      description="Fast and transparent cfDNA command-line analysis"
    >
      <main className={clsx('heroSection')}>
        <div className="heroInner">
          <img
            src="/img/cfdnalab_logo_257x285_250dpi.png"
            alt="cfDNAlab logo"
            className="heroLogo"
          />
          <h1>cfDNAlab</h1>
          <p>
            A command-line toolkit for cfDNA analytics with clear contracts,
            reproducible workflows, and release-grade validation.
          </p>
          <HeroButtons />
        </div>
      </main>
    </Layout>
  );
}
