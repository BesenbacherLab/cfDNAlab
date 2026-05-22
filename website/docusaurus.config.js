// @ts-check

/** @type {import('@docusaurus/types').Config} */
const config = {
  title: 'cfDNAlab',
  tagline: 'Fast and transparent cfDNA command-line analysis',
  favicon: 'img/cfdnalab_logo_little_guy_172x200_144dpi.png',

  url: 'https://cfDNAlab.tools',
  baseUrl: '/',

  organizationName: 'BesenbacherLab',
  projectName: 'cfdnalab',

  onBrokenLinks: 'throw',
  onBrokenAnchors: 'warn',

  i18n: {
    defaultLocale: 'en',
    locales: ['en']
  },

  presets: [
    [
      'classic',
      /** @type {import('@docusaurus/preset-classic').Options} */
      ({
        docs: {
          routeBasePath: 'docs',
          sidebarPath: require.resolve('./sidebars.js')
        },
        blog: false,
        theme: {
          customCss: require.resolve('./src/css/custom.css')
        }
      })
    ]
  ],

  plugins: [
    [
      require.resolve('@easyops-cn/docusaurus-search-local'),
      {
        indexDocs: true,
        docsRouteBasePath: '/docs',
        indexBlog: false,
        indexPages: true,
        language: ['en'],
        hashed: true
      }
    ]
  ],

  themeConfig:
    /** @type {import('@docusaurus/preset-classic').ThemeConfig} */
    ({
      prism: {
        additionalLanguages: ['r', 'bash']
      },
      navbar: {
        title: 'cfDNAlab',
        logo: {
          alt: 'cfDNAlab logo',
          src: 'img/cfdnalab_logo_little_guy_172x200_144dpi.png'
        },
        items: [
          {
            type: 'docSidebar',
            sidebarId: 'docsSidebar',
            position: 'left',
            label: 'Docs'
          },
          {
            href: 'https://github.com/BesenbacherLab/cfdnalab',
            label: 'GitHub',
            position: 'right'
          }
        ]
      },
      footer: {
        style: 'dark',
        links: [
          {
            title: 'Docs',
            items: [
              {
                label: 'Get Started',
                to: '/docs/get-started/installation'
              },
              {
                label: 'CLI Reference',
                to: '/docs/generated/cli/overview'
              }
            ]
          }
        ],
        copyright: `Copyright © ${new Date().getFullYear()}<br /><a href="https://ludvigolsen.dk" target="_blank" rel="noopener">Ludvig Renbo Olsen</a>, <a href="https://pure.au.dk/portal/da/persons/besenbacher@clin.au.dk/" target="_blank" rel="noopener">Søren Besenbacher</a>, <a href="https://www.moma.dk/bioinformatics/computational-genomics" target="_blank" rel="noopener">BesenbacherLab</a>`
      }
    }),

  customFields: {
    generatedDirNotice: 'Generated files live in docs/generated and must not be edited manually'
  },

  staticDirectories: ['static'],
  trailingSlash: false,
  future: {
    v4: true
  },

  themes: [],
  markdown: {
    mermaid: false,
    hooks: {
      onBrokenMarkdownLinks: 'throw'
    }
  },

  scripts: [
    // Site analytics with no tracking etc.
    {
      src: 'https://gc.zgo.at/count.js',
      async: true,
      'data-goatcounter': 'https://cfdnalab.goatcounter.com/count'
    }
  ],

  headTags: [
    {
      tagName: 'meta',
      attributes: {
        name: 'description',
        content: 'cfDNAlab command reference and user documentation'
      }
    }
  ],

  clientModules: [],
  titleDelimiter: '·'
};

module.exports = config;
