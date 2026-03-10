// @ts-check

/** @type {import('@docusaurus/types').Config} */
const config = {
  title: 'cfDNAlab',
  tagline: 'Fast and transparent cfDNA command-line analysis',
  favicon: 'img/cfdnalab_logo_little_guy_200x181_144dpi.png',

  url: 'https://ludvigolsen.github.io',
  baseUrl: '/cfdnalab/',

  organizationName: 'ludvigolsen',
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
      image: 'img/cfdnalab_logo_257x285_250dpi.png',
      navbar: {
        title: 'cfDNAlab',
        logo: {
          alt: 'cfDNAlab logo',
          src: 'img/cfdnalab_logo_little_guy_200x181_144dpi.png'
        },
        items: [
          {
            type: 'docSidebar',
            sidebarId: 'docsSidebar',
            position: 'left',
            label: 'Docs'
          },
          {
            href: 'https://github.com/ludvigolsen/cfdnalab',
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
        copyright: `Copyright © ${new Date().getFullYear()} cfDNAlab`
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
