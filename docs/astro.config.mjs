// @ts-check
import { defineConfig } from 'astro/config';
import starlight from '@astrojs/starlight';

// https://astro.build/config
export default defineConfig({
  site: 'https://combine-lab.github.io',
  base: '/gravlax',
  integrations: [
    starlight({
      logo: {
        light: './src/assets/logo-light.svg',
        dark: './src/assets/logo-dark.svg',
        replacesTitle: true,
      },
      title: 'Gravlax',
      description:
        'Align once and query forever — a compact molecular-evidence index for annotation replay in single-cell RNA-seq.',
      social: [
        {
          icon: 'github',
          label: 'GitHub',
          href: 'https://github.com/COMBINE-lab/gravlax',
        },
      ],
      editLink: {
        baseUrl: 'https://github.com/COMBINE-lab/gravlax/edit/main/docs/',
      },
      sidebar: [
        {
          label: 'Getting started',
          items: [
            { label: 'Overview', link: '/' },
            { label: 'Installation', link: '/installation/' },
            { label: 'Distribution and integrity', link: '/distribution/' },
            { label: 'Quick start', link: '/quickstart/' },
          ],
        },
        {
          label: 'Guides',
          items: [
            { label: 'Workflow and interfaces', link: '/workflow/' },
            { label: 'Capabilities', link: '/capabilities/' },
            { label: 'Python and AnnData', link: '/python/' },
            { label: 'The .aie format', link: '/format/' },
          ],
        },
        {
          label: 'CLI reference',
          items: [
            { label: 'Overview', link: '/cli/' },
            { label: 'Projects and plans', link: '/cli/projects/' },
            { label: 'aie doctor', link: '/cli/doctor/' },
            { label: 'aie explore', link: '/cli/explore/' },
            { label: 'Shell completions', link: '/cli/completions/' },
            { label: 'aie resolve', link: '/cli/resolve/' },
            { label: 'Ingest setup and preflight', link: '/cli/ingest/' },
            { label: 'aie ingest-archive', link: '/cli/ingest-archive/' },
            { label: 'aie compile-annotation', link: '/cli/compile-annotation/' },
            { label: 'aie export-molecule-bam', link: '/cli/export-molecule-bam/' },
            { label: 'aie replay-rows', link: '/cli/replay-rows/' },
            { label: 'Compare annotations', link: '/cli/compare-annotations/' },
            { label: 'aie query', link: '/cli/query/' },
            { label: 'aie collection', link: '/cli/collection/' },
            { label: 'Transcript equivalence classes', link: '/cli/transcript-ecs/' },
            { label: 'aie federate', link: '/cli/federate/' },
            { label: 'aie cohort', link: '/cli/cohort/' },
            { label: 'aie extend', link: '/cli/extend/' },
            { label: 'aie stamp-genome', link: '/cli/stamp-genome/' },
            { label: 'Archive identity and sealing', link: '/cli/archive-integrity/' },
            { label: 'aie dev em', link: '/cli/em/' },
            { label: 'Development commands', link: '/cli/development/' },
          ],
        },
      ],
    }),
  ],
});
