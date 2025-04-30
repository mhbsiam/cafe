// @ts-check
import { defineConfig } from "astro/config";
import starlight from "@astrojs/starlight";

// https://astro.build/config
export default defineConfig({
  site: "https://mhbsiam.github.io",
  base: "/cafe",

  integrations: [
    starlight({
      title: "CAFE Docs",
      description:
        "Documentation for CAFE, an open-source, free, no-code, web-app platform for high-dimensional spectral flow cytometry data (SFCM) analysis.",
      social: [
        {
          icon: "github",
          label: "GitHub",
          href: "https://github.com/mhbsiam/cafe",
        },
      ],
      sidebar: [
        {
          label: "Getting Started",
          autogenerate: { directory: "getting-started" },
        },
        {
          label: "Setting Up",
          items: [
            "installation/system_requirements",
            "installation/download",
            {
              label: "Installation",
              items: ["installation/pixi_way", "installation/conda_way"],
            },
            "installation/run_cafe",
          ],
        },
        {
          label: "Workflow",
          items: [
            "workflow/diagram",
            "workflow/preparation",
            "workflow/scaled_csv_file_structure",
          ],
        },

        "others/faq",
        "others/data_processing_hpc",
        "others/tool_dependencies",
        "others/known_issues",
        "others/citation",
      ],
    }),
  ],
});
