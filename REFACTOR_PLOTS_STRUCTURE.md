REFACTOR_PLOTS_STRUCTURE
========================

Purpose
-------
Map plotting functions into a clearer file/folder layout and record where each function is called so we can split safely.

New folders
-----------
- `R/plots/` — plotting modules (small focused files)
- `R/utils/` — shared helpers and small utilities

Planned file mapping and roles
------------------------------
- `R/plots/heatmaps.R`
  - Functions: `heatmap.3`, `plotCNA`, `plotSubclones`
  - Role: heatmap rendering and high-level CNA heatmap plots used by pipelines and classifiers

- `R/plots/segmentation.R`
  - Functions: `segmPos`, `segmPosSpec`, `plotCNAline`, `plotCNAlineOnlyTumor`, `plotConsensusCNA`
  - Role: segmentation helpers and line-based CNA visualizations

- `R/plots/tsne_and_misc.R`
  - Functions: `plotTSNE`, `plotCNsubclones`, `plotCN`
  - Role: t-SNE and miscellaneous per-sample visualizations

- `R/plots/onco_heat.R`
  - Functions: `plotOncoHeat`, `plotOncoHeatSubclones`, `plotOncoHeat2`
  - Role: oncogenic event heatmaps and subclone-centric heatmaps

- `R/plots/clone_tree.R`
  - Functions: `plotCloneTree`, `plotCloneTreeNew`
  - Role: clone tree visualizations across samples or within a sample

- `R/utils/plot_helpers.R`
  - Functions: `modifyNameCNV`, color palettes, small helpers (e.g. name/pos helpers, segm helpers)
  - Role: shared utilities used across plotting modules (internal)

Current callers (found via grep)
--------------------------------
- `plotCNA`
  - Called in: `R/classifyTumor.R` (lines ~269, 315), `R/pipelineCNA.R` (line ~134)
  - Role: used by classification and pipeline scripts to produce sample heatmaps

- `plotSubclones`
  - Called in: `R/testFunc.R` (~line 192)

- `plotTSNE`
  - Called in: `R/pipelineCNA.R` (~line 236)

- `plotOncoHeatSubclones` / `plotOncoHeat`
  - Called in: `R/multiSampleComparisonClonalCN.R` (~lines 238, 242), `R/pipelineCNA.R` (~line 234)

- `plotCloneTree` / `plotCloneTreeNew`
  - Called in: `R/multiSampleComparisonClonalCN.R` (~line 242), `R/pipelineCNA.R` (~line 238, conditional)

- `heatmap.3`
  - Called in: `R/plotSegmentation.R` (~line 225)

Notes & guidance
-----------------
- Keep public functions (those called from other `R/*.R` files) exported using roxygen (`@export`) after moving them.
- Mark low-level helpers as internal (`@keywords internal` or `@noRd`).
- After splitting, run the following in R at the package root to update docs/NAMESPACE and validate the package:

```r
install.packages("devtools")
devtools::document()
devtools::check()
```

Next steps
----------
- Complete the split: create `R/plots/*.R` and `R/utils/plot_helpers.R` with implementations copied from `R/plots/heatmaps.R`.
- Add roxygen headers and `@export` for public functions.
- Run `devtools::document()` and `devtools::check()` and iterate on issues (missing imports, namespace problems).

If you want, I can now generate the concrete new R files with roxygen skeletons and move functions into them as patches you can apply.

Function Reference
------------------
The table below lists plotting functions, their main inputs/outputs, a one-line description, and the planned file path where the function will live after the refactor.

| Function | Inputs | Outputs | Description | Planned Path |
|---|---|---|---|---|
| heatmap.3 | numeric matrix, annotation (data.frame), plotting params | draws heatmap / returns ggplot/gtable | Low-level heatmap renderer used by segmentation and heatmap plots | R/plots/heatmaps.R |
| plotCNA | CNA matrix or processed matrix, sample name, options | plot (printed) and optionally saved file | High-level CNA heatmap for a sample (calls `heatmap.3`) | R/plots/heatmaps.R |
| plotSubclones | subclone assignment list, CNA matrix, sample | plot of subclone heatmaps | Visualise subclone-specific CN patterns | R/plots/heatmaps.R |
| segmPos | segmentation object, chromosome map | numeric positions / segments table | Helper to compute segment positions for plotting | R/plots/segmentation.R |
| segmPosSpec | segmentation object, filter params | segments subset | Specialized segmentation position selector | R/plots/segmentation.R |
| plotCNAline | CNA profile matrix, sample, plotting params | line plot of CN across genome | Line-based genomic CNA visualization per sample | R/plots/segmentation.R |
| plotCNAlineOnlyTumor | CNA matrix restricted to tumor cells, sample | line plot | CNA line plot showing only tumour cells (no normals) | R/plots/segmentation.R |
| plotConsensusCNA | list of CNA profiles (cells/subclones) | consensus plot/table | Summarise consensus CNAs across groups/subclones | R/plots/segmentation.R |
| plotTSNE | expression or reduced-dim matrix, annotation | t-SNE scatter plot | t-SNE visualization colored by annotations (tumor/normal/subclone) | R/plots/tsne_and_misc.R |
| plotCNsubclones | CNA matrix, subclone labels | clustered heatmap / plots | Per-subclone CNA visualisation (compact) | R/plots/tsne_and_misc.R |
| plotCN | CNA matrix, gene/region selection | heatmap or matrix plot | Generic CNA matrix plotting helper | R/plots/tsne_and_misc.R |
| plotOncoHeat | event matrix (binary or counts), sample/meta | oncoprint-style heatmap | Oncogenic event heatmap across cells or subclones | R/plots/onco_heat.R |
| plotOncoHeatSubclones | event matrix, subclone assignments | oncoprint by subclone | Onco-heat focused on subclone membership | R/plots/onco_heat.R |
| plotOncoHeat2 | event matrix, custom ordering | alternative oncoprint variant | Variant oncoprint with different layout/annotations | R/plots/onco_heat.R |
| plotCloneTree | phylo / tree object, sample metadata | clone-tree plot (ggtree) | Visualise inferred clone relationships as a tree | R/plots/clone_tree.R |
| plotCloneTreeNew | enhanced tree + annotations | enhanced clone-tree plot | Newer clone-tree visualization with annotations | R/plots/clone_tree.R |
| modifyNameCNV | string vector (names), mapping rules | cleaned name vector | Utility: normalize/shorten CNV region/gene names | R/utils/plot_helpers.R |
| (palette helpers) | none or parameters | color palette object | Shared color palettes for consistent plotting | R/utils/plot_helpers.R |

Notes:

- After moving functions, add `@export` to public functions and `@keywords internal` to helper utilities.
- If you want, I can create the skeleton files in `R/plots/` and `R/utils/` and copy the implementations over, preserving existing code comments and roxygen blocks.

