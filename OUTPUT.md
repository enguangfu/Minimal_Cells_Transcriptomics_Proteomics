# Conventions for Output

## Tabular data

- When output tabular data, prefer TSV over CSV and Excel file because TSV can be viewed in VS Code better.

## Scripts

- Write important information including background, algorithm and summary of outcome into docstring in the head of Python and Jupyter Notebook files.
- Always generate a txt file to hold the output info from the Python script. 

## Figure generation 


Follow these conventions for every matplotlib figure.

### Canvas size
- Use inches as the unit throughout. Do not use pixels or millimeters.
- Full figure: `figsize=(7, 7)`.
- Half figure: `figsize=(7/2, 7/2)` (i.e., 3.5 × 3.5 in).
- Third figure: `figsize=(7/3, 7/3)` (i.e., ~2.33 × 2.33 in).
- Other sizes: specify explicitly as needed.
- 7 inches matches Nature's 2-column width (180 mm ≈ 7.09 in), rounded down
  for a safety margin.
- Never resize figures after generation — they should be born at final print size.

### Fonts
- Family: Arial, sans-serif.
- Size range: 5 pt minimum, 7 pt maximum (Nature requirement, strictly enforced).
- Default rcParams:
  - `font.size`: 7
  - `axes.titlesize`: 7
  - `axes.labelsize`: 7
  - `xtick.labelsize`: 6
  - `ytick.labelsize`: 6
  - `legend.fontsize`: 6
- Inside-plot annotations (text drawn inside axes via `ax.text()`,
  `ax.annotate()`, panel labels, etc.): use 5 pt.
- Never set any text above 7 pt or below 5 pt.

### PDF export for Illustrator compatibility
- Set `pdf.fonttype: 42` and `ps.fonttype: 42` in rcParams. This keeps text as
  editable text in Illustrator rather than outlined paths.
- Export as PDF: `fig.savefig('figure_name.pdf', dpi=300)`.
- DPI only affects embedded rasters; pure vector content is resolution-independent.
- Avoid `bbox_inches='tight'` when exact figure dimensions matter — it can trim
  the figure below the specified `figsize`. Use `constrained_layout=True` or
  manual `subplots_adjust()` instead.

### Color
- Use RGB color mode (Nature requires RGB for original research content).
- Matplotlib defaults to RGB, so no special action needed — just do not convert
  to CMYK anywhere in the pipeline.

### What to avoid
- No drop shadows, 3D effects, bevels, or glows. Nature exports these as low-res
  bitmaps and may reject the figure.
- No rasterization of vector elements. Keep `rasterized=False` (the default)
  unless a specific element genuinely needs rasterization (e.g., a dense
  scatter plot with millions of points).
- Do not embed bitmapped images in plots for this project — all figures are
  pure vector.
- By default, do NOT set title.

### Standard preamble
Every plotting script should begin with this configuration block:

​```python
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({
    'font.size': 7,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial'],
    'axes.titlesize': 7,
    'axes.labelsize': 7,
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,
    'legend.fontsize': 6,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})
​```

For inside-plot annotations, pass `fontsize=5` explicitly:

​```python
ax.text(x, y, 'label', fontsize=5)
ax.annotate('note', xy=(x, y), fontsize=5)
​```

### File output
- Save all figures as PDF in the project's figure output directory.
- Name files descriptively (e.g., `fig1_de_volcano.pdf`, `fig2_gsea_dotplot.pdf`).
- Each figure should be a single PDF file at the appropriate canvas size.

## Operon plot with RNA isoforms

Per-operon visualization produced by `Syn1_Operon/Operon_Visualization.py`
(`plot_one_operon`), driven by `Operon_Visualization.ipynb`. Three stacked tracks in
**transcript space** (5'->3' left->right; the operon strand flips genomic
coordinates). Born at print size, OUTPUT.md fonts (Arial 5-7 pt, `pdf.fonttype 42`).

Tracks (top -> bottom):
- **Gene track** -- every gene **touching** the padded operon span `[s0−200, e0+200]`
  (pad = `200 bp`). The plot window is then **expanded** so every touched gene is shown
  with its FULL body (`plot_s0 = min(s0−200, min touched start)`, `plot_e0 = max(e0+200,
  max touched end)`) -- no flanking gene is clipped at the window edge.
  Polygon-arrows pointing by gene strand. **Gray (`#7a7a7a`) for both strands;
  pseudogenes purple (`#b0a0c8`).** Genes SENSE to the operon (same strand) = solid
  outline, opaque (alpha 0.95); genes ANTISENSE to the operon = **dotted outline,
  transparent fill (alpha 0.38)**, lighter label. **EVERY gene is labelled** (5 pt):
  `<gene name>/locusNum` using the **syn3A ortholog's gene name** (from
  `Genomes_Input/syn3a_genome.gff3`, matched by the preserved locus number), falling back
  to the syn1 name, then `pseudo/locusNum`, then the bare locus number. If that full label
  would be wider than the gene's box (estimated cheaply with no renderer), the label
  **collapses to the bare locus number** (e.g. `0671`) so it never overflows the box.
- **Isoform track** -- PacBio FLNC isoforms on the operon strand overlapping the operon
  (`n_reads >= isoform_reads_threshold`, top 100 by reads), as arrows. **Coloured by
  5' end (TSS group), cycled Okabe-Ito palette.** Line width `= clip(0.3 + 0.7*log10(reads), 0.5, 3.0)`;
  alpha `= min(1, 0.45 + 0.5*reads/max_reads)`. Greedy row packing.
- **Depth track** (optional) -- strand-correct PacBio depth (plus/minus bedGraph chosen
  by operon strand), light-blue fill (`#9ecae1`) + line (`#3182bd`); top y-tick = max
  rounded to 1 significant figure. Loaded on demand via an `awk` slice of the bedGraph.

Other conventions:
- **syn3A deletion**: by default, syn3A-deleted regions
  (`Genome_Reduction/aln/raw/syn1_deleted_regions.bed`) are shaded on the gene track
  (`#e8736a`, alpha 0.17) with a small note. Toggle with `MARK_SYN3A_DELETION`.
- **Axes**: bottom x = `Transcript coordinate (nt)` (0 = operon TSS); a secondary top
  axis on the gene track = `Genome position (kb)`. No title (the operon id is the filename).
- **Operon boundary**: marked on all panels by **dashed vertical guide lines only**
  (at tx 0 and tx = operon length); no `TSS`/`TTS` text labels.
- **Canvas**: default `figsize=(7, 7/2)`; panel `height_ratios` genes:isoforms:depth =
  `1.0:2.6:1.1` (or `1.2:2.8` without depth). Tunable per call via `fig_w`/`fig_h`.
  Use `constrained_layout=True` (never `bbox_inches='tight'`).
- **Output**: two PDFs per operon in `Syn1_Operon/operon_plots/` -- `<operon_id>.pdf`
  (genes + isoforms) and `<operon_id>_wdepth.pdf` (+ depth). These exact names are
  consumed by `Genome_Reduction/09,10` and `Operon_Annotation.py`'s long-UTR copies.
- **Shared API**: `OperonCoord` and `draw_isoforms(ax, oc, iso_df, plot_s0, plot_e0)`
  are imported by `Operon_Annotation.py` (its terminator-section plots reuse the isoform
  drawer, so the isoform style stays consistent across the two scripts).