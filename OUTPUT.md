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