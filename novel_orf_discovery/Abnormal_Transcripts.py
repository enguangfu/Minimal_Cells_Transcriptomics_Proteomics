
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

ISOFORM_XLSX  = '../isoform_annotation/isoform_clusters_annotated.xlsx'
MIN_READS     = 10


PERCENTILES = [0, 25, 50, 75, 90, 95, 96, 97, 98, 99, 100]


def _weighted_median(vals: np.ndarray, w: np.ndarray) -> float:
    """Weighted median: smallest value v such that cumulative weight >= 0.5 * total."""
    order = np.argsort(vals)
    vals_s, w_s = vals[order], w[order]
    cumw = np.cumsum(w_s)
    return float(vals_s[np.searchsorted(cumw, 0.5 * cumw[-1])])


def _weighted_percentile(vals: np.ndarray, w: np.ndarray, pct: float) -> float:
    """Weighted percentile (0–100) via sorted cumulative weight."""
    order = np.argsort(vals)
    vals_s, w_s = vals[order], w[order]
    cumw = np.cumsum(w_s)
    return float(vals_s[np.searchsorted(cumw, (pct / 100) * cumw[-1])])


def _print_percentiles(vals: np.ndarray, w, field: str, subtitle: str):
    """Print a percentile table to stdout."""
    print(f'\n  {field}  [{subtitle}]')
    header = '  ' + ''.join(f'  p{p:<5}' for p in PERCENTILES)
    print(header)
    if w is None:
        row = '  ' + ''.join(f'  {np.percentile(vals, p):<7.4f}' for p in PERCENTILES)
    else:
        row = '  ' + ''.join(
            f'  {_weighted_percentile(vals, w, p):<7.4f}' for p in PERCENTILES)
    print(row)


def _plot_one_frac(vals, weights, field, color, ylabel, subtitle,
                   n_total, min_reads, out_pdf):
    """
    Draw and save a single figure for one (field, weighting) combination.

    vals    : Series of frac values for all isoforms (including zeros).
    weights : None for isoform-count histogram; Series of n_reads for weighted.

    Markers drawn on the plot:
      - Unweighted median (dashed black line, always shown).
      - Weighted median   (dashed gray  line, only when weights are provided).
    Percentile table is printed to stdout for every call.
    """
    bins = np.linspace(0, 1, 51)   # 50 equal bins; first bin [0, 0.02] captures zeros

    v = vals.values.astype(float)
    w = weights.loc[vals.index].astype(float).values if weights is not None else None

    med_unweighted = float(np.median(v))

    fig, ax = plt.subplots(figsize=(6, 3))
    # fig.suptitle(
    #     f'{field}  —  {subtitle}\n'
    #     f'(n_reads > {min_reads},  n = {n_total:,})',
    #     fontsize=10, fontweight='bold')

    ax.hist(v, bins=bins, weights=w, color=color, alpha=0.75,
            edgecolor='white', lw=0.3)
    ax.set_yscale('log')
    ax.set_xlim(-0.01, 1.01)
    ax.set_xlabel(field, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=9)

    # Always show the unweighted median
    ax.axvline(med_unweighted, color='black', lw=1.2, ls='--',
               label=f'median (Isoform kinds) = {med_unweighted:.3f}')

    # Weighted median only when weights are provided
    if w is not None:
        med_weighted = _weighted_median(v, w)
        ax.axvline(med_weighted, color='dimgray', lw=1.2, ls=':',
                   label=f'median (weighted by count)   = {med_weighted:.3f}')

    ax.legend(fontsize=8, framealpha=0.8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f'{int(x):,}'))

    plt.tight_layout()
    fig.savefig(out_pdf, dpi=300, bbox_inches='tight')
    # plt.show()
    print(f'Saved: {out_pdf}')

    # Percentile report
    _print_percentiles(v, w, field, subtitle)


def plot_frac_distributions(
        xlsx_path: str = ISOFORM_XLSX,
        min_reads: int = MIN_READS):
    """
    Produce four separate figures (each saved to its own PDF):
      1. frac_antisense   — isoform count
      2. frac_intergenic  — isoform count
      3. frac_antisense   — read-count weighted
      4. frac_intergenic  — read-count weighted

    Zeros are included in the first histogram bin [0, 0.02].
    The median (over all isoforms including zeros) is marked.
    """
    df  = pd.read_excel(xlsx_path)
    iso = df[df['n_reads'] > min_reads].copy()
    n_total = len(iso)
    
    print(f'Total isoforms with n_reads > {min_reads}: {n_total:,}')
    
    specs = [
        # (field,             color,     weights,          ylabel,                        subtitle,              out_pdf)
        ('frac_antisense',  '#d62728', None,              'Unique isoforms (log)',     'Isoform count',       'frac_antisense_count.pdf'),
        ('frac_intergenic', '#1f77b4', None,              'Unique isoforms (log)',     'Isoform count',       'frac_intergenic_count.pdf'),
        ('frac_antisense',  '#d62728', iso['n_reads'],    'Isoform read count (log)',      'Read-count weighted', 'frac_antisense_weighted.pdf'),
        ('frac_intergenic', '#1f77b4', iso['n_reads'],    'Isoform read count (log)',      'Read-count weighted', 'frac_intergenic_weighted.pdf'),
    ]

    for field, color, weights, ylabel, subtitle, out_pdf in specs:
        vals = iso[field].dropna()
        _plot_one_frac(vals, weights, field, color, ylabel,
                       subtitle, n_total, min_reads, out_pdf)


if __name__ == '__main__':
    plot_frac_distributions()
