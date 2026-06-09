#!/usr/bin/env python
"""
build_S4_qc.py
==============
Assemble Supplementary Data S4 (RNA-sample quality control): zip the four
facility QC reports together with a small README index. The PDFs are kept as
separate files inside the archive (not merged), renamed to descriptive names.

    python Manuscript/SI/build_S4_qc.py

Output (next to this script, Manuscript/SI/):
    Supplementary_Data_S4_QC.zip   the four QC PDFs + README.txt
"""
import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent          # Manuscript/SI
PROJECT = HERE.parents[1]                        # project root
ZIP_OUT = HERE / "Supplementary_Data_S4_QC.zip"

# (organism, platform label, archive name, source path relative to project root)
REPORTS = [
    ("Syn1",  "PacBio (cDNA / library QC)",  "S4_Syn1_PacBio_cDNA_library_QC.pdf",
     "Syn1_Transcriptomics/PacBio/PacBio_Raw/LutheySchulten cDNA and libraries.pdf"),
    ("Syn1",  "Illumina RNA QC (2023-08)",   "S4_Syn1_Illumina_RNA_QC_2023-08.pdf",
     "Syn1_Transcriptomics/Illumina/Illumina_Raw/Syn1_Illumina_2023_08/Mehta_1RNA.pdf"),
    ("Syn1",  "Illumina RNA QC (2023-09)",   "S4_Syn1_Illumina_RNA_QC_2023-09.pdf",
     "Syn1_Transcriptomics/Illumina/Illumina_Raw/Syn1_Illumina_2023_09/QC_Mehta_9_6_23.pdf"),
    ("Syn3A", "ONT direct-RNA QC",           "S4_Syn3A_ONT_directRNA_QC.pdf",
     "Syn3A_Transcriptomics/ONT/ONT_Raw/ont_qc.pdf"),
]


def build_readme(rows):
    w_org = max(8, max(len(r[0]) for r in rows))
    w_plat = max(8, max(len(r[1]) for r in rows))
    header = f"{'Organism':<{w_org}}  {'Platform':<{w_plat}}  File"
    lines = [
        "Supplementary Data S4 -- RNA-sample quality control",
        "=" * 51,
        "",
        "Raw RNA-sample quality-control reports for the sequencing libraries used",
        "in this study, as provided by the sequencing facilities. They report RNA",
        "concentration and fragment-size integrity (Qubit / Agilent TapeStation /",
        "Bioanalyzer) and, where applicable, ribosomal-RNA depletion (loss of the",
        "16S and 23S peaks). Each report is a separate PDF in this archive:",
        "",
        header,
        "-" * len(header),
    ]
    for organism, platform, arcname, _ in rows:
        lines.append(f"{organism:<{w_org}}  {platform:<{w_plat}}  {arcname}")
    lines += [
        "",
        "Note: Syn3A Illumina reads came from a public dataset (Sandberg et al.)",
        "and therefore have no in-house RNA-sample QC report.",
        "",
    ]
    return "\n".join(lines)


def main():
    missing = [rel for *_, rel in REPORTS if not (PROJECT / rel).is_file()]
    if missing:
        sys.exit("ERROR: missing source PDF(s):\n  " + "\n  ".join(missing))

    readme = build_readme(REPORTS)
    with zipfile.ZipFile(ZIP_OUT, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("README.txt", readme)
        for _, _, arcname, rel in REPORTS:
            z.write(PROJECT / rel, arcname)

    print(f"wrote {ZIP_OUT.relative_to(PROJECT)}  ({len(REPORTS)} PDFs + README.txt)")
    for _, platform, arcname, _ in REPORTS:
        print(f"  {arcname:42s} <- {platform}")


if __name__ == "__main__":
    main()
