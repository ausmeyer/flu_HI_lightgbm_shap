# Manuscript snapshot (created 18 September; updated 19 September 2026)

This folder is a **dated copy** of the H3N2 LightGBM+SHAP HI manuscript and supplement for the public analysis repository.

## Which version is authoritative?

**Overleaf is the writing workspace.** Edit there:

- Main text: `69b778fd3e7b181fe1c2943b/Research_report_ve.tex`
- Supplement: `69b7797c762f515edcff3ad6/research_report_supplement.tex`

This snapshot is a GitHub-facing copy of those files plus our generated figures and tables. If the Overleaf sources and this folder disagree, trust Overleaf.

## What is not included

Copyrighted publisher PDFs and other full texts in `all_citations/` (and `relevant_literature/*.pdf`) are **gitignored** and must never be staged, committed, or pushed. They are not part of this snapshot.

Compiled manuscript and supplement PDFs are omitted from this snapshot (gitignored at the snapshot root) so they are not committed. Figure PDFs in `figures/` remain. Overleaf is the compile workspace.

GISAID FASTA files, raw downloads, and rich internal metadata exports are excluded from the current public package. A separate identifier-only WIC inventory is in `H3N2-WIC/publication/` at the repository root. Complete internal inputs remain local; earlier Git history still contains previously tracked metadata. A study-specific GISAID acknowledgment artifact remains a submission item.

Tables S4/S5 are included as accompanying TSV files in `tables/`. The five additional exploratory TSVs retained in the SI Overleaf project's `supplementary_data/` directory are not part of this public snapshot and are not needed for compilation.

The 19 September update applies the Cursor audit corrections to captions, citation scope, manual GISAID retrieval wording, Table S2 layout, and package documentation. S18's source-PDF page bounds were expanded without changing its drawing content. S19 was synchronized from the existing recorded final-run figure whose input summary matches the current result table; no analysis was rerun.

## Compile

From this folder, after placing `figures/` and `tables/` next to the `.tex` files as in the Overleaf layouts:

- Main: `latexmk -pdf Research_report_ve.tex`
- Supplement: `latexmk -pdf research_report_supplement.tex`

`latexmk` runs BibTeX and repeats PDFLaTeX as needed to resolve citations and cross-references.

The supplement `\input{tables/...}` and `\includegraphics{figures/...}` paths match the Overleaf SI layout.
