# Manuscript snapshot (18 September 2026)

This folder is a **dated copy** of the H3N2 LightGBM+SHAP HI manuscript and supplement for the public analysis repository.

## Which version is authoritative?

**Overleaf is the writing workspace.** Edit there:

- Main text: `69b778fd3e7b181fe1c2943b/Research_report_ve.tex`
- Supplement: `69b7797c762f515edcff3ad6/research_report_supplement.tex`

This snapshot is a GitHub-facing copy of those files plus our generated figures and tables. If the Overleaf sources and this folder disagree, trust Overleaf.

## What is not included

Copyrighted publisher PDFs and other full texts in `all_citations/` (and `relevant_literature/*.pdf`) are **gitignored** and must never be staged, committed, or pushed. They are not part of this snapshot.

Compiled manuscript and supplement PDFs are omitted from this snapshot (gitignored at the snapshot root) so they are not committed. Figure PDFs in `figures/` remain. Overleaf is the compile workspace.

GISAID FASTA files and raw GISAID metadata downloads are excluded from this snapshot. The raw metadata directory is ignored by the analysis repository; earlier Git history may still contain previously tracked metadata.

## Compile

From this folder, after placing `figures/` and `tables/` next to the `.tex` files as in the Overleaf layouts:

- Main: `latexmk -pdf Research_report_ve.tex`
- Supplement: `latexmk -pdf research_report_supplement.tex`

`latexmk` runs BibTeX and repeats PDFLaTeX as needed to resolve citations and cross-references.

The supplement `\input{tables/...}` and `\includegraphics{figures/...}` paths match the Overleaf SI layout.
