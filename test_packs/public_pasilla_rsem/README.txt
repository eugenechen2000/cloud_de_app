Public RSEM-mode test pack (derived from pasilla counts)

Purpose:
- Internal app testing for the RSEM upload path.
- Files are synthetic .genes.results wrappers around public pasilla gene counts.

Source matrix:
https://raw.githubusercontent.com/bioc/pasilla/master/inst/extdata/pasilla_gene_counts.tsv

How generated:
- For each sample column, wrote a tab-delimited file with columns:
  gene_id, expected_count
- Named files: <sample_id>.genes.results

Use in app:
1) Create run
2) Input type: rsem_expected_counts
3) Upload sample_sheet.csv + all *.genes.results files in this folder
4) Start run with group_a=untreated, group_b=treated (or leave blank)

Note:
- These are not full native RSEM outputs; they include only columns required by this app parser.
