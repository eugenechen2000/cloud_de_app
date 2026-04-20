# Cloud Upload-to-DESeq2 (MVP)

This app lets users upload RNA-seq quantification outputs plus a sample sheet,
then run two-group DESeq2 and download results.

## Supported input modes

1. **Raw count matrix** (`raw_counts_matrix`)  
   - Count matrix CSV: first column `gene_id`, remaining columns are sample IDs (integer/count-like values).  
   - Sample sheet columns: `sample_id,group`.

2. **RSEM outputs** (`rsem_expected_counts`)  
   - Upload many `*.genes.results` files (uses `expected_count`).  
   - Sample sheet columns: `sample_id,file_name,group`.

3. **Salmon/Kallisto estimated counts** (`salmon_kallisto_estimated_counts`)  
   - Upload a **gene-level** estimated-count matrix CSV (same shape as raw count matrix mode).  
   - Sample sheet columns: `sample_id,group`.

4. **TPM/FPKM-only** (`tpm_fpkm_only`)  
   - **Blocked by design** for DESeq2 because these are normalized abundance values, not count-like inputs.

## API endpoints

- `POST /projects/{id}/runs`
- `POST /runs/{id}/upload-rsem`
- `POST /runs/{id}/upload-count-matrix`
- `POST /runs/{id}/start`
- `GET /runs/{id}`
- `GET /runs/{id}/artifacts`
- `GET /runs/{id}/artifacts/{artifact}`

## Quickstart (local)

### 1) Clone and start API

```bash
cd "/Users/eugenechen/Downloads/RNAErnie TP53 plots/cloud_de_app"
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
export RUN_INLINE=true
uvicorn backend.app.main:app --reload --port 8000
```

Then open [http://localhost:8000](http://localhost:8000).

> `RUN_INLINE=true` runs jobs in-process for local smoke tests.

### 2) Run a real public human dataset test

Use the GEO test pack included in this repo:

- `test_packs/public_gse164073_human/count_matrix.csv`
- `test_packs/public_gse164073_human/sample_sheet.csv`

In the UI:

1. Create run (any project ID)
2. Input type: **Raw count matrix**
3. Upload those two files
4. Start run with:
   - `group_a=mock`
   - `group_b=infected`

Expected result: run completes and generates DESeq2 + GSEA + CIBERSORTx artifacts.

## Docker compose (API + worker + Redis)

```bash
docker compose up --build
```

## Outputs per run

Each run writes artifacts under `data/runs/<run_id>/artifacts/`.
All filenames are prefixed with sanitized `project_id` (for example, `soft_demo_...`):

- `<project>_deseq2_all_genes.csv`
- `<project>_deseq2_significant_genes.csv`
- `<project>_top_genes_up.csv`
- `<project>_top_genes_down.csv`
- `<project>_qc_sample_counts.csv`
- `<project>_run_parameters.json`
- `<project>_volcano.png`
- `<project>_volcano.pdf`
- `<project>_ma_plot.png`
- `<project>_ma_plot.pdf`
- `<project>_summary.txt`
- `<project>_results_bundle.zip` (when `zip` binary is available)
- `<project>_gsea_hallmark.csv`, `<project>_gsea_kegg.csv`, `<project>_gsea_reactome.csv`
- `<project>_gsea_*_top.png/.pdf` and `<project>_gsea_summary_top10.csv`
- `<project>_cibersortx_mixture_normalized_counts.csv`
- `<project>_cibersortx_sample_metadata.csv`
- `<project>_cibersortx_instructions.txt`


## Object storage (S3-compatible)

For hosted deployments (separate API + worker), configure object storage so both services share metadata/uploads/artifacts:

- `S3_BUCKET`
- `S3_REGION`
- `S3_ENDPOINT_URL` (optional; required for R2/MinIO)
- `S3_ACCESS_KEY_ID`
- `S3_SECRET_ACCESS_KEY`
- `S3_PREFIX` (optional, default `deapp`)
- `S3_PRESIGN_EXPIRES` (seconds, default 3600)

If `S3_BUCKET` is unset, the app uses local filesystem mode.
