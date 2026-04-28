# BIOF3004 Annotation Pipeline

This repository provides a **local variant annotation pipeline** for `.tsv` and `.vcf` files from NGS data.

**Supported annotations**

- PanelApp (Mendeliome, Incidentalome, Paediatric)
- REVEL
- AlphaMissense
- ClinVar
- gnomAD
- SpliceAI

---

## Folder Structure

```text
CAPNGSETB24 for BIOF/
├── sample/                   # Input folder (example) 
├── result/                   # Output files (auto-created)
├── annotation/               # Reference annotation files
├── automation/
│   ├── auto_annotate_generated.py
│   └── download_gnomad.py
├── run_annotation_windows.bat
└── README.md
```

---

## Files Not Included in GitHub

Some large annotation data files are intentionally not pushed to GitHub (size/storage reasons).

You must download/provide them manually in `annotation/` before running the full pipeline.

Expected large files:

- `annotation/revel_with_transcript_ids`
- `annotation/AlphaMissense_hg19.tsv.gz`
- `annotation/clinvar.vcf.gz`
- `annotation/gnomad*.vcf.gz` or `annotation/gnomad*.vcf.bgz`
- `annotation/hg19.fa` (reference FASTA)

Small PanelApp files are usually included:

- `annotation/Mendeliome.tsv`
- `annotation/Incidentalome.tsv`
- `annotation/Additional findings_Paediatric.tsv`

---

## How To Download Annotation Files

Large annotation files are intentionally **not included** in GitHub.

Required files in `annotation/`:


| Required file(s)            | How to get                                                                           |
| --------------------------- | ------------------------------------------------------------------------------------ |
| `gnomad*.vcf.bgz` + `.tbi`  | `python3 automation/download_gnomad.py --mode legacy_grch37 --output-dir annotation` |
| `revel_with_transcript_ids` | Download from REVEL source or copy from lab storage                                  |
| `AlphaMissense_hg19.tsv.gz` | Download from `dm_alphamissense` storage                                             |
| `clinvar.vcf.gz`            | Download from NCBI ClinVar FTP                                                       |


Official sources:

- gnomAD releases: [gnomAD Downloads](https://gnomad.broadinstitute.org/downloads)
- ClinVar VCF: [NCBI ClinVar FTP](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/)
- AlphaMissense: [Google Cloud Storage Browser (dm_alphamissense)](https://console.cloud.google.com/storage/browser/dm_alphamissense;tab=objects?pli=1&prefix=&forceOnObjectsSortingFiltering=false)
- REVEL: [REVEL Scores](https://sites.google.com/site/revelgenomics/downloads)

---

## Quick Start

### 1) Environment setup (first time only)

```bash
conda create -n biof_annotation python=3.9
conda activate biof_annotation

conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda

conda install tensorflow-cpu keras
pip install spliceai --no-deps
pip install pandas numpy requests pysam pyfaidx
```

### 2) Prepare annotation files

Download gnomAD GRCh37-compatible file:

```bash
python3 automation/download_gnomad.py --mode legacy_grch37 --output-dir annotation
```

### 3) Run the pipeline

All TSV files in a folder:

```bash
python3 automation/auto_annotate_generated.py --all-tsv --input-dir "sample"
```

Single file:

```bash
python3 automation/auto_annotate_generated.py "sample/your_file.tsv"
```

---

## Windows (Anaconda Prompt)

Default launcher behavior (`run_annotation_windows.bat`):

- Activates conda environment: `biof_annotation`
- Runs in **batch mode** on the `sample` folder
- Executes this command:
  - `python3 automation/auto_annotate_generated.py --all-tsv --input-dir "sample" --skip-gnomad`

1. Open **Anaconda Prompt**
2. Go to repository root:

```bat
cd /d "C:\path\to\CAPNGSETB24 for BIOF"
```

1. Run launcher:

```bat
run_annotation_windows.bat
```

The launcher runs immediately with the default command above (no prompt for input directory).

---

## Outputs

For each input file `<name>.tsv`, pipeline generates:

1. Main annotated table:
  - `result/<name>.tsv`
2. Separate HGMD upload file:
  - `result/<name>_hgmd.txt`

`HGMD_input` is also included as the first column in main output.  
The separate HGMD file is written early, so you can upload to HGMD while the remaining steps continue.

---

## Optional Flags

```bash
--skip-hgmd
--skip-panelapp
--skip-revel
--skip-alphamissense
--skip-clinvar
--skip-gnomad
--skip-spliceai
```

Batch mode:

```bash
--all-tsv --input-dir input
```

---

## Troubleshooting

- `SpliceAI executable was not found on PATH: spliceai`
  - Activate correct conda environment and verify with `which spliceai`.
- `ModuleNotFoundError: No module named 'pysam'`
  - Install missing runtime dependencies: `pip install pysam pyfaidx`.
- Very slow REVEL step
  - Expected for large files; use stronger machine or adjust chunk sizes.
- No output file generated
  - Confirm input path and run from repository root.

