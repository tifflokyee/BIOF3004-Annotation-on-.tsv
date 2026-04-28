# BIOF3004 Annotation Pipeline

This repository provides a **local variant annotation pipeline** for `.tsv` and `.vcf` files from NGS data.

**Supported annotations**

- PanelApp (Mendeliome, Incidentalome, Paediatric)
- REVEL
- AlphaMissense
- ClinVar
- gnomAD
- Local SpliceAI

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
| `hg19.fa` + `hg19.fa.fai`   | Download hg19 FASTA from UCSC, then run `samtools faidx annotation/hg19.fa`          |


Official sources:

- gnomAD releases: [gnomAD Downloads](https://gnomad.broadinstitute.org/downloads)
- ClinVar VCF: [NCBI ClinVar FTP](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/)
- AlphaMissense: [Google Cloud Storage Browser (dm_alphamissense)](https://console.cloud.google.com/storage/browser/dm_alphamissense;tab=objects?pli=1&prefix=&forceOnObjectsSortingFiltering=false)
- REVEL: [REVEL Scores](https://sites.google.com/site/revelgenomics/downloads)
- hg19 FASTA: [UCSC hg19 FASTA](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz)

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

Download local hg19 FASTA for SpliceAI CLI mode:

```bash
wget -O annotation/hg19.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip -f annotation/hg19.fa.gz
samtools faidx annotation/hg19.fa
```

### 3) Run the pipeline

All TSV files in a folder:

```bash
python3 automation/auto_annotate_generated.py --all-tsv --input-dir "sample" --local-spliceai-reference "annotation/hg19.fa"
```

Single file:

```bash
python3 automation/auto_annotate_generated.py "sample/your_file.tsv" --local-spliceai-reference "annotation/hg19.fa"
```

---

## Windows (Anaconda Prompt)

Default launcher behavior (`run_annotation_windows.bat`):

- Activates conda environment: `biof_annotation`
- Uses local SpliceAI FASTA path: `annotation\hg19.fa` (must exist)
- Runs in **batch format** (`--all-tsv`) on the selected input directory
- Includes by default:
  - `--all-tsv`
  - `--input-dir "<your_folder>"`
  - `--skip-gnomad`
  - `--local-spliceai-reference "annotation\hg19.fa"`



1. Open **Anaconda Prompt**
2. Go to repository root:

```bat
cd /d "C:\path\to\CAPNGSETB24 for BIOF"
```

1. Run launcher:

```bat
run_annotation_windows.bat
```

1. At prompt:

`Enter input directory containing TSV files [default: sample]:`

Type your folder name (example: `dragon`) or press Enter for `sample`.

---

## SpliceAI (2 Plans)

### Plan A (recommended): CLI mode with local FASTA

Use this for stable, local reference-based execution.

```bash
python3 automation/auto_annotate_generated.py --all-tsv --input-dir "sample" --local-spliceai-reference "annotation/hg19.fa"
```

### Plan B: UCSC mode (no local FASTA argument)

If `--local-spliceai-reference` is omitted, pipeline uses UCSC runner mode.

```bash
python3 automation/auto_annotate_generated.py --all-tsv --input-dir "sample"
```

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

- `FileNotFoundError: SPLICEAI_REFERENCE does not exist`
  - Check `--local-spliceai-reference` path.
  - If you only have `annotation/hg19.fa.gz`, decompress to `annotation/hg19.fa`.
- `SpliceAI executable was not found on PATH: spliceai`
  - Activate correct conda environment and verify with `which spliceai`.
- `ModuleNotFoundError: No module named 'pysam'`
  - Install missing runtime dependencies: `pip install pysam pyfaidx`.
- Very slow REVEL step
  - Expected for large files; use stronger machine or adjust chunk sizes.
- No output file generated
  - Confirm input path and run from repository root.

