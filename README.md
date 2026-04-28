

# BIOF3004 Annotation Pipeline

This repository provides a local annotation pipeline for `.tsv` / `.vcf` variant files.

Main script:

- `automation/auto_annotate_generated.py`

It supports:

- single-file mode
- batch mode (`--all-tsv`)
- HGMD output generation
- PanelApp / REVEL / AlphaMissense / ClinVar / gnomAD / local SpliceAI annotations

## Folder Structure

 Project layout:

```text
CAPNGSETB24 for BIOF/
├── annotation/
│   ├── Mendeliome.tsv
│   ├── Incidentalome.tsv
│   ├── Additional findings_Paediatric.tsv
│   ├── revel_with_transcript_ids
│   ├── AlphaMissense_hg19.tsv.gz
│   ├── clinvar.vcf.gz
│   ├── gnomad*.vcf.gz
│   └── hg19_spliceai_windows_cache.json (auto-created/updated)
├── automation/
│   └── auto_annotate_generated.py
├── sample/                  # your input TSV files (example)
└── result/                  # output files
```

You can keep `automation/` and your sample folders separate.  
Run from the project root and point to inputs with `--input-dir`.

## Input Requirements

For TSV input, include either:

- `Chr`, `Coordinate`, `ref`, `alt`

or:

- `Chr`, `Coordinate`, `Variant` (for example `C>T` or `C>C/T`)

VCF input is also supported (`.vcf`, `.vcf.gz`, `.vcf.bgz`).

## Environment Setup

Recommended: Python 3.9+ in a conda environment.

```bash
conda create -n biof_annotation python=3.9
conda activate biof_annotation
```

Install core dependencies:

```bash
conda config –-add channels defaults 
conda config –-add channels conda-forge 
conda config –-add channels bioconda 
conda install tensorflow-cpu keras 
pip install spliceai --no-deps
```

Notes:

- `spliceai`, `tensorflow`/`keras` are required for the local SpliceAI step.
- If you skip SpliceAI (`--skip-spliceai`), those heavy dependencies are not needed for that run.

## How To Run

Run commands from repository root.

### 1) Batch mode (all TSV files in a folder)

```bash
python3 automation/auto_annotate_generated.py --all-tsv --input-dir "sample"
```

### 2) Single file mode

```bash
python3 automation/auto_annotate_generated.py "sample/your_file.tsv"
```

## Outputs

For each input file `<name>.tsv`, the pipeline now generates 2 outputs:

1. Main annotated table:
  - `result/<name>.tsv`
2. Separate HGMD upload file:
  - `result/<name>_hgmd.txt`
The separate HGMD file is written early (after HGMD step), so you can upload it to the HGMD website while the remaining annotation steps are still running.

## Optional Flags

Skip individual steps if needed:

```bash
--skip-hgmd
--skip-panelapp
--skip-revel
--skip-alphamissense
--skip-clinvar
--skip-gnomad
--skip-spliceai
```

Useful performance options:

```bash
--revel-chunk-size 500000
--alphamissense-chunk-size 500000
```

SpliceAI options:

```bash
--spliceai-distance 50
--spliceai-mask 0
--local-spliceai-cache "annotation/hg19_spliceai_windows_cache.json"
```

## Example

Process all TSV files in `sample/` and skip gnomAD:

```bash
python3 automation/auto_annotate_generated.py --all-tsv --input-dir "sample" --skip-gnomad
```

## Troubleshooting

- `ModuleNotFoundError: No module named 'pandas'`
  - Install dependencies in the active environment (`pip install pandas`).
- Very slow REVEL step
  - This is expected on large files; tune chunk size or run on a stronger machine.
- No output file generated
  - Confirm input folder path and TSV files exist.
  - Ensure you run command from repo root.

