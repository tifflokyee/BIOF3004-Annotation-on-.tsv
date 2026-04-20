# Batch Annotation System - Quick Start Guide

## What Was Created

I've created a complete **batch annotation management system** that lets you:

1. ✅ **Track multiple samples** through the annotation pipeline
2. ✅ **Automate batch processing** with a single command
3. ✅ **Monitor progress** with status indicators  
4. ✅ **Resume from checkpoints** if something fails
5. ✅ **Handle the manual VEP step** gracefully

## Files Created

| File | Purpose |
|------|---------|
| `sample_list.tsv` | Manifest file tracking all samples and their status |
| `automation/batch_annotate.py` | Main Python orchestrator |
| `automation/batch_annotate.sh` | Linux/Mac launcher script |
| `automation/batch_annotate.bat` | Windows launcher script |
| `SAMPLE_LIST_README.md` | Full documentation (read this!) |

## Quick Start (3 Steps)

### Step 1: Add Your VCF Files

Place your VCF files in the project root:
```
CAPNGSETB24 for BIOF/
├── your_sample_001.vcf.gz
├── your_sample_002.vcf.gz
├── sample_list.tsv
└── ...
```

### Step 2: Update the Sample List

Edit `sample_list.tsv` and add your samples:

```
sample_id	vcf_file	status
my_sample_001	your_sample_001.vcf.gz	pending
my_sample_002	your_sample_002.vcf.gz	pending
```

### Step 3: Run Batch Annotation

From the `automation/` folder:

**Linux/Mac:**
```bash
./batch_annotate.sh ../sample_list.tsv
```

**Windows:**
```batch
batch_annotate.bat ..\sample_list.tsv
```

**Or directly with Python:**
```bash
python3 batch_annotate.py ../sample_list.tsv
```

## Key Features

### 1. View Sample List
```bash
python3 batch_annotate.py ../sample_list.tsv --list
```

Shows:
- All samples and their current status
- Progress through each annotation step
- Any error messages or notes

### 2. Process Single Sample
```bash
python3 batch_annotate.py ../sample_list.tsv --sample-id sample_001
```

### 3. Resume from a Specific Step
If annotation fails at step X, fix the issue and resume:
```bash
python3 batch_annotate.py ../sample_list.tsv --sample-id sample_001 --start-from clinvar
```

### 4. Reset and Retry
```bash
python3 batch_annotate.py ../sample_list.tsv --reset
```

### 5. Process Failed Samples
```bash
python3 batch_annotate.py ../sample_list.tsv --status failed
```

## Annotation Pipeline

The system runs these steps automatically (in order):

1. **PanelApp** - Panel app annotation
2. **REVEL** - REVEL score annotation
3. **AlphaMissense** - AlphaMissense predictions
4. **ClinVar** - ClinVar annotations
5. **gnomAD** - Population frequency data
6. **VEP** ⚠️ **MANUAL STEP** - Requires user interaction
7. **SpliceAI** - SpliceAI score annotations

### About the VEP Step

Since VEP requires manual upload to Ensembl, the script:
- Pauses at VEP with clear instructions
- Creates the small VCF file for you
- Tells you exactly where to upload it
- Resumes when you provide the VEP output

To include VEP (after you've uploaded and downloaded results):
```bash
python3 batch_annotate.py ../sample_list.tsv --sample-id sample_001 --start-from vep --include-vep
```

## Output Files

Results go into the `result/` folder:

```
result/
├── sample_001_step1_panelapp.tsv
├── sample_001_step2_revel.tsv
├── sample_001_step3_alphamissense.tsv
├── sample_001_step4_clinvar.tsv
├── sample_001_step5_gnomad.tsv           ← Main output before VEP
└── sample_001_step6_final_with_spliceai.tsv  ← Final output with SpliceAI
```

## Status Indicators in sample_list.tsv

The TSV file shows progress with symbols:

- **✓** = Step completed successfully
- **✗** = Step failed
- **→** = Step in progress
- **(empty)** = Not started yet

## Workflow Example

```bash
# 1. Add samples to sample_list.tsv
#    (edit file with your VCF filenames)

# 2. Start batch processing (auto, skips VEP)
python3 batch_annotate.py ../sample_list.tsv

# 3. When it pauses at VEP:
#    - Create small VCF
#    - Upload to https://grch37.ensembl.org/Tools/VEP
#    - Enable SpliceAI
#    - Download and save result

# 4. Check progress
python3 batch_annotate.py ../sample_list.tsv --list

# 5. Resume from VEP after downloading results
python3 batch_annotate.py ../sample_list.tsv --start-from vep --include-vep

# 6. Your final annotated files are in result/
```

## Troubleshooting

| Problem | Solution |
|---------|----------|
| "VCF file not found" | Check file path in sample_list.tsv and verify file exists |
| Step fails | Check error in `notes` column, fix issue, re-run with `--start-from <step>` |
| VEP step stuck | It's waiting for manual interaction - follow instructions in console output |
| Need to retry | Use `--reset` to clear all status marks, then re-run |

## Full Documentation

For detailed documentation including advanced usage, see: **`SAMPLE_LIST_README.md`**

## Command Reference

```bash
# View all options
python3 batch_annotate.py --help

# Process all pending samples
python3 batch_annotate.py ../sample_list.tsv

# List samples
python3 batch_annotate.py ../sample_list.tsv --list

# Reset all samples
python3 batch_annotate.py ../sample_list.tsv --reset

# Process one sample
python3 batch_annotate.py ../sample_list.tsv --sample-id sample_001

# Resume from specific step
python3 batch_annotate.py ../sample_list.tsv --sample-id sample_001 --start-from gnomad

# Reprocess failed samples
python3 batch_annotate.py ../sample_list.tsv --status failed

# Include VEP step (requires manual upload first)
python3 batch_annotate.py ../sample_list.tsv --include-vep
```

---

**Enjoy smooth, automated annotation!** 🧬

For questions or issues, check `SAMPLE_LIST_README.md` or run `--help`.
