# Sample List Management & Batch Annotation Guide

## Overview

The sample list system allows you to:
- **Track multiple samples** through the annotation pipeline
- **Automate batch processing** of VCF files
- **Monitor progress** with status tracking
- **Resume from checkpoints** if errors occur
- **Handle the manual VEP step** gracefully

## Files Included

1. **`sample_list.tsv`** - Tab-separated manifest of all samples and their annotation status
2. **`batch_annotate.py`** - Python script to orchestrate batch annotation
3. **`sample_list_template.txt`** - This documentation file

## Quick Start

### Step 1: Add Your Samples to the Sample List

Edit `sample_list.tsv` and add your VCF files. For example:

```
sample_id	vcf_file	status	panelapp_result	revel_result	alphamissense_result	clinvar_result	gnomad_result	vep_input	vep_output	spliceai_result	final_result	notes
sample_001	your_sample_001.vcf.gz	pending						
sample_002	your_sample_002.vcf.gz	pending						
```

**Important:** Place your VCF files in the project root directory (same level as `sample_list.tsv`).

### Step 2: Run Batch Annotation

From the `automation/` folder, run:

```bash
# Process all pending samples (automated, skips manual VEP)
python batch_annotate.py ../sample_list.tsv

# Process single sample
python batch_annotate.py ../sample_list.tsv --sample-id sample_001

# Show current progress
python batch_annotate.py ../sample_list.tsv --list
```

### Step 3: Handle the Manual VEP Step

When the script reaches the VEP step, it will pause and give you instructions:

1. Create small VCF for VEP
2. Upload to Ensembl GRCh37 VEP (https://grch37.ensembl.org/Tools/VEP)
   - Enable **SpliceAI** during submission
3. Download the VEP output file
4. Save it with your sample ID in the project root: `sample_001_vep_output.vcf`
5. Re-run the batch script with `--include-vep` to complete

Or, continue with other samples and come back to VEP later.

## Sample List Columns

| Column | Purpose | Example |
|--------|---------|---------|
| `sample_id` | Unique identifier for the sample | `sample_001`, `patient_A` |
| `vcf_file` | Path to input VCF (relative to project root) | `your_sample.vcf.gz` |
| `status` | Current processing status | `pending`, `in_progress`, `completed`, `failed` |
| `panelapp_result` | Status of PanelApp annotation | `✓`, `✗`, `→`, or blank |
| `revel_result` | Status of REVEL annotation | `✓`, `✗`, `→`, or blank |
| `alphamissense_result` | Status of AlphaMissense annotation | `✓`, `✗`, `→`, or blank |
| `clinvar_result` | Status of ClinVar annotation | `✓`, `✗`, `→`, or blank |
| `gnomad_result` | Status of gnomAD annotation | `✓`, `✗`, `→`, or blank |
| `vep_input` | Status of VEP input preparation | `✓`, `✗`, `→`, or blank |
| `vep_output` | Path to VEP output file (if provided) | `sample_001_vep_output.vcf` |
| `spliceai_result` | Status of SpliceAI annotation | `✓`, `✗`, `→`, or blank |
| `final_result` | Path to final annotated result | `result/sample_001_step6_...` |
| `notes` | Additional information or error messages | Error details, timestamps |

## Usage Examples

### Example 1: Basic Batch Annotation

```bash
cd automation
python batch_annotate.py ../sample_list.tsv
```

This will:
- Process all samples with status = "pending"
- Run: PanelApp → REVEL → AlphaMissense → ClinVar → gnomAD
- Skip VEP and SpliceAI (manual steps)
- Update `sample_list.tsv` with progress

### Example 2: Process Single Sample

```bash
python batch_annotate.py ../sample_list.tsv --sample-id sample_001
```

### Example 3: Restart from a Specific Step

If a sample failed at the ClinVar step, fix the issue and resume:

```bash
python batch_annotate.py ../sample_list.tsv --sample-id sample_001 --start-from clinvar
```

**Available steps:**
- `panelapp`
- `revel`
- `alphamissense`
- `clinvar`
- `gnomad`
- `vep` (requires `--include-vep` flag)
- `spliceai`

### Example 4: View Sample List

```bash
python batch_annotate.py ../sample_list.tsv --list
```

Shows a formatted table of all samples and their current status.

### Example 5: Reset and Retry

If something goes wrong, reset all samples to pending:

```bash
python batch_annotate.py ../sample_list.tsv --reset
```

Then re-run the batch annotation.

### Example 6: Process Failed Samples

```bash
python batch_annotate.py ../sample_list.tsv --status failed
```

This will re-process only samples that previously failed.

### Example 7: Include VEP and SpliceAI

After you've uploaded your files to VEP and downloaded the results:

```bash
python batch_annotate.py ../sample_list.tsv --sample-id sample_001 --start-from vep --include-vep
```

**Important:** Save your VEP output with the naming pattern:
- `sample_001_vep_output.vcf`
- The script will auto-detect this file

## Output Files

After processing, check the `result/` folder for:

```
result/
├── sample_001_step1_panelapp.tsv
├── sample_001_step2_revel.tsv
├── sample_001_step3_alphamissense.tsv
├── sample_001_step4_clinvar.tsv
├── sample_001_step5_gnomad.tsv
├── sample_001_vep_input.vcf              # For manual VEP submission
├── sample_001_step6_final_with_spliceai.tsv    # Final result
└── ...
```

Your main output files are:
- **After gnomAD:** `sample_001_step5_gnomad.tsv`
- **Final (with SpliceAI):** `sample_001_step6_final_with_spliceai.tsv`

## Troubleshooting

### Problem: "VCF file not found"

**Solution:** Make sure your VCF file path in `sample_list.tsv` is correct and the file is in the project root folder.

```bash
ls -la *.vcf*    # Check if your VCF files are here
```

### Problem: Script fails at a specific step

1. Check the error message in the `notes` column of `sample_list.tsv`
2. Verify the required annotation files exist in `annotation/` folder
3. Fix the issue
4. Re-run from the failed step:

```bash
python batch_annotate.py ../sample_list.tsv --sample-id sample_001 --start-from <failed_step>
```

### Problem: VEP step is stuck

The VEP step is **manual**. The script will show you instructions. You must:

1. Create the small VCF file
2. Upload to https://grch37.ensembl.org/Tools/VEP manually
3. Enable SpliceAI
4. Download and save the result
5. Re-run the script with `--include-vep`

### Problem: Memory errors or timeouts

Some steps may timeout on large VCF files. If this happens:

1. Re-run from the failed step
2. Or, process one sample at a time instead of batch

## Advanced Usage

### Modify the Script

Edit `batch_annotate.py` to:
- Change timeout duration (default: 1 hour)
- Add custom filtering logic
- Integrate with external databases
- Send notifications on completion

The script uses the existing Python scripts from `automation/`, so all customization options of those scripts are available.

### Integration with CI/CD

You can run batch annotation in a scheduled job:

```bash
# crontab example: run every day at 9 AM
0 9 * * * cd /path/to/CAPNGSETB24\ for\ BIOF/automation && python batch_annotate.py ../sample_list.tsv >> annotation.log 2>&1
```

## Tips for Smooth Workflow

1. **Organize your VCF files:** Keep all VCF files in the project root
2. **Name consistently:** Use consistent sample IDs and filenames
3. **Monitor progress:** Regularly check `sample_list.tsv` with `--list` flag
4. **Backup results:** Copy the `result/` folder periodically
5. **Document notes:** Add notes to the `notes` column for your reference
6. **Batch VEP submissions:** Group VEP submissions to save time

## Summary

| Task | Command |
|------|---------|
| Add samples | Edit `sample_list.tsv` |
| Start processing | `python batch_annotate.py ../sample_list.tsv` |
| Check progress | `python batch_annotate.py ../sample_list.tsv --list` |
| Process one sample | `python batch_annotate.py ../sample_list.tsv --sample-id sample_001` |
| Resume from step | `python batch_annotate.py ../sample_list.tsv --sample-id sample_001 --start-from <step>` |
| Process failed | `python batch_annotate.py ../sample_list.tsv --status failed` |
| Reset all | `python batch_annotate.py ../sample_list.tsv --reset` |

Enjoy automated annotation! 🧬
