# Local SpliceAI Model Run

This project now has a local SpliceAI wrapper:

```powershell
automation\run_local_spliceai.bat result\step5_gnomad.tsv --reference C:\path\to\hg19.fa -o result\step6_final_with_local_spliceai.tsv
```

The script does three things:

1. Creates a small VCF from clinically relevant TSV rows.
2. Runs the local `spliceai` command.
3. Merges the SpliceAI result back into the TSV as `SpliceAI_max` and `SpliceAI_Impact`.

## Requirements

You need a separate environment with the SpliceAI command available. The default Python in this workspace is Python 3.14, which is not a good target for TensorFlow-based SpliceAI.

Local SpliceAI model inference does not run on Windows XP. The batch launchers avoid newer command-line conveniences where possible, but the model step requires a modern Python/TensorFlow runtime.

On Windows XP, `run_stepwise_annotation.bat` can still finish if an existing local SpliceAI VCF is already present at:

```text
result\local_spliceai_model_output.vcf
```

If the VCF is missing on XP, the automation creates this handoff file:

```text
result\local_spliceai_input_for_external_run.vcf
```

Copy that input VCF to Windows 10/11 or Linux, run local SpliceAI there, and save the output as:

```text
local_spliceai_model_output.vcf
```

Then copy it back to the XP computer as:

```text
result\local_spliceai_model_output.vcf
```

You can then rerun `run_stepwise_annotation.bat`, or merge only the SpliceAI step without rerunning steps 1-5:

```bat
automation\merge_local_spliceai_xp.bat result\step5_gnomad.tsv result\local_spliceai_model_output.vcf result\step6_final_with_local_spliceai.tsv
```

On a supported system, you can also set `SPLICEAI_PYTHON` to a compatible Python command:

```bat
set SPLICEAI_PYTHON=C:\Path\To\Python313\python.exe
automation\run_stepwise_annotation.bat input.tsv
```

Recommended environment setup:

```bash
conda create -n spliceai -c bioconda -c conda-forge spliceai
conda activate spliceai
spliceai -h
```

You also need an uncompressed GRCh37/hg19 FASTA file, for example:

```text
C:\references\hg19.fa
```

This pipeline is hg19/GRCh37-style, so the wrapper defaults to:

```text
--annotation grch37
```

## Run With Conda Prefix

If `spliceai` is only available inside a conda environment, run:

```powershell
automation\run_local_spliceai.bat result\step5_gnomad.tsv --reference C:\references\hg19.fa --runner-prefix "conda run -n spliceai" -o result\step6_final_with_local_spliceai.tsv
```

## Direct Python Form

```powershell
python automation\run_local_spliceai.py result\step5_gnomad.tsv --reference C:\references\hg19.fa -o result\step6_final_with_local_spliceai.tsv
```

Useful options:

```text
--distance 50
--mask 0
--spliceai-exe spliceai
--spliceai-vcf result\local_spliceai_output.vcf
--keep-input-vcf
```

The local SpliceAI VCF is preserved by default. The temporary input VCF is deleted unless `--keep-input-vcf` is used.

## Actual Windows Method Used Here

This machine did not have conda, WSL, Docker, or a local hg19 FASTA. The actual completed run used:

```powershell
automation\run_local_spliceai_ucsc.bat result\variants_with_AlphaMissense_and_ClinVar_and_gnomAD.tsv -o result\step6_final_with_local_spliceai.tsv --spliceai-vcf result\local_spliceai_model_output.vcf
```

That runner:

1. Uses Python 3.13.
2. Loads the local SpliceAI model files installed by the `spliceai` package.
3. Fetches only the required hg19 sequence windows from UCSC and caches them in `annotation\hg19_spliceai_windows_cache.json`.
4. Runs model inference locally.
5. Writes native SpliceAI output to `result\local_spliceai_model_output.vcf`.
6. Merges `SpliceAI_max` and `SpliceAI_Impact` into `result\step6_final_with_local_spliceai.tsv`.
