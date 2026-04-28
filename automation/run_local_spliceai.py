import argparse
import shlex
import shutil
import subprocess
import tempfile
from pathlib import Path

from add_spliceai_from_vep import add_spliceai_to_tsv, build_spliceai_lookup
from create_small_vcf_for_vep import create_small_vcf


ROOT = Path(__file__).resolve().parent.parent


def default_output_paths(input_tsv):
    input_tsv = Path(input_tsv)
    name = input_tsv.name
    for suffix in (".tsv.gz", ".tsv"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    work_dir = ROOT / "result"
    return {
        "input_vcf": work_dir / f"{name}.local_spliceai.input.vcf",
        "spliceai_vcf": work_dir / f"{name}.local_spliceai.output.vcf",
        "output_tsv": work_dir / f"{name}_with_local_SpliceAI.tsv",
    }


def split_prefix(value):
    if not value:
        return []
    return shlex.split(value, posix=False)


def build_command(args, input_vcf, output_vcf):
    command = split_prefix(args.runner_prefix)
    command.extend(
        [
            args.spliceai_exe,
            "-I",
            str(input_vcf),
            "-O",
            str(output_vcf),
            "-R",
            str(args.reference),
            "-A",
            args.annotation,
            "-D",
            str(args.distance),
            "-M",
            str(args.mask),
        ]
    )
    return command


def validate_inputs(args):
    input_tsv = Path(args.input_tsv)
    if not input_tsv.exists():
        raise FileNotFoundError(f"Input TSV not found: {input_tsv}")

    reference = Path(args.reference)
    if not reference.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {reference}")

    if not args.runner_prefix and shutil.which(args.spliceai_exe) is None:
        raise FileNotFoundError(
            f"SpliceAI executable was not found on PATH: {args.spliceai_exe}. "
            "Install SpliceAI in a separate environment or pass --spliceai-exe."
        )


def detect_reference_chrom_style(reference_path):
    reference_path = Path(reference_path)
    fai_path = reference_path.with_suffix(reference_path.suffix + ".fai")
    contigs = set()

    if fai_path.exists():
        with open(fai_path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if parts and parts[0]:
                    contigs.add(parts[0].strip())

    if not contigs:
        return None

    has_chr = any(name == "chr1" for name in contigs)
    has_plain = any(name == "1" for name in contigs)
    if has_chr and not has_plain:
        return "chr"
    if has_plain and not has_chr:
        return "plain"
    return None


def normalize_chrom_for_reference(chrom, style):
    chrom = str(chrom).strip()
    chrom_no_chr = chrom[3:] if chrom.lower().startswith("chr") else chrom
    if chrom_no_chr == "M":
        chrom_no_chr = "MT"
    if chrom_no_chr == "MT" and style == "chr":
        return "chrM"
    if style == "chr":
        return f"chr{chrom_no_chr}"
    if style == "plain":
        return "MT" if chrom_no_chr == "M" else chrom_no_chr
    return chrom


def rewrite_vcf_chromosomes(input_vcf, reference_path):
    style = detect_reference_chrom_style(reference_path)
    if style is None:
        return input_vcf, None

    input_vcf = Path(input_vcf)
    changed = False
    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".vcf",
        prefix=f".{input_vcf.stem}.contig_fix.",
        dir=input_vcf.parent,
        delete=False,
        encoding="utf-8",
    ) as temp_handle:
        rewritten_vcf = Path(temp_handle.name)
        with open(input_vcf, "rt", encoding="utf-8", errors="replace") as src:
            for line in src:
                if line.startswith("#"):
                    temp_handle.write(line)
                    continue
                fields = line.rstrip("\n").split("\t")
                if not fields:
                    temp_handle.write(line)
                    continue
                updated = normalize_chrom_for_reference(fields[0], style)
                if updated != fields[0]:
                    fields[0] = updated
                    changed = True
                temp_handle.write("\t".join(fields) + "\n")

    if not changed:
        rewritten_vcf.unlink(missing_ok=True)
        return input_vcf, None

    print(f"Adjusted VCF chromosome naming to match reference style: {style}")
    return rewritten_vcf, style


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Run the local SpliceAI model on clinically relevant variants from a TSV, "
            "then add SpliceAI_max and SpliceAI_Impact columns to the TSV."
        )
    )
    parser.add_argument("input_tsv", help="Annotated TSV to score with local SpliceAI")
    parser.add_argument(
        "-r",
        "--reference",
        required=True,
        help="GRCh37/hg19 reference FASTA used by SpliceAI, e.g. hg19.fa",
    )
    parser.add_argument("-o", "--output", help="Output TSV path")
    parser.add_argument("--input-vcf", help="VCF path to create for the SpliceAI input")
    parser.add_argument("--spliceai-vcf", help="VCF path written by SpliceAI")
    parser.add_argument(
        "--annotation",
        default="grch37",
        help="SpliceAI annotation argument. Use grch37 for this hg19-style pipeline.",
    )
    parser.add_argument("--distance", type=int, default=50, help="SpliceAI -D value")
    parser.add_argument("--mask", type=int, default=0, choices=[0, 1], help="SpliceAI -M value")
    parser.add_argument("--spliceai-exe", default="spliceai", help="SpliceAI executable name or full path")
    parser.add_argument(
        "--runner-prefix",
        default="",
        help='Optional command prefix, for example: "conda run -n spliceai"',
    )
    parser.add_argument("--keep-input-vcf", action="store_true", help="Keep the generated SpliceAI input VCF")
    args = parser.parse_args()

    defaults = default_output_paths(args.input_tsv)
    input_vcf = Path(args.input_vcf) if args.input_vcf else defaults["input_vcf"]
    spliceai_vcf = Path(args.spliceai_vcf) if args.spliceai_vcf else defaults["spliceai_vcf"]
    output_tsv = Path(args.output) if args.output else defaults["output_tsv"]

    input_vcf.parent.mkdir(parents=True, exist_ok=True)
    spliceai_vcf.parent.mkdir(parents=True, exist_ok=True)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    validate_inputs(args)

    print("Creating VCF for local SpliceAI...")
    create_small_vcf(args.input_tsv, input_vcf)

    spliceai_input_vcf, style = rewrite_vcf_chromosomes(input_vcf, args.reference)
    if style:
        print(f"Using chromosome style '{style}' from reference index.")

    command = build_command(args, spliceai_input_vcf, spliceai_vcf)
    print("Running local SpliceAI:")
    print(" ".join(str(part) for part in command))
    temp_spliceai_input = spliceai_input_vcf if spliceai_input_vcf != input_vcf else None
    try:
        subprocess.run(command, cwd=ROOT, check=True)
    finally:
        if temp_spliceai_input is not None:
            temp_spliceai_input.unlink(missing_ok=True)

    print("Merging local SpliceAI scores into TSV...")
    lookup = build_spliceai_lookup(spliceai_vcf)
    add_spliceai_to_tsv(lookup, args.input_tsv, output_tsv)

    if not args.keep_input_vcf:
        input_vcf.unlink(missing_ok=True)

    print(f"Local SpliceAI VCF: {spliceai_vcf}")
    print(f"Annotated TSV: {output_tsv}")


if __name__ == "__main__":
    main()
