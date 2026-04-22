import argparse
from pathlib import Path

from auto_annotate_generated import default_output_path, run_local_spliceai


def main():
    parser = argparse.ArgumentParser(description="Step 6: add local SpliceAI annotations.")
    parser.add_argument("input", help="Input annotated TSV file")
    parser.add_argument("-o", "--output", help="Output TSV file")
    parser.add_argument("--local-spliceai-vcf", help="Output VCF containing local SpliceAI annotations")
    parser.add_argument("--local-spliceai-cache", help="JSON cache for fetched hg19 sequence windows")
    parser.add_argument("--spliceai-distance", type=int, default=50)
    parser.add_argument("--spliceai-mask", type=int, default=0, choices=[0, 1])
    args = parser.parse_args()

    output = Path(args.output) if args.output else default_output_path(args.input)
    output.parent.mkdir(parents=True, exist_ok=True)

    print("=== Step: Local SpliceAI annotation ===")
    print(f"Input: {args.input}")
    run_local_spliceai(
        args.input,
        output,
        args.local_spliceai_vcf,
        args.local_spliceai_cache,
        args.spliceai_distance,
        args.spliceai_mask,
    )
    print(f"Output: {output}")


if __name__ == "__main__":
    main()
