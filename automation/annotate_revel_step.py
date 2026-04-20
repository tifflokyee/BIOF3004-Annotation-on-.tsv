import argparse
from pathlib import Path

from auto_annotate_generated import annotate_revel, default_output_path, load_variant_table


def main():
    parser = argparse.ArgumentParser(description="Step 2: add REVEL annotations.")
    parser.add_argument("input", help="Input TSV/VCF/VCF.GZ file")
    parser.add_argument("-o", "--output", help="Output TSV file")
    parser.add_argument("--annotation-dir", default="annotation")
    parser.add_argument("--chunk-size", type=int, default=500000)
    args = parser.parse_args()

    output = Path(args.output) if args.output else default_output_path(args.input)
    output.parent.mkdir(parents=True, exist_ok=True)

    print("=== Step: REVEL annotation ===")
    print(f"Input: {args.input}")
    df = load_variant_table(args.input)
    before_columns = set(df.columns)
    df = annotate_revel(df, Path(args.annotation_dir), args.chunk_size)
    added_columns = [column for column in df.columns if column not in before_columns]
    matched = (df["REVEL_score"] != ".").sum()
    df.to_csv(output, sep="\t", index=False)
    print(f"Rows annotated by REVEL: {matched:,} / {len(df):,}")
    print(f"Columns added: {', '.join(added_columns) if added_columns else 'none'}")
    print(f"Output: {output}")


if __name__ == "__main__":
    main()
