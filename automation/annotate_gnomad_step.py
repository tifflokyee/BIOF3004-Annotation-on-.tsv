import argparse
from pathlib import Path

from auto_annotate_generated import annotate_gnomad, default_output_path, load_variant_table


def main():
    parser = argparse.ArgumentParser(description="Step 5: add gnomAD annotations.")
    parser.add_argument("input", help="Input TSV/VCF/VCF.GZ file")
    parser.add_argument("-o", "--output", help="Output TSV file")
    parser.add_argument("--annotation-dir", default="annotation")
    args = parser.parse_args()

    output = Path(args.output) if args.output else default_output_path(args.input)
    output.parent.mkdir(parents=True, exist_ok=True)

    print("=== Step: gnomAD annotation ===")
    print(f"Input: {args.input}")
    df = load_variant_table(args.input)
    before_columns = set(df.columns)
    df = annotate_gnomad(df, Path(args.annotation_dir))
    added_columns = [column for column in df.columns if column not in before_columns]
    matched = (df["gnomad_af"] != ".").sum()
    df.to_csv(output, sep="\t", index=False)
    print(f"Rows annotated by gnomAD: {matched:,} / {len(df):,}")
    print(f"Columns added: {', '.join(added_columns) if added_columns else 'none'}")
    print(f"Output: {output}")


if __name__ == "__main__":
    main()
