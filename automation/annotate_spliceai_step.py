import argparse
from pathlib import Path

from auto_annotate_generated import annotate_spliceai_from_vep, default_output_path, load_variant_table


def main():
    parser = argparse.ArgumentParser(description="Step 6: add real SpliceAI annotations from VEP output.")
    parser.add_argument("input", help="Input TSV/VCF/VCF.GZ file")
    parser.add_argument("-o", "--output", help="Output TSV file")
    parser.add_argument("--vep-spliceai", required=True, help="VEP output VCF/TSV containing SpliceAI DS scores")
    args = parser.parse_args()

    output = Path(args.output) if args.output else default_output_path(args.input)
    output.parent.mkdir(parents=True, exist_ok=True)

    print("=== Step: SpliceAI annotation from VEP ===")
    print(f"Input: {args.input}")
    print(f"VEP SpliceAI file: {args.vep_spliceai}")
    df = load_variant_table(args.input)
    before_columns = set(df.columns)
    df = annotate_spliceai_from_vep(df, args.vep_spliceai)
    added_columns = [column for column in df.columns if column not in before_columns]
    matched = (df["SpliceAI_max"] != ".").sum()
    df.to_csv(output, sep="\t", index=False)
    print(f"Rows annotated by SpliceAI: {matched:,} / {len(df):,}")
    print(f"Columns added: {', '.join(added_columns) if added_columns else 'none'}")
    print(f"Output: {output}")


if __name__ == "__main__":
    main()
