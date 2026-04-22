from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
DEFAULT_INPUTS = [
    ROOT / "CAPNGSETB24-modified proband-leukoencephalopathy_freq_annotated.tsv",
    ROOT / "CAPNGSETB24-modified proband-leukoencephalopathy_annotated.tsv",
    ROOT / "result" / "variants_with_AlphaMissense_and_ClinVar_and_gnomAD.tsv",
    ROOT / "result" / "variants_with_AlphaMissense_and_ClinVar.tsv",
    ROOT / "result" / "freq_with_PanelApp_columns.tsv",
]
OUTPUT_FILE = ROOT / "result" / "local_spliceai_input_for_external_run.vcf"


PANEL_CATEGORY_VALUES = {
    "Mendeliome (primary diagnostic panel)",
    "Paediatric Additional Findings (childhood-onset actionable)",
    "Incidentalome (adult-onset actionable)",
}

PANEL_COLUMNS = [
    "Mode of Inheritance (Mendeliome)",
    "Incidentalome Status",
    "Paediatric Additional Status",
]


def normalize_chrom(value):
    chrom = str(value).strip().replace("chr", "")
    return f"chr{chrom}" if chrom else ""


def parse_ref_alt(row):
    ref = str(row.get("ref", "")).strip().upper()
    alt = str(row.get("alt", "")).strip().upper()
    if ref and alt and ref not in {"NAN", "<NA>"} and alt not in {"NAN", "<NA>"}:
        return ref, alt

    variant = str(row.get("Variant", "")).strip()
    if ">" not in variant:
        return "", ""
    ref_part, alt_part = variant.split(">", 1)
    ref = ref_part.strip().upper()
    alt_options = [item.strip().upper() for item in alt_part.split("/") if item.strip()]
    if len(alt_options) > 1 and alt_options[0] == ref:
        return ref, alt_options[1]
    return ref, alt_options[0] if alt_options else ""


def is_clinically_relevant(row):
    if "Panel_Category" in row.index:
        return str(row.get("Panel_Category", "")).strip() in PANEL_CATEGORY_VALUES

    found_any_panel_col = any(column in row.index for column in PANEL_COLUMNS)
    if not found_any_panel_col:
        return True

    for column in PANEL_COLUMNS:
        value = str(row.get(column, "")).strip()
        if value and value not in {".", "nan", "NaN", "<NA>"}:
            return True
    return False


def load_default_input():
    for path in DEFAULT_INPUTS:
        if path.exists():
            print(f"Using input TSV: {path}")
            return pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    raise FileNotFoundError("No known annotated TSV file was found. Pass a TSV path as the first argument.")


def create_small_vcf(input_tsv=None, output_file=OUTPUT_FILE):
    if input_tsv:
        input_tsv = Path(input_tsv)
        print(f"Using input TSV: {input_tsv}")
        df = pd.read_csv(input_tsv, sep="\t", dtype=str, low_memory=False)
    else:
        df = load_default_input()

    df_relevant = df[df.apply(is_clinically_relevant, axis=1)].copy()
    print(f"Clinically relevant rows selected: {len(df_relevant):,} / {len(df):,}")

    seen = set()
    vcf_lines = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]

    for _, row in df_relevant.iterrows():
        chrom = normalize_chrom(row.get("Chr", ""))
        pos = str(row.get("Coordinate", "")).strip()
        ref, alt = parse_ref_alt(row)
        if not chrom or not pos or not ref or not alt:
            continue
        key = (chrom, pos, ref, alt)
        if key in seen:
            continue
        seen.add(key)
        vcf_lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.")

    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.write_text("\n".join(vcf_lines) + "\n", encoding="utf-8")
    print(f"Small VCF created: {output_file}")
    print(f"Variants written: {len(vcf_lines) - 2:,}")
    print("Use this VCF as the input for local SpliceAI on a newer Windows or Linux computer.")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Create a small VCF for local SpliceAI scoring on another computer.")
    parser.add_argument("input_tsv", nargs="?", help="Annotated TSV to filter. Defaults to known local outputs.")
    parser.add_argument("-o", "--output", default=str(OUTPUT_FILE), help="Output VCF path.")
    args = parser.parse_args()
    create_small_vcf(args.input_tsv, args.output)
