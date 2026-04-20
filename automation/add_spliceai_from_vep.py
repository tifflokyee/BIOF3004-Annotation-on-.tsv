import argparse
import gzip
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
MAX_COL = "SpliceAI_max"
IMPACT_COL = "SpliceAI_Impact"
DELTA_SCORE_NAMES = ("DS_AG", "DS_AL", "DS_DG", "DS_DL")
SPLICEAI_DELTA_SCORE_NAMES = (
    "SpliceAI_pred_DS_AG",
    "SpliceAI_pred_DS_AL",
    "SpliceAI_pred_DS_DG",
    "SpliceAI_pred_DS_DL",
)


def open_text(path):
    path = Path(path)
    if path.suffix.lower() in {".gz", ".bgz"}:
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def normalize_chrom(value):
    text = str(value).strip()
    text = re.sub(r"^chr", "", text, flags=re.IGNORECASE)
    return "MT" if text == "M" else text


def parse_info(info_text):
    parsed = {}
    for item in str(info_text).split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            parsed[key] = value
        elif item:
            parsed[item] = True
    return parsed


def parse_csq_fields_from_header(vep_file):
    with open_text(vep_file) as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=CSQ"):
                match = re.search(r"Format: ([^\">]+)", line)
                if match:
                    return [field.strip() for field in match.group(1).split("|")]
            if line.startswith("#CHROM"):
                break
    return []


def numeric_delta_scores(value):
    if value is None or pd.isna(value):
        return []
    text = str(value)
    if not text or text == ".":
        return []

    values = []
    for number in re.findall(r"(?<!\d)(?:0(?:\.\d+)?|1(?:\.0+)?)(?!\d)", text):
        try:
            parsed = float(number)
        except ValueError:
            continue
        if 0.0 <= parsed <= 1.0:
            values.append(parsed)
    return values


def extract_max_spliceai_from_values(values):
    scores = []
    for value in values:
        scores.extend(numeric_delta_scores(value))
    return round(max(scores), 4) if scores else 0.0


def add_lookup_score(lookup, key, score):
    if score is None:
        return
    current = lookup.get(key)
    if current is None or score > current:
        lookup[key] = score


def build_lookup_from_table(vep_file):
    df = pd.read_csv(vep_file, sep="\t", comment="#", dtype=str, low_memory=False)
    splice_cols = [
        col
        for col in df.columns
        if col in SPLICEAI_DELTA_SCORE_NAMES or col.upper() in DELTA_SCORE_NAMES
    ]
    if not splice_cols:
        return {}

    chrom_col = next((col for col in ("#CHROM", "CHROM", "Chr", "chr", "Location") if col in df.columns), None)
    pos_col = next((col for col in ("POS", "Coordinate", "Start", "Uploaded_variation") if col in df.columns), None)
    ref_col = next((col for col in ("REF", "ref") if col in df.columns), None)
    alt_col = next((col for col in ("ALT", "Allele", "alt") if col in df.columns), None)

    if not chrom_col or not pos_col or not ref_col or not alt_col:
        return {}

    lookup = {}
    for _, row in df.iterrows():
        chrom = normalize_chrom(row.get(chrom_col, ""))
        pos = str(row.get(pos_col, "")).strip()
        ref = str(row.get(ref_col, "")).strip().upper()
        alt = str(row.get(alt_col, "")).strip().upper()
        if not chrom or not pos or not ref or not alt or pos == ".":
            continue
        score = extract_max_spliceai_from_values([row[col] for col in splice_cols])
        add_lookup_score(lookup, (chrom, pos, ref, alt), score)
    return lookup


def build_lookup_from_vcf(vep_file):
    csq_fields = parse_csq_fields_from_header(vep_file)
    csq_splice_indexes = [
        index
        for index, field in enumerate(csq_fields)
        if field in SPLICEAI_DELTA_SCORE_NAMES or field.upper() in DELTA_SCORE_NAMES
    ]
    lookup = {}

    with open_text(vep_file) as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            chrom, pos, _record_id, ref, alts, _qual, _filt, info_text = fields[:8]
            chrom = normalize_chrom(chrom)
            pos = str(pos).strip()
            ref = ref.strip().upper()
            info = parse_info(info_text)
            alt_list = [alt.strip().upper() for alt in alts.split(",")]

            direct_values = []
            for key, value in info.items():
                if key in SPLICEAI_DELTA_SCORE_NAMES or key.upper() in DELTA_SCORE_NAMES:
                    direct_values.append(value)

            for alt in alt_list:
                score_values = list(direct_values)
                if "CSQ" in info and csq_splice_indexes:
                    for record in str(info["CSQ"]).split(","):
                        parts = record.split("|")
                        if parts and parts[0].strip().upper() not in {"", alt}:
                            continue
                        for index in csq_splice_indexes:
                            if index < len(parts):
                                score_values.append(parts[index])

                score = extract_max_spliceai_from_values(score_values)
                add_lookup_score(lookup, (chrom, pos, ref, alt), score)

    return lookup


def build_spliceai_lookup(vep_file):
    vep_file = Path(vep_file)
    print(f"Reading VEP output file: {vep_file}")

    table_lookup = build_lookup_from_table(vep_file)
    if table_lookup:
        print(f"Found tabular SpliceAI annotations for {len(table_lookup):,} variants.")
        return table_lookup

    vcf_lookup = build_lookup_from_vcf(vep_file)
    if vcf_lookup:
        print(f"Found VCF SpliceAI annotations for {len(vcf_lookup):,} variants.")
        return vcf_lookup

    raise ValueError(
        "Could not find SpliceAI annotations. Expected SpliceAI_pred, SpliceAI, "
        "DS_AG, DS_AL, DS_DG, or DS_DL in the VEP output."
    )


def parse_tsv_ref_alt(row):
    ref = str(row.get("ref", "")).strip().upper()
    alt = str(row.get("alt", "")).strip().upper()
    if ref and alt and ref != "<NA>" and alt != "<NA>":
        return ref, alt

    variant = str(row.get("Variant", "")).strip()
    if ">" not in variant:
        return "", ""
    ref_part, alt_part = variant.split(">", 1)
    ref = ref_part.strip().upper()
    alts = [item.strip().upper() for item in alt_part.split("/") if item.strip()]
    if len(alts) > 1 and alts[0] == ref:
        alt = alts[1]
    elif alts:
        alt = alts[0]
    else:
        alt = ""
    return ref, alt


def add_spliceai_to_tsv(lookup, tsv_path, output_path=None):
    tsv_path = Path(tsv_path)
    df = pd.read_csv(tsv_path, sep="\t", dtype=str, low_memory=False)

    def get_score(row):
        chrom = normalize_chrom(row.get("Chr", ""))
        pos = str(row.get("Coordinate", row.get("POS", ""))).strip()
        ref, alt = parse_tsv_ref_alt(row)
        if not chrom or not pos or not ref or not alt:
            return pd.NA
        return lookup.get((chrom, pos, ref, alt), pd.NA)

    df[MAX_COL] = df.apply(get_score, axis=1)
    df[IMPACT_COL] = df[MAX_COL].apply(
        lambda score: (
            "High"
            if pd.notna(score) and float(score) >= 0.2
            else "Medium"
            if pd.notna(score) and float(score) >= 0.05
            else "Low"
            if pd.notna(score)
            else "."
        )
    )
    df[MAX_COL] = df[MAX_COL].apply(lambda score: "." if pd.isna(score) else round(float(score), 4))

    if output_path is None:
        output_path = tsv_path.with_name(tsv_path.stem + "_with_real_SpliceAI.tsv")
    else:
        output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)

    matched = (df[MAX_COL] != ".").sum()
    print(f"Saved {output_path} ({matched:,}/{len(df):,} rows matched SpliceAI).")
    return output_path


def default_tsv_files():
    candidates = [
        ROOT / "CAPNGSETB24-modified proband-leukoencephalopathy_freq_annotated.tsv",
        ROOT / "CAPNGSETB24-modified proband-leukoencephalopathy_annotated.tsv",
        ROOT / "result" / "variants_with_AlphaMissense_and_ClinVar_and_gnomAD.tsv",
        ROOT / "result" / "variants_with_AlphaMissense_and_ClinVar.tsv",
    ]
    return [path for path in candidates if path.exists()]


def build_parser():
    parser = argparse.ArgumentParser(
        description="Add real max-delta SpliceAI scores from VEP output to annotated TSV file(s)."
    )
    parser.add_argument(
        "--vep",
        default=str(ROOT / "VEP_output_file.vcf"),
        help="VEP output file containing SpliceAI annotations.",
    )
    parser.add_argument(
        "tsv",
        nargs="*",
        help="TSV file(s) to annotate. Defaults to known annotated TSV outputs in this folder.",
    )
    return parser


def main():
    args = build_parser().parse_args()
    vep_file = Path(args.vep)
    if not vep_file.exists():
        print(f"ERROR: VEP output not found: {vep_file}")
        print("Upload small_clinical_variants_for_vep.vcf to GRCh37 VEP with SpliceAI enabled,")
        print("then save the downloaded output as VEP_output_file.vcf in this folder.")
        print("Or run: add_spliceai_from_vep.bat --vep path\\to\\your_vep_output.vcf")
        return 1

    tsv_files = [Path(path) for path in args.tsv] if args.tsv else default_tsv_files()
    if not tsv_files:
        print("ERROR: No TSV files found to annotate. Pass one or more TSV paths.")
        return 1

    try:
        lookup = build_spliceai_lookup(vep_file)
    except Exception as error:
        print(f"ERROR: {error}")
        return 1

    for tsv_path in tsv_files:
        add_spliceai_to_tsv(lookup, tsv_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
