import argparse
import gzip
from pathlib import Path

import pandas as pd

from add_spliceai_from_vep import normalize_chrom, parse_info, parse_tsv_ref_alt


DELTA_COLUMNS = {
    "SpliceAI_DS_AG": "acceptor gain",
    "SpliceAI_DS_AL": "acceptor loss",
    "SpliceAI_DS_DG": "donor gain",
    "SpliceAI_DS_DL": "donor loss",
}


def open_text(path):
    path = Path(path)
    if path.suffix.lower() in {".gz", ".bgz"}:
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def parse_float(value):
    if value in {None, "", "."} or pd.isna(value):
        return pd.NA
    try:
        return float(value)
    except ValueError:
        return pd.NA


def better_scores(current, candidate):
    if current is None:
        return candidate
    current_max = max(value for value in current if pd.notna(value))
    candidate_max = max(value for value in candidate if pd.notna(value))
    return candidate if candidate_max > current_max else current


def build_delta_lookup(spliceai_vcf):
    lookup = {}
    with open_text(spliceai_vcf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            chrom, pos, _record_id, ref, alts, _qual, _filt, info_text = fields[:8]
            chrom = normalize_chrom(chrom)
            ref = ref.upper()
            info = parse_info(info_text)
            if "SpliceAI" not in info:
                continue

            for alt in [item.upper() for item in alts.split(",")]:
                best = None
                for record in str(info["SpliceAI"]).split(","):
                    parts = record.split("|")
                    if len(parts) < 6 or parts[0].upper() != alt:
                        continue
                    scores = tuple(parse_float(value) for value in parts[2:6])
                    if any(pd.isna(value) for value in scores):
                        continue
                    best = better_scores(best, scores)
                if best is not None:
                    lookup[(chrom, pos, ref, alt)] = best
    return lookup


def add_delta_scores(input_tsv, spliceai_vcf, output_tsv):
    df = pd.read_csv(input_tsv, sep="\t", dtype=str, low_memory=False)
    lookup = build_delta_lookup(spliceai_vcf)

    def get_scores(row):
        chrom = normalize_chrom(row.get("Chr", ""))
        pos = str(row.get("Coordinate", row.get("POS", ""))).strip()
        ref, alt = parse_tsv_ref_alt(row)
        return lookup.get((chrom, pos, ref, alt), (pd.NA, pd.NA, pd.NA, pd.NA))

    scores = df.apply(get_scores, axis=1, result_type="expand")
    scores.columns = list(DELTA_COLUMNS)
    for column in DELTA_COLUMNS:
        df[column] = scores[column].apply(lambda value: "." if pd.isna(value) else round(float(value), 4))

    output_tsv = Path(output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_tsv, sep="\t", index=False)
    matched = (df["SpliceAI_DS_AG"] != ".").sum()
    print(f"Added SpliceAI delta-score columns to {output_tsv} ({matched:,}/{len(df):,} rows matched).")
    print("Columns added:")
    for column, label in DELTA_COLUMNS.items():
        print(f"  {column}: {label}")


def main():
    parser = argparse.ArgumentParser(description="Add the four SpliceAI delta-score columns to a TSV.")
    parser.add_argument("input_tsv")
    parser.add_argument("--spliceai-vcf", required=True, help="VCF containing native SpliceAI INFO annotations")
    parser.add_argument("-o", "--output", required=True, help="Output TSV path")
    args = parser.parse_args()
    add_delta_scores(args.input_tsv, args.spliceai_vcf, args.output)


if __name__ == "__main__":
    main()
