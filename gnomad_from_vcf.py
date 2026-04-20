import argparse
import gzip
from glob import glob
from pathlib import Path

import pandas as pd


DEFAULT_INPUT_FILE = "./result/variants_with_AlphaMissense_and_ClinVar.tsv"
DEFAULT_OUTPUT_FILE = "./result/variants_with_AlphaMissense_and_ClinVar_and_gnomAD.tsv"
DEFAULT_GNOMAD_PATTERNS = ["./annotation/gnomad*.vcf.gz", "./annotation/gnomad*.vcf.bgz"]

INFO_FIELDS = [
    ("gnomad_filter", "__FILTER__"),
    ("gnomad_ac", "AC"),
    ("gnomad_an", "AN"),
    ("gnomad_af", "AF"),
    ("gnomad_nhomalt", "nhomalt"),
    ("gnomad_ac_afr", "AC_AFR"),
    ("gnomad_an_afr", "AN_AFR"),
    ("gnomad_af_afr", "AF_AFR"),
    ("gnomad_ac_amr", "AC_AMR"),
    ("gnomad_an_amr", "AN_AMR"),
    ("gnomad_af_amr", "AF_AMR"),
    ("gnomad_ac_asj", "AC_ASJ"),
    ("gnomad_an_asj", "AN_ASJ"),
    ("gnomad_af_asj", "AF_ASJ"),
    ("gnomad_ac_eas", "AC_EAS"),
    ("gnomad_an_eas", "AN_EAS"),
    ("gnomad_af_eas", "AF_EAS"),
    ("gnomad_ac_fin", "AC_FIN"),
    ("gnomad_an_fin", "AN_FIN"),
    ("gnomad_af_fin", "AF_FIN"),
    ("gnomad_ac_nfe", "AC_NFE"),
    ("gnomad_an_nfe", "AN_NFE"),
    ("gnomad_af_nfe", "AF_NFE"),
    ("gnomad_ac_oth", "AC_OTH"),
    ("gnomad_an_oth", "AN_OTH"),
    ("gnomad_af_oth", "AF_OTH"),
    ("gnomad_ac_sas", "AC_SAS"),
    ("gnomad_an_sas", "AN_SAS"),
    ("gnomad_af_sas", "AF_SAS"),
    ("gnomad_popmax", "popmax"),
    ("gnomad_faf95_popmax", "faf95_popmax"),
    ("gnomad_source_file", "__SOURCE__"),
    ("gnomad_info", "__INFO__"),
]


def clean_text(value):
    if value is None or pd.isna(value):
        return "."
    text = str(value).strip()
    return text if text else "."


def normalize_chromosome(value):
    chrom = str(value).strip()
    chrom = chrom.removeprefix("chr").strip()
    if chrom == "M":
        chrom = "MT"
    return chrom


def parse_info_field(info_text):
    info = {}
    for item in info_text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        elif item:
            info[item] = True
    return info


def get_info_value(info, info_key, allele_index, alt_count, source_name, info_text, filt):
    if info_key == "__FILTER__":
        return filt
    if info_key == "__INFO__":
        return info_text
    if info_key == "__SOURCE__":
        return source_name

    raw_value = info.get(info_key)
    if raw_value is None:
        return None
    if raw_value is True:
        return "TRUE"

    text = str(raw_value).strip()
    if "," not in text:
        return text

    parts = [part.strip() for part in text.split(",")]
    if len(parts) == alt_count:
        return parts[allele_index]
    return text


def pick_gnomad_files(patterns, candidate_chromosomes):
    resolved = []
    for pattern in patterns:
        matches = sorted(glob(pattern))
        if matches:
            resolved.extend(matches)
        elif Path(pattern).exists():
            resolved.append(pattern)

    unique_files = []
    seen = set()
    for file_path in resolved:
        if file_path not in seen:
            unique_files.append(file_path)
            seen.add(file_path)

    if not unique_files:
        raise FileNotFoundError(
            "No gnomAD VCF files were found. Put them under ./annotation or pass --gnomad."
        )

    selected = []
    wanted_tokens = {chrom.lower() for chrom in candidate_chromosomes}
    for file_path in unique_files:
        name = Path(file_path).name.lower()
        if ".chr" not in name:
            selected.append(file_path)
            continue

        matched = False
        for chrom in wanted_tokens:
            if f".chr{chrom}." in name:
                matched = True
                break
        if matched:
            selected.append(file_path)

    return selected or unique_files


def build_argument_parser():
    parser = argparse.ArgumentParser(
        description="Append gnomAD annotations to a TSV by exact Chr/pos/ref/alt matching."
    )
    parser.add_argument("--input", default=DEFAULT_INPUT_FILE, help="Input TSV file")
    parser.add_argument("--output", default=DEFAULT_OUTPUT_FILE, help="Output TSV file")
    parser.add_argument(
        "--gnomad",
        nargs="+",
        default=DEFAULT_GNOMAD_PATTERNS,
        help="One or more VCF(.gz/.bgz) paths or glob patterns",
    )
    parser.add_argument(
        "--coordinate-column",
        default="Coordinate",
        help="Input column to use as genomic position. Use pos/pos38 if needed.",
    )
    return parser


def main():
    args = build_argument_parser().parse_args()

    print("\nLoading variant table...")
    df = pd.read_csv(args.input, sep="\t", low_memory=False)
    print(f"Variants loaded: {len(df):,} rows")

    required_columns = ["Chr", args.coordinate_column, "ref", "alt"]
    missing_columns = [column for column in required_columns if column not in df.columns]
    if missing_columns:
        raise KeyError(f"Missing required column(s): {', '.join(missing_columns)}")

    df["Chr"] = df["Chr"].astype("string").map(normalize_chromosome)
    df[args.coordinate_column] = pd.to_numeric(df[args.coordinate_column], errors="coerce").astype("Int64")
    df["ref"] = df["ref"].astype("string").str.upper().str.strip()
    df["alt"] = df["alt"].astype("string").str.upper().str.strip()

    candidate_variants = df.loc[
        df[args.coordinate_column].notna() & df["ref"].notna() & df["alt"].notna(),
        ["Chr", args.coordinate_column, "ref", "alt"],
    ].drop_duplicates()

    candidate_keys = set(
        zip(
            candidate_variants["Chr"].astype(str),
            candidate_variants[args.coordinate_column].astype(int),
            candidate_variants["ref"].astype(str),
            candidate_variants["alt"].astype(str),
        )
    )

    candidate_chromosomes = sorted(candidate_variants["Chr"].dropna().astype(str).unique().tolist())
    print(f"Candidate variant keys to look up: {len(candidate_keys):,}")
    print(f"Chromosomes requested from gnomAD: {candidate_chromosomes}")

    gnomad_files = pick_gnomad_files(args.gnomad, candidate_chromosomes)
    print("\nUsing gnomAD files:")
    for file_path in gnomad_files:
        print(f"  - {file_path}")

    print("\nScanning gnomAD VCF...")
    matches = []
    records_seen = 0
    records_variant_matched = 0

    for file_path in gnomad_files:
        source_name = Path(file_path).name
        with gzip.open(file_path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue

                records_seen += 1
                chrom, pos, record_id, ref, alt, qual, filt, info_text = line.rstrip("\n").split("\t", 8)
                chrom = normalize_chromosome(chrom)
                if chrom not in candidate_chromosomes:
                    continue

                pos = int(pos)
                ref = ref.strip().upper()
                alt_alleles = [alt_allele.strip().upper() for alt_allele in alt.split(",")]
                info = None

                for allele_index, alt_allele in enumerate(alt_alleles):
                    key = (chrom, pos, ref, alt_allele)
                    if key not in candidate_keys:
                        continue

                    records_variant_matched += 1
                    if info is None:
                        info = parse_info_field(info_text)

                    row = {
                        "Chr": chrom,
                        args.coordinate_column: pos,
                        "ref": ref,
                        "alt": alt_allele,
                    }

                    for output_column, info_key in INFO_FIELDS:
                        value = get_info_value(
                            info=info,
                            info_key=info_key,
                            allele_index=allele_index,
                            alt_count=len(alt_alleles),
                            source_name=source_name,
                            info_text=info_text,
                            filt=filt,
                        )
                        row[output_column] = clean_text(value)

                    matches.append(row)

    print(f"gnomAD records scanned: {records_seen:,}")
    print(f"Exact variant matches found: {records_variant_matched:,}")

    gnomad_df = pd.DataFrame(matches)
    value_columns = [name for name, _ in INFO_FIELDS]
    if gnomad_df.empty:
        gnomad_df = pd.DataFrame(columns=["Chr", args.coordinate_column, "ref", "alt", *value_columns])
    else:
        gnomad_df = gnomad_df.drop_duplicates(subset=["Chr", args.coordinate_column, "ref", "alt"])

    print(f"Rows retained for merge: {len(gnomad_df):,}")

    print("\nMerging gnomAD data into variant table...")
    df = df.merge(
        gnomad_df,
        on=["Chr", args.coordinate_column, "ref", "alt"],
        how="left",
    )

    for output_column in value_columns:
        df[output_column] = df[output_column].fillna(".")

    matched_rows = (df["gnomad_af"] != ".").sum()
    print(
        f"\nDone. gnomAD data added to {matched_rows:,} out of {len(df):,} rows "
        f"({matched_rows / len(df):.0%} match rate)"
    )

    df.to_csv(args.output, sep="\t", index=False)
    print(f"\nAnnotated file saved to: {args.output}")


if __name__ == "__main__":
    main()
