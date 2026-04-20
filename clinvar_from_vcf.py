import gzip
from collections import defaultdict

import pandas as pd

# ====================== CONFIG ======================
INPUT_FILE = "./result/variants_with_AlphaMissense.tsv"
CLINVAR_VCF = "./annotation/clinvar.vcf.gz"
OUTPUT_FILE = "./result/variants_with_AlphaMissense_and_ClinVar.tsv"

OUTPUT_COLUMNS = [
    "clinvar_vcf_id",
    "clinvar_vcf_alleleid",
    "clinvar_vcf_clnsig",
    "clinvar_vcf_clnrevstat",
    "clinvar_vcf_clndn",
    "clinvar_vcf_clndisdb",
    "clinvar_vcf_clnhgvs",
    "clinvar_vcf_clnvc",
    "clinvar_vcf_mc",
    "clinvar_vcf_origin",
    "clinvar_vcf_geneinfo",
    "clinvar_vcf_rs",
    "clinvar_vcf_info",
]


def clean_gene_symbol(value):
    if pd.isna(value):
        return ""
    return str(value).strip().upper()


def clean_text(value):
    if value is None or pd.isna(value):
        return "."
    text = str(value).strip()
    return text if text else "."


def parse_info_field(info_text):
    info = {}
    for item in info_text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        elif item:
            info[item] = True
    return info


def parse_gene_symbols(geneinfo_text):
    if not geneinfo_text:
        return set()

    symbols = set()
    for item in str(geneinfo_text).split("|"):
        if not item:
            continue
        symbol = item.split(":", 1)[0].strip().upper()
        if symbol:
            symbols.add(symbol)
    return symbols


print("\nLoading AlphaMissense-annotated variants...")
df = pd.read_csv(INPUT_FILE, sep="\t", low_memory=False)
print(f"Variants loaded: {len(df):,} rows")

required_columns = ["Gene", "Chr", "Coordinate", "ref", "alt"]
missing_columns = [column for column in required_columns if column not in df.columns]
if missing_columns:
    raise KeyError(f"Missing required column(s): {', '.join(missing_columns)}")

df["Gene"] = df["Gene"].astype("string").map(clean_gene_symbol)
df["Chr"] = df["Chr"].astype("string").str.replace(r"^chr", "", regex=True).str.strip()
df["Coordinate"] = pd.to_numeric(df["Coordinate"], errors="coerce").astype("Int32")
df["ref"] = df["ref"].astype("string").str.upper().str.strip()
df["alt"] = df["alt"].astype("string").str.upper().str.strip()

candidate_rows = df.loc[
    df["Gene"].ne("")
    & df["Coordinate"].notna()
    & df["ref"].notna()
    & df["alt"].notna(),
    ["Gene", "Chr", "Coordinate", "ref", "alt"],
].drop_duplicates()

candidate_keys = set(
    zip(
        candidate_rows["Gene"].astype(str),
        candidate_rows["Chr"].astype(str),
        candidate_rows["Coordinate"].astype(int),
        candidate_rows["ref"].astype(str),
        candidate_rows["alt"].astype(str),
    )
)

candidate_variant_keys = set(
    zip(
        candidate_rows["Chr"].astype(str),
        candidate_rows["Coordinate"].astype(int),
        candidate_rows["ref"].astype(str),
        candidate_rows["alt"].astype(str),
    )
)

print(f"Candidate gene + variant keys to look up: {len(candidate_keys):,}")

print("\nScanning ClinVar VCF...")
matches_by_key = defaultdict(list)
records_seen = 0
records_variant_matched = 0
records_gene_matched = 0

with gzip.open(CLINVAR_VCF, "rt", encoding="utf-8", errors="replace") as handle:
    for line in handle:
        if line.startswith("#"):
            continue

        records_seen += 1
        chrom, pos, record_id, ref, alt, qual, filt, info_text = line.rstrip("\n").split("\t", 8)
        chrom = chrom.replace("chr", "").strip()
        pos = int(pos)
        ref = ref.strip().upper()

        for alt_allele in alt.split(","):
            alt_allele = alt_allele.strip().upper()
            variant_key = (chrom, pos, ref, alt_allele)
            if variant_key not in candidate_variant_keys:
                continue

            records_variant_matched += 1
            info = parse_info_field(info_text)
            geneinfo = info.get("GENEINFO", "")
            gene_symbols = parse_gene_symbols(geneinfo)
            if not gene_symbols:
                continue

            record = {
                "clinvar_vcf_id": clean_text(record_id),
                "clinvar_vcf_alleleid": clean_text(info.get("ALLELEID")),
                "clinvar_vcf_clnsig": clean_text(info.get("CLNSIG")),
                "clinvar_vcf_clnrevstat": clean_text(info.get("CLNREVSTAT")),
                "clinvar_vcf_clndn": clean_text(info.get("CLNDN")),
                "clinvar_vcf_clndisdb": clean_text(info.get("CLNDISDB")),
                "clinvar_vcf_clnhgvs": clean_text(info.get("CLNHGVS")),
                "clinvar_vcf_clnvc": clean_text(info.get("CLNVC")),
                "clinvar_vcf_mc": clean_text(info.get("MC")),
                "clinvar_vcf_origin": clean_text(info.get("ORIGIN")),
                "clinvar_vcf_geneinfo": clean_text(geneinfo),
                "clinvar_vcf_rs": clean_text(info.get("RS")),
                "clinvar_vcf_info": clean_text(info_text),
            }

            for gene_symbol in gene_symbols:
                key = (gene_symbol, chrom, pos, ref, alt_allele)
                if key in candidate_keys:
                    matches_by_key[key].append(record)
                    records_gene_matched += 1

print(f"ClinVar records scanned: {records_seen:,}")
print(f"Variant-level matches found: {records_variant_matched:,}")
print(f"Gene-filtered matches found: {records_gene_matched:,}")

annotated_records = []
for key, records in matches_by_key.items():
    merged = {}
    for column in OUTPUT_COLUMNS:
        values = []
        seen = set()
        for record in records:
            value = record.get(column, ".")
            if value == "." or value in seen:
                continue
            values.append(value)
            seen.add(value)
        merged[column] = " | ".join(values) if values else "."

    gene_symbol, chrom, pos, ref, alt = key
    merged["Gene"] = gene_symbol
    merged["Chr"] = chrom
    merged["Coordinate"] = pos
    merged["ref"] = ref
    merged["alt"] = alt
    annotated_records.append(merged)

clinvar_df = pd.DataFrame(annotated_records)
if clinvar_df.empty:
    clinvar_df = pd.DataFrame(columns=["Gene", "Chr", "Coordinate", "ref", "alt", *OUTPUT_COLUMNS])

print(f"Rows retained for merge: {len(clinvar_df):,}")

print("\nMerging ClinVar data into variant table...")
df = df.merge(
    clinvar_df,
    on=["Gene", "Chr", "Coordinate", "ref", "alt"],
    how="left",
)

for column in OUTPUT_COLUMNS:
    df[column] = df[column].fillna(".")

matched_rows = (df["clinvar_vcf_id"] != ".").sum()
print(
    f"\nDone. ClinVar VCF data added to {matched_rows:,} out of {len(df):,} rows "
    f"({matched_rows / len(df):.0%} match rate)"
)

df.to_csv(OUTPUT_FILE, sep="\t", index=False)
print(f"\nAnnotated file saved to: {OUTPUT_FILE}")
