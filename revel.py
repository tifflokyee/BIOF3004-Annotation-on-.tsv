import pandas as pd

# ====================== CONFIG ======================
REVEL_FILE = "./annotation/revel_with_transcript_ids"
VARIANT_FILE = "./result/freq_with_PanelApp_columns.tsv"
GENOME_BUILD = "hg19"  # using hg19_pos column
REVEL_CHUNK_SIZE = 500_000

# ====================== LOAD YOUR VARIANTS ======================
print("\nLoading variant data...")
df = pd.read_csv(VARIANT_FILE, sep="\t", low_memory=False)

print(f"Variants loaded: {len(df):,} rows")

# Clean chromosome names consistently
df["Chr"] = df["Chr"].astype("string").str.replace(r"^chr", "", regex=True).str.strip()

# Use the correct position column (confirmed: 'Coordinate')
df["pos"] = pd.to_numeric(df["Coordinate"], errors="coerce").astype("Int32")

# Drop rows with invalid / missing positions (should be 0 in your case)
invalid_pos = df["pos"].isna().sum()
if invalid_pos > 0:
    print(f"-> Dropped {invalid_pos} rows with invalid position")
df = df[df["pos"].notna()].copy()

# ====================== PARSE REF / ALT ======================
print("Parsing ref and alt alleles...")
if "Variant" not in df.columns:
    raise KeyError("Expected a 'Variant' column containing allele strings such as C>C/T")

variant_parts = df["Variant"].astype("string").str.split(">", n=1, expand=True)
if variant_parts.shape[1] < 2:
    raise ValueError("Could not parse allele strings from 'Variant' column")

df["ref"] = variant_parts[0].str.strip().str.upper().astype("string")
alt_options = variant_parts[1].fillna("").astype("string").str.split("/", expand=True)


def pick_alt(row):
    ref = row["ref"]
    for value in row.drop(labels=["ref"]).tolist():
        if pd.isna(value):
            continue
        allele = str(value).strip().upper()
        if allele and allele != ref:
            return allele
    first = row.iloc[1] if len(row) > 1 else pd.NA
    if pd.isna(first):
        return pd.NA
    first = str(first).strip().upper()
    return first or pd.NA


alt_candidates = alt_options.copy()
alt_candidates.insert(0, "ref", df["ref"])
df["alt"] = alt_candidates.apply(pick_alt, axis=1).astype("string")

# Quick check: how many look like clean SNVs (REVEL only scores missense SNVs)
is_snv = (
    (df["ref"].str.len() == 1)
    & (df["alt"].str.len() == 1)
    & df["ref"].notna()
    & df["alt"].notna()
    & df["pos"].notna()
)
print(f"SNVs eligible for REVEL: {is_snv.sum():,} / {len(df):,} ({is_snv.mean():.1%})")

candidate_variants = df.loc[is_snv, ["Chr", "pos", "ref", "alt"]].drop_duplicates().copy()
candidate_variants["pos"] = candidate_variants["pos"].astype("Int32")
candidate_keys = set(
    zip(
        candidate_variants["Chr"].astype(str),
        candidate_variants["pos"].astype(int),
        candidate_variants["ref"].astype(str),
        candidate_variants["alt"].astype(str),
    )
)
candidate_chromosomes = sorted(candidate_variants["Chr"].dropna().astype(str).unique().tolist())

print(f"Candidate SNV keys to look up: {len(candidate_keys):,}")
print(f"Chromosomes requested from REVEL: {candidate_chromosomes}")

# Diagnostic: show what we're trying to match
print("\nFirst 10 variants to be matched:")
print(df[["Gene", "Chr", "pos", "Variant", "ref", "alt", "Consequence"]].head(10))

# ====================== LOAD REVEL SCORES ======================
print("\nScanning REVEL scores in chunks...")
matched_chunks = []
seen_chromosomes = set()
chunk_count = 0

for chunk in pd.read_csv(
    REVEL_FILE,
    sep=",",
    usecols=["chr", "hg19_pos", "grch38_pos", "ref", "alt", "REVEL"],
    dtype={
        "chr": "string",
        "hg19_pos": "Int32",
        "grch38_pos": "Int32",
        "ref": "string",
        "alt": "string",
        "REVEL": "float32",
    },
    na_values=["."],
    low_memory=False,
    chunksize=REVEL_CHUNK_SIZE,
):
    chunk_count += 1
    chunk["chr"] = chunk["chr"].str.replace(r"^chr", "", regex=True).str.strip()
    seen_chromosomes.update(chunk["chr"].dropna().astype(str).unique().tolist())
    chunk = chunk[chunk["chr"].isin(candidate_chromosomes)].copy()
    if chunk.empty:
        continue

    chunk["key"] = list(
        zip(
            chunk["chr"].astype(str),
            chunk["hg19_pos"].astype("Int64"),
            chunk["ref"].astype(str),
            chunk["alt"].astype(str),
        )
    )
    matched = chunk[chunk["key"].isin(candidate_keys)].drop(columns=["key"])
    if not matched.empty:
        matched_chunks.append(matched)

print(f"Processed REVEL chunks: {chunk_count}")
chromosomes = sorted(seen_chromosomes)
print(f"Chromosomes present in REVEL: {chromosomes}")
if len(chromosomes) < 22:
    print(
        "WARNING: REVEL source appears incomplete. "
        "Expected broad chromosome coverage, but loaded only "
        f"{len(chromosomes)} chromosome(s)."
    )

if matched_chunks:
    revel_df = pd.concat(matched_chunks, ignore_index=True).drop_duplicates(
        subset=["chr", "hg19_pos", "ref", "alt"]
    )
else:
    revel_df = pd.DataFrame(columns=["chr", "hg19_pos", "ref", "alt", "REVEL"])

print(f"Relevant REVEL rows retained for merge: {len(revel_df):,}")

# ====================== MERGE WITH REVEL ======================
print("\nMatching variants to REVEL scores...")
df = df.merge(
    revel_df[["chr", "hg19_pos", "ref", "alt", "REVEL"]],
    left_on=["Chr", "pos", "ref", "alt"],
    right_on=["chr", "hg19_pos", "ref", "alt"],
    how="left",
)

# Clean up duplicate columns
df = df.drop(columns=["chr", "hg19_pos"], errors="ignore")

# Rename for clarity
df = df.rename(columns={"REVEL": "REVEL_score"})
df["REVEL_score"] = df["REVEL_score"].fillna(".")

# ====================== SUMMARY ======================
matched = (df["REVEL_score"] != ".").sum()
print(
    f"\nDone. REVEL scores added to {matched:,} out of {len(df):,} variants "
    f"({matched / len(df):.0%} match rate)"
)

unmatched = len(df) - matched
if unmatched > 0:
    print(
        f"-> {unmatched:,} variants did not receive a score "
        "(likely indels, non-missense, or no match in REVEL subset)"
    )

# Show results
print("\nFirst 12 rows with REVEL scores:")
print(df[["Gene", "Chr", "pos", "Variant", "Consequence", "REVEL_score"]].head(12))

# Optional: save the annotated file
output_file = "./result/variants_with_REVEL.tsv"
df.to_csv(output_file, sep="\t", index=False)
print(f"\nAnnotated file saved to: {output_file}")
