import pandas as pd

# ====================== CONFIG ======================
INPUT_FILE = "./result/variants_with_REVEL.tsv"
ALPHAMISSENSE_FILE = "./annotation/AlphaMissense_hg19.tsv.gz"
OUTPUT_FILE = "./result/variants_with_REVEL_and_AlphaMissense.tsv"
ALPHAMISSENSE_CHUNK_SIZE = 500_000

# ====================== LOAD VARIANTS ======================
print("\nLoading REVEL-annotated variants...")
df = pd.read_csv(INPUT_FILE, sep="\t", low_memory=False)
print(f"Variants loaded: {len(df):,} rows")

required_columns = ["Chr", "Coordinate", "ref", "alt"]
missing_columns = [column for column in required_columns if column not in df.columns]
if missing_columns:
    raise KeyError(f"Missing required column(s): {', '.join(missing_columns)}")

df["Chr"] = df["Chr"].astype("string").str.replace(r"^chr", "", regex=True).str.strip()
df["Coordinate"] = pd.to_numeric(df["Coordinate"], errors="coerce").astype("Int32")
df["ref"] = df["ref"].astype("string").str.upper().str.strip()
df["alt"] = df["alt"].astype("string").str.upper().str.strip()

candidate_variants = df.loc[
    df["Coordinate"].notna() & df["ref"].notna() & df["alt"].notna(),
    ["Chr", "Coordinate", "ref", "alt"],
].drop_duplicates()
candidate_keys = set(
    zip(
        candidate_variants["Chr"].astype(str),
        candidate_variants["Coordinate"].astype(int),
        candidate_variants["ref"].astype(str),
        candidate_variants["alt"].astype(str),
    )
)
candidate_chromosomes = sorted(candidate_variants["Chr"].dropna().astype(str).unique().tolist())

print(f"Candidate variant keys to look up: {len(candidate_keys):,}")
print(f"Chromosomes requested from AlphaMissense: {candidate_chromosomes}")

# ====================== LOAD ALPHAMISSENSE ======================
print("\nScanning AlphaMissense in chunks...")
matched_chunks = []
seen_chromosomes = set()
chunk_count = 0

for chunk in pd.read_csv(
    ALPHAMISSENSE_FILE,
    sep="\t",
    skiprows=3,
    usecols=["#CHROM", "POS", "REF", "ALT", "am_pathogenicity", "am_class"],
    dtype={
        "#CHROM": "string",
        "POS": "Int32",
        "REF": "string",
        "ALT": "string",
        "am_pathogenicity": "float32",
        "am_class": "string",
    },
    low_memory=False,
    chunksize=ALPHAMISSENSE_CHUNK_SIZE,
):
    chunk_count += 1
    chunk = chunk.rename(columns={"#CHROM": "Chr"})
    chunk["Chr"] = chunk["Chr"].str.replace(r"^chr", "", regex=True).str.strip()
    chunk["REF"] = chunk["REF"].str.upper().str.strip()
    chunk["ALT"] = chunk["ALT"].str.upper().str.strip()
    seen_chromosomes.update(chunk["Chr"].dropna().astype(str).unique().tolist())
    chunk = chunk[chunk["Chr"].isin(candidate_chromosomes)].copy()
    if chunk.empty:
        continue

    chunk["key"] = list(
        zip(
            chunk["Chr"].astype(str),
            chunk["POS"].astype("Int64"),
            chunk["REF"].astype(str),
            chunk["ALT"].astype(str),
        )
    )
    matched = chunk[chunk["key"].isin(candidate_keys)].drop(columns=["key"])
    if not matched.empty:
        matched_chunks.append(matched)

print(f"Processed AlphaMissense chunks: {chunk_count}")
print(f"Chromosomes present in AlphaMissense: {sorted(seen_chromosomes)}")

if matched_chunks:
    am_df = pd.concat(matched_chunks, ignore_index=True).drop_duplicates(
        subset=["Chr", "POS", "REF", "ALT"]
    )
else:
    am_df = pd.DataFrame(columns=["Chr", "POS", "REF", "ALT", "am_pathogenicity", "am_class"])

print(f"Relevant AlphaMissense rows retained for merge: {len(am_df):,}")

# ====================== MERGE ======================
print("\nMatching variants to AlphaMissense...")
df = df.merge(
    am_df[["Chr", "POS", "REF", "ALT", "am_pathogenicity", "am_class"]],
    left_on=["Chr", "Coordinate", "ref", "alt"],
    right_on=["Chr", "POS", "REF", "ALT"],
    how="left",
)

df = df.drop(columns=["POS", "REF", "ALT"], errors="ignore")
df["am_class"] = df["am_class"].fillna(".")
df["am_pathogenicity"] = df["am_pathogenicity"].fillna(".")

# ====================== SUMMARY ======================
matched = (df["am_class"] != ".").sum()
print(
    f"\nDone. AlphaMissense added to {matched:,} out of {len(df):,} variants "
    f"({matched / len(df):.0%} match rate)"
)

unmatched = len(df) - matched
if unmatched > 0:
    print(f"-> {unmatched:,} variants did not receive an AlphaMissense annotation")

# ====================== SAVE ======================
df.to_csv(OUTPUT_FILE, sep="\t", index=False)
print(f"\nAnnotated file saved to: {OUTPUT_FILE}")
