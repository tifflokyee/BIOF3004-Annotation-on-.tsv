import pandas as pd
import sys
from pathlib import Path

print("HGMD Batch File Generator")
print("=" * 50)

# ================================================================
# Get input file from user (command line or prompt)
# ================================================================

if len(sys.argv) > 1:
    # Use filename passed as argument: python generate_hgmd.py myfile.tsv
    input_tsv = sys.argv[1]
    print(f"Using file from argument: {input_tsv}")
else:
    # Ask user to type the filename
    print("Enter the path to your .tsv file (or just the filename if it's in the current folder):")
    input_tsv = input("> ").strip().strip('"\'')   # remove quotes if user copies with them

# ================================================================
# Configuration
# ================================================================
output_hgmd = f"result/{Path(input_tsv).stem}_HGMD_batch_simple.txt"

# Column names in your TSV (change only if different)
chr_column = "Chr"
pos_column = "Coordinate"

# ================================================================
# Main logic
# ================================================================
print(f"\nInput file : {input_tsv}")
print(f"Output file: {output_hgmd}\n")

# Check if file exists
if not Path(input_tsv).exists():
    print(f"ERROR: File not found: {input_tsv}")
    print("Please check the path and try again.")
    input("\nPress Enter to exit...")
    sys.exit(1)

# Load the TSV
print("Loading TSV file...")
try:
    df = pd.read_csv(input_tsv, sep='\t')
except Exception as e:
    print(f"ERROR: Failed to read the file: {e}")
    input("\nPress Enter to exit...")
    sys.exit(1)

print(f"Loaded {len(df)} variants.")

# Create HGMD format: chr1:45975088
df['Chr_with_prefix'] = 'chr' + df[chr_column].astype(str).str.replace(r'^chr', '', regex=True, flags=0)
df['HGMD_input'] = df['Chr_with_prefix'] + ':' + df[pos_column].astype(str)

# Save to file
output_path = Path(output_hgmd)
output_path.parent.mkdir(parents=True, exist_ok=True)   # create 'result' folder if not exists

df['HGMD_input'].to_csv(
    output_hgmd,
    index=False,
    header=False,
    quoting=3,
    lineterminator='\n'
)

# Verification
print(f"\n✅ HGMD batch file successfully created:")
print(f"   {output_hgmd}")
print(f"   Total lines written: {len(df)}")

print("\nFirst 12 lines for verification:")
print(df['HGMD_input'].head(12).to_string(index=False))

print("\nDone! You can now upload this file to HGMD Professional.")

input("\nPress Enter to exit...")