import pandas as pd
from pathlib import Path

print("HGMD Batch File Generator - Multi-File Mode")
print("=" * 65)

# ================================================================
# Configuration
# ================================================================

input_folder = Path(".")                    # Current folder
output_folder = Path("result")              # Output folder

# Column names in your TSV files
chr_column = "Chr"
pos_column = "Coordinate"

# ================================================================
# Create output folder if it doesn't exist
# ================================================================
output_folder.mkdir(parents=True, exist_ok=True)

# ================================================================
# Find all .tsv files (excluding already processed ones)
# ================================================================
tsv_files = list(input_folder.glob("*.tsv"))

# Optional: Skip files that already end with _hgmd
tsv_files = [f for f in tsv_files if not f.stem.lower().endswith("_hgmd")]

if not tsv_files:
    print("No .tsv files found in the current folder.")
    print("Please put your input .tsv files here and run again.")
    input("\nPress Enter to exit...")
    exit()

print(f"Found {len(tsv_files)} .tsv file(s) to process:\n")

for i, file in enumerate(tsv_files, 1):
    print(f"  {i:2d}. {file.name}")

print("\n" + "=" * 65)

# ================================================================
# Process each file
# ================================================================
successful = 0

for input_tsv in tsv_files:
    print(f"\nProcessing → {input_tsv.name}")
    
    try:
        # Load data
        df = pd.read_csv(input_tsv, sep='\t')
        print(f"   Loaded {len(df)} variants")
        
        # Create HGMD format: chr1:45975088
        df['Chr_with_prefix'] = 'chr' + df[chr_column].astype(str).str.replace(r'^chr', '', regex=True)
        df['HGMD_input'] = df['Chr_with_prefix'] + ':' + df[pos_column].astype(str)
        
        # New naming: same name + _hgmd.txt
        output_filename = f"{input_tsv.stem}_hgmd.txt"
        output_path = output_folder / output_filename
        
        # Save the file
        df['HGMD_input'].to_csv(
            output_path,
            index=False,
            header=False,
            quoting=3,
            lineterminator='\n'
        )
        
        print(f"   ✅ Success → {output_filename}")
        print(f"      {len(df)} lines written")
        
        # Show first few lines for verification
        print("      Preview (first 5):")
        print(df['HGMD_input'].head(5).to_string(index=False))
        
        successful += 1
        
    except Exception as e:
        print(f"   ❌ Error processing {input_tsv.name}: {e}")

# ================================================================
# Summary
# ================================================================
print("\n" + "=" * 65)
print(f"Batch processing finished!")
print(f"Successfully created: {successful} / {len(tsv_files)} files")
print(f"Output folder: {output_folder.resolve()}")

print("\n✅ All done! You can now upload these *_hgmd.txt files to HGMD.")