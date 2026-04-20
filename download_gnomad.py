import argparse
import csv
import os
import sys
import urllib.request
from pathlib import Path


LATEST_RELEASE = "4.1"
LEGACY_RELEASE = "2.1.1"
GCP_BASE = "https://storage.googleapis.com/gcp-public-data--gnomad/release"


def infer_chromosomes_from_tsv(tsv_path):
    chromosomes = []
    seen = set()

    with open(tsv_path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "Chr" not in reader.fieldnames:
            raise KeyError("Input TSV must contain a 'Chr' column.")

        for row in reader:
            chrom = str(row["Chr"]).strip()
            chrom = chrom.removeprefix("chr").strip()
            if chrom == "M":
                chrom = "MT"
            if chrom and chrom not in seen:
                chromosomes.append(chrom)
                seen.add(chrom)

    return chromosomes


def latest_joint_urls(chromosomes):
    urls = []
    for chrom in chromosomes:
        name = f"gnomad.joint.v{LATEST_RELEASE}.sites.chr{chrom}.vcf.bgz"
        base = f"{GCP_BASE}/{LATEST_RELEASE}/vcf/joint/{name}"
        urls.append(base)
        urls.append(f"{base}.tbi")
    return urls


def legacy_grch37_urls():
    name = f"gnomad.exomes.r{LEGACY_RELEASE}.sites.vcf.bgz"
    base = f"{GCP_BASE}/{LEGACY_RELEASE}/vcf/exomes/{name}"
    return [base, f"{base}.tbi"]


def download_file(url, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {url}")
    print(f"-> {destination}")
    with urllib.request.urlopen(url) as response, open(destination, "wb") as output_handle:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            output_handle.write(chunk)


def build_argument_parser():
    parser = argparse.ArgumentParser(
        description="Download official gnomAD files into the annotation directory."
    )
    parser.add_argument(
        "--mode",
        choices=["latest_joint", "legacy_grch37"],
        default="latest_joint",
        help="latest_joint downloads current joint GRCh38 VCFs legacy_grch37 downloads v2.1.1 exomes.",
    )
    parser.add_argument(
        "--input-tsv",
        default="./result/variants_with_AlphaMissense_and_ClinVar.tsv",
        help="TSV used to infer chromosomes for latest_joint mode.",
    )
    parser.add_argument(
        "--scope",
        choices=["matched", "all"],
        default="matched",
        help="matched downloads only chromosomes present in the TSV; all downloads chr1-22,X,Y.",
    )
    parser.add_argument(
        "--output-dir",
        default="./annotation",
        help="Directory where downloaded files will be stored.",
    )
    return parser


def main():
    args = build_argument_parser().parse_args()
    output_dir = Path(args.output_dir)

    if args.mode == "latest_joint":
        if args.scope == "all":
            chromosomes = [str(i) for i in range(1, 23)] + ["X", "Y"]
        else:
            chromosomes = infer_chromosomes_from_tsv(args.input_tsv)
        urls = latest_joint_urls(chromosomes)
        print(
            "Selected latest gnomAD joint release. "
            "This is GRCh38 and will only match a GRCh38 TSV."
        )
    else:
        urls = legacy_grch37_urls()
        print(
            "Selected gnomAD v2.1.1 exomes. "
            "This is the compatible GRCh37 option for your current hg19-style pipeline."
        )

    for url in urls:
        destination = output_dir / Path(url).name
        if destination.exists() and destination.stat().st_size > 0:
            print(f"Skipping existing file: {destination}")
            continue
        download_file(url, destination)

    print("\nDownload complete.")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(130)
    except Exception as exc:
        print(f"Download failed: {exc}", file=sys.stderr)
        sys.exit(1)
