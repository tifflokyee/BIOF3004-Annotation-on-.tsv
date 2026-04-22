import argparse
import json
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from keras.models import load_model
from pkg_resources import resource_filename

from add_spliceai_from_vep import add_spliceai_to_tsv, build_spliceai_lookup, normalize_chrom, parse_tsv_ref_alt
from create_small_vcf_for_vep import is_clinically_relevant


ROOT = Path(__file__).resolve().parent.parent
SPLICEAI_INFO_HEADER = (
    '##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAIv1.3.1 local model annotation. '
    "Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL\">"
)


def one_hot_encode(seq):
    mapping = np.asarray(
        [
            [0, 0, 0, 0],
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ]
    )
    seq = seq.upper().replace("A", "\x01").replace("C", "\x02")
    seq = seq.replace("G", "\x03").replace("T", "\x04").replace("N", "\x00")
    return mapping[np.frombuffer(seq.encode("latin1"), dtype=np.int8) % 5]


class WindowAnnotator:
    def __init__(self, annotation="grch37"):
        if annotation == "grch37":
            annotation = resource_filename("spliceai", "annotations/grch37.txt")
        elif annotation == "grch38":
            annotation = resource_filename("spliceai", "annotations/grch38.txt")

        df = pd.read_csv(annotation, sep="\t", dtype={"CHROM": object})
        self.genes = df["#NAME"].to_numpy()
        self.chroms = df["CHROM"].to_numpy()
        self.strands = df["STRAND"].to_numpy()
        self.tx_starts = df["TX_START"].to_numpy() + 1
        self.tx_ends = df["TX_END"].to_numpy()
        self.exon_starts = [np.asarray(list(map(int, re.split(",", c)[:-1]))) + 1 for c in df["EXON_START"].to_numpy()]
        self.exon_ends = [np.asarray(list(map(int, re.split(",", c)[:-1]))) for c in df["EXON_END"].to_numpy()]

        paths = (f"models/spliceai{x}.h5" for x in range(1, 6))
        self.models = [load_model(resource_filename("spliceai", path)) for path in paths]

    def get_name_and_strand(self, chrom, pos):
        chrom = normalize_chrom(chrom)
        idxs = np.intersect1d(
            np.nonzero(self.chroms == chrom)[0],
            np.intersect1d(np.nonzero(self.tx_starts <= pos)[0], np.nonzero(pos <= self.tx_ends)[0]),
        )
        if len(idxs) >= 1:
            return self.genes[idxs], self.strands[idxs], idxs
        return [], [], []

    def get_pos_data(self, idx, pos):
        dist_tx_start = self.tx_starts[idx] - pos
        dist_tx_end = self.tx_ends[idx] - pos
        dist_exon_bdry = min(np.union1d(self.exon_starts[idx], self.exon_ends[idx]) - pos, key=abs)
        return dist_tx_start, dist_tx_end, dist_exon_bdry


def default_paths(input_tsv):
    input_tsv = Path(input_tsv)
    name = input_tsv.name
    for suffix in (".tsv.gz", ".tsv"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    return {
        "vcf": ROOT / "result" / f"{name}.local_spliceai_ucsc.vcf",
        "tsv": ROOT / "result" / f"{name}_with_local_SpliceAI.tsv",
        "cache": ROOT / "annotation" / "hg19_spliceai_windows_cache.json",
    }


def load_cache(path):
    path = Path(path)
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def save_cache(path, cache):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cache, indent=2, sort_keys=True), encoding="utf-8")


def fetch_hg19_sequence(chrom, start, end, cache):
    ucsc_chrom = chrom if str(chrom).startswith("chr") else f"chr{chrom}"
    key = f"hg19:{ucsc_chrom}:{start}-{end}"
    if key in cache:
        return cache[key]

    url = "https://api.genome.ucsc.edu/getData/sequence"
    response = requests.get(
        url,
        params={"genome": "hg19", "chrom": ucsc_chrom, "start": start, "end": end},
        timeout=60,
    )
    response.raise_for_status()
    payload = response.json()
    if "dna" not in payload:
        raise ValueError(f"UCSC did not return sequence for {key}: {payload}")
    seq = payload["dna"].upper()
    cache[key] = seq
    return seq


def spliceai_scores_for_variant(ann, chrom, pos, ref, alt, seq, dist_var, mask):
    cov = 2 * dist_var + 1
    wid = 10000 + cov
    delta_scores = []

    genes, strands, idxs = ann.get_name_and_strand(chrom, pos)
    if len(idxs) == 0:
        return delta_scores

    if seq[wid // 2 : wid // 2 + len(ref)].upper() != ref:
        logging.warning("Skipping %s:%s %s>%s (reference mismatch)", chrom, pos, ref, alt)
        return delta_scores

    if len(seq) != wid:
        logging.warning("Skipping %s:%s %s>%s (sequence window length issue)", chrom, pos, ref, alt)
        return delta_scores

    if len(ref) > 2 * dist_var:
        logging.warning("Skipping %s:%s %s>%s (ref too long)", chrom, pos, ref, alt)
        return delta_scores

    for i in range(len(idxs)):
        if any(token in alt for token in [".", "-", "*", "<", ">"]):
            continue

        if len(ref) > 1 and len(alt) > 1:
            delta_scores.append(f"{alt}|{genes[i]}|.|.|.|.|.|.|.|.")
            continue

        dist_ann = ann.get_pos_data(idxs[i], pos)
        pad_size = [max(wid // 2 + dist_ann[0], 0), max(wid // 2 - dist_ann[1], 0)]
        ref_len = len(ref)
        alt_len = len(alt)
        del_len = max(ref_len - alt_len, 0)

        x_ref = "N" * pad_size[0] + seq[pad_size[0] : wid - pad_size[1]] + "N" * pad_size[1]
        x_alt = x_ref[: wid // 2] + alt + x_ref[wid // 2 + ref_len :]

        x_ref = one_hot_encode(x_ref)[None, :]
        x_alt = one_hot_encode(x_alt)[None, :]

        if strands[i] == "-":
            x_ref = x_ref[:, ::-1, ::-1]
            x_alt = x_alt[:, ::-1, ::-1]

        y_ref = np.mean([model.predict(x_ref, verbose=0) for model in ann.models], axis=0)
        y_alt = np.mean([model.predict(x_alt, verbose=0) for model in ann.models], axis=0)

        if strands[i] == "-":
            y_ref = y_ref[:, ::-1]
            y_alt = y_alt[:, ::-1]

        if ref_len > 1 and alt_len == 1:
            y_alt = np.concatenate(
                [y_alt[:, : cov // 2 + alt_len], np.zeros((1, del_len, 3)), y_alt[:, cov // 2 + alt_len :]],
                axis=1,
            )
        elif ref_len == 1 and alt_len > 1:
            y_alt = np.concatenate(
                [
                    y_alt[:, : cov // 2],
                    np.max(y_alt[:, cov // 2 : cov // 2 + alt_len], axis=1)[:, None, :],
                    y_alt[:, cov // 2 + alt_len :],
                ],
                axis=1,
            )

        y = np.concatenate([y_ref, y_alt])

        idx_pa = (y[1, :, 1] - y[0, :, 1]).argmax()
        idx_na = (y[0, :, 1] - y[1, :, 1]).argmax()
        idx_pd = (y[1, :, 2] - y[0, :, 2]).argmax()
        idx_nd = (y[0, :, 2] - y[1, :, 2]).argmax()

        mask_pa = np.logical_and((idx_pa - cov // 2 == dist_ann[2]), mask)
        mask_na = np.logical_and((idx_na - cov // 2 != dist_ann[2]), mask)
        mask_pd = np.logical_and((idx_pd - cov // 2 == dist_ann[2]), mask)
        mask_nd = np.logical_and((idx_nd - cov // 2 != dist_ann[2]), mask)

        delta_scores.append(
            "{}|{}|{:.2f}|{:.2f}|{:.2f}|{:.2f}|{}|{}|{}|{}".format(
                alt,
                genes[i],
                (y[1, idx_pa, 1] - y[0, idx_pa, 1]) * (1 - mask_pa),
                (y[0, idx_na, 1] - y[1, idx_na, 1]) * (1 - mask_na),
                (y[1, idx_pd, 2] - y[0, idx_pd, 2]) * (1 - mask_pd),
                (y[0, idx_nd, 2] - y[1, idx_nd, 2]) * (1 - mask_nd),
                idx_pa - cov // 2,
                idx_na - cov // 2,
                idx_pd - cov // 2,
                idx_nd - cov // 2,
            )
        )

    return delta_scores


def variants_from_tsv(input_tsv):
    df = pd.read_csv(input_tsv, sep="\t", dtype=str, low_memory=False)
    df = df[df.apply(is_clinically_relevant, axis=1)].copy()
    seen = set()
    variants = []
    for _, row in df.iterrows():
        chrom = str(row.get("Chr", "")).strip()
        pos_text = str(row.get("Coordinate", row.get("POS", ""))).strip()
        ref, alt = parse_tsv_ref_alt(row)
        if not chrom or not pos_text or not ref or not alt:
            continue
        pos = int(float(pos_text))
        key = (normalize_chrom(chrom), pos, ref, alt)
        if key not in seen:
            seen.add(key)
            variants.append(key)
    return variants


def write_vcf(records, output_vcf):
    output_vcf = Path(output_vcf)
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    with output_vcf.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write(SPLICEAI_INFO_HEADER + "\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, ref, alt, scores in records:
            info = "." if not scores else "SpliceAI=" + ",".join(scores)
            handle.write(f"chr{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t{info}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Run local SpliceAI models using hg19 sequence windows fetched from UCSC."
    )
    parser.add_argument("input_tsv", help="Annotated TSV to score")
    parser.add_argument("-o", "--output", help="Output TSV path")
    parser.add_argument("--spliceai-vcf", help="Output VCF containing local SpliceAI INFO annotations")
    parser.add_argument("--cache", help="JSON cache for fetched hg19 windows")
    parser.add_argument("--distance", type=int, default=50)
    parser.add_argument("--mask", type=int, default=0, choices=[0, 1])
    args = parser.parse_args()

    defaults = default_paths(args.input_tsv)
    output_tsv = Path(args.output) if args.output else defaults["tsv"]
    output_vcf = Path(args.spliceai_vcf) if args.spliceai_vcf else defaults["vcf"]
    cache_path = Path(args.cache) if args.cache else defaults["cache"]

    variants = variants_from_tsv(args.input_tsv)
    print(f"Variants to score: {len(variants):,}")
    print("Loading local SpliceAI models...")
    ann = WindowAnnotator("grch37")
    cache = load_cache(cache_path)

    cov = 2 * args.distance + 1
    wid = 10000 + cov
    records = []
    for index, (chrom, pos, ref, alt) in enumerate(variants, 1):
        start = pos - wid // 2 - 1
        end = pos + wid // 2
        if start < 0:
            logging.warning("Skipping %s:%s %s>%s near chromosome start", chrom, pos, ref, alt)
            records.append((chrom, pos, ref, alt, []))
            continue
        print(f"[{index}/{len(variants)}] chr{chrom}:{pos} {ref}>{alt}")
        seq = fetch_hg19_sequence(chrom, start, end, cache)
        scores = spliceai_scores_for_variant(ann, chrom, pos, ref, alt, seq, args.distance, args.mask)
        records.append((chrom, pos, ref, alt, scores))

    save_cache(cache_path, cache)
    write_vcf(records, output_vcf)

    lookup = build_spliceai_lookup(output_vcf)
    add_spliceai_to_tsv(lookup, args.input_tsv, output_tsv)
    print(f"SpliceAI VCF saved to: {output_vcf}")
    print(f"Annotated TSV saved to: {output_tsv}")


if __name__ == "__main__":
    main()
