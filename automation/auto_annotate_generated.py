#!/usr/bin/env python3
import argparse
import gzip
import os
import platform
import re
import shlex
import subprocess
from collections import defaultdict
from pathlib import Path

import pandas as pd


PANEL_COLUMNS = [
    "Mode of Inheritance (Mendeliome)",
    "Incidentalome Status",
    "Paediatric Additional Status",
    "Mendeliome Evidence",
    "Incidentalome Evidence",
    "Paediatric Evidence",
]

CLINVAR_COLUMNS = [
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

GNOMAD_FIELDS = [
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


def open_text(path):
    path = Path(path)
    if path.suffix in {".gz", ".bgz"}:
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def clean_text(value):
    if value is None or pd.isna(value):
        return "."
    text = str(value).strip()
    return text if text else "."


def normalize_chromosome(value):
    chrom = str(value).strip()
    chrom = re.sub(r"^chr", "", chrom, flags=re.IGNORECASE)
    return "MT" if chrom == "M" else chrom


def parse_info_field(info_text):
    info = {}
    for item in str(info_text).split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        elif item:
            info[item] = True
    return info


def pick_first_existing(df, names):
    for name in names:
        if name in df.columns:
            return name
    return None


def infer_variant_type(ref, alt):
    ref = "" if pd.isna(ref) else str(ref)
    alt = "" if pd.isna(alt) else str(alt)
    if len(ref) == len(alt) == 1:
        return "SNV"
    if len(ref) < len(alt):
        return "insertion"
    if len(ref) > len(alt):
        return "deletion"
    return "substitution"


def extract_gene_from_info(info, alt_allele):
    for key in ("Gene", "GENE", "SYMBOL", "HGNC"):
        if info.get(key):
            return str(info[key]).split(",")[0].strip()

    geneinfo = info.get("GENEINFO")
    if geneinfo:
        return ",".join(
            sorted({part.split(":", 1)[0].strip() for part in str(geneinfo).split("|") if part})
        )

    ann = info.get("ANN")
    if ann:
        genes = []
        for record in str(ann).split(","):
            fields = record.split("|")
            if len(fields) > 3 and (not alt_allele or fields[0] == alt_allele):
                genes.append(fields[3])
        genes = sorted({gene for gene in genes if gene})
        if genes:
            return ",".join(genes)

    csq = info.get("CSQ")
    if csq:
        genes = []
        for record in str(csq).split(","):
            fields = record.split("|")
            if len(fields) > 3 and (not alt_allele or fields[0] == alt_allele):
                genes.append(fields[3])
        genes = sorted({gene for gene in genes if gene})
        if genes:
            return ",".join(genes)

    return "."


def vcf_to_dataframe(input_path):
    rows = []
    sample_names = []
    with open_text(input_path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                sample_names = fields[9:]
                continue
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            chrom, pos, record_id, ref, alts, qual, filt, info_text = fields[:8]
            fmt = fields[8] if len(fields) > 8 else ""
            sample_values = fields[9:] if len(fields) > 9 else []
            info = parse_info_field(info_text)
            format_keys = fmt.split(":") if fmt else []

            for alt in alts.split(","):
                alt = alt.strip().upper()
                row = {
                    "Gene": extract_gene_from_info(info, alt),
                    "Chr": normalize_chromosome(chrom),
                    "Coordinate": int(pos),
                    "Variant": f"{ref}>{alt}",
                    "ref": ref.strip().upper(),
                    "alt": alt,
                    "Type": infer_variant_type(ref, alt),
                    "Filters": filt,
                    "Quality": qual,
                    "dbSNP ID": clean_text(record_id),
                    "vcf_info": info_text,
                }

                if sample_names and sample_values:
                    genotypes = []
                    for sample_name, sample_value in zip(sample_names, sample_values):
                        values = sample_value.split(":")
                        sample_map = dict(zip(format_keys, values))
                        if "GT" in sample_map:
                            genotypes.append(f"{sample_name}:{sample_map['GT']}")
                    row["Genotype"] = ";".join(genotypes) if genotypes else "."

                rows.append(row)

    return pd.DataFrame(rows)


def load_variant_table(input_path):
    input_path = Path(input_path)
    name = input_path.name.lower()
    if name.endswith(".vcf") or name.endswith(".vcf.gz") or name.endswith(".vcf.bgz"):
        print("Input detected as VCF; converting records to an annotation table.")
        df = vcf_to_dataframe(input_path)
    else:
        print("Input detected as TSV.")
        df = pd.read_csv(input_path, sep="\t", low_memory=False)

    if "Chr" not in df.columns:
        chrom_col = pick_first_existing(df, ["#CHROM", "CHROM", "Chromosome", "chrom", "chr"])
        if chrom_col:
            df["Chr"] = df[chrom_col]
    if "Coordinate" not in df.columns:
        pos_col = pick_first_existing(df, ["POS", "Position", "pos", "Start"])
        if pos_col:
            df["Coordinate"] = df[pos_col]
    if "ref" not in df.columns:
        ref_col = pick_first_existing(df, ["REF", "Reference", "reference"])
        if ref_col:
            df["ref"] = df[ref_col]
    if "alt" not in df.columns:
        alt_col = pick_first_existing(df, ["ALT", "Alternate", "alternate"])
        if alt_col:
            df["alt"] = df[alt_col]

    if ("ref" not in df.columns or "alt" not in df.columns) and "Variant" in df.columns:
        parts = df["Variant"].astype("string").str.split(">", n=1, expand=True)
        if parts.shape[1] == 2:
            df["ref"] = parts[0].str.strip().str.upper()
            alt_options = parts[1].fillna("").astype("string").str.split("/", expand=True)

            def pick_alt(row):
                ref = row["ref"]
                for value in row.drop(labels=["ref"]).tolist():
                    if pd.isna(value):
                        continue
                    allele = str(value).strip().upper()
                    if allele and allele != ref:
                        return allele
                return pd.NA

            alt_candidates = alt_options.copy()
            alt_candidates.insert(0, "ref", df["ref"])
            df["alt"] = alt_candidates.apply(pick_alt, axis=1).astype("string")

    required = ["Chr", "Coordinate", "ref", "alt"]
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise KeyError(
            "Input must contain or allow inference of these columns: "
            + ", ".join(missing)
            + ". For TSV input, provide Chr/Coordinate/ref/alt or a Variant column like C>C/T."
        )

    if "Gene" not in df.columns:
        df["Gene"] = "."
    if "Variant" not in df.columns:
        df["Variant"] = df["ref"].astype(str) + ">" + df["alt"].astype(str)

    df["Chr"] = df["Chr"].astype("string").map(normalize_chromosome)
    df["Coordinate"] = pd.to_numeric(df["Coordinate"], errors="coerce").astype("Int64")
    df["ref"] = df["ref"].astype("string").str.upper().str.strip()
    df["alt"] = df["alt"].astype("string").str.upper().str.strip()
    df["Gene"] = df["Gene"].astype("string").str.strip()
    return df


def create_panel_dict(panel_file):
    df_panel = pd.read_csv(panel_file, sep="\t", low_memory=False)

    def map_status(value):
        if pd.isna(value):
            return "Unknown"
        try:
            value = int(value)
            if value == 3:
                return "Green"
            if value == 2:
                return "Amber"
            return "Red/Low"
        except Exception:
            return str(value)

    df_panel = df_panel.dropna(subset=["Gene Symbol"]).copy()
    df_panel["label"] = (
        df_panel["Model_Of_Inheritance"].fillna(".").astype(str)
        + " ("
        + df_panel["GEL_Status"].apply(map_status)
        + ")"
    )
    return dict(zip(df_panel["Gene Symbol"].astype(str).str.upper(), df_panel["label"]))


def get_panel_status(gene_str, panel_dict):
    if pd.isna(gene_str) or str(gene_str).strip() in {"", "."}:
        return "."
    results = []
    for gene in re.split(r"[,;|]", str(gene_str)):
        gene = gene.strip().upper()
        if gene:
            results.append(panel_dict.get(gene, "."))
    return "; ".join(results) if results else "."


def annotate_panelapp(df, annotation_dir):
    panel_files = {
        "Mode of Inheritance (Mendeliome)": annotation_dir / "Mendeliome.tsv",
        "Incidentalome Status": annotation_dir / "Incidentalome.tsv",
        "Paediatric Additional Status": annotation_dir / "Additional findings_Paediatric.tsv",
    }
    if not all(path.exists() for path in panel_files.values()):
        print("PanelApp files missing; PanelApp columns will be filled with '.'.")
        for column in PANEL_COLUMNS:
            df[column] = "."
        return df

    print("Annotating PanelApp gene panels.")
    for output_column, panel_file in panel_files.items():
        panel_dict = create_panel_dict(panel_file)
        df[output_column] = df["Gene"].apply(lambda gene: get_panel_status(gene, panel_dict))

    evidence_map = {
        "Mendeliome Evidence": "Mode of Inheritance (Mendeliome)",
        "Incidentalome Evidence": "Incidentalome Status",
        "Paediatric Evidence": "Paediatric Additional Status",
    }
    for evidence_col, source_col in evidence_map.items():
        df[evidence_col] = df[source_col].astype("string").str.extract(r"\(([^)]+)\)", expand=False).fillna(".")
    return df


def candidate_variants(df, snv_only=False):
    mask = df["Coordinate"].notna() & df["ref"].notna() & df["alt"].notna()
    if snv_only:
        mask &= df["ref"].str.len().eq(1) & df["alt"].str.len().eq(1)
    variants = df.loc[mask, ["Chr", "Coordinate", "ref", "alt"]].drop_duplicates().copy()
    keys = set(
        zip(
            variants["Chr"].astype(str),
            variants["Coordinate"].astype(int),
            variants["ref"].astype(str),
            variants["alt"].astype(str),
        )
    )
    chromosomes = sorted(variants["Chr"].dropna().astype(str).unique().tolist())
    return variants, keys, chromosomes


def annotate_revel(df, annotation_dir, chunk_size):
    revel_file = annotation_dir / "revel_with_transcript_ids"
    if not revel_file.exists():
        print("REVEL file missing; REVEL_score will be filled with '.'.")
        df["REVEL_score"] = "."
        return df

    _, keys, chromosomes = candidate_variants(df, snv_only=True)
    print(f"Annotating REVEL for {len(keys):,} candidate SNVs.")
    matched_chunks = []
    for chunk in pd.read_csv(
        revel_file,
        sep=",",
        usecols=["chr", "hg19_pos", "ref", "alt", "REVEL"],
        dtype={"chr": "string", "hg19_pos": "Int64", "ref": "string", "alt": "string", "REVEL": "float32"},
        na_values=["."],
        low_memory=False,
        chunksize=chunk_size,
    ):
        chunk["chr"] = chunk["chr"].map(normalize_chromosome)
        chunk = chunk[chunk["chr"].isin(chromosomes)].copy()
        if chunk.empty:
            continue
        chunk["key"] = list(
            zip(chunk["chr"].astype(str), chunk["hg19_pos"].astype("Int64"), chunk["ref"].astype(str), chunk["alt"].astype(str))
        )
        matched = chunk[chunk["key"].isin(keys)].drop(columns=["key"])
        if not matched.empty:
            matched_chunks.append(matched)

    if matched_chunks:
        revel_df = pd.concat(matched_chunks, ignore_index=True).drop_duplicates(["chr", "hg19_pos", "ref", "alt"])
    else:
        revel_df = pd.DataFrame(columns=["chr", "hg19_pos", "ref", "alt", "REVEL"])

    df = df.merge(
        revel_df[["chr", "hg19_pos", "ref", "alt", "REVEL"]],
        left_on=["Chr", "Coordinate", "ref", "alt"],
        right_on=["chr", "hg19_pos", "ref", "alt"],
        how="left",
    ).drop(columns=["chr", "hg19_pos"], errors="ignore")
    df = df.rename(columns={"REVEL": "REVEL_score"})
    df["REVEL_score"] = df["REVEL_score"].fillna(".")
    print(f"REVEL matches: {(df['REVEL_score'] != '.').sum():,} / {len(df):,}")
    return df


def annotate_alphamissense(df, annotation_dir, chunk_size):
    am_file = annotation_dir / "AlphaMissense_hg19.tsv.gz"
    if not am_file.exists():
        print("AlphaMissense file missing; AlphaMissense columns will be filled with '.'.")
        df["am_pathogenicity"] = "."
        df["am_class"] = "."
        return df

    _, keys, chromosomes = candidate_variants(df)
    print(f"Annotating AlphaMissense for {len(keys):,} candidate variants.")
    matched_chunks = []
    for chunk in pd.read_csv(
        am_file,
        sep="\t",
        skiprows=3,
        usecols=["#CHROM", "POS", "REF", "ALT", "am_pathogenicity", "am_class"],
        dtype={"#CHROM": "string", "POS": "Int64", "REF": "string", "ALT": "string", "am_pathogenicity": "float32", "am_class": "string"},
        low_memory=False,
        chunksize=chunk_size,
    ):
        chunk = chunk.rename(columns={"#CHROM": "Chr"})
        chunk["Chr"] = chunk["Chr"].map(normalize_chromosome)
        chunk = chunk[chunk["Chr"].isin(chromosomes)].copy()
        if chunk.empty:
            continue
        chunk["key"] = list(
            zip(chunk["Chr"].astype(str), chunk["POS"].astype("Int64"), chunk["REF"].astype(str), chunk["ALT"].astype(str))
        )
        matched = chunk[chunk["key"].isin(keys)].drop(columns=["key"])
        if not matched.empty:
            matched_chunks.append(matched)

    if matched_chunks:
        am_df = pd.concat(matched_chunks, ignore_index=True).drop_duplicates(["Chr", "POS", "REF", "ALT"])
    else:
        am_df = pd.DataFrame(columns=["Chr", "POS", "REF", "ALT", "am_pathogenicity", "am_class"])

    df = df.merge(
        am_df[["Chr", "POS", "REF", "ALT", "am_pathogenicity", "am_class"]],
        left_on=["Chr", "Coordinate", "ref", "alt"],
        right_on=["Chr", "POS", "REF", "ALT"],
        how="left",
    ).drop(columns=["POS", "REF", "ALT"], errors="ignore")
    df["am_pathogenicity"] = df["am_pathogenicity"].fillna(".")
    df["am_class"] = df["am_class"].fillna(".")
    print(f"AlphaMissense matches: {(df['am_class'] != '.').sum():,} / {len(df):,}")
    return df


def parse_gene_symbols(geneinfo_text):
    if not geneinfo_text:
        return set()
    return {item.split(":", 1)[0].strip().upper() for item in str(geneinfo_text).split("|") if item}


def annotate_clinvar(df, annotation_dir):
    clinvar_file = annotation_dir / "clinvar.vcf.gz"
    if not clinvar_file.exists():
        print("ClinVar VCF missing; ClinVar columns will be filled with '.'.")
        for column in CLINVAR_COLUMNS:
            df[column] = "."
        return df

    candidate_rows = df.loc[
        df["Gene"].notna() & df["Gene"].ne(".") & df["Coordinate"].notna() & df["ref"].notna() & df["alt"].notna(),
        ["Gene", "Chr", "Coordinate", "ref", "alt"],
    ].drop_duplicates()
    candidate_keys = set()
    for _, row in candidate_rows.iterrows():
        for gene in re.split(r"[,;|]", str(row["Gene"])):
            gene = gene.strip().upper()
            if gene:
                candidate_keys.add((gene, str(row["Chr"]), int(row["Coordinate"]), str(row["ref"]), str(row["alt"])))
    variant_keys = {(chrom, pos, ref, alt) for _, chrom, pos, ref, alt in candidate_keys}
    print(f"Annotating ClinVar for {len(variant_keys):,} candidate variants.")

    matches_by_key = defaultdict(list)
    with open_text(clinvar_file) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            chrom, pos, record_id, ref, alt, _qual, _filt, info_text = line.rstrip("\n").split("\t", 8)
            chrom = normalize_chromosome(chrom)
            pos = int(pos)
            ref = ref.strip().upper()
            for alt_allele in alt.split(","):
                alt_allele = alt_allele.strip().upper()
                if (chrom, pos, ref, alt_allele) not in variant_keys:
                    continue
                info = parse_info_field(info_text)
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
                    "clinvar_vcf_geneinfo": clean_text(info.get("GENEINFO")),
                    "clinvar_vcf_rs": clean_text(info.get("RS")),
                    "clinvar_vcf_info": clean_text(info_text),
                }
                for gene in parse_gene_symbols(info.get("GENEINFO")):
                    key = (gene, chrom, pos, ref, alt_allele)
                    if key in candidate_keys:
                        matches_by_key[key].append(record)

    annotated_records = []
    for key, records in matches_by_key.items():
        merged = {}
        for column in CLINVAR_COLUMNS:
            values = []
            seen = set()
            for record in records:
                value = record.get(column, ".")
                if value == "." or value in seen:
                    continue
                values.append(value)
                seen.add(value)
            merged[column] = " | ".join(values) if values else "."
        gene, chrom, pos, ref, alt = key
        merged.update({"Gene": gene, "Chr": chrom, "Coordinate": pos, "ref": ref, "alt": alt})
        annotated_records.append(merged)

    clinvar_df = pd.DataFrame(annotated_records)
    if clinvar_df.empty:
        clinvar_df = pd.DataFrame(columns=["Gene", "Chr", "Coordinate", "ref", "alt", *CLINVAR_COLUMNS])

    merge_rows = []
    for _, row in df.iterrows():
        genes = [gene.strip().upper() for gene in re.split(r"[,;|]", str(row["Gene"])) if gene.strip()]
        if not genes:
            genes = [str(row["Gene"]).strip().upper()]
        for gene in genes:
            new_row = row.copy()
            new_row["Gene"] = gene
            merge_rows.append(new_row)
    expanded = pd.DataFrame(merge_rows)
    expanded = expanded.merge(clinvar_df, on=["Gene", "Chr", "Coordinate", "ref", "alt"], how="left")
    for column in CLINVAR_COLUMNS:
        expanded[column] = expanded[column].fillna(".")
    print(f"ClinVar matches: {(expanded['clinvar_vcf_id'] != '.').sum():,} / {len(expanded):,}")
    return expanded


def gnomad_info_value(info, info_key, allele_index, alt_count, source_name, info_text, filt):
    if info_key == "__FILTER__":
        return filt
    if info_key == "__INFO__":
        return info_text
    if info_key == "__SOURCE__":
        return source_name
    value = info.get(info_key)
    if value is None:
        return None
    if value is True:
        return "TRUE"
    text = str(value).strip()
    if "," not in text:
        return text
    parts = [part.strip() for part in text.split(",")]
    return parts[allele_index] if len(parts) == alt_count else text


def pick_gnomad_files(annotation_dir, chromosomes):
    files = sorted(annotation_dir.glob("gnomad*.vcf.gz")) + sorted(annotation_dir.glob("gnomad*.vcf.bgz"))
    selected = []
    wanted = {chrom.lower() for chrom in chromosomes}
    for file_path in files:
        name = file_path.name.lower()
        if ".chr" not in name:
            selected.append(file_path)
        elif any(f".chr{chrom}." in name for chrom in wanted):
            selected.append(file_path)
    return selected or files


def annotate_gnomad(df, annotation_dir):
    _, keys, chromosomes = candidate_variants(df)
    files = pick_gnomad_files(annotation_dir, chromosomes)
    value_columns = [name for name, _ in GNOMAD_FIELDS]
    if not files:
        print("gnomAD files missing; gnomAD columns will be filled with '.'.")
        for column in value_columns:
            df[column] = "."
        return df

    print(f"Annotating gnomAD from {len(files)} local file(s).")
    matches = []
    for file_path in files:
        source_name = file_path.name
        with open_text(file_path) as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                chrom, pos, _record_id, ref, alt, _qual, filt, info_text = line.rstrip("\n").split("\t", 8)
                chrom = normalize_chromosome(chrom)
                if chrom not in chromosomes:
                    continue
                pos = int(pos)
                ref = ref.strip().upper()
                alt_alleles = [item.strip().upper() for item in alt.split(",")]
                info = None
                for allele_index, alt_allele in enumerate(alt_alleles):
                    key = (chrom, pos, ref, alt_allele)
                    if key not in keys:
                        continue
                    if info is None:
                        info = parse_info_field(info_text)
                    row = {"Chr": chrom, "Coordinate": pos, "ref": ref, "alt": alt_allele}
                    for output_column, info_key in GNOMAD_FIELDS:
                        row[output_column] = clean_text(
                            gnomad_info_value(info, info_key, allele_index, len(alt_alleles), source_name, info_text, filt)
                        )
                    matches.append(row)

    if matches:
        gnomad_df = pd.DataFrame(matches).drop_duplicates(["Chr", "Coordinate", "ref", "alt"])
    else:
        gnomad_df = pd.DataFrame(columns=["Chr", "Coordinate", "ref", "alt", *value_columns])

    df = df.merge(gnomad_df, on=["Chr", "Coordinate", "ref", "alt"], how="left")
    for column in value_columns:
        df[column] = df[column].fillna(".")
    print(f"gnomAD matches: {(df['gnomad_af'] != '.').sum():,} / {len(df):,}")
    return df


def run_local_spliceai(input_tsv, output_tsv, spliceai_vcf=None, cache=None, distance=50, mask=0):
    script = Path(__file__).resolve().parent / "run_local_spliceai_ucsc.py"
    if spliceai_vcf is None:
        spliceai_vcf = Path(output_tsv).with_suffix(".local_spliceai.vcf")

    if platform.system() == "Windows" and platform.release() == "XP":
        spliceai_vcf = Path(spliceai_vcf)
        if spliceai_vcf.exists():
            from add_spliceai_from_vep import add_spliceai_to_tsv, build_spliceai_lookup

            print(f"Windows XP detected; using existing local SpliceAI VCF: {spliceai_vcf}")
            lookup = build_spliceai_lookup(spliceai_vcf)
            add_spliceai_to_tsv(lookup, input_tsv, output_tsv)
            return
        raise RuntimeError(
            "Local SpliceAI model inference cannot run on Windows XP because it requires a modern "
            "Python/TensorFlow environment. Provide --local-spliceai-vcf with an existing local "
            "SpliceAI VCF, or run the SpliceAI step on Windows 10/11 or Linux."
        )

    python_command = os.environ.get("SPLICEAI_PYTHON", "py -3.13")
    command = shlex.split(python_command) + [
        str(script),
        str(input_tsv),
        "-o",
        str(output_tsv),
        "--spliceai-vcf",
        str(spliceai_vcf),
        "--distance",
        str(distance),
        "--mask",
        str(mask),
    ]
    if cache:
        command.extend(["--cache", str(cache)])

    print("Running local SpliceAI model annotation...")
    print(" ".join(command))
    subprocess.run(command, check=True)


def build_parser():
    parser = argparse.ArgumentParser(
        description="Automatically annotate a TSV or VCF(.gz/.bgz) using local PanelApp, REVEL, AlphaMissense, ClinVar, gnomAD, and local SpliceAI."
    )
    parser.add_argument("input", help="Input .tsv, .vcf, .vcf.gz, or .vcf.bgz file")
    parser.add_argument("-o", "--output", help="Output TSV path. Default: result/<input>.all_annotations.tsv")
    parser.add_argument("--annotation-dir", default="annotation", help="Directory containing local annotation files")
    parser.add_argument("--local-spliceai-vcf", help="Output VCF containing local SpliceAI annotations")
    parser.add_argument("--local-spliceai-cache", help="JSON cache for fetched hg19 sequence windows")
    parser.add_argument("--spliceai-distance", type=int, default=50)
    parser.add_argument("--spliceai-mask", type=int, default=0, choices=[0, 1])
    parser.add_argument("--revel-chunk-size", type=int, default=500000)
    parser.add_argument("--alphamissense-chunk-size", type=int, default=500000)
    parser.add_argument("--skip-panelapp", action="store_true")
    parser.add_argument("--skip-revel", action="store_true")
    parser.add_argument("--skip-alphamissense", action="store_true")
    parser.add_argument("--skip-clinvar", action="store_true")
    parser.add_argument("--skip-gnomad", action="store_true")
    parser.add_argument("--skip-spliceai", action="store_true")
    return parser


def default_output_path(input_path):
    input_path = Path(input_path)
    name = input_path.name
    for suffix in (".vcf.gz", ".vcf.bgz", ".tsv.gz", ".tsv", ".vcf"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    return Path("result") / f"{name}.all_annotations.tsv"


def main():
    args = build_parser().parse_args()
    annotation_dir = Path(args.annotation_dir)
    output = Path(args.output) if args.output else default_output_path(args.input)
    output.parent.mkdir(parents=True, exist_ok=True)

    print(f"Loading input: {args.input}")
    df = load_variant_table(args.input)
    print(f"Variants loaded: {len(df):,}")

    if not args.skip_panelapp:
        df = annotate_panelapp(df, annotation_dir)
    if not args.skip_revel:
        df = annotate_revel(df, annotation_dir, args.revel_chunk_size)
    if not args.skip_alphamissense:
        df = annotate_alphamissense(df, annotation_dir, args.alphamissense_chunk_size)
    if not args.skip_clinvar:
        df = annotate_clinvar(df, annotation_dir)
    if not args.skip_gnomad:
        df = annotate_gnomad(df, annotation_dir)
    if args.skip_spliceai:
        df.to_csv(output, sep="\t", index=False)
        print(f"Annotated file saved to: {output}")
    else:
        pre_spliceai_output = output.with_name(output.stem + ".pre_spliceai.tsv")
        df.to_csv(pre_spliceai_output, sep="\t", index=False)
        print(f"Pre-SpliceAI annotations saved to: {pre_spliceai_output}")
        run_local_spliceai(
            pre_spliceai_output,
            output,
            args.local_spliceai_vcf,
            args.local_spliceai_cache,
            args.spliceai_distance,
            args.spliceai_mask,
        )
        print(f"Annotated file saved to: {output}")


if __name__ == "__main__":
    main()
