#!/usr/bin/env python3
"""Resolve as many WIC H3N2 strains as possible to public NCBI accessions."""

from __future__ import annotations

import argparse
import json
import os
import time
import urllib.request
from pathlib import Path

import pandas as pd
from Bio import Entrez

from common import load_config, normalize_strain_name
from wic_name_utils import canonical_strain, is_searchable_lookup_strain, normalize_key, search_variants


NON_HA_KEYWORDS = (
    "neuraminidase",
    "matrix",
    "polymerase",
    "nucleoprotein",
    "nonstructural",
    "ns1",
    "ns2",
    "pa-x",
    "pb1-f2",
)

HA_HINTS = (
    "hemagglutinin",
    "ha ",
    " ha",
    "ha0",
    "ha1",
    "segment 4",
    "ha gene",
)

HA_KEYWORDS = ("hemagglutinin", "ha", "ha0", "ha1", "hemaglutinin")
BV_BRC_API = "https://www.bv-brc.org/api/genome_feature/"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        default="H3N2-WIC/config/h3n2_wic_config.json",
        help="Path to WIC config JSON.",
    )
    parser.add_argument(
        "--manifest",
        default="H3N2-WIC/manifests/intermediate/H3N2_WIC_seq_manifest_repofill.tsv",
        help="Input manifest to resolve after repo-FASTA matching.",
    )
    parser.add_argument(
        "--output-manifest",
        default="H3N2-WIC/manifests/H3N2_WIC_seq_manifest_publicfill.tsv",
        help="Updated manifest after public accession resolution.",
    )
    parser.add_argument(
        "--cache",
        default="H3N2-WIC/output/H3N2_WIC_public_resolution_cache.tsv",
        help="Cache of public lookup results by canonical lookup strain.",
    )
    parser.add_argument(
        "--remaining",
        default="H3N2-WIC/metadata/intermediate/H3N2_WIC_public_unresolved.tsv",
        help="Remaining unresolved manifest rows after public resolution.",
    )
    parser.add_argument(
        "--gisaid-query",
        default="H3N2-WIC/metadata/H3N2_WIC_gisaid_remaining_query.txt",
        help="Remaining canonical strain list for batch GISAID search.",
    )
    parser.add_argument("--email", required=True, help="NCBI Entrez email.")
    parser.add_argument("--api-key", default=None, help="NCBI Entrez API key.")
    parser.add_argument(
        "--retmax",
        type=int,
        default=20,
        help="Maximum NCBI search hits to inspect per query.",
    )
    parser.add_argument(
        "--top-candidates",
        type=int,
        default=5,
        help="Top candidate accessions to validate per search term.",
    )
    parser.add_argument(
        "--sleep-seconds",
        type=float,
        default=0.34,
        help="Delay between NCBI requests (ignored if --api-key and lower than 0.11).",
    )
    parser.add_argument(
        "--max-strains",
        type=int,
        default=None,
        help="Optional limit on unresolved lookup strains to process this run.",
    )
    parser.add_argument(
        "--refresh-cache",
        action="store_true",
        help="Ignore any existing cached public lookup results.",
    )
    return parser.parse_args()


def choose_ha_record(records):
    preferred = []
    fallback = []
    for rec in records:
        desc = rec.description.lower()
        if any(key in desc for key in HA_KEYWORDS):
            preferred.append(rec)
        else:
            fallback.append(rec)
    if preferred:
        return max(preferred, key=lambda r: len(r.seq))
    if fallback:
        return max(fallback, key=lambda r: len(r.seq))
    return None

def safe_year(value) -> str:
    if pd.isna(value):
        return ""
    s = str(value).strip()
    if s in {"", "<NA>", "nan", "None"}:
        return ""
    try:
        return str(int(float(s)))
    except ValueError:
        return s


def is_blank(value) -> bool:
    if pd.isna(value):
        return True
    return str(value).strip() in {"", "<NA>", "nan", "None"}

def search_terms(lookup_strain: str, year: str) -> list[str]:
    year_term = f' AND "{year}"[All Fields]' if year else ""
    terms = []
    for variant in search_variants(lookup_strain):
        quoted = f'"{variant}"[All Fields]'
        candidate_terms = [
            f"{quoted}{year_term} AND influenza a virus[Organism] AND (hemagglutinin[Title] OR \"segment 4\"[Title])",
            f"{quoted} AND influenza a virus[Organism] AND (hemagglutinin[Title] OR \"segment 4\"[Title])",
            f"{quoted}{year_term} AND influenza a virus[Organism] AND H3N2[All Fields]",
            f"{quoted} AND influenza a virus[Organism]",
        ]
        for term in candidate_terms:
            if term not in terms:
                terms.append(term)
    return terms


def fetch_search_ids(term: str, retmax: int) -> list[str]:
    handle = Entrez.esearch(db="nucleotide", term=term, retmax=retmax, sort="relevance")
    try:
        record = Entrez.read(handle)
    finally:
        handle.close()
    return list(record.get("IdList", []))


def fetch_summaries(id_list: list[str]) -> list[dict]:
    if not id_list:
        return []
    handle = Entrez.esummary(db="nucleotide", id=",".join(id_list), retmode="xml")
    try:
        record = Entrez.read(handle)
    finally:
        handle.close()
    if isinstance(record, list):
        return list(record)
    if isinstance(record, dict) and "DocumentSummarySet" in record:
        return list(record["DocumentSummarySet"]["DocumentSummary"])
    return []


def retry_entrez_call(func, *args, retries: int = 3, sleep_seconds: float = 1.0, **kwargs):
    last_exc = None
    for attempt in range(1, retries + 1):
        try:
            return func(*args, **kwargs)
        except Exception as exc:  # noqa: BLE001
            last_exc = exc
            if attempt == retries:
                raise
            time.sleep(sleep_seconds * attempt)
    raise last_exc


def candidate_accession(summary: dict) -> str:
    for key in ("AccessionVersion", "Caption"):
        value = str(summary.get(key, "")).strip()
        if value:
            return value.split(".")[0].upper()
    return ""


def candidate_title(summary: dict) -> str:
    for key in ("Title",):
        value = str(summary.get(key, "")).strip()
        if value:
            return value
    return ""


def score_candidate(lookup_strain: str, year: str, title: str) -> int:
    title_upper = title.upper()
    normalized_title = normalize_strain_name(title_upper)
    lookup_key = normalize_strain_name(lookup_strain)
    score = 0

    if lookup_key and lookup_key in normalized_title:
        score += 200
    if lookup_strain.upper() in title_upper:
        score += 100
    if any(hint in title.lower() for hint in HA_HINTS):
        score += 40
    if "WHOLE GENOME" in title_upper or "COMPLETE GENOME" in title_upper:
        score += 15
    if "H3N2" in title_upper or "H3" in title_upper:
        score += 10
    if year and year in title_upper:
        score += 5
    if any(bad in title.lower() for bad in NON_HA_KEYWORDS):
        score -= 80
    return score


def fetch_bvbrc_candidates(lookup_strain: str, limit: int = 10) -> list[dict]:
    query = (
        f"and(eq(product,hemagglutinin),keyword({lookup_strain}))"
        f"&select(accession,genome_name,product,annotation)&limit({limit})"
        "&http_accept=application/json"
    )
    url = BV_BRC_API + "?" + query
    with urllib.request.urlopen(url) as response:  # noqa: S310
        payload = response.read().decode("utf-8")
    rows = json.loads(payload)
    if not isinstance(rows, list):
        return []
    return rows


def score_bvbrc_candidate(lookup_strain: str, year: str, row: dict) -> int:
    genome_name = str(row.get("genome_name", ""))
    annotation = str(row.get("annotation", ""))
    score = score_candidate(lookup_strain, year, genome_name)
    if annotation.upper() == "REFSEQ":
        score += 10
    if annotation.upper() == "PATRIC":
        score += 5
    return score


def validate_accession(accession: str) -> tuple[bool, str, int]:
    try:
        handle = Entrez.efetch(
            db="nucleotide",
            id=accession,
            rettype="fasta_cds_aa",
            retmode="text",
        )
        try:
            from Bio import SeqIO

            records = list(SeqIO.parse(handle, "fasta"))
        finally:
            handle.close()
        rec = choose_ha_record(records)
    except Exception:  # noqa: BLE001
        return False, "", 0
    if rec is None:
        return False, "", 0
    seq = str(rec.seq).upper().replace("*", "")
    if len(seq) < 500:
        return False, rec.description, len(seq)
    return True, rec.description, len(seq)


def resolve_lookup_strain(
    lookup_strain: str,
    year: str,
    retmax: int,
    top_candidates: int,
    sleep_seconds: float,
) -> dict:
    attempted_terms = []
    seen_accessions = set()

    if not is_searchable_lookup_strain(lookup_strain):
        return {
            "lookup_strain": lookup_strain,
            "status": "skipped_malformed_lookup",
            "accession": "",
            "database": "",
            "matched_title": "",
            "matched_ha_record": "",
            "aa_length": 0,
            "matched_term_index": pd.NA,
            "candidate_rank": pd.NA,
            "note": "lookup strain not suitable for automated public search",
        }

    try:
        bvbrc_rows = retry_entrez_call(
            fetch_bvbrc_candidates,
            lookup_strain,
            sleep_seconds=sleep_seconds,
        )
        ranked_bvbrc = []
        for row in bvbrc_rows:
            accession = str(row.get("accession", "")).split(".")[0].upper()
            if not accession or accession in seen_accessions:
                continue
            ranked_bvbrc.append(
                {
                    "accession": accession,
                    "title": str(row.get("genome_name", "")),
                    "score": score_bvbrc_candidate(lookup_strain, year, row),
                    "annotation": str(row.get("annotation", "")),
                }
            )
            seen_accessions.add(accession)

        ranked_bvbrc = sorted(ranked_bvbrc, key=lambda row: row["score"], reverse=True)
        for rank_idx, candidate in enumerate(ranked_bvbrc[:top_candidates], start=1):
            ok, ha_desc, aa_len = retry_entrez_call(
                validate_accession,
                candidate["accession"],
                sleep_seconds=sleep_seconds,
            )
            time.sleep(sleep_seconds)
            if ok:
                return {
                    "lookup_strain": lookup_strain,
                    "status": "ird_match",
                    "accession": candidate["accession"],
                    "database": "IRD",
                    "matched_title": candidate["title"],
                    "matched_ha_record": ha_desc,
                    "aa_length": aa_len,
                    "matched_term_index": 0,
                    "candidate_rank": rank_idx,
                    "note": f"bvbrc_annotation={candidate['annotation']}; validated_ha_cds_aa",
                }
    except Exception as exc:  # noqa: BLE001
        attempted_terms.append(f"BV-BRC error: {type(exc).__name__}: {exc}")

    for term_idx, term in enumerate(search_terms(lookup_strain, year), start=1):
        attempted_terms.append(term)
        ids = retry_entrez_call(fetch_search_ids, term, retmax, sleep_seconds=sleep_seconds)
        time.sleep(sleep_seconds)
        if not ids:
            continue

        summaries = retry_entrez_call(fetch_summaries, ids, sleep_seconds=sleep_seconds)
        time.sleep(sleep_seconds)
        ranked = []
        for summary in summaries:
            accession = candidate_accession(summary)
            title = candidate_title(summary)
            if not accession or accession in seen_accessions:
                continue
            ranked.append(
                {
                    "accession": accession,
                    "title": title,
                    "score": score_candidate(lookup_strain, year, title),
                }
            )
            seen_accessions.add(accession)

        ranked = sorted(ranked, key=lambda row: row["score"], reverse=True)
        for rank_idx, candidate in enumerate(ranked[:top_candidates], start=1):
            ok, ha_desc, aa_len = retry_entrez_call(
                validate_accession,
                candidate["accession"],
                sleep_seconds=sleep_seconds,
            )
            time.sleep(sleep_seconds)
            if ok:
                return {
                    "lookup_strain": lookup_strain,
                    "status": "public_match",
                    "accession": candidate["accession"],
                    "database": "NCBI",
                    "matched_title": candidate["title"],
                    "matched_ha_record": ha_desc,
                    "aa_length": aa_len,
                    "matched_term_index": term_idx,
                    "candidate_rank": rank_idx,
                    "note": "validated_ha_cds_aa",
                }

    return {
        "lookup_strain": lookup_strain,
        "status": "no_public_match",
        "accession": "",
        "database": "",
        "matched_title": "",
        "matched_ha_record": "",
        "aa_length": 0,
        "matched_term_index": pd.NA,
        "candidate_rank": pd.NA,
        "note": " | ".join(attempted_terms),
    }


def load_cache(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(
            columns=[
                "lookup_strain",
                "status",
                "accession",
                "database",
                "matched_title",
                "matched_ha_record",
                "aa_length",
                "matched_term_index",
                "candidate_rank",
                "note",
            ]
        )
    return pd.read_csv(path, sep="\t", dtype={"lookup_strain": str})


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)
    manifest_path = Path(args.manifest or cfg["seq_path"])
    output_manifest_path = Path(args.output_manifest)
    cache_path = Path(args.cache)
    remaining_path = Path(args.remaining)
    gisaid_query_path = Path(args.gisaid_query)

    manifest = pd.read_csv(manifest_path, sep="\t")
    if "lookup_strain" not in manifest.columns:
        manifest["lookup_strain"] = manifest.get("base_strain", manifest["strain"]).map(canonical_strain)
    else:
        manifest["lookup_strain"] = manifest["lookup_strain"].fillna(manifest.get("base_strain", manifest["strain"])).map(canonical_strain)

    manifest["year"] = pd.to_numeric(manifest["year"], errors="coerce").astype("Int64")
    manifest["total_measurement_count"] = (
        pd.to_numeric(manifest.get("virus_measurement_count", 0), errors="coerce").fillna(0).astype(int)
        + pd.to_numeric(manifest.get("reference_measurement_count", 0), errors="coerce").fillna(0).astype(int)
    )

    unresolved = manifest[manifest["accession"].map(is_blank)].copy()
    lookup_order = (
        unresolved.groupby("lookup_strain", dropna=False)
        .agg(total_measurement_count=("total_measurement_count", "max"), year=("year", "min"))
        .reset_index()
        .sort_values(["total_measurement_count", "lookup_strain"], ascending=[False, True])
    )

    cache = load_cache(cache_path)
    cache["lookup_strain"] = cache["lookup_strain"].astype(str)
    if args.refresh_cache:
        cache = cache.iloc[0:0].copy()

    api_key = args.api_key or os.environ.get("NCBI_API_KEY") or os.environ.get("ENTREZ_API_KEY")

    Entrez.email = args.email
    if api_key:
        Entrez.api_key = api_key
    min_sleep = 0.11 if api_key else 0.34
    sleep_seconds = max(args.sleep_seconds, min_sleep)

    cached_lookup = set(cache["lookup_strain"]) if not cache.empty else set()
    to_process = lookup_order[~lookup_order["lookup_strain"].isin(cached_lookup)].copy()
    if args.max_strains is not None:
        to_process = to_process.head(args.max_strains).copy()

    if to_process.empty:
        print("No new unresolved public lookup strains remain to process.")
        print(f"Cached lookup strains already processed: {len(cached_lookup)}/{total_lookup}")

    new_rows = []
    total = len(to_process)
    for idx, row in enumerate(to_process.itertuples(index=False), start=1):
        lookup_strain = str(row.lookup_strain)
        year = safe_year(row.year)
        try:
            result = resolve_lookup_strain(
                lookup_strain=lookup_strain,
                year=year,
                retmax=args.retmax,
                top_candidates=args.top_candidates,
                sleep_seconds=sleep_seconds,
            )
        except Exception as exc:  # noqa: BLE001
            result = {
                "lookup_strain": lookup_strain,
                "status": "lookup_error",
                "accession": "",
                "database": "",
                "matched_title": "",
                "matched_ha_record": "",
                "aa_length": 0,
                "matched_term_index": pd.NA,
                "candidate_rank": pd.NA,
                "note": f"{type(exc).__name__}: {exc}",
            }
        new_rows.append(result)
        print(
            f"[{idx}/{total}] {lookup_strain} -> {result['status']}"
            + (f" ({result['accession']})" if result["accession"] else "")
        )
        if idx % 10 == 0:
            partial = pd.concat([cache, pd.DataFrame(new_rows)], ignore_index=True)
            partial = partial.drop_duplicates(subset=["lookup_strain"], keep="last")
            partial = partial.sort_values(["status", "lookup_strain"]).reset_index(drop=True)
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            partial.to_csv(cache_path, sep="\t", index=False)

    if new_rows:
        cache = pd.concat([cache, pd.DataFrame(new_rows)], ignore_index=True)
        cache = cache.drop_duplicates(subset=["lookup_strain"], keep="last")

    cache = cache.sort_values(["status", "lookup_strain"]).reset_index(drop=True)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache.to_csv(cache_path, sep="\t", index=False)

    success_map = (
        cache[cache["status"] == "public_match"]
        .set_index("lookup_strain")[["accession", "database", "matched_title", "matched_ha_record", "aa_length", "note"]]
        .to_dict("index")
    )
    status_map = cache.set_index("lookup_strain")[["status"]].to_dict("index")

    enriched_rows = []
    for row in manifest.to_dict("records"):
        lookup_strain = row["lookup_strain"]
        public = success_map.get(lookup_strain)
        status = status_map.get(lookup_strain, {}).get("status", "")

        row["public_lookup_strain"] = lookup_strain
        row["public_status"] = status
        row["public_matched_title"] = ""
        row["public_matched_ha_record"] = ""
        row["public_note"] = ""
        row["public_aa_length"] = pd.NA

        if public:
            row["public_matched_title"] = public["matched_title"]
            row["public_matched_ha_record"] = public["matched_ha_record"]
            row["public_note"] = public["note"]
            row["public_aa_length"] = public["aa_length"]
            if is_blank(row.get("accession", "")):
                row["accession"] = public["accession"]
                row["database"] = public["database"]
        enriched_rows.append(row)

    enriched = pd.DataFrame(enriched_rows)
    output_manifest_path.parent.mkdir(parents=True, exist_ok=True)
    enriched.to_csv(output_manifest_path, sep="\t", index=False)

    remaining = enriched[enriched["accession"].map(is_blank)].copy()
    remaining = remaining.sort_values(
        ["reference_measurement_count", "virus_measurement_count", "strain"],
        ascending=[False, False, True],
    )
    remaining.to_csv(remaining_path, sep="\t", index=False)

    gisaid_queries = (
        remaining["lookup_strain"]
        .dropna()
        .astype(str)
        .drop_duplicates()
        .sort_values()
    )
    gisaid_query_path.parent.mkdir(parents=True, exist_ok=True)
    gisaid_queries.to_csv(gisaid_query_path, index=False, header=False)

    total_lookup = lookup_order.shape[0]
    cache_successes = int((cache["status"] == "public_match").sum())
    total_accession_rows = int((~enriched["accession"].map(is_blank)).sum())

    print(f"Input manifest: {manifest_path}")
    print(f"Resolved lookup strains in cache: {cache_successes}/{total_lookup}")
    print(f"Manifest rows with accessions after public fill: {total_accession_rows}/{len(enriched)}")
    print(f"Updated manifest: {output_manifest_path}")
    print(f"Public resolution cache: {cache_path}")
    print(f"Remaining unresolved rows: {remaining_path}")
    print(f"Remaining GISAID query list: {gisaid_query_path}")
    if total == 0:
        print("Public lookup queue is exhausted. Use the remaining GISAID query list for the unresolved strains.")


if __name__ == "__main__":
    main()
