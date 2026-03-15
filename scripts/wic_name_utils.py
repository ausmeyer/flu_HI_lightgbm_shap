#!/usr/bin/env python3
"""Shared WIC strain-name normalization helpers."""

from __future__ import annotations

import re

SEARCH_TOKEN_ALIASES = {
    "HONGKONG": ["HONG KONG"],
    "NEWYORK": ["NEW YORK"],
    "SOUTHAUSTRALIA": ["SOUTH AUSTRALIA"],
    "NEWCALEDONIA": ["NEW CALEDONIA"],
    "SAUDIARABIA": ["SAUDI ARABIA"],
    "NORTHCAROLINA": ["NORTH CAROLINA"],
    "NEWHAMPSHIRE": ["NEW HAMPSHIRE"],
    "SOUTHAFRICA": ["SOUTH AFRICA"],
    "NORTHCAROLINA": ["NORTH CAROLINA"],
    "SOUTHCAROLINA": ["SOUTH CAROLINA"],
    "NORTHCAROLINA": ["NORTH CAROLINA"],
    "BADENWURTTEMBERG": ["BADEN WURTTEMBERG"],
    "BOSNIAANDHERZEGOVINA": ["BOSNIA AND HERZEGOVINA"],
    "COTEDIVOIRE": ["COTE DIVOIRE"],
    "DOMREPUBLIC": ["DOMINICAN REPUBLIC"],
    "STPETERSBURG": ["ST PETERSBURG"],
    "LAREUNION": ["LA REUNION"],
    "SRILANKA": ["SRI LANKA"],
    "LARIOJA": ["LA RIOJA"],
    "CASTILLALAMANCHA": ["CASTILLA LA MANCHA"],
    "PAISVASCO": ["PAIS VASCO"],
    "BUENOSAIRES": ["BUENOS AIRES"],
    "SANSEBASTIAN": ["SAN SEBASTIAN"],
    "SOUTHAUCKLAND": ["SOUTH AUCKLAND"],
}

SEARCH_TOKEN_CORRECTIONS = {
    "GUANDONG": ["GUANGDONG"],
    "GWANGIU": ["GWANGJU"],
}

STRAIN_PATTERN = re.compile(r"^[AB]/[^/]+/[^/]+/\d{4}$")
HEADER_STRAIN_PATTERN = re.compile(r"[AB]/[A-Z0-9 .'\-]+/[A-Z0-9._\-]+/\d{4}")


def strip_decorator(strain: str) -> str:
    return re.split(r"[@#]", str(strain).strip(), maxsplit=1)[0]


def canonical_strain(strain: str) -> str:
    base = strip_decorator(strain).strip().upper()
    parts = [p for p in base.split("/") if p]
    canon_parts = []
    for part in parts:
        if part.isdigit():
            canon_parts.append(str(int(part)))
        else:
            canon_parts.append(part)
    return "/".join(canon_parts)


def normalize_key(strain: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", canonical_strain(strain))


def is_searchable_lookup_strain(strain: str) -> bool:
    return bool(STRAIN_PATTERN.fullmatch(canonical_strain(strain)))


def _token_variants(token: str) -> list[str]:
    token = token.upper()
    variants = [token]
    for alt in SEARCH_TOKEN_CORRECTIONS.get(token, []):
        if alt not in variants:
            variants.append(alt)
    for alt in SEARCH_TOKEN_ALIASES.get(token, []):
        if alt not in variants:
            variants.append(alt)
    for corrected in SEARCH_TOKEN_CORRECTIONS.get(token, []):
        for spaced in SEARCH_TOKEN_ALIASES.get(corrected, []):
            if spaced not in variants:
                variants.append(spaced)
    return variants


def search_variants(strain: str) -> list[str]:
    canonical = canonical_strain(strain)
    parts = [p for p in canonical.split("/") if p]
    if len(parts) != 4:
        return [canonical]

    country_variants = _token_variants(parts[1])
    variants = []
    for country in country_variants:
        candidate = "/".join([parts[0], country, parts[2], parts[3]])
        if candidate not in variants:
            variants.append(candidate)
    if canonical not in variants:
        variants.insert(0, canonical)
    return variants


def preferred_query_strain(strain: str) -> str:
    variants = search_variants(strain)
    for variant in variants:
        if " " in variant:
            return variant
    return variants[0]


def extract_header_strain_candidates(header: str) -> list[str]:
    upper = str(header).upper().replace("_", "/")
    candidates = []
    for match in HEADER_STRAIN_PATTERN.findall(upper):
        cleaned = canonical_strain(match)
        if cleaned not in candidates:
            candidates.append(cleaned)
    return candidates
