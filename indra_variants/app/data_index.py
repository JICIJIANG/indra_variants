"""TSV discovery and endpoint indexing.  Computed once at import time."""
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

from indra_variants.app.config import DATA_DIR
from indra_variants.app.utils import _sort_text


TSV_RE = re.compile(r"^(?P<prot>.+)_variant_effects_with_clinvar_with_domains\.tsv$", re.I)


TSV_FILES = {
    TSV_RE.match(p.name).group("prot"): p
    for p in Path(DATA_DIR).iterdir()
    if TSV_RE.match(p.name)
}
PROTS = sorted(TSV_FILES, key=_sort_text)
PROT_OPTIONS = [{'label': p, 'value': p} for p in PROTS]


# Upstream INDRA chains occasionally place gene-family / receptor-family
# symbols into the biological_process/disease column.  Exclude these so they
# do not appear in the phenotype-centric browser as if they were endpoints.
_GENE_LIKE_RE = re.compile(r"^[A-Z][A-Z0-9]{1,9}$")
_GENE_LIKE_BLACKLIST = {
    "ADRA", "ADRA2", "CAPN", "CHRN", "DRD", "EIF4EBP", "FGF", "FGFR",
    "GABR", "GAP", "P2RX", "PDE1", "PDGF", "PPP1", "S1PR", "SLC2A",
    "UBE2", "VEGFR",
}


def _looks_like_gene(name: str) -> bool:
    return name in _GENE_LIKE_BLACKLIST or bool(_GENE_LIKE_RE.match(name))


def _build_endpoint_index() -> tuple[list[str], dict[str, dict[str, dict[str, int]]]]:
    endpoint_index: dict[str, dict[str, dict[str, int]]] = defaultdict(dict)

    for prot, tsv_path in TSV_FILES.items():
        try:
            df = pd.read_csv(
                tsv_path,
                sep="\t",
                usecols=["biological_process/disease", "variant_info"]
            ).fillna('')
        except ValueError:
            df = pd.read_csv(tsv_path, sep="\t").fillna('')
            if "biological_process/disease" not in df.columns:
                continue

        work = pd.DataFrame({
            "endpoint": df["biological_process/disease"].astype(str).str.strip(),
            "variant": (
                df["variant_info"].astype(str).str.strip()
                if "variant_info" in df.columns
                else pd.Series("", index=df.index)
            ),
        })
        work = work[work["endpoint"].ne("")]
        work = work[~work["endpoint"].map(_looks_like_gene)]
        if work.empty:
            continue

        grouped = work.groupby("endpoint", sort=False).agg(
            row_count=("endpoint", "size"),
            variant_count=("variant", lambda s: s[s.ne("")].nunique()),
        )
        for endpoint_name, stats in grouped.iterrows():
            endpoint_index[endpoint_name][prot] = {
                "row_count": int(stats["row_count"]),
                "variant_count": int(stats["variant_count"]),
            }

    endpoint_names = sorted(endpoint_index, key=_sort_text)
    return endpoint_names, dict(endpoint_index)


ENDPOINTS, ENDPOINT_INDEX = _build_endpoint_index()
ENDPOINT_OPTIONS = [{'label': e, 'value': e} for e in ENDPOINTS]
