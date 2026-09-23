# ===== Notebook code cell 1 =====
# Section 1: Install/import required packages
# openpyxl, requests, pandas, pdfplumber are used. All are common, well-maintained libraries.
import sys, subprocess

def _ensure(pkg):
    try:
        __import__(pkg)
    except ImportError:
        subprocess.run([sys.executable, "-m", "pip", "install", "--quiet", pkg], check=False)

for _pkg in ["requests", "pandas", "openpyxl", "pdfplumber", "biopython"]:
    _ensure(_pkg)

import re
import io
import gzip
import time
import warnings
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from collections import defaultdict

import requests
import pandas as pd
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.utils import get_column_letter
import xml.etree.ElementTree as ET
import pdfplumber
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

warnings.filterwarnings("ignore")
print("All packages imported successfully.")


# ===== Notebook code cell 2 =====
# Section 2: User input -------------------------------------------------------
protein_name = "ITK"          # <-- CHANGE THIS to any gene/protein name, e.g. "EGFR", "ABL1"

# Optional: restrict the UniProt search to one organism to avoid cross-species ambiguity.
# NCBI Taxonomy ID, e.g. 9606 = Homo sapiens. Set to None to search all organisms.
ORGANISM_ID = "9606"

# Optional manual override: if Section 3 reports an ambiguous UniProt match, paste the
# correct accession here (e.g. "P06241") and re-run the notebook.
UNIPROT_ACCESSION_OVERRIDE = None

# Contact email used in the HTTP User-Agent header, per UniProt/RCSB API etiquette.
CONTACT_EMAIL = "your_email@example.com"

# Networking behaviour
REQUEST_TIMEOUT = 30      # seconds
MAX_RETRIES = 3           # retried automatically for transient HTTP errors (429/500/502/503/504)
RETRY_BACKOFF = 1.5       # seconds, exponential backoff factor
RATE_LIMIT_DELAY = 0.15   # seconds between outgoing requests (politeness / rate-limit avoidance)

OUTPUT_FILENAME = f"Protein_PDB_Data_{protein_name}.xlsx"

print(f"Protein to search: {protein_name}")
print(f"Organism filter (NCBI taxid): {ORGANISM_ID}")
print(f"Output file will be: {OUTPUT_FILENAME}")


# ===== Notebook code cell 3 =====
# Section 3: UniProt search and validation -------------------------------------

def make_session():
    s = requests.Session()
    retries = Retry(
        total=MAX_RETRIES,
        backoff_factor=RETRY_BACKOFF,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
    )
    s.mount("https://", HTTPAdapter(max_retries=retries))
    s.headers.update({"User-Agent": f"ProteinPDBDataCollector/1.0 (mailto:{CONTACT_EMAIL})"})
    return s

SESSION = make_session()

def safe_get_json(url, params=None, timeout=REQUEST_TIMEOUT):
    """GET a URL and return (json_or_None, error_message_or_None). Never raises."""
    time.sleep(RATE_LIMIT_DELAY)
    try:
        r = SESSION.get(url, params=params, timeout=timeout)
        if r.status_code == 404:
            return None, "404 Not Found"
        r.raise_for_status()
        return r.json(), None
    except requests.exceptions.RequestException as e:
        return None, str(e)
    except ValueError as e:
        return None, f"Malformed JSON response: {e}"

def safe_get_bytes(url, timeout=REQUEST_TIMEOUT):
    """GET raw bytes from a URL. Returns (bytes_or_None, error_message_or_None)."""
    time.sleep(RATE_LIMIT_DELAY)
    try:
        r = SESSION.get(url, timeout=timeout)
        if r.status_code == 404:
            return None, "404 Not Found"
        r.raise_for_status()
        return r.content, None
    except requests.exceptions.RequestException as e:
        return None, str(e)


UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_ENTRY_URL_TMPL = "https://rest.uniprot.org/uniprotkb/{accession}.json"


def search_uniprot_candidates(gene_name, reviewed_only=True, organism_id=None):
    query_parts = [f"gene:{gene_name}"]
    if reviewed_only:
        query_parts.append("reviewed:true")
    if organism_id:
        query_parts.append(f"organism_id:{organism_id}")
    query = " AND ".join(query_parts)
    params = {
        "query": query,
        "format": "json",
        "size": 25,
        "fields": "accession,id,protein_name,organism_name,gene_names,reviewed",
    }
    data, err = safe_get_json(UNIPROT_SEARCH_URL, params=params)
    if err or not data:
        return []
    return data.get("results", [])


def gene_matches(entry, gene_name):
    for g in entry.get("genes", []) or []:
        gn = (g.get("geneName", {}) or {}).get("value", "")
        if gn.lower() == gene_name.lower():
            return True
        for syn in g.get("synonyms", []) or []:
            if syn.get("value", "").lower() == gene_name.lower():
                return True
    return False


def resolve_uniprot_entry(protein_name, organism_id=None, manual_override=None):
    if manual_override:
        print(f"Using manually specified UniProt accession override: {manual_override}")
        return manual_override

    stages = [
        dict(reviewed_only=True,  organism_id=organism_id, label="reviewed + organism + exact gene"),
        dict(reviewed_only=True,  organism_id=None,        label="reviewed + exact gene (any organism)"),
        dict(reviewed_only=False, organism_id=organism_id, label="any status + organism + exact gene"),
        dict(reviewed_only=False, organism_id=None,        label="any status + exact gene (any organism)"),
    ]
    exact, stage_label = [], None
    for stage in stages:
        candidates = search_uniprot_candidates(
            protein_name, reviewed_only=stage["reviewed_only"], organism_id=stage["organism_id"]
        )
        exact = [c for c in candidates if gene_matches(c, protein_name)]
        stage_label = stage["label"]
        if exact:
            break

    if not exact:
        raise ValueError(
            f"No UniProt entry found with an exact gene-name match for '{protein_name}'. "
            f"Verify the spelling, or set UNIPROT_ACCESSION_OVERRIDE manually in Section 2."
        )

    if len(exact) > 1:
        print("=" * 78)
        print(f"AMBIGUITY WARNING: {len(exact)} UniProt entries match gene "
              f"'{protein_name}' at resolution stage [{stage_label}]:")
        for c in exact:
            print(f"  - {c.get('primaryAccession')} | {c.get('uniProtkbId')} | "
                  f"{(c.get('organism') or {}).get('scientificName')} | {c.get('entryType')}")
        print("This pipeline will NOT silently pick one. Set UNIPROT_ACCESSION_OVERRIDE in Section 2")
        print("to the correct accession from the list above, then re-run the notebook.")
        print("=" * 78)
        raise ValueError("Ambiguous UniProt entry — manual disambiguation required (see list above).")

    chosen = exact[0]
    print(f"Selected UniProt entry (stage: {stage_label}):")
    print(f"  Accession : {chosen.get('primaryAccession')}")
    print(f"  Entry ID  : {chosen.get('uniProtkbId')}")
    print(f"  Organism  : {(chosen.get('organism') or {}).get('scientificName')}")
    print(f"  Status    : {chosen.get('entryType')}")
    return chosen.get("primaryAccession")


UNIPROT_ACCESSION = resolve_uniprot_entry(
    protein_name, organism_id=ORGANISM_ID, manual_override=UNIPROT_ACCESSION_OVERRIDE
)
PROTEIN_NAME_COLUMN_VALUE = f"{protein_name} ({UNIPROT_ACCESSION})"
print(f"\nColumn 1 value for every row will be: {PROTEIN_NAME_COLUMN_VALUE}")


# ===== Notebook code cell 4 =====
# Section 4: Retrieve associated PDB IDs from UniProt ---------------------------

def get_uniprot_full_entry(accession):
    url = UNIPROT_ENTRY_URL_TMPL.format(accession=accession)
    data, err = safe_get_json(url)
    if err:
        raise RuntimeError(f"Failed to fetch UniProt entry {accession}: {err}")
    return data


def is_xray_method(method_text):
    """Return True only for structures determined by X-ray diffraction."""
    text = (method_text or "").strip().lower()
    return "x-ray" in text or "xray" in text


def extract_pdb_ids_from_uniprot(entry_json):
    """Return unique X-ray diffraction PDB records from UniProt cross-references."""
    rows, seen, duplicates = [], set(), 0
    for xref in entry_json.get("uniProtKBCrossReferences", []) or []:
        if xref.get("database") != "PDB":
            continue
        props = {p["key"]: p["value"] for p in xref.get("properties", []) or []}
        method = props.get("Method", "N/A")
        if not is_xray_method(method):
            continue
        pdb_id = (xref.get("id") or "").upper()
        if not pdb_id:
            continue
        if pdb_id in seen:
            duplicates += 1
            continue
        seen.add(pdb_id)
        rows.append({
            "pdb_id": pdb_id,
            "uniprot_method": method,
            "uniprot_resolution": props.get("Resolution", "N/A"),
            "uniprot_chains": props.get("Chains", "N/A"),
        })
    return rows, duplicates


uniprot_entry_json = get_uniprot_full_entry(UNIPROT_ACCESSION)
pdb_id_records, duplicate_pdb_count = extract_pdb_ids_from_uniprot(uniprot_entry_json)

print(f"Total distinct PDB structures found for {protein_name} ({UNIPROT_ACCESSION}): {len(pdb_id_records)}")
if duplicate_pdb_count:
    print(f"Note: {duplicate_pdb_count} duplicate PDB ID entries in UniProt's cross-references were removed.")
if not pdb_id_records:
    print(
        "No X-ray diffraction crystal structures were found for this UniProt entry. "
        "Consider using an AlphaFold-predicted structure (or another suitable predicted "
        "structure) for the study. The rest of the notebook will produce an empty "
        "(header-only) workbook."
    )


# ===== Notebook code cell 5 =====
# Section 5: Retrieve RCSB PDB entry information --------------------------------

RCSB_ENTRY_URL_TMPL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"


def get_rcsb_entry(pdb_id):
    url = RCSB_ENTRY_URL_TMPL.format(pdb_id=pdb_id)
    return safe_get_json(url)


def format_method(rcsb_entry):
    methods = [e.get("method") for e in (rcsb_entry.get("exptl") or []) if e.get("method")]
    return "; ".join(methods) if methods else "N/A"


def format_resolution(rcsb_entry):
    """rcsb_entry_info.resolution_combined is a list of floats. RCSB always reports resolution
    to 2 decimal places, matching the 'PDB-resolution' field also seen in the validation XML
    (e.g. 1.70, 2.00), so we format to exactly 2 decimals rather than using Python's default
    float string (which would drop a trailing zero, e.g. 2.3 instead of 2.30)."""
    res_list = (rcsb_entry.get("rcsb_entry_info") or {}).get("resolution_combined")
    if not res_list:
        return "N/A"
    try:
        return ", ".join(f"{float(r):.2f}" for r in res_list)
    except (TypeError, ValueError):
        return "N/A"


print("Section 5 helper functions defined (used inside the main processing loop in Section 10).")


# R-values for the "R-Value Free" and "R-Value Work" columns are collected only from
# the depositor-reported mmCIF fields:
# _refine.ls_R_factor_R_free and _refine.ls_R_factor_R_work.
# DCC values are not used for these two columns, because DCC values may be absent for
# some PDB entries. The separate "R-free" validation-graph column is left unchanged.
MMCIF_URL_TMPL = "https://files.rcsb.org/download/{pdb_id}.cif"

def get_depositor_r_values(pdb_id):
    """Return raw deposited R-free and R-work strings from mmCIF, or None."""
    url = MMCIF_URL_TMPL.format(pdb_id=pdb_id.upper())
    raw, err = safe_get_bytes(url)
    if err or raw is None:
        return None, None, f"mmCIF unavailable ({err})"
    try:
        mmcif_text = raw.decode("utf-8", errors="replace")
        cif_path = io.StringIO(mmcif_text)
        cif = MMCIF2Dict(cif_path)
        def first_value(key):
            vals = cif.get(key)
            if isinstance(vals, list):
                vals = next((v for v in vals if str(v).strip() not in ("", ".", "?")), None)
            if vals is None or str(vals).strip() in ("", ".", "?"):
                return None
            return str(vals).strip()
        return (first_value("_refine.ls_R_factor_R_free"),
                first_value("_refine.ls_R_factor_R_work"), None)
    except Exception as e:
        return None, None, f"mmCIF parse failed ({e})"


# ===== Notebook code cell 6 =====
# Section 6: Deposited / modeled residue counts ---------------------------------

def get_residue_counts(rcsb_entry):
    info = rcsb_entry.get("rcsb_entry_info") or {}
    deposited = info.get("deposited_polymer_monomer_count")
    modeled = info.get("deposited_modeled_polymer_monomer_count")
    deposited = int(deposited) if deposited is not None else None
    modeled = int(modeled) if modeled is not None else None
    diff = (deposited - modeled) if (deposited is not None and modeled is not None) else None
    return deposited, modeled, diff


print("Section 6 helper function defined.")


# ===== Notebook code cell 7 =====
# Section 7: Exact report-display formatting
# ------------------------------------------------------------

from decimal import Decimal, InvalidOperation, ROUND_HALF_EVEN


def decimal_from_text(value):
    """
    Parse a numeric string using Decimal without introducing
    binary floating-point rounding.
    """
    if value is None:
        return None

    s = str(value).strip()

    if s == "" or s.lower() in {
        "n/a",
        "na",
        "notavailable",
        "not available",
        "."
    }:
        return None

    try:
        return Decimal(s)
    except (InvalidOperation, ValueError, TypeError):
        return None


def format_decimal(value, places):
    """
    Format a numeric value to the requested number of decimal places.

    Decimal is used instead of float so that Python does not introduce
    binary floating-point precision artifacts.

    The returned value is intentionally a STRING so that Excel does not
    reinterpret the number as a binary float.
    """
    d = decimal_from_text(value)

    if d is None:
        return "N/A"

    quantum = Decimal("1").scaleb(-places)

    rounded = d.quantize(
        quantum,
        rounding=ROUND_HALF_EVEN
    )

    return format(rounded, f".{places}f")


def format_report_integer(value):
    """
    Format report metrics such as Clashscore as an integer-looking value.
    """
    d = decimal_from_text(value)

    if d is None:
        return "N/A"

    rounded = d.quantize(
        Decimal("1"),
        rounding=ROUND_HALF_EVEN
    )

    return format(rounded, "f")


def format_report_percent(value, places=1):
    """
    Format validation percentages using the requested report precision.
    
    Example:
        0.87  -> 0.9
        4.92  -> 4.9
    """
    return format_decimal(value, places)

def format_ramachandran_outliers(value):
    """Display a reported zero as '0' to match the validation graph's 0% label.
    Nonzero values retain the existing one-decimal-place display precision.
    """
    d = decimal_from_text(value)
    if d is None:
        return "N/A"
    if d == 0:
        return "0"
    return format_decimal(value, 1)



def compute_displayed_r_difference(rfree_display, rwork_display):
    """
    Calculate R-value difference from the already displayed R-free
    and R-work values.

    Example:
        0.219 - 0.182 = 0.037
    """
    a = decimal_from_text(rfree_display)
    b = decimal_from_text(rwork_display)

    if a is None or b is None:
        return "N/A"

    difference = a - b

    return format(
        difference.quantize(
            Decimal("0.001"),
            rounding=ROUND_HALF_EVEN
        ),
        ".3f"
    )


def na_if_missing(value):
    """
    Normalize missing values from XML/API responses to a consistent marker.
    """
    if value is None:
        return "N/A"

    s = str(value).strip()

    if s == "" or s.lower() in {
        "n/a",
        "na",
        "notavailable",
        "not available",
        "."
    }:
        return "N/A"

    return s


# ------------------------------------------------------------
# 9NWX regression tests
# ------------------------------------------------------------

assert format_decimal("0.2186", 3) == "0.219"

assert format_decimal("0.1825", 3) == "0.182"

assert compute_displayed_r_difference(
    "0.219",
    "0.182"
) == "0.037"

assert format_report_integer("3.88") == "4"

assert format_report_percent(
    "0.87",
    1
) == "0.9"

assert format_report_percent(
    "4.92",
    1
) == "4.9"


print(
    "Report-display precision helpers loaded successfully; "
    "9NWX regression tests passed."
)


# ===== Notebook code cell 8 =====
# Section 8: wwPDB validation metrics (raw XML -> report precision)
# -----------------------------------------------------------------

import gzip
import xml.etree.ElementTree as ET


VALIDATION_XML_URL_TMPL = (
    "https://files.rcsb.org/pub/pdb/validation_reports/"
    "{mid}/{pid}/{pid}_validation.xml.gz"
)


def get_validation_xml_entry_attrs(pdb_id):
    """
    Download and parse the wwPDB validation XML.

    Returns the raw XML attributes as STRINGS.

    Keeping the raw values as strings is important because we apply
    the report's display precision ourselves later.
    """

    pid = str(pdb_id).strip().lower()
    mid = pid[1:3]

    url = VALIDATION_XML_URL_TMPL.format(
        mid=mid,
        pid=pid
    )

    raw, err = safe_get_bytes(url)

    if err or raw is None:
        return {}, f"validation XML unavailable ({err})"

    try:
        xml_bytes = gzip.decompress(raw)
        root = ET.fromstring(xml_bytes)

        entry_el = root.find("Entry")

        if entry_el is None:
            return {}, "validation XML had no <Entry> element"

        return dict(entry_el.attrib), None

    except (OSError, ET.ParseError) as e:
        return {}, f"failed to parse validation XML ({e})"


print("Validation XML helper defined.")


# ===== Notebook code cell 9 =====
# Section 9: Average B, all atoms (Å²) -------------------------------------------

FULL_VALIDATION_PDF_URL_TMPL = "https://files.rcsb.org/pub/pdb/validation_reports/{mid}/{pid}/{pid}_full_validation.pdf.gz"

# Matches "Average B, all atoms (Å2) 37.0" or "Average B, all atoms (Å²) 28.0", tolerating
# whatever the Å² superscript unit renders as inside the parentheses in the extracted text.
AVERAGE_B_PATTERN = re.compile(
    r"Average\s*B,\s*all\s*atoms\s*\([^)]*\)\s*(-?\d+\.\d+)", re.IGNORECASE
)
# Fallback: if the exact pattern above doesn't match (e.g. line-wrapped text extraction),
# take the last decimal number on any line that mentions both "Average B" and "atoms".
_NUMBER_PATTERN = re.compile(r"-?\d+\.\d+")


def get_average_b_all_atoms(pdb_id):
    """Return (value_as_string_or_None, error_message_or_None)."""
    pid = pdb_id.lower()
    mid = pid[1:3]
    url = FULL_VALIDATION_PDF_URL_TMPL.format(mid=mid, pid=pid)
    raw, err = safe_get_bytes(url)
    if err or raw is None:
        return None, f"full validation PDF unavailable ({err})"
    try:
        pdf_bytes = gzip.decompress(raw)
        with pdfplumber.open(io.BytesIO(pdf_bytes)) as pdf:
            for page in pdf.pages:
                text = page.extract_text() or ""
                if "Average B" not in text:
                    continue
                m = AVERAGE_B_PATTERN.search(text)
                if m:
                    return m.group(1), None
                # Fallback: scan line-by-line for the label, take the last number on that line.
                for line in text.splitlines():
                    if "Average B" in line and "atom" in line.lower():
                        nums = _NUMBER_PATTERN.findall(line)
                        if nums:
                            return nums[-1], None
        return None, "'Average B, all atoms' label not found in report text (method may not report it, e.g. some NMR/EM entries)"
    except Exception as e:
        return None, f"failed to parse full validation PDF ({e})"


print("Section 9 helper function defined.")


# ===== Notebook code cell 10 =====
# Section 10a: Ligand information + activity helpers -----------------------------

NONPOLY_ENTITY_URL_TMPL = "https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{pdb_id}/{entity_id}"

# Documented exclusion list: water, common ions, and common crystallization/cryo additives.
# This is a transparency heuristic, not a claim of biological irrelevance — the raw component
# list is always preserved separately so nothing is silently discarded.
COMMON_NON_LIGAND_COMPONENTS = {
    "HOH", "H2O", "WAT",
    "NA", "CL", "K", "MG", "CA", "ZN", "MN", "CD", "CO", "NI", "CU", "FE", "FE2",
    "BR", "IOD", "CS", "LI", "SR", "BA", "RB", "AL", "AG", "AU", "HG",
    "SO4", "PO4", "GOL", "EDO", "PEG", "PG4", "MPD", "ACT", "TRS", "IPA", "DMS",
    "BME", "MES", "HEPES", "CIT", "FMT", "ACY", "1PE", "PGE", "P6G", "EPE", "IMD",
    "NH4", "OXL", "UNX", "UNL", "DTT", "TAM", "BOG", "CO3", "NO3", "AZI",
}


def get_ligand_info(rcsb_entry, pdb_id):
    """Return (relevant_ligands, all_raw_components) — each a list of {'comp_id','name'} dicts."""
    container_ids = rcsb_entry.get("rcsb_entry_container_identifiers") or {}
    nonpoly_entity_ids = container_ids.get("non_polymer_entity_ids") or []
    raw_components = []
    for ent_id in nonpoly_entity_ids:
        url = NONPOLY_ENTITY_URL_TMPL.format(pdb_id=pdb_id, entity_id=ent_id)
        data, err = safe_get_json(url)
        if err or not data:
            continue
        comp_info = data.get("pdbx_entity_nonpoly") or {}
        comp_id = comp_info.get("comp_id", "")
        comp_name = comp_info.get("name", comp_id)
        if comp_id:
            raw_components.append({"comp_id": comp_id, "name": comp_name})
    relevant = [c for c in raw_components if c["comp_id"].upper() not in COMMON_NON_LIGAND_COMPONENTS]
    return relevant, raw_components


def get_ligand_activity(rcsb_entry, relevant_ligand_ids):
    """Format entry.rcsb_binding_affinity records for the given ligand IDs, in the style:
    '[LIG] Ki: 0.7 (nM) from 1 assay(s); IC50: min: 45, max: 55 (nM) from 2 assay(s)'.
    Returns 'N/A' if no binding-affinity data exists for these ligands."""
    affinities = rcsb_entry.get("rcsb_binding_affinity") or []
    if not affinities or not relevant_ligand_ids:
        return "N/A"
    by_ligand = defaultdict(lambda: defaultdict(list))
    for aff in affinities:
        comp_id = aff.get("comp_id", "")
        if comp_id not in relevant_ligand_ids:
            continue
        by_ligand[comp_id][aff.get("type", "?")].append((aff.get("value"), aff.get("unit", "")))
    if not by_ligand:
        return "N/A"
    ligand_parts = []
    for comp_id, types in by_ligand.items():
        type_parts = []
        for atype, vals in types.items():
            numeric_vals = [v for v, _u in vals if v is not None]
            unit = vals[0][1] if vals else ""
            n_assays = len(vals)
            if len(numeric_vals) == 1:
                type_parts.append(f"{atype}: {numeric_vals[0]} ({unit}) from {n_assays} assay(s)")
            elif len(numeric_vals) > 1:
                type_parts.append(
                    f"{atype}: min: {min(numeric_vals)}, max: {max(numeric_vals)} ({unit}) from {n_assays} assay(s)"
                )
        if type_parts:
            ligand_parts.append(f"[{comp_id}] " + "; ".join(type_parts))
    return " | ".join(ligand_parts) if ligand_parts else "N/A"


# NOTE: na_if_missing is defined once, in Section 7, and used from there. (A second,
# weaker redefinition used to live here and silently overrode the Section 7 version,
# which stopped it from recognizing placeholder text like "NotAvailable".)


print("Section 10a helper functions defined.")


# ===== Notebook code cell 11 =====
# Section 10b: Main per-PDB processing loop --------------------------------------

REQUIRED_COLUMNS = [
    "Protein Name", "PDB ID", "PDB Link", "Method", "Resolution (Å)",
    "Deposited Residue Count", "Modeled Residue Count", "Residue Count Difference",
    "R-Value Free", "R-Value Work", "R-Value Difference",
    "R-free", "Clashscore", "Ramachandran Outliers (%)", "Sidechain Outliers (%)",
    "RSRZ Outliers (%)", "Average B, all atoms (Å²)", "Ligand ID", "Ligand Name", "Ligand Activity",
]
HELPER_COLUMNS = ["All HET Components (raw)", "Notes / Errors"]

processing_log = {"success": [], "failed": [], "missing_validation": [], "missing_ligand": [], "missing_activity": []}
all_rows = []

print(f"Processing {len(pdb_id_records)} PDB structures for {protein_name} ({UNIPROT_ACCESSION})...\n")

for i, pdb_row in enumerate(pdb_id_records, 1):
    pdb_id = pdb_row["pdb_id"]
    print(f"[{i}/{len(pdb_id_records)}] {pdb_id} ...", end=" ")

    row = {c: "N/A" for c in REQUIRED_COLUMNS + HELPER_COLUMNS}
    row["Protein Name"] = PROTEIN_NAME_COLUMN_VALUE
    row["PDB ID"] = pdb_id
    row["PDB Link"] = pdb_id  # display text; hyperlink target is attached during Excel export
    row["Notes / Errors"] = ""

    try:
        rcsb_entry, err = get_rcsb_entry(pdb_id)
        if err or not rcsb_entry:
            row["Notes / Errors"] = f"RCSB entry fetch failed: {err}"
            processing_log["failed"].append((pdb_id, row["Notes / Errors"]))
            all_rows.append(row)
            print("FAILED (RCSB entry)")
            continue

        # --- Section 5: method + resolution ---
        row["Method"] = format_method(rcsb_entry)
        row["Resolution (Å)"] = format_resolution(rcsb_entry)

        # --- Section 6: residue counts ---
        deposited, modeled, diff = get_residue_counts(rcsb_entry)
        row["Deposited Residue Count"] = deposited if deposited is not None else "N/A"
        row["Modeled Residue Count"] = modeled if modeled is not None else "N/A"
        row["Residue Count Difference"] = diff if diff is not None else "N/A"

        # --- Sections 7 & 8: wwPDB validation XML -------------------------------
        # Keep XML values as strings and apply report-display rounding ourselves
        # (Section 7's format_decimal) rather than using whatever precision the
        # source happens to store -- that mismatch (e.g. a stored "0.26437" vs. the
        # "0.264" the PDB page actually displays) was the cause of the "extra digit"
        # values in the workbook.
        val_attrs, val_err = get_validation_xml_entry_attrs(pdb_id)
        dcc_rfree_raw = val_attrs.get("DCC_Rfree") if val_attrs else None
        if val_err:
            processing_log["missing_validation"].append((pdb_id, val_err))
            row["Notes / Errors"] += f"Validation XML: {val_err}. "

        # --- "R-Value Free" / "R-Value Work": depositor values only ---------------
        # DCC values are deliberately not used for these columns.
        depositor_rfree_raw, depositor_rwork_raw, cif_err = get_depositor_r_values(pdb_id)
        row["R-Value Free"] = format_decimal(depositor_rfree_raw, 3)
        row["R-Value Work"] = format_decimal(depositor_rwork_raw, 3)
        row["R-Value Difference"] = compute_displayed_r_difference(
            row["R-Value Free"], row["R-Value Work"]
        )

        # --- "R-free": the wwPDB validation-report graph value (DCC only) -------
        # This column mirrors the "wwPDB Validation" percentile-bar graph on the PDB
        # page, which plots the DCC-recalculated Rfree specifically. It NEVER falls
        # back to the depositor value: if the graph doesn't show an Rfree bar for an
        # entry, DCC_Rfree is absent/placeholder in the XML, and this column must be
        # "N/A" rather than showing an unrelated number.
        row["R-free"] = format_decimal(dcc_rfree_raw, 3)

        if cif_err:
            row["Notes / Errors"] += f"Depositor R-values: {cif_err}. "
        if row["R-Value Free"] == "N/A" or row["R-Value Work"] == "N/A":
            row["Notes / Errors"] += "No depositor R-value reported for this entry. "
        if row["R-free"] == "N/A" and not val_err:
            row["Notes / Errors"] += "R-free not shown in wwPDB validation graph for this entry. "

        if not val_err:
            row["Clashscore"] = format_report_integer(val_attrs.get("clashscore"))
            row["Ramachandran Outliers (%)"] = format_ramachandran_outliers(
                val_attrs.get("percent-rama-outliers")
            )
            row["Sidechain Outliers (%)"] = format_report_percent(
                val_attrs.get("percent-rota-outliers"), 1
            )
            row["RSRZ Outliers (%)"] = format_report_percent(
                val_attrs.get("percent-RSRZ-outliers"), 1
            )

        # --- Section 9: Average B, all atoms ---
        avg_b, avg_b_err = get_average_b_all_atoms(pdb_id)
        if avg_b is not None:
            row["Average B, all atoms (Å²)"] = avg_b
        else:
            row["Notes / Errors"] += f"Average B: {avg_b_err}. "

        # --- Section 10a: ligands + activity ---
        relevant_ligands, raw_components = get_ligand_info(rcsb_entry, pdb_id)
        if raw_components:
            row["All HET Components (raw)"] = ", ".join(c["comp_id"] for c in raw_components)
        if relevant_ligands:
            row["Ligand ID"] = ", ".join(c["comp_id"] for c in relevant_ligands)
            row["Ligand Name"] = ", ".join(c["name"] for c in relevant_ligands)
            activity = get_ligand_activity(rcsb_entry, {c["comp_id"] for c in relevant_ligands})
            row["Ligand Activity"] = activity
            if activity == "N/A":
                processing_log["missing_activity"].append(pdb_id)
        else:
            processing_log["missing_ligand"].append(pdb_id)

        processing_log["success"].append(pdb_id)
        print("done")

    except Exception as e:  # noqa: BLE001 - intentionally broad: one bad PDB must not halt the run
        row["Notes / Errors"] += f"Unexpected error: {e}. "
        processing_log["failed"].append((pdb_id, str(e)))
        print(f"FAILED ({e})")

    all_rows.append(row)

df = pd.DataFrame(all_rows, columns=REQUIRED_COLUMNS + HELPER_COLUMNS)
print(f"\nFinished. {len(df)} rows built.")


# ===== Notebook code cell 12 =====
# Section 11: Data cleaning and validation ---------------------------------------

issues = []

if df.empty:
    issues.append("The results table is empty (no PDB structures were found or processed).")
else:
    if (df["Protein Name"].isna() | (df["Protein Name"] == "")).any():
        issues.append("Some rows are missing a Protein Name.")
    if (df["PDB ID"].isna() | (df["PDB ID"] == "")).any():
        issues.append("Some rows are missing a PDB ID.")

    bad_ids = [pid for pid in df["PDB ID"] if not re.match(r"^[A-Za-z0-9]{4,8}$", str(pid))]
    if bad_ids:
        issues.append(f"{len(bad_ids)} PDB ID(s) do not look like valid PDB identifiers: {bad_ids}")

    dup_ids = df["PDB ID"][df["PDB ID"].duplicated()].tolist()
    if dup_ids:
        issues.append(f"Duplicate PDB ID rows present (kept, but flagged): {dup_ids}")

    for count_col in ["Deposited Residue Count", "Modeled Residue Count"]:
        non_numeric = [v for v in df[count_col] if v != "N/A" and not isinstance(v, (int, float))]
        if non_numeric:
            issues.append(f"Non-numeric values found in '{count_col}': {non_numeric}")

    def _diff_ok(row):
        d, m, diff = row["Deposited Residue Count"], row["Modeled Residue Count"], row["Residue Count Difference"]
        if d == "N/A" or m == "N/A" or diff == "N/A":
            return True
        return diff == (d - m)
    bad_diff_rows = df[~df.apply(_diff_ok, axis=1)]["PDB ID"].tolist()
    if bad_diff_rows:
        issues.append(f"Residue Count Difference does not match Deposited-Modeled for: {bad_diff_rows}")

print("Pre-export validation results")
print("-" * 40)
if not issues:
    print("All checks passed — no issues found.")
else:
    for i in issues:
        print(f"  ! {i}")
print("-" * 40)


# ===== Notebook code cell 13 =====
# Section 12: Create the Excel workbook -------------------------------------------

final_columns = REQUIRED_COLUMNS + HELPER_COLUMNS
df = df[final_columns]

wb = Workbook()
ws = wb.active
ws.title = "PDB_Data"

ws.append(final_columns)
header_font = Font(bold=True, color="FFFFFF")
header_fill = PatternFill(start_color="305496", end_color="305496", fill_type="solid")
for col_idx in range(1, len(final_columns) + 1):
    cell = ws.cell(row=1, column=col_idx)
    cell.font = header_font
    cell.fill = header_fill
    cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)

pdb_link_col_idx = final_columns.index("PDB Link") + 1
numeric_cols = {"Deposited Residue Count", "Modeled Residue Count", "Residue Count Difference"}

for r_idx, row_dict in enumerate(df.to_dict("records"), start=2):
    pdb_id_val = row_dict["PDB ID"]
    for c_idx, col_name in enumerate(final_columns, start=1):
        value = row_dict[col_name]
        cell = ws.cell(row=r_idx, column=c_idx, value=value)
        if col_name == "PDB Link":
            cell.value = pdb_id_val
            cell.hyperlink = f"https://www.rcsb.org/structure/{pdb_id_val}"
            cell.font = Font(color="0563C1", underline="single")
        elif col_name in numeric_cols and isinstance(value, (int, float)):
            cell.number_format = "0"
        else:
            cell.alignment = Alignment(wrap_text=True, vertical="top")

ws.freeze_panes = "A2"
ws.auto_filter.ref = ws.dimensions

for col_idx, col_name in enumerate(final_columns, start=1):
    values_in_col = [str(col_name)] + [str(v) for v in df[col_name].tolist()]
    max_len = max(len(v) for v in values_in_col)
    ws.column_dimensions[get_column_letter(col_idx)].width = min(max(max_len + 2, 10), 55)

# --- Sources sheet -----------------------------------------------------------------
ws2 = wb.create_sheet("Sources")
ws2.append(["Category", "Detail"])
for c_idx in (1, 2):
    cell = ws2.cell(row=1, column=c_idx)
    cell.font = header_font
    cell.fill = header_fill

retrieval_time_utc = datetime.now(timezone.utc).isoformat()
sources_info = [
    ("Protein searched", protein_name),
    ("UniProt accession", UNIPROT_ACCESSION),
    ("UniProt entry (API)", UNIPROT_ENTRY_URL_TMPL.format(accession=UNIPROT_ACCESSION)),
    ("UniProt entry (browser)", f"https://www.uniprot.org/uniprotkb/{UNIPROT_ACCESSION}"),
    ("RCSB entry API pattern", RCSB_ENTRY_URL_TMPL),
    ("RCSB structure page pattern", "https://www.rcsb.org/structure/{pdb_id}"),
    ("wwPDB validation XML pattern", VALIDATION_XML_URL_TMPL),
    ("wwPDB full validation PDF pattern", FULL_VALIDATION_PDF_URL_TMPL),
    ("Non-polymer / ligand entity API pattern", NONPOLY_ENTITY_URL_TMPL),
    ("Binding-affinity source field", "entry.rcsb_binding_affinity (RCSB-aggregated from BindingDB / PDBBind / Binding MOAD)"),
    ("Retrieval date/time (UTC)", retrieval_time_utc),
    ("Total PDB IDs found (UniProt)", len(pdb_id_records)),
    ("Duplicate PDB IDs removed at UniProt stage", duplicate_pdb_count),
]
for k, v in sources_info:
    ws2.append([k, v])
ws2.column_dimensions["A"].width = 42
ws2.column_dimensions["B"].width = 95
ws2.freeze_panes = "A2"

wb.save(OUTPUT_FILENAME)
print(f"Workbook saved: {OUTPUT_FILENAME}")


# ===== Notebook code cell 14 =====
# Section 12b: Validation-metric audit -------------------------------------------

DISPLAY_COLUMNS = [
    "R-Value Free", "R-Value Work", "R-Value Difference",
    "R-free", "Clashscore", "Ramachandran Outliers (%)",
    "Sidechain Outliers (%)", "RSRZ Outliers (%)",
]

def audit_display_precision(df):
    problems = []
    for _, r in df.iterrows():
        pid = r["PDB ID"]

        for col in ["R-Value Free", "R-Value Work", "R-free"]:
            v = r[col]
            if v != "N/A" and decimal_from_text(v) is None:
                problems.append((pid, col, v, "expected 3 decimal places"))

        v = r["R-Value Difference"]
        if v != "N/A" and decimal_from_text(v) is None:
            problems.append((pid, "R-Value Difference", v, "expected 3 decimal places"))

        v = r["Clashscore"]
        if v != "N/A" and not re.fullmatch(r"-?\d+", str(v)):
            problems.append((pid, "Clashscore", v, "expected integer display"))

        for col in ["Ramachandran Outliers (%)", "Sidechain Outliers (%)", "RSRZ Outliers (%)"]:
            v = r[col]
            valid_percent_display = (
                re.fullmatch(r"\d+\.\d", str(v)) is not None
                or (col == "Ramachandran Outliers (%)" and str(v) == "0")
            )
            if v != "N/A" and not valid_percent_display:
                problems.append((pid, col, v, "expected 1 decimal place, except 0% Ramachandran outliers displayed as 0"))

    if problems:
        print("DISPLAY-PRECISION AUDIT FAILED")
        for p in problems[:50]:
            print("  ", p)
        if len(problems) > 50:
            print(f"  ... and {len(problems)-50} more")
    else:
        print("DISPLAY-PRECISION AUDIT PASSED for all populated validation metrics.")

audit_display_precision(df)


# ===== Notebook code cell 15 =====
# Section 13: Generate processing summary ----------------------------------------

print("=" * 60)
print("PROCESSING SUMMARY")
print("=" * 60)
print(f"Protein searched: {protein_name}")
print(f"UniProt accession: {UNIPROT_ACCESSION}")
print(f"Total PDB IDs found: {len(pdb_id_records)}")
print(f"Successfully processed: {len(processing_log['success'])}")
print(f"Failed: {len(processing_log['failed'])}")
print(f"Missing validation data: {len(processing_log['missing_validation'])}")
print(f"Missing ligand data: {len(processing_log['missing_ligand'])}")
print(f"Missing activity data: {len(processing_log['missing_activity'])}")

if processing_log["failed"]:
    print("\nFailed PDB IDs and reasons:")
    for pid, reason in processing_log["failed"]:
        print(f"  - {pid}: {reason}")

if processing_log["missing_validation"]:
    print("\nPDB IDs with missing/partial validation data:")
    for pid, reason in processing_log["missing_validation"]:
        print(f"  - {pid}: {reason}")

print("=" * 60)
print(f"Output file: {OUTPUT_FILENAME}")
print("=" * 60)

