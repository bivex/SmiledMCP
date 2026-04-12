"""Clinical evidence tools — PubMed, ClinicalTrials.gov, OpenFDA."""

import os
import json

import httpx
from xml.etree import ElementTree

from server import mcp

# ---------------------------------------------------------------------------
# Shared HTTP helper — fresh client per request survives asyncio.run() in tests
# ---------------------------------------------------------------------------
_TIMEOUT = 15.0


def _client():
    return httpx.AsyncClient(timeout=_TIMEOUT)


# ---------------------------------------------------------------------------
# PubMed (NCBI E-utilities)
# ---------------------------------------------------------------------------
_PUBMED_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_NCBI_KEY = os.getenv("NCBI_API_KEY", "")


def _pm_params(**kw):
    p = {k: v for k, v in kw.items() if v}
    if _NCBI_KEY:
        p["api_key"] = _NCBI_KEY
    return p


async def _pubmed_search(query: str, max_results: int = 20, article_types: str = "") -> list[str]:
    term = query
    if article_types:
        term += f" AND {article_types}[pt]"
    async with _client() as c:
        resp = await c.get(f"{_PUBMED_BASE}/esearch.fcgi", params=_pm_params(
            db="pubmed", term=term, retmode="json", retmax=max_results, sort="relevance",
        ))
        resp.raise_for_status()
        return resp.json().get("esearchresult", {}).get("idlist", [])


async def _pubmed_summaries(pmids: list[str]) -> list[dict]:
    if not pmids:
        return []
    async with _client() as c:
        resp = await c.get(f"{_PUBMED_BASE}/esummary.fcgi", params=_pm_params(
            db="pubmed", id=",".join(pmids), retmode="json",
        ))
        resp.raise_for_status()
        result = resp.json().get("result", {})
    articles = []
    for pmid in pmids:
        info = result.get(pmid, {})
        if not info or "error" in info:
            continue
        pmc_id = next((e["value"] for e in info.get("articleids", []) if e.get("idtype") == "pmc"), "")
        articles.append({
            "pmid": pmid,
            "title": info.get("title", ""),
            "authors": ", ".join(a.get("name", "") for a in info.get("authors", [])[:5]),
            "journal": info.get("fulljournalname", info.get("source", "")),
            "pub_date": info.get("pubdate", ""),
            "doi": next((e["value"] for e in info.get("articleids", []) if e.get("idtype") == "doi"), ""),
            "pmc_id": pmc_id,
            "pmc_url": f"https://pmc.ncbi.nlm.nih.gov/articles/{pmc_id}/" if pmc_id else "",
            "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
        })
    return articles


async def _pubmed_abstracts(pmids: list[str]) -> dict[str, str]:
    if not pmids:
        return {}
    async with _client() as c:
        resp = await c.get(f"{_PUBMED_BASE}/efetch.fcgi", params=_pm_params(
            db="pubmed", id=",".join(pmids), rettype="abstract", retmode="xml",
        ))
        resp.raise_for_status()
        text = resp.text
    abstracts = {}
    root = ElementTree.fromstring(text)
    for article in root.findall(".//PubmedArticle"):
        pmid_el = article.find(".//PMID")
        abstract_el = article.find(".//AbstractText")
        if pmid_el is not None and abstract_el is not None:
            abstracts[pmid_el.text] = abstract_el.text or ""
    return abstracts


async def _pubmed_search_with_details(query: str, max_results: int = 10, article_types: str = "") -> list[dict]:
    pmids = await _pubmed_search(query, max_results, article_types)
    if not pmids:
        return []
    summaries = await _pubmed_summaries(pmids)
    abstracts = await _pubmed_abstracts(pmids)
    for a in summaries:
        a["abstract"] = abstracts.get(a["pmid"], "")
    return summaries


# ---------------------------------------------------------------------------
# ClinicalTrials.gov API v2
# ---------------------------------------------------------------------------
_CT_BASE = "https://clinicaltrials.gov/api/v2"


async def _ct_search(condition: str = "", intervention: str = "", status: str = "", max_results: int = 10) -> list[dict]:
    params: dict = {"pageSize": min(max_results, 100)}
    if condition:
        params["query.cond"] = condition
    if intervention:
        params["query.intr"] = intervention
    if status:
        params["filter.overallStatus"] = status
    async with _client() as c:
        resp = await c.get(f"{_CT_BASE}/studies", params=params)
        resp.raise_for_status()
        data = resp.json()
    trials = []
    for study in data.get("studies", []):
        proto = study.get("protocolSection", {})
        ident = proto.get("identificationModule", {})
        status_mod = proto.get("statusModule", {})
        desc = proto.get("descriptionModule", {})
        design = proto.get("designModule", {})
        arms = proto.get("armsInterventionsModule", {})
        nct = ident.get("nctId", "")
        trials.append({
            "nct_id": nct,
            "title": ident.get("briefTitle", ""),
            "status": status_mod.get("overallStatus", ""),
            "brief_summary": desc.get("briefSummary", ""),
            "phases": design.get("phases", []),
            "interventions": [i.get("name", "") for i in arms.get("interventions", [])],
            "url": f"https://clinicaltrials.gov/study/{nct}",
        })
    return trials


async def _ct_detail(nct_id: str) -> dict | None:
    async with _client() as c:
        resp = await c.get(f"{_CT_BASE}/studies/{nct_id}")
        if resp.status_code in (404, 400):
            return None
        resp.raise_for_status()
        data = resp.json()
    proto = data.get("protocolSection", {})
    ident = proto.get("identificationModule", {})
    desc = proto.get("descriptionModule", {})
    eligibility = proto.get("eligibilityModule", {})
    contacts = proto.get("contactsLocationsModule", {})
    nct = ident.get("nctId", "")
    return {
        "nct_id": nct,
        "title": ident.get("officialTitle", ident.get("briefTitle", "")),
        "status": proto.get("statusModule", {}).get("overallStatus", ""),
        "brief_summary": desc.get("briefSummary", ""),
        "detailed_description": desc.get("detailedDescription", ""),
        "eligibility_criteria": eligibility.get("eligibilityCriteria", ""),
        "locations": [
            {"facility": l.get("facility", ""), "city": l.get("city", ""),
             "state": l.get("state", ""), "country": l.get("country", "")}
            for l in contacts.get("locations", [])[:10]
        ],
        "url": f"https://clinicaltrials.gov/study/{nct}",
    }


# ---------------------------------------------------------------------------
# OpenFDA
# ---------------------------------------------------------------------------
_FDA_BASE = "https://api.fda.gov"
_FDA_KEY = os.getenv("OPENFDA_API_KEY", "")


def _fda_params(**kw):
    p = {k: v for k, v in kw.items() if v}
    if _FDA_KEY:
        p["api_key"] = _FDA_KEY
    return p


async def _fda_drug_labels(query: str, limit: int = 5) -> list[dict]:
    async with _client() as c:
        resp = await c.get(f"{_FDA_BASE}/drug/label.json", params=_fda_params(
            search=f'openfda.brand_name:"{query}" OR openfda.generic_name:"{query}"',
            limit=min(limit, 20),
        ))
        if resp.status_code == 404:
            return []
        resp.raise_for_status()
        raw = resp.json()
    results = []
    for r in raw.get("results", []):
        fda = r.get("openfda", {})
        results.append({
            "brand_name": ", ".join(fda.get("brand_name", ["Unknown"])),
            "generic_name": ", ".join(fda.get("generic_name", ["Unknown"])),
            "manufacturer": ", ".join(fda.get("manufacturer_name", [])),
            "route": ", ".join(fda.get("route", [])),
            "indications_and_usage": (r.get("indications_and_usage", [""])[0] or "")[:500],
            "warnings": (r.get("warnings", [""])[0] or "")[:500],
            "contraindications": (r.get("contraindications", [""])[0] or "")[:500],
            "adverse_reactions": (r.get("adverse_reactions", [""])[0] or "")[:500],
        })
    return results


async def _fda_adverse_events(drug_name: str, limit: int = 10) -> dict:
    search = (
        f'patient.drug.medicinalproduct:"{drug_name}"'
        f' OR patient.drug.openfda.brand_name:"{drug_name}"'
        f' OR patient.drug.openfda.generic_name:"{drug_name}"'
    )
    async with _client() as c:
        # Get total report count
        r_total = await c.get(f"{_FDA_BASE}/drug/event.json", params=_fda_params(
            search=search, limit=1,
        ))
        if r_total.status_code == 404:
            return {"drug": drug_name, "total_reports": 0, "top_reactions": []}
        r_total.raise_for_status()
        total = r_total.json().get("meta", {}).get("results", {}).get("total", 0)

        # Get top reactions
        r_rxn = await c.get(f"{_FDA_BASE}/drug/event.json", params=_fda_params(
            search=search,
            count="patient.reaction.reactionmeddrapt.exact",
            limit=min(limit, 25),
        ))
        reactions = []
        if r_rxn.status_code == 200:
            reactions = [{"reaction": r["term"], "count": r["count"]} for r in r_rxn.json().get("results", [])]

    return {"drug": drug_name, "total_reports": total, "top_reactions": reactions}


async def _fda_drug_approvals(drug_name: str, limit: int = 5) -> list[dict]:
    async with _client() as c:
        resp = await c.get(f"{_FDA_BASE}/drug/drugsfda.json", params=_fda_params(
            search=f'openfda.brand_name:"{drug_name}" OR openfda.generic_name:"{drug_name}"',
            limit=min(limit, 10),
        ))
        if resp.status_code == 404:
            return []
        resp.raise_for_status()
        raw = resp.json()
    results = []
    for r in raw.get("results", []):
        fda = r.get("openfda", {})
        results.append({
            "brand_name": ", ".join(fda.get("brand_name", [])),
            "generic_name": ", ".join(fda.get("generic_name", [])),
            "application_number": r.get("application_number", ""),
            "sponsor_name": r.get("sponsor_name", ""),
            "submissions": [
                {"type": s.get("submission_type", ""), "number": s.get("submission_number", ""),
                 "status": s.get("submission_status", ""), "status_date": s.get("submission_status_date", "")}
                for s in r.get("submissions", [])[:5]
            ],
        })
    return results


# ===========================================================================
# MCP Tools
# ===========================================================================

@mcp.tool()
async def pubmed_search(
    query: str,
    max_results: int = 10,
    article_type: str = "",
) -> str:
    """Search PubMed for biomedical literature. Returns titles, authors, journals, and abstracts.

    Args:
        query: Search terms (e.g., "CPAP therapy obstructive sleep apnea efficacy")
        max_results: Number of results (default 10, max 20)
        article_type: Filter by type — "systematic review", "randomized controlled trial", "meta-analysis", "clinical trial", or "" for all
    """
    articles = await _pubmed_search_with_details(query, min(max_results, 20), article_type)
    if not articles:
        return json.dumps({"error": f"No PubMed articles found for: {query}"})

    results = []
    for a in articles:
        entry = {
            "title": a["title"],
            "authors": a["authors"],
            "journal": a["journal"],
            "pub_date": a["pub_date"],
            "pmid": a["pmid"],
            "doi": a["doi"] or None,
            "url": a["url"],
        }
        if a.get("pmc_url"):
            entry["pmc_url"] = a["pmc_url"]
        if a.get("abstract"):
            entry["abstract"] = a["abstract"][:600] + ("..." if len(a["abstract"]) > 600 else "")
        results.append(entry)
    return json.dumps(results, indent=2, ensure_ascii=False)


@mcp.tool()
async def clinical_trial_search(
    condition: str = "",
    intervention: str = "",
    status: str = "",
    max_results: int = 10,
) -> str:
    """Search ClinicalTrials.gov for clinical trials by condition and/or intervention.

    Args:
        condition: Disease or condition (e.g., "major depressive disorder")
        intervention: Treatment or drug (e.g., "ketamine", "cognitive behavioral therapy")
        status: Filter by status — "RECRUITING", "COMPLETED", "ACTIVE_NOT_RECRUITING", or "" for all
        max_results: Number of results (default 10, max 20)
    """
    trials = await _ct_search(condition, intervention, status, min(max_results, 20))
    if not trials:
        return json.dumps({"error": f"No clinical trials found for condition='{condition}', intervention='{intervention}'"})

    results = []
    for t in trials:
        entry = {
            "nct_id": t["nct_id"],
            "title": t["title"],
            "status": t["status"],
            "phases": t["phases"],
            "interventions": t["interventions"][:3],
            "url": t["url"],
        }
        if t.get("brief_summary"):
            entry["brief_summary"] = t["brief_summary"][:400] + ("..." if len(t["brief_summary"]) > 400 else "")
        results.append(entry)
    return json.dumps(results, indent=2, ensure_ascii=False)


@mcp.tool()
async def clinical_trial_detail(nct_id: str) -> str:
    """Get full details for a specific clinical trial by NCT ID.

    Args:
        nct_id: The ClinicalTrials.gov identifier (e.g., "NCT03872505")
    """
    trial = await _ct_detail(nct_id)
    if not trial:
        return json.dumps({"error": f"Trial {nct_id} not found."})
    return json.dumps(trial, indent=2, ensure_ascii=False)


@mcp.tool()
async def drug_label(drug_name: str) -> str:
    """Look up FDA-approved drug labeling — indications, warnings, contraindications, adverse reactions.

    Args:
        drug_name: Brand or generic drug name (e.g., "Ozempic", "metformin", "duloxetine")
    """
    labels = await _fda_drug_labels(drug_name, limit=3)
    if not labels:
        return json.dumps({"error": f"No FDA drug labels found for: {drug_name}"})
    return json.dumps(labels, indent=2, ensure_ascii=False)


@mcp.tool()
async def drug_adverse_events(drug_name: str, limit: int = 15) -> str:
    """Search FDA adverse event reports (FAERS) for a drug — shows most commonly reported reactions.

    Args:
        drug_name: Brand or generic drug name
        limit: Number of top reactions to return (default 15)
    """
    result = await _fda_adverse_events(drug_name, min(limit, 25))
    if result["total_reports"] == 0:
        return json.dumps({"error": f"No adverse event reports found for: {drug_name}"})
    return json.dumps(result, indent=2, ensure_ascii=False)


@mcp.tool()
async def drug_approvals(drug_name: str) -> str:
    """Search FDA drug approval history — application numbers, sponsors, submission timeline.

    Args:
        drug_name: Brand or generic drug name
    """
    approvals = await _fda_drug_approvals(drug_name)
    if not approvals:
        return json.dumps({"error": f"No FDA approval records found for: {drug_name}"})
    return json.dumps(approvals, indent=2, ensure_ascii=False)


@mcp.tool()
async def medical_evidence_search(
    condition: str,
    treatment: str = "",
    author: str = "",
    evidence_type: str = "all",
) -> str:
    """Comprehensive evidence search across PubMed, ClinicalTrials.gov, and FDA.

    Args:
        condition: The medical condition (e.g., "treatment-resistant depression")
        treatment: The specific treatment in question (e.g., "transcranial magnetic stimulation")
        author: Optional author name to filter by
        evidence_type: "all", "pubmed", "trials", or "fda"
    """
    query = f"{condition} {treatment}".strip()
    if author:
        query += f" {author}[Author]"
    sections = {}

    if evidence_type in ("all", "pubmed"):
        articles = await _pubmed_search_with_details(query, max_results=5, article_types="systematic review")
        if not articles:
            articles = await _pubmed_search_with_details(query, max_results=5, article_types="")
        if articles:
            sections["pubmed"] = [
                {"title": a["title"], "journal": a["journal"], "pub_date": a["pub_date"],
                 "pmid": a["pmid"], "url": a["url"],
                 "pmc_url": a.get("pmc_url") or None,
                 "abstract": (a.get("abstract") or "")[:600] + ("..." if len(a.get("abstract") or "") > 600 else "")}
                for a in articles
            ]

    if evidence_type in ("all", "trials") and treatment:
        trials = await _ct_search(condition, treatment, max_results=5)
        if trials:
            sections["trials"] = [
                {"nct_id": t["nct_id"], "title": t["title"], "status": t["status"],
                 "phases": t["phases"], "url": t["url"]}
                for t in trials
            ]

    if evidence_type in ("all", "fda") and treatment:
        labels = await _fda_drug_labels(treatment, limit=2)
        if labels:
            sections["fda"] = [
                {"brand_name": l["brand_name"], "generic_name": l["generic_name"],
                 "manufacturer": l["manufacturer"],
                 "indications": l.get("indications_and_usage", "")[:200]}
                for l in labels
            ]

    if not sections:
        return json.dumps({"error": f"No evidence found for condition='{condition}', treatment='{treatment}'"})

    return json.dumps({"condition": condition, "treatment": treatment or None, "evidence": sections}, indent=2, ensure_ascii=False)
