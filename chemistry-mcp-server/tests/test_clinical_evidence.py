"""Tests for clinical evidence tools — PubMed, ClinicalTrials.gov, OpenFDA.

These tests make real API calls and require network access.
"""

import json
import sys
import os
import asyncio

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tools.clinical_evidence import (
    pubmed_search,
    clinical_trial_search,
    clinical_trial_detail,
    drug_label,
    drug_adverse_events,
    drug_approvals,
    medical_evidence_search,
)


def _json(result: str) -> dict | list:
    return json.loads(result)


def _run(coro):
    return asyncio.run(coro)


# ════════════════════════════════════════════════════════════
#  PubMed
# ════════════════════════════════════════════════════════════

class TestPubMedSearch:
    def test_basic_search(self):
        result = _json(_run(pubmed_search("aspirin headache", max_results=3)))
        assert isinstance(result, list)
        assert len(result) >= 1

    def test_result_has_expected_keys(self):
        result = _json(_run(pubmed_search("aspirin", max_results=1)))
        article = result[0]
        for key in ("title", "authors", "journal", "pub_date", "pmid", "url"):
            assert key in article, f"Missing key: {key}"

    def test_pmid_is_string(self):
        result = _json(_run(pubmed_search("caffeine", max_results=1)))
        assert isinstance(result[0]["pmid"], str)

    def test_url_is_valid_pubmed_link(self):
        result = _json(_run(pubmed_search("caffeine", max_results=1)))
        url = result[0]["url"]
        assert "pubmed.ncbi.nlm.nih.gov" in url

    def test_max_results_respected(self):
        result = _json(_run(pubmed_search("aspirin", max_results=3)))
        assert len(result) <= 3

    def test_max_results_capped_at_20(self):
        result = _json(_run(pubmed_search("cancer", max_results=100)))
        assert len(result) <= 20

    def test_article_type_filter(self):
        result = _json(_run(pubmed_search(
            "CPAP sleep apnea", max_results=3,
            article_type="randomized controlled trial",
        )))
        assert isinstance(result, list)

    def test_abstract_truncated_at_600(self):
        result = _json(_run(pubmed_search("depression treatment", max_results=5)))
        for article in result:
            if "abstract" in article and article["abstract"]:
                assert len(article["abstract"]) <= 603  # 600 + "..."

    def test_no_results_returns_error(self):
        result = _json(_run(pubmed_search("xyzzyqwert12345nonexistent", max_results=1)))
        assert isinstance(result, dict)
        assert "error" in result

    def test_doi_field_present(self):
        result = _json(_run(pubmed_search("aspirin", max_results=3)))
        for article in result:
            assert "doi" in article

    def test_pmc_url_field(self):
        result = _json(_run(pubmed_search("aspirin", max_results=5)))
        for article in result:
            if "pmc_url" in article and article["pmc_url"]:
                assert "pmc.ncbi.nlm.nih.gov" in article["pmc_url"]


# ════════════════════════════════════════════════════════════
#  ClinicalTrials.gov search
# ════════════════════════════════════════════════════════════

class TestClinicalTrialSearch:
    def test_search_by_condition(self):
        result = _json(_run(clinical_trial_search(condition="diabetes", max_results=3)))
        assert isinstance(result, list)
        assert len(result) >= 1

    def test_search_by_intervention(self):
        result = _json(_run(clinical_trial_search(intervention="metformin", max_results=3)))
        assert isinstance(result, list)
        assert len(result) >= 1

    def test_search_condition_and_intervention(self):
        result = _json(_run(clinical_trial_search(
            condition="depression", intervention="ketamine", max_results=3,
        )))
        assert isinstance(result, list)

    def test_result_has_expected_keys(self):
        result = _json(_run(clinical_trial_search(condition="asthma", max_results=1)))
        trial = result[0]
        for key in ("nct_id", "title", "status", "phases", "interventions", "url"):
            assert key in trial, f"Missing key: {key}"

    def test_nct_id_format(self):
        result = _json(_run(clinical_trial_search(condition="asthma", max_results=1)))
        nct = result[0]["nct_id"]
        assert nct.startswith("NCT")

    def test_status_filter(self):
        result = _json(_run(clinical_trial_search(
            condition="diabetes", status="RECRUITING", max_results=3,
        )))
        assert isinstance(result, list)
        for trial in result:
            assert trial["status"] == "RECRUITING"

    def test_max_results_capped_at_20(self):
        result = _json(_run(clinical_trial_search(condition="cancer", max_results=100)))
        assert len(result) <= 20

    def test_no_results_returns_error(self):
        result = _json(_run(clinical_trial_search(
            condition="xyzzyqwertnonexistent12345", max_results=1,
        )))
        assert isinstance(result, dict)
        assert "error" in result

    def test_interventions_list_truncated(self):
        result = _json(_run(clinical_trial_search(condition="cancer", max_results=3)))
        for trial in result:
            assert len(trial["interventions"]) <= 3

    def test_brief_summary_truncated(self):
        result = _json(_run(clinical_trial_search(condition="diabetes", max_results=3)))
        for trial in result:
            if "brief_summary" in trial and trial["brief_summary"]:
                assert len(trial["brief_summary"]) <= 403

    def test_url_links_to_clinicaltrials(self):
        result = _json(_run(clinical_trial_search(condition="asthma", max_results=1)))
        assert "clinicaltrials.gov" in result[0]["url"]


# ════════════════════════════════════════════════════════════
#  Clinical trial detail
# ════════════════════════════════════════════════════════════

class TestClinicalTrialDetail:
    def test_valid_nct_id(self):
        result = _json(_run(clinical_trial_detail("NCT00217802")))
        assert isinstance(result, dict)
        assert "nct_id" in result
        assert "title" in result
        assert "status" in result

    def test_brief_summary_present(self):
        result = _json(_run(clinical_trial_detail("NCT00217802")))
        assert "brief_summary" in result

    def test_locations_is_list(self):
        result = _json(_run(clinical_trial_detail("NCT00217802")))
        assert isinstance(result.get("locations", []), list)

    def test_not_found_returns_error(self):
        result = _json(_run(clinical_trial_detail("NCT00000000")))
        assert isinstance(result, dict)
        assert "error" in result

    def test_url_field_present(self):
        result = _json(_run(clinical_trial_detail("NCT00217802")))
        assert "url" in result
        assert "clinicaltrials.gov" in result["url"]


# ════════════════════════════════════════════════════════════
#  FDA Drug Label
# ════════════════════════════════════════════════════════════

class TestDrugLabel:
    def test_search_by_brand_name(self):
        result = _json(_run(drug_label("Ozempic")))
        assert isinstance(result, list)
        assert len(result) >= 1

    def test_search_by_generic_name(self):
        result = _json(_run(drug_label("metformin")))
        assert isinstance(result, list)
        assert len(result) >= 1

    def test_result_has_expected_keys(self):
        result = _json(_run(drug_label("aspirin")))
        label = result[0]
        for key in ("brand_name", "generic_name", "manufacturer", "route"):
            assert key in label, f"Missing key: {key}"

    def test_indications_field(self):
        result = _json(_run(drug_label("Ozempic")))
        assert any(l.get("indications_and_usage") for l in result)

    def test_no_results_returns_error(self):
        result = _json(_run(drug_label("XYZNOTAREALDRUG12345")))
        assert isinstance(result, dict)
        assert "error" in result

    def test_text_fields_truncated_at_500(self):
        result = _json(_run(drug_label("Ozempic")))
        for label in result:
            for field in ("indications_and_usage", "warnings", "contraindications", "adverse_reactions"):
                if label.get(field):
                    assert len(label[field]) <= 500


# ════════════════════════════════════════════════════════════
#  FDA Adverse Events
# ════════════════════════════════════════════════════════════

class TestDrugAdverseEvents:
    def test_search_known_drug(self):
        result = _json(_run(drug_adverse_events("aspirin")))
        assert isinstance(result, dict)
        assert "total_reports" in result
        assert "top_reactions" in result

    def test_total_reports_positive(self):
        result = _json(_run(drug_adverse_events("aspirin")))
        assert result["total_reports"] > 0

    def test_reactions_have_reaction_and_count(self):
        result = _json(_run(drug_adverse_events("aspirin", limit=5)))
        for r in result["top_reactions"]:
            assert "reaction" in r
            assert "count" in r
            assert isinstance(r["count"], int)

    def test_limit_respected(self):
        result = _json(_run(drug_adverse_events("aspirin", limit=5)))
        assert len(result["top_reactions"]) <= 5

    def test_limit_capped_at_25(self):
        result = _json(_run(drug_adverse_events("aspirin", limit=100)))
        assert len(result["top_reactions"]) <= 25

    def test_no_results_returns_error(self):
        result = _json(_run(drug_adverse_events("XYZNOTAREALDRUG12345")))
        assert isinstance(result, dict)
        assert "error" in result


# ════════════════════════════════════════════════════════════
#  FDA Drug Approvals
# ════════════════════════════════════════════════════════════

class TestDrugApprovals:
    def test_search_known_drug(self):
        result = _json(_run(drug_approvals("metformin")))
        assert isinstance(result, list)
        assert len(result) >= 1

    def test_result_has_expected_keys(self):
        result = _json(_run(drug_approvals("metformin")))
        approval = result[0]
        for key in ("brand_name", "generic_name", "application_number", "sponsor_name"):
            assert key in approval, f"Missing key: {key}"

    def test_submissions_structure(self):
        result = _json(_run(drug_approvals("metformin")))
        approval = result[0]
        if approval.get("submissions"):
            sub = approval["submissions"][0]
            for key in ("type", "number", "status", "status_date"):
                assert key in sub, f"Missing submission key: {key}"

    def test_no_results_returns_error(self):
        result = _json(_run(drug_approvals("XYZNOTAREALDRUG12345")))
        assert isinstance(result, dict)
        assert "error" in result


# ════════════════════════════════════════════════════════════
#  Comprehensive Medical Evidence Search
# ════════════════════════════════════════════════════════════

class TestMedicalEvidenceSearch:
    def test_all_sources(self):
        result = _json(_run(medical_evidence_search(
            condition="depression", treatment="ketamine",
        )))
        assert isinstance(result, dict)
        assert "evidence" in result
        assert "condition" in result

    def test_pubmed_only(self):
        result = _json(_run(medical_evidence_search(
            condition="headache", evidence_type="pubmed",
        )))
        assert isinstance(result, dict)
        evidence = result.get("evidence", {})
        assert "pubmed" in evidence
        assert "trials" not in evidence
        assert "fda" not in evidence

    def test_trials_only(self):
        result = _json(_run(medical_evidence_search(
            condition="diabetes", treatment="metformin", evidence_type="trials",
        )))
        assert isinstance(result, dict)
        evidence = result.get("evidence", {})
        assert "trials" in evidence
        assert "pubmed" not in evidence

    def test_fda_only(self):
        result = _json(_run(medical_evidence_search(
            condition="diabetes", treatment="metformin", evidence_type="fda",
        )))
        assert isinstance(result, dict)
        evidence = result.get("evidence", {})
        assert "fda" in evidence

    def test_author_filter(self):
        result = _json(_run(medical_evidence_search(
            condition="sleep apnea", author="Smith",
        )))
        assert isinstance(result, dict)

    def test_condition_without_treatment(self):
        result = _json(_run(medical_evidence_search(condition="diabetes")))
        assert isinstance(result, dict)
        # Without treatment, trials and fda sections should not be present
        evidence = result.get("evidence", {})
        assert "pubmed" in evidence

    def test_no_results_returns_error(self):
        result = _json(_run(medical_evidence_search(
            condition="xyzzyqwertnonexistent12345", treatment="abcnotreal67890",
        )))
        assert isinstance(result, dict)
        assert "error" in result

    def test_treatment_reflected_in_output(self):
        result = _json(_run(medical_evidence_search(
            condition="depression", treatment="ketamine",
        )))
        assert result.get("treatment") == "ketamine"

    def test_condition_reflected_in_output(self):
        result = _json(_run(medical_evidence_search(
            condition="diabetes", treatment="metformin",
        )))
        assert result.get("condition") == "diabetes"


# ════════════════════════════════════════════════════════════
#  JSON output consistency
# ════════════════════════════════════════════════════════════

class TestClinicalEvidenceJsonOutput:
    """All tools must return valid JSON strings."""

    def test_pubmed_search_valid_json(self):
        result = _run(pubmed_search("aspirin", max_results=1))
        parsed = json.loads(result)
        assert isinstance(parsed, (dict, list))

    def test_clinical_trial_search_valid_json(self):
        result = _run(clinical_trial_search(condition="asthma", max_results=1))
        parsed = json.loads(result)
        assert isinstance(parsed, (dict, list))

    def test_drug_label_valid_json(self):
        result = _run(drug_label("aspirin"))
        parsed = json.loads(result)
        assert isinstance(parsed, (dict, list))

    def test_drug_adverse_events_valid_json(self):
        result = _run(drug_adverse_events("aspirin", limit=3))
        parsed = json.loads(result)
        assert isinstance(parsed, (dict, list))

    def test_drug_approvals_valid_json(self):
        result = _run(drug_approvals("aspirin"))
        parsed = json.loads(result)
        assert isinstance(parsed, (dict, list))

    def test_medical_evidence_search_valid_json(self):
        result = _run(medical_evidence_search(condition="headache", treatment="aspirin"))
        parsed = json.loads(result)
        assert isinstance(parsed, (dict, list))
