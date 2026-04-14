"""Shared test fixtures for HIRN Literature Retrieve."""
import pytest


SAMPLE_HIRN_HTML = """
<div class="row">
    <div class="col-sm-2">12/03/2025</div>
    <div class="col-sm-2">CBDS</div>
    <div class="col-sm-4">
        <p class="pub-title">
            <a href="https://pubmed.ncbi.nlm.nih.gov/39630627" target="_blank">
                Islet autoimmunity and HLA markers of presymptomatic and clinical type 1 diabetes.
            </a>
        </p>
        <p>
            Jacobsen LM, Bocchino LE, Evans-Molina C, et al.<br />
            Diabetologia. 2024 Dec;67(12):2611-2637.<br />
            doi: 10.1007/s00125-024-06244-0.<br />
            PMID:39630627
        </p>
    </div>
</div>
<div class="row">
    <div class="col-sm-2">11/15/2025</div>
    <div class="col-sm-2">HPAC, PanKbase</div>
    <div class="col-sm-4">
        <p class="pub-title">
            <a href="https://pubmed.ncbi.nlm.nih.gov/38012345" target="_blank">
                Single-cell transcriptomics of human pancreatic islets.
            </a>
        </p>
        <p>
            Smith J, Doe A, Johnson B.<br />
            Nature. 2024 Nov;615(7950):123-130.<br />
            doi: 10.1038/s41586-024-00001-1.<br />
            PMID:38012345
        </p>
    </div>
</div>
""".strip()


SAMPLE_AJAX_RESPONSE = {
    "html": SAMPLE_HIRN_HTML,
    "meta": {
        "postcount": 2,
        "totalposts": 2,
        "debug": False,
    },
}


SAMPLE_BIOC_JSON = {
    "source": "Auto-CuratedFull",
    "date": "20240101",
    "key": "autocuratedfull.key",
    "documents": [
        {
            "id": "39630627",
            "passages": [
                {
                    "infons": {"section_type": "TITLE", "type": "front"},
                    "offset": 0,
                    "text": "Islet autoimmunity and HLA markers of presymptomatic and clinical type 1 diabetes.",
                },
                {
                    "infons": {"section_type": "ABSTRACT", "type": "abstract"},
                    "offset": 83,
                    "text": "Type 1 diabetes is an autoimmune disease that destroys pancreatic beta cells. "
                    "Islet autoantibodies are the primary biomarkers for disease prediction. "
                    "HLA genotyping provides additional risk stratification.",
                },
                {
                    "infons": {"section_type": "INTRO", "type": "paragraph"},
                    "offset": 300,
                    "text": "The natural history of type 1 diabetes begins with genetic susceptibility, "
                    "progresses through islet autoimmunity, and culminates in clinical disease. "
                    "Understanding this progression is critical for intervention strategies.",
                },
                {
                    "infons": {"section_type": "METHODS", "type": "paragraph"},
                    "offset": 550,
                    "text": "We analyzed data from the TrialNet Pathway to Prevention study cohort. "
                    "Participants were screened for islet autoantibodies and HLA genotyped. "
                    "Follow-up extended to 15 years from initial autoantibody detection.",
                },
                {
                    "infons": {"section_type": "RESULTS", "type": "paragraph"},
                    "offset": 800,
                    "text": "Among 1,500 autoantibody-positive relatives, 45% progressed to clinical diabetes. "
                    "HLA-DR3/DR4-DQ8 carriers had the highest risk. "
                    "Multiple autoantibodies predicted faster progression.",
                },
                {
                    "infons": {"section_type": "DISCUSS", "type": "paragraph"},
                    "offset": 1050,
                    "text": "Our findings confirm that combined autoantibody and HLA profiling improves prediction "
                    "of type 1 diabetes progression. This has implications for clinical trial enrollment "
                    "and early intervention strategies.",
                },
            ],
        }
    ],
}


SAMPLE_ID_CONVERTER_RESPONSE = {
    "status": "ok",
    "response-date": "2024-01-01 00:00:00",
    "request": {"ids": ["39630627"], "format": "json"},
    "records": [
        {
            "pmid": 39630627,
            "pmcid": "PMC11615173",
            "doi": "10.1007/s00125-024-06244-0",
            "requested-id": "39630627",
        }
    ],
}


@pytest.fixture
def sample_hirn_html():
    return SAMPLE_HIRN_HTML


@pytest.fixture
def sample_ajax_response():
    return SAMPLE_AJAX_RESPONSE


@pytest.fixture
def sample_bioc_json():
    return SAMPLE_BIOC_JSON


@pytest.fixture
def sample_id_converter_response():
    return SAMPLE_ID_CONVERTER_RESPONSE


@pytest.fixture
def tmp_cache_dir(tmp_path):
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    return cache_dir
