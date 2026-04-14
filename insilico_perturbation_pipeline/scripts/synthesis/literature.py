"""Phase 6: Literature integration via PubMed E-utilities API.

Searches PubMed for relevant publications per perturbation gene and top pathways.
Uses properly parenthesized, field-tagged queries with species constraints.
"""

import time
import xml.etree.ElementTree as ET

import requests

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
RATE_LIMIT_SECONDS = 0.35  # ~3 requests/sec without API key


def search_pubmed(query, max_results=5, api_key=None):
    """Search PubMed and return article metadata.

    Args:
        query: PubMed search query string
        max_results: maximum number of results
        api_key: NCBI API key (optional, increases rate limit)

    Returns:
        list of dicts with pmid, title, authors, journal, year, snippet
    """
    params = {
        "db": "pubmed",
        "term": query,
        "retmax": max_results,
        "sort": "relevance",
        "retmode": "json",
    }
    if api_key:
        params["api_key"] = api_key

    try:
        time.sleep(RATE_LIMIT_SECONDS)
        resp = requests.get(f"{EUTILS_BASE}/esearch.fcgi", params=params, timeout=10)
        resp.raise_for_status()
        data = resp.json()
    except Exception as e:
        return [{"error": str(e)}]

    pmids = data.get("esearchresult", {}).get("idlist", [])
    if not pmids:
        return []

    return _fetch_summaries(pmids, api_key)


def _fetch_summaries(pmids, api_key=None):
    """Fetch article summaries for a list of PMIDs."""
    params = {
        "db": "pubmed",
        "id": ",".join(pmids),
        "retmode": "xml",
    }
    if api_key:
        params["api_key"] = api_key

    try:
        time.sleep(RATE_LIMIT_SECONDS)
        resp = requests.get(f"{EUTILS_BASE}/efetch.fcgi", params=params, timeout=15)
        resp.raise_for_status()
    except Exception as e:
        return [{"error": str(e)}]

    articles = []
    try:
        root = ET.fromstring(resp.text)
        for article in root.findall(".//PubmedArticle"):
            citation = article.find(".//MedlineCitation")
            art = citation.find(".//Article") if citation is not None else None
            if art is None:
                continue

            title_el = art.find("ArticleTitle")
            title = title_el.text if title_el is not None and title_el.text else "No title"

            # Authors
            author_list = art.find("AuthorList")
            authors = ""
            if author_list is not None:
                names = []
                for author in author_list.findall("Author")[:3]:
                    last = author.find("LastName")
                    if last is not None and last.text:
                        names.append(last.text)
                authors = ", ".join(names)
                if len(author_list.findall("Author")) > 3:
                    authors += " et al."

            # Year
            pub_date = art.find(".//PubDate")
            year = ""
            if pub_date is not None:
                year_el = pub_date.find("Year")
                year = year_el.text if year_el is not None else ""

            # Journal
            journal_el = art.find(".//Journal/Title")
            journal = journal_el.text if journal_el is not None and journal_el.text else ""

            # Abstract snippet
            abstract_el = art.find(".//Abstract/AbstractText")
            snippet = ""
            if abstract_el is not None and abstract_el.text:
                snippet = abstract_el.text[:300]
                if len(abstract_el.text) > 300:
                    snippet += "..."

            # PMID
            pmid_el = citation.find("PMID")
            pmid = pmid_el.text if pmid_el is not None else ""

            articles.append({
                "pmid": pmid,
                "title": title,
                "authors": authors,
                "journal": journal,
                "year": year,
                "snippet": snippet,
            })
    except ET.ParseError:
        pass

    return articles


def search_for_perturbation(gene, tissue, cell_type, top_pathways=None,
                            species="human", api_key=None):
    """Run a set of PubMed queries for one perturbation gene.

    Args:
        gene: perturbed gene symbol
        tissue: tissue context (e.g., "gut", "islet")
        cell_type: cell type (e.g., "Stem cells")
        top_pathways: list of top enriched pathway names to search
        species: "human" or "mouse"
        api_key: NCBI API key

    Returns:
        dict with query → list of article dicts
    """
    # Use simpler species filter — Mesh with quotes can cause issues in URL encoding
    species_filter = "human[tiab]" if species == "human" else "mouse[tiab]"

    queries = {
        "gene_in_context": (
            f"({gene}[tiab]) AND ({cell_type}[tiab] OR {tissue}[tiab])"
        ),
        "gene_perturbation": (
            f"({gene}[tiab]) AND (knockdown[tiab] OR knockout[tiab] OR perturbation[tiab] OR KO[tiab])"
        ),
        "gene_function": (
            f"({gene}[tiab]) AND ({species_filter})"
        ),
    }

    # Add pathway-specific queries for top 2 pathways
    if top_pathways:
        for i, pathway in enumerate(top_pathways[:2]):
            pw_clean = pathway.replace("_", " ").replace("HALLMARK ", "")
            queries[f"gene_pathway_{i}"] = (
                f"({gene}[tiab]) AND ({pw_clean}[tiab])"
            )

    results = {}
    for query_name, query in queries.items():
        articles = search_pubmed(query, max_results=5, api_key=api_key)
        results[query_name] = {
            "query": query,
            "n_results": len([a for a in articles if "error" not in a]),
            "articles": articles,
        }

    return results


def get_gene_summary(gene, api_key=None):
    """Get a brief gene summary from NCBI Gene database.

    Returns:
        dict with gene_id, description, summary, aliases
    """
    params = {
        "db": "gene",
        "term": f"{gene}[Gene Name] AND Homo sapiens[Organism]",
        "retmax": 1,
        "retmode": "json",
    }
    if api_key:
        params["api_key"] = api_key

    try:
        time.sleep(RATE_LIMIT_SECONDS)
        resp = requests.get(f"{EUTILS_BASE}/esearch.fcgi", params=params, timeout=10)
        data = resp.json()
        gene_ids = data.get("esearchresult", {}).get("idlist", [])

        if not gene_ids:
            return {"gene": gene, "description": "Not found in NCBI Gene"}

        gene_id = gene_ids[0]
        time.sleep(RATE_LIMIT_SECONDS)
        summary_resp = requests.get(
            f"{EUTILS_BASE}/esummary.fcgi",
            params={"db": "gene", "id": gene_id, "retmode": "json",
                    **({"api_key": api_key} if api_key else {})},
            timeout=10,
        )
        summary_data = summary_resp.json()
        result = summary_data.get("result", {}).get(gene_id, {})

        return {
            "gene": gene,
            "gene_id": gene_id,
            "description": result.get("description", ""),
            "summary": result.get("summary", "")[:500],
            "aliases": result.get("otheraliases", ""),
        }
    except Exception as e:
        return {"gene": gene, "error": str(e)}
