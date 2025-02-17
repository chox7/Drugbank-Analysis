import pytest
import os
from unittest.mock import patch, MagicMock
from drugbank import *

test_namespace = {"ns": "http://www.drugbank.ca"}
test_file = "test_drugbank.xml"

def create_test_xml(file_path):
    """Tworzy testowy plik XML z minimalnym zestawem danych."""
    root = etree.Element("{http://www.drugbank.ca}drugbank")

    drug = etree.SubElement(root, "{http://www.drugbank.ca}drug", type="small molecule")
    drug_id = etree.SubElement(drug, "{http://www.drugbank.ca}drugbank-id", primary="true")
    drug_id.text = "DB00001"
    name = etree.SubElement(drug, "{http://www.drugbank.ca}name")
    name.text = "Test Drug"

    synonyms = etree.SubElement(drug, "{http://www.drugbank.ca}synonyms")
    synonym = etree.SubElement(synonyms, "{http://www.drugbank.ca}synonym")
    synonym.text = "Test Synonym"

    etree.ElementTree(root).write(file_path, encoding="utf-8", xml_declaration=True, pretty_print=True)


@pytest.fixture(scope="module", autouse=True)
def setup_test_file():
    create_test_xml(test_file)
    yield
    os.remove(test_file)


@pytest.mark.parametrize("extractor, expected_columns", [
    (extract_drug_info, ["drugbank-id", "name", "type", "description", "state", "indication", "mechanism-of-action",
                         "food-interactions"]),
    (extract_synonyms, ["drugbank-id", "synonym"]),
    (extract_drug_status, ["drugbank-id", "approved", "withdrawn", "experimental", "veterinary"]),
])
def test_extractors(extractor, expected_columns):
    df = extractor(test_file, test_namespace)
    assert isinstance(df, pd.DataFrame)
    for column in expected_columns:
        assert column in df.columns
    assert not df.empty


def test_extract_targets():
    df = extract_targets(test_file, test_namespace)
    assert isinstance(df, pd.DataFrame)
    assert df.empty  # Brak docelowych białek w teście


def test_extract_interactions():
    df = extract_interactions(test_file, test_namespace)
    assert isinstance(df, pd.DataFrame)
    assert df.empty  # Brak interakcji w teście


def test_create_pathway_count_df():
    pathways_data = pd.DataFrame({
        "drugbank-id": ["DB00001", "DB00002", "DB00001"],
        "pathway-name": ["Pathway1", "Pathway2", "Pathway3"]
    })
    drugs_data = pd.DataFrame({
        "drugbank-id": ["DB00001", "DB00002", "DB00003"],
        "name": ["Drug1", "Drug2", "Drug3"]
    })
    result = create_pathway_count_df(pathways_data, drugs_data)
    assert isinstance(result, pd.DataFrame)
    assert "pathway-count" in result.columns
    assert result.loc[result["drugbank-id"] == "DB00001", "pathway-count"].values[0] == 2
    assert result.loc[result["drugbank-id"] == "DB00002", "pathway-count"].values[0] == 1
    assert result.loc[result["drugbank-id"] == "DB00003", "pathway-count"].values[0] == 0

def test_extract_drug_info_empty_file():
    empty_test_file = "empty_drugbank.xml"
    with open(empty_test_file, "w") as f:
        f.write("<drugbank></drugbank>")
    df = extract_drug_info(empty_test_file, test_namespace)
    assert df.empty
    os.remove(empty_test_file)


def test_extract_synonyms_invalid_namespace():
    df = extract_synonyms(test_file, {"ns": "http://invalid.namespace"})
    assert df.empty


def test_create_pathway_count_df_large():
    pathways_data = pd.DataFrame({
        "drugbank-id": [f"DB{i:05d}" for i in range(1000)],
        "pathway-name": [f"Pathway{i%10}" for i in range(1000)]
    })
    drugs_data = pd.DataFrame({
        "drugbank-id": [f"DB{i:05d}" for i in range(1000)],
        "name": [f"Drug{i}" for i in range(1000)]
    })
    result = create_pathway_count_df(pathways_data, drugs_data)
    assert isinstance(result, pd.DataFrame)
    assert len(result) == 1000
    assert "pathway-count" in result.columns


@patch("drugbank.requests.get")
def test_fetch_protein_data(mock_get):
    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.content = b"""
        <uniprot xmlns="http://uniprot.org/uniprot">
            <entry>
                <protein>
                    <recommendedName>
                        <fullName>Example Name</fullName>
                    </recommendedName>
                </protein>
            </entry>
        </uniprot>
    """
    mock_get.return_value = mock_response
    result = fetch_protein_data("P12345")
    assert isinstance(result, etree._Element)

    result = get_name("P12345")
    assert result == "Example Name"
