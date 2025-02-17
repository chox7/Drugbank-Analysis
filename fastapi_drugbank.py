from fastapi import FastAPI, HTTPException
from drugbank import extract_pathways_drugs_interaction, extract_drug_info, create_pathway_count_df

# Inicjalizacja aplikacji FastAPI
app = FastAPI()

xml_file = "drugbank_partial.xml"
namespace = {"ns": "http://www.drugbank.ca"}

# Wczytanie danych
pathways_interaction_df = extract_pathways_drugs_interaction(xml_file, namespace)
drugs_df = extract_drug_info(xml_file, namespace)

drug_pathway_df = create_pathway_count_df(pathways_interaction_df, drugs_df)
drug_pathway_dict = drug_pathway_df.set_index("drugbank-id")["pathway-count"].to_dict()

@app.post("/get_pathway_count/")
def get_pathway_count(drug_id: str):
    """Zwraca liczbę szlaków, z którymi dany lek wchodzi w interakcje."""
    if drug_id not in drug_pathway_dict:
        raise HTTPException(status_code=404, detail="Lek o podanym ID nie został znaleziony.")
    return {"drugbank-id": drug_id, "pathway-count": drug_pathway_dict[drug_id]}