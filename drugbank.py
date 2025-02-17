from lxml import etree
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import requests
import networkx as nx
from typing import Optional, Callable, Iterator, Union, Generator

def iter_rows(
    xml_file: str,
    namespace: dict,
    row_extractor: Callable[[etree._Element, dict], list[dict]]
) -> Iterator[dict]:
    """
    Generator iterujący po elementach <drug> w XML i zwracający zwracający słowniki z danymi.

    :param xml_file: Ścieżka do pliku XML.
    :param namespace: Słownik namespace'ów do wyszukiwania elementów.
    :param row_extractor: Funkcja przetwarzająca element XML na listę słowników.
    :return: Iterator zwracający kolejne rekordy w postaci słowników.
    """
    tag =f"{{{namespace['ns']}}}drug"
    parent_tag = f"{{{namespace['ns']}}}drugbank"
    context = etree.iterparse(xml_file, events=("end",), tag=tag)

    for event, elem in context:
        if elem.tag == tag and elem.getparent().tag == parent_tag:
            yield from row_extractor(elem, namespace)

            # Czyszczenie pamięci po przetworzonym elemencie
            elem.clear()
            while elem.getprevious() is not None:
                del elem.getparent()[0]


def extract_drug_info_row(drug: etree._Element, namespace: dict) -> list[dict]:
    """
    Pobiera informacje o leku z pojedynczego elementu <drug> w XML.

    :param drug: Element XML reprezentujący lek.
    :param namespace: Słownik namespace'ów do wyszukiwania elementów.
    :return: Lista słowników zawierających informacje o leku.
    """
    id_element = drug.find("ns:drugbank-id[@primary='true']", namespace)
    id_txt = id_element.text if id_element is not None else None

    name_element = drug.find("ns:name", namespace)
    name_txt = name_element.text if name_element is not None else None

    type_txt = str(drug.get("type"))

    description_element = drug.find("ns:description", namespace)
    description_txt = description_element.text if description_element is not None else None

    state_element = drug.find("ns:state", namespace)
    state_txt = state_element.text if state_element is not None else None

    indication_element = drug.find("ns:indication", namespace)
    indication_txt = indication_element.text if indication_element is not None else None

    mechanism_of_action_element = drug.find("ns:mechanism-of-action", namespace)
    mechanism_of_action_txt = mechanism_of_action_element.text if mechanism_of_action_element is not None else None

    food_interactions_elements = drug.findall("ns:food-interactions/ns:food-interaction", namespace)
    food_interactions_list = [fi.text for fi in food_interactions_elements if fi.text]

    return [{
        'drugbank-id': id_txt,
        'name': name_txt,
        'type': type_txt,
        'description': description_txt,
        'state': state_txt,
        'indication': indication_txt,
        'mechanism-of-action': mechanism_of_action_txt,
        'food-interactions': '\n\n'.join(food_interactions_list) if food_interactions_list else None
    }]


def extract_drug_info(xml_file: str, namespace: dict) -> pd.DataFrame:
    """
    Pobiera informacje o lekach z XML i zwraca je w postaci DataFrame.

    :param xml_file: Ścieżka do pliku XML.
    :param namespace: Słownik namespace'ów do wyszukiwania elementów.
    :return: DataFrame z informacjami o lekach.
    """
    rows_iter = iter_rows(xml_file, namespace, extract_drug_info_row)
    rows = [row for row in rows_iter if row]
    return pd.DataFrame(rows)


def extract_synonyms_row(drug: etree._Element, namespace: dict) -> list[dict]:
    """
    Pobiera synonimy dla pojedynczego elementu <drug>.

    :param drug: Element XML reprezentujący lek.
    :param namespace: Słownik namespace'ów do wyszukiwania elementów.
    :return: Lista słowników zawierających ID leku i jego synonimy.
    """
    id_element = drug.find("ns:drugbank-id[@primary='true']", namespace)
    id_txt = id_element.text if id_element is not None else None

    synonyms_elements = drug.findall("ns:synonyms/ns:synonym", namespace)
    return [
        {"drugbank-id": id_txt, "synonym": synonym.text}
        for synonym in synonyms_elements if synonym.text
    ]


def extract_synonyms(xml_file: str, namespace: dict) -> pd.DataFrame:
    """
    Pobiera synonimy leków z XML i zwraca je w postaci DataFrame.

    :param xml_file: Ścieżka do pliku XML.
    :param namespace: Słownik namespace'ów do wyszukiwania elementów.
    :return: DataFrame z synonimami dla leków.
    """
    rows_iter = iter_rows(xml_file, namespace, extract_synonyms_row)
    rows = [row for row in rows_iter if row]
    return pd.DataFrame(rows)


def create_synonym_graph(drugbank_id: str, drugs_df: pd.DataFrame, synonym_df: pd.DataFrame) -> None:
    """
    Tworzy i wyświetla graf synonimów dla danego leku na podstawie DrugBank ID.

    :param drugbank_id: Unikalny identyfikator leku (DrugBank ID).
    :param drugs_df: DataFrame z informacjami o lekach.
    :param synonym_df: DataFrame z synonimami leków.
    """
    if drugbank_id not in drugs_df['drugbank-id'].values:
        print(f"Lek o DrugBank ID {drugbank_id} nie został znaleziony.")
        return

    drug_name = drugs_df.loc[drugs_df['drugbank-id'] == drugbank_id, 'name'].values[0]
    drug_synonyms = synonym_df.loc[synonym_df['drugbank-id'] == drugbank_id, 'synonym']
    if drug_synonyms.empty:
        print(f"Nie znaleziono synonimów dla leku o DrugBank ID: {drugbank_id}")
        return

    G = nx.Graph()
    G.add_node(drug_name)

    for synonym in drug_synonyms:
        if synonym and synonym != drug_name:
            G.add_node(synonym)
            G.add_edge(drug_name, synonym)

    # Rysowanie grafu
    plt.figure(figsize=(10, 8))
    pos = nx.spring_layout(G, seed=42, k=0.1)

    nx.draw(
        G, pos, with_labels=False,
        node_size=5000, node_color="skyblue",
        font_size=10, font_weight="bold", width=2
    )

    for node, (x, y) in pos.items():
        plt.text(
            x, y, node, ha='center', va='center',
            bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="white")
        )

    plt.title(f"Graf synonimów dla leku o DrugBank ID: {drugbank_id}")
    plt.axis('off')
    plt.show()


def extract_drug_products_row(drug: etree._Element, namespace: dict) -> list[dict]:
    """
    Pobiera informacje o produktach farmaceutycznych powiązanych z lekiem.

    :param drug: Element XML reprezentujący lek.
    :param namespace: Słownik namespace'ów do wyszukiwania elementów.
    :return: Lista słowników z informacjami o produktach farmaceutycznych.
    """
    products = []
    id_element = drug.find("ns:drugbank-id[@primary='true']", namespace)
    drugbank_id = id_element.text if id_element is not None else None

    drug_name_element = drug.find("ns:name", namespace)
    drug_name = drug_name_element.text if drug_name_element is not None else None

    products_element = drug.find("ns:products", namespace)
    if products_element is None:
        return []

    ndc_code = None
    for product in products_element.findall("ns:product", namespace):
        source_element = product.find("ns:source", namespace)
        source = source_element.text if source_element is not None else None

        if source == "FDA NDC":
            ndc_code_element = product.find("ns:ndc-product-code", namespace)
            ndc_code = ndc_code_element.text if ndc_code_element is not None else None
            break

    for product in products_element.findall("ns:product", namespace):
        name_element = product.find("ns:name", namespace)
        name = name_element.text if name_element is not None else None

        producer_element = product.find("ns:labeller", namespace)
        producer = producer_element.text if producer_element is not None else None

        dosage_form_element = product.find("ns:dosage-form", namespace)
        dosage_form = dosage_form_element.text if dosage_form_element is not None else None

        route_element = product.find("ns:route", namespace)
        route = route_element.text if route_element is not None else None

        strength_element = product.find("ns:strength", namespace)
        strength = strength_element.text if strength_element is not None else None

        country_element = product.find("ns:country", namespace)
        country = country_element.text if country_element is not None else None

        agency_element = product.find("ns:source", namespace)
        agency = agency_element.text if agency_element is not None else None

        products.append({
            "drugbank-id": drugbank_id,
            "drug-name": drug_name,
            "product-name": name,
            "producer": producer,
            "ndc-code": ndc_code,
            "form": dosage_form,
            "route": route,
            "strength": strength,
            "country": country,
            "agency": agency
        })

    return products


def extract_drug_products(xml_file: str, namespace: dict) -> pd.DataFrame:
    """
    Pobiera informacje o produktach farmaceutycznych z XML i zwraca je w formie DataFrame.

    :param xml_file: Ścieżka do pliku XML.
    :param namespace: Słownik namespace'ów do wyszukiwania elementów.
    :return: DataFrame zawierający dane o produktach farmaceutycznych.
    """
    rows_iter = iter_rows(xml_file, namespace, extract_drug_products_row)
    rows = [row for row in rows_iter if row]
    return pd.DataFrame(rows)


def extract_pathways_row(drug: etree._Element, namespace: dict) -> list[dict]:
    """
    Pobiera informacje o szlakach (pathways) z pojedynczego elementu <drug>.

    :param drug: Element XML reprezentujący lek.
    :param namespace: Słownik namespace'ów do wyszukiwania elementów.
    :return: Lista słowników z informacjami o szlakach.
    """
    pathways_element = drug.find("ns:pathways", namespace)
    if pathways_element is None:
        return []

    pathways = []
    for pathway in pathways_element.findall("ns:pathway", namespace):
        smpdb_id_element = pathway.find("ns:smpdb-id", namespace)
        smpdb_id = smpdb_id_element.text if smpdb_id_element is not None else None

        pathway_name_element = pathway.find("ns:name", namespace)
        pathway_name = pathway_name_element.text if pathway_name_element is not None else None

        category_element = pathway.find("ns:category", namespace)
        category = category_element.text if category_element is not None else None

        pathways.append({
            "smpdb-id": smpdb_id,
            "pathway-name": pathway_name,
            "category": category,
        })

    return pathways


def extract_pathways(xml_file: str, namespace: dict) -> pd.DataFrame:
    """
    Pobiera szlaki sygnałowe i metaboliczne z XML i zwraca je w formie DataFrame.

    :param xml_file: Ścieżka do pliku XML.
    :param namespace: Słownik namespace'ów do wyszukiwania elementów.
    :return: DataFrame zawierający informacje o szlakach.
    """
    rows_iter = iter_rows(xml_file, namespace, extract_pathways_row)
    rows = [row for row in rows_iter if row]
    return pd.DataFrame(rows)


def extract_pathways_drugs_interaction_row(drug: etree._Element, namespace: dict) -> list[dict]:
    """
    Pobiera informację o interakcjach leków ze szlakami metabolicznymi.

    :param drug: Element XML reprezentujący lek.
    :param namespace: Słownik namespace'ów do wyszukiwania elementów.
    :return: Lista słowników zawierających interakcje leków ze szlakami.
    """
    pathways_element = drug.find("ns:pathways", namespace)
    if pathways_element is None:
        return []

    interactions = []
    for pathway in pathways_element.findall("ns:pathway", namespace):
        smpdb_id_element = pathway.find("ns:smpdb-id", namespace)
        smpdb_id = smpdb_id_element.text if smpdb_id_element is not None else None

        pathway_name_element = pathway.find("ns:name", namespace)
        pathway_name = pathway_name_element.text if pathway_name_element is not None else None

        drugs_interaction_element = pathway.find("ns:drugs", namespace)
        if drugs_interaction_element is None:
            return []

        for drug_interaction in drugs_interaction_element.findall("ns:drug", namespace):
            id_element = drug_interaction.find("ns:drugbank-id", namespace)
            id = id_element.text if id_element is not None else None

            name_element = drug_interaction.find("ns:name", namespace)
            name = name_element.text if name_element is not None else None

            interactions.append({
                "smpdb-id": smpdb_id,
                "pathway-name": pathway_name,
                "drugbank-id": id,
                "drug-name": name
            })

    return interactions


def extract_pathways_drugs_interaction(xml_file: str, namespace: dict) -> pd.DataFrame:
    """
    Tworzy DataFrame zawierający interakcje leków ze szlakami metabolicznymi.

    :param xml_file: Ścieżka do pliku XML.
    :param namespace: Słownik namespace'ów do wyszukiwania elementów.
    :return: DataFrame z informacjami o interakcjach leków ze szlakami.
    """
    rows_iter = iter_rows(xml_file, namespace, extract_pathways_drugs_interaction_row)
    rows = [row for row in rows_iter if row]
    return pd.DataFrame(rows)


def create_pathways_graph(pathways_interaction_df: pd.DataFrame) -> None:
    """
    Tworzy i wizualizuje graf przedstawiający interakcje leków i szlaków.

    :param pathways_interaction_df: DataFrame zawierający interakcje leków ze szlakami metabolicznymi.
    """
    G = nx.Graph()

    # Wierzchołki
    drugs = pathways_interaction_df["drug-name"].dropna().unique()
    pathways = pathways_interaction_df["pathway-name"].dropna().unique()

    # Dodanie wierzchołków
    G.add_nodes_from(drugs, bipartite=0)
    G.add_nodes_from(pathways, bipartite=1)

    # Dodanie krawędzi
    edges = list(zip(pathways_interaction_df["drug-name"], pathways_interaction_df["pathway-name"]))
    G.add_edges_from(edges)

    # Rysowanie grafu
    plt.figure(figsize=(12, 8))

    pos = {}
    pos.update((node, (0, i)) for i, node in enumerate(drugs))
    pos.update((node, (1, i)) for i, node in enumerate(pathways))

    nx.draw(G, pos, node_size=1000, font_size=8, edge_color='gray')

    for node, (x, y) in pos.items():
        plt.text(
            x, y, node, ha='center', va='center',
            bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="white")
        )

    plt.title("Graf interakcji leków i szlaków")
    plt.show()

def create_pathway_count_df(pathways_interaction_df: pd.DataFrame, drugs_df: pd.DataFrame) -> pd.DataFrame:
    """
    Tworzy DataFrame zawierający liczbę szlaków dla każdego leku.

    :param pathways_interaction_df: DataFrame zawierający informacje o interakcjach leków ze szlakami.
    :param drugs_df: DataFrame zawierający informacje o lekach.
    :return: DataFrame zawierający identyfikatory leków, ich nazwy oraz liczbę szlaków, z którymi wchodzą w interakcję.
    """
    pathway_counts = (
        pathways_interaction_df.groupby("drugbank-id")["pathway-name"]
        .count()
        .reset_index()
        .rename(columns={"pathway-name": "pathway-count"})
    )
    merged_df = drugs_df.merge(pathway_counts, on="drugbank-id", how="left")
    merged_df["pathway-count"] = merged_df["pathway-count"].fillna(0).astype(int)
    return merged_df[["drugbank-id", "name", "pathway-count"]]

def plot_pathway_count(pathways_count_df: pd.DataFrame) -> None:
    """
    Tworzy histogram przedstawiający rozkład liczby szlaków, z którymi wchodzą w interakcję poszczególne leki.

    :param pathways_count_df: DataFrame zawierający identyfikatory leków, ich nazwy oraz liczbę szlaków,
        z którymi wchodzą w interakcję.
    """
    # Tworzenie histogramu
    plt.figure(figsize=(10, 6))

    min_val = pathways_count_df["pathway-count"].min()
    max_val = pathways_count_df["pathway-count"].max()
    bins = np.arange(min_val, max_val + 2) - 0.5

    # Rysowanie histogramu
    plt.hist(pathways_count_df["pathway-count"], bins=bins, color="cornflowerblue", edgecolor="black", alpha=0.8)

    # Ustawienie etykiet na osi X
    xticks = np.arange(min_val, max_val + 1)
    step = max(1, len(xticks) // 15)  # Maksymalnie 15 etykiet na osi X
    plt.xticks(xticks[::step])

    plt.xlabel("Liczba szlaków, z którymi lek wchodzi w interakcję")
    plt.ylabel("Liczba leków")
    plt.title("Histogram liczby szlaków z którymi dany lek wchodzi w interakcje")

    plt.grid(axis="y", linestyle="--", alpha=0.7)

    plt.show()


def extract_targets_row(drug: etree._Element, namespace: dict) -> list[dict]:
    """
    Pobiera informacje o białkach docelowych dla danego leku.

    :param drug: Element XML reprezentujący lek.
    :param namespace: Słownik namespace'ów do wyszukiwania elementów.
    :return: Lista słowników zawierających informacje o białkach docelowych.
    """
    drug_id_element = drug.find("ns:drugbank-id[@primary='true']", namespace)
    drug_id = drug_id_element.text if drug_id_element is not None else None

    drug_name_element = drug.find("ns:name", namespace)
    drug_name = drug_name_element.text if drug_name_element is not None else None


    targets_element = drug.find("ns:targets", namespace)
    if targets_element is None:
        return []

    targets = []
    for target in targets_element.findall("ns:target", namespace):
        drugbank_id_element = target.find("ns:id", namespace)
        drugbank_id = drugbank_id_element.text if drugbank_id_element is not None else None

        polypeptide_element = target.find("ns:polypeptide", namespace)
        if polypeptide_element is not None:
            source = str(polypeptide_element.get("source"))
            external_id = str(polypeptide_element.get("id"))

            name_element = polypeptide_element.find("ns:name", namespace)
            name = name_element.text if name_element is not None else None

            gene_name_element = polypeptide_element.find("ns:gene-name", namespace)
            gene_name = gene_name_element.text if gene_name_element is not None else None

            identifiers_element = polypeptide_element.find("ns:external-identifiers", namespace)

            genatlas_id = None
            if identifiers_element is not None:
                for identifier in identifiers_element.findall("ns:external-identifier", namespace):
                    resource_element = identifier.find("ns:resource", namespace)
                    resource = resource_element.text if resource_element is not None else None
                    if resource == "GenAtlas":
                        genatlas_id_element = identifier.find("ns:identifier", namespace)
                        genatlas_id = genatlas_id_element.text if genatlas_id_element is not None else None
                        break

            chromosome_location_element = polypeptide_element.find("ns:chromosome-location", namespace)
            chromosome_location = chromosome_location_element.text if chromosome_location_element is not None else None

            cellular_location_element = polypeptide_element.find("ns:cellular-location", namespace)
            cellular_location = cellular_location_element.text if cellular_location_element is not None else None

            targets.append({
                "drug-drugbank-id": drug_id,
                "drug-name": drug_name,
                "target-drugbank-id": drugbank_id,
                "source": source,
                "external-id": external_id,
                "polypepetide-name": name,
                "gene-name": gene_name,
                "gene-genatlas-id": genatlas_id,
                "chromosome-location": chromosome_location,
                "cellular-location": cellular_location
            })

    return targets


def extract_targets(xml_file: str, namespace: dict) -> pd.DataFrame:
    """
    Tworzy ramke danych zawierajaca informacje o bialkach, z ktorymi poszczegolne leki wchodza w interakcje.

    :param xml_file: Ścieżka do pliku XML.
    :param namespace: Słownik namespace'ów do wyszukiwania elementów.
    :return: DataFrame zawierający informacje o białkach docelowych.
    """
    rows_iter = iter_rows(xml_file, namespace, extract_targets_row)
    rows = [row for row in rows_iter if row]
    return pd.DataFrame(rows)


def plot_cellural_location(targets_df: pd.DataFrame) -> None:
    """
    Tworzy wykres kołowy przedstawiający procentowe występowanie targetów w różnych częściach komórki.

    :param targets_df: DataFrame z informacjami o lokalizacji komórkowej.
    """
    location_counts = targets_df["cellular-location"].value_counts()

    fig, ax = plt.subplots(figsize=(10, 8))
    wedges, texts, autotexts = ax.pie(
        location_counts,
        autopct=lambda pct: f"{pct:.1f}%" if pct > 5 else "",
        wedgeprops={"edgecolor": "white"},
        colors=plt.cm.Paired.colors,
        textprops={"fontsize": 15, "color": "black"}
    )

    # Dodanie legendy obok wykresu
    ax.legend(wedges, location_counts.index, title="Lokalizacja w komórce",
              loc="center left", bbox_to_anchor=(1, 0.5), fontsize=10)

    plt.title("Procentowe występowanie targetów w różnych częściach komórki", fontsize=14, pad=20)
    plt.tight_layout()
    plt.show()


def extract_drug_status_row(drug: etree._Element, namespace: dict) -> list[dict]:
    """
    Pobiera informacje o statusie leku z pojedynczego elementu <drug>.

    :param drug: Element XML reprezentujący lek.
    :param namespace: Słownik namespace'ów do wyszukiwania elementów.
    :return: Lista słowników zawierających informacje o statusie leku.
    """
    id_element = drug.find("ns:drugbank-id[@primary='true']", namespace)
    id_txt = id_element.text if id_element is not None else None

    groups = [g.text for g in drug.findall("ns:groups/ns:group", namespace)]

    return [{
        "drugbank-id": id_txt,
        "approved": "approved" in groups,
        "withdrawn": "withdrawn" in groups,
        "experimental": any(x in groups for x in ["experimental", "investigational"]),
        "veterinary": "vet_approved" in groups
    }]


def extract_drug_status(xml_file: str, namespace: dict) -> pd.DataFrame:
    """
    Pobiera dane o statusie leków z pliku XML i zwraca DataFrame.

    :param xml_file: Ścieżka do pliku XML.
    :param namespace: Słownik namespace'ów do wyszukiwania elementów.
    :return:  DataFrame zawierający informacje o statusie leków.
    """
    rows_iter = iter_rows(xml_file, namespace, extract_drug_status_row)
    rows = [row for row in rows_iter if row]
    return pd.DataFrame(rows)


def plot_drug_status(drug_status_df: pd.DataFrame) -> None:
    """
    Tworzy wykres kołowy przedstawiający liczbę leków o danym statusie.

    :param drug_status_df: DataFrame zawierający kolumny z informacjami o statusie leków:
                           - "approved"
                           - "withdrawn"
                           - "experimental"
                           - "veterinary"
    """
    summary = drug_status_df[["approved", "withdrawn", "experimental", "veterinary"]].sum()

    fig, ax = plt.subplots(figsize=(10, 8))
    wedges, texts, autotexts = ax.pie(
        summary,
        autopct=lambda pct: f"{pct * sum(summary) / 100:.0f}",
        wedgeprops={"edgecolor": "white"},
        colors=plt.cm.Paired.colors,
        textprops={"fontsize": 15, "color": "black",  "weight": "bold"}
    )

    ax.legend(wedges, summary.index, loc="center left", bbox_to_anchor=(1, 0.5), fontsize=12, title="Status leku")

    plt.title("Liczba leków o danym statusie", fontsize=16, pad=20)
    plt.tight_layout()
    plt.show()


def extract_interactions_row(drug: etree._Element, namespace: dict) -> list[dict]:
    """
    Pobiera informacje o interakcjach lekowych z pojedynczego elementu <drug>.

    :param drug: Element XML reprezentujący lek.
    :param namespace: Słownik namespace'ów do wyszukiwania elementów.
    :return: Lista słowników zawierających informacje o interakcjach lekowych.
    """
    data = []
    id_element = drug.find("ns:drugbank-id[@primary='true']", namespace)
    id_txt = id_element.text if id_element is not None else None

    for drug_interaction in drug.findall("ns:drug-interactions/ns:drug-interaction", namespace):
        interaction_id_element = drug_interaction.find("ns:drugbank-id", namespace)
        interaction_id_txt = interaction_id_element.text if interaction_id_element is not None else None

        data.append({
            "drugbank-id": id_txt,
            "interaction-drugbank-id": interaction_id_txt,
        })

    return data


def extract_interactions(xml_file: str, namespace: dict) -> pd.DataFrame:
    """
    Dla kazdego szlaku w bazie funckja zwraca leki, ktore wchodza z nim w interakcje.

    :param xml_file: Ścieżka do pliku XML.
    :param namespace: Słownik namespace'ów do wyszukiwania elementów.
    :return: DataFrame zawierający informacje o interakcjach lekowych.
    """
    rows_iter = iter_rows(xml_file, namespace, extract_interactions_row)
    rows = [row for row in rows_iter if row]
    return pd.DataFrame(rows)


def draw_gene_interaction_network(targets_df: pd.DataFrame, products_info: pd.DataFrame, gene_name: str) -> None:
    """
    Tworzy i wizualizuje sieć interakcji dla danego genu, pokazując powiązania między genem, lekami i produktami leczniczymi.

    :param targets_df: DataFrame zawierający informacje o interakcjach między lekami a genami.
    :param products_info: DataFrame zawierający informacje o produktach.
    :param gene_name: Nazwa genu, dla którego ma zostać wygenerowana sieć interakcji.
    """
    # Filtrowanie interakcji dla danego genu
    gene_targets = targets_df[targets_df["gene-genatlas-id"] == gene_name]

    # Pobranie powiązanych substancji leczniczych
    drugs = gene_targets.drop_duplicates(subset="drug-drugbank-id")[["drug-drugbank-id", "drug-name"]]

    # Pobranie produktów zawierających te substancje
    filtered_products = products_info[products_info["drugbank-id"].isin(drugs["drug-drugbank-id"])][
        ["drugbank-id", "drug-name", "product-name"]]
    filtered_products = filtered_products.drop_duplicates(subset=["drugbank-id", "product-name"])
    # Tworzenie grafu
    G = nx.Graph()
    G.add_node(gene_name + " (gene)", category="gene")
    G.add_nodes_from(drugs['drug-name'] + '[' + drugs['drug-drugbank-id'] + ']' + ' (drug)', category="drug")
    G.add_nodes_from(filtered_products["product-name"].unique() + " (product)", category="product")

    # Dodawanie krawędzi
    for index, drug in drugs.iterrows():
        G.add_edge(gene_name + " (gene)", f"{drug['drug-name']}[{drug['drug-drugbank-id']}] (drug)")
    for _, row in filtered_products.iterrows():
        G.add_edge(f"{row['drug-name']}[{row['drugbank-id']}] (drug)", row["product-name"] + " (product)")

    # Ustalanie pozycji wierzchołków
    pos = {}
    height = max(len(drugs["drug-name"]), len(filtered_products["product-name"].unique()))
    pos[gene_name + " (gene)"] = (0, height // 2)
    pos.update((drug, (1, int(i * (height - 1) / (len(drugs["drug-name"]) - 1)))) for i, drug in
               enumerate(drugs['drug-name'] + '[' + drugs['drug-drugbank-id'] + ']' + ' (drug)'))
    pos.update(
        (product, (2, int(i * (height - 1) / (len(filtered_products["product-name"].unique()) - 1)))) for i, product in
        enumerate(filtered_products["product-name"].unique() + " (product)"))

    # Rysowanie grafu
    plt.figure(figsize=(10, 8))

    color_map = {
        'gene': 'blue',
        'drug': 'green',
        'product': 'red'
    }
    node_colors = [color_map[G.nodes[n]['category']] for n in G.nodes]

    nx.draw(G, pos, node_color=node_colors, with_labels=False, node_size=500, edge_color='gray')

    # Dodanie etykiet z obramowaniem
    for node, (x, y) in pos.items():
        plt.text(
            x, y, node, ha='center', va='center',
            bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="white"),
            fontsize=8
        )

    plt.title(f"Interakcje dla genu: {gene_name}")
    plt.show()


def fetch_protein_data(protein_id: str) -> etree._Element | None:
    """
    Pobiera dane o białku z bazy UniProt w formacie XML na podstawie podanego identyfikatora białka.

    :param protein_id: Identyfikator białka w bazie UniProt.
    :return: Drzewo XML z danymi o białku lub None, jeśli nie udało się pobrać danych.
    """
    url = f"https://www.uniprot.org/uniprot/{protein_id}.xml"
    response = requests.get(url)
    if response.status_code == 200:
        tree = etree.fromstring(response.content)
        return tree
    else:
        return None


def get_name(protein_id: str) -> str | None:
    """
    Pobiera pełną nazwę białka na podstawie jego identyfikatora w bazie UniProt.

    :param protein_id: Identyfikator białka w bazie UniProt.
    :return: Pełna nazwa białka lub None, jeśli nie udało się jej uzyskać.
    """
    tree = fetch_protein_data(protein_id)
    if tree is not None:
        namespace = {"ns": "http://uniprot.org/uniprot"}
        entry = tree.find("ns:entry", namespace)
        if entry is not None:
            protein = entry.find("ns:protein/ns:recommendedName/ns:fullName", namespace)
            if protein is not None:
                return protein.text
            else:
                return None
        else:
            return None
    else:
        return None


def draw_protein_interactions(protein_id: str, targets_df: pd.DataFrame) -> None:
    """
    Tworzy i wizualizuje sieć interakcji dla danego białka, pokazując powiązania z innymi białkami oraz lekami.

    :param protein_id: Identyfikator białka w bazie UniProt.
    :param targets_df: DataFrame zawierający informacje o interakcjach między lekami a białkami.
    """
    protein_data = fetch_protein_data(protein_id)
    protein_name = get_name(protein_id)
    if protein_name is None or protein_data is None:
        return

    namespace = {"ns": "http://uniprot.org/uniprot"}
    interactions = protein_data.findall("ns:entry/ns:comment/ns:interactant/ns:id", namespace)
    interactions = set([get_name(int.text) for int in interactions])
    interactions = [int for int in interactions if int is not None and int != protein_name]

    drugs = targets_df[targets_df['external-id'] == protein_id]

    G = nx.Graph()

    G.add_node(protein_name, category="main-protein")

    for interaction in interactions:
        G.add_node(interaction, category="protein")
        G.add_edge(protein_name, interaction)

    if not drugs.empty:
        for _, row in drugs.iterrows():
            G.add_node(row['drug-name'] + '[' + row['drug-drugbank-id'] + ']', category="drug")
            G.add_edge(protein_name, row['drug-name'] + '[' + row['drug-drugbank-id'] + ']')

    plt.figure(figsize=(10, 8))

    pos = {}
    height = max(len(drugs["drug-name"]), len(interactions))
    pos[protein_name] = (1, height // 2)
    pos.update((drug, (0, int(i * (height - 1) / (len(drugs["drug-name"]) - 1)))) for i, drug in
               enumerate(drugs["drug-name"] + '[' + drugs['drug-drugbank-id'] + ']'))
    pos.update((interaction, (2, int(i * (height - 1) / (len(interactions) - 1)))) for i, interaction in
               enumerate(interactions))

    color_map = {
        'main-protein': 'red',
        'drug': 'green',
        'protein': 'blue'
    }
    node_colors = [color_map[G.nodes[n]['category']] for n in G.nodes]

    nx.draw(G, pos, node_color=node_colors, with_labels=False, node_size=1000, edge_color='gray')
    for node, (x, y) in pos.items():
        plt.text(
            x, y, node, ha='center', va='center',
            bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="white")
        )

    plt.title(f"Interakcje białka {protein_name} z lekami i innymi białkami")
    plt.axis('off')
    plt.show()
