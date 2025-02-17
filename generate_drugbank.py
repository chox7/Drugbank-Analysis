import random
import copy
from lxml import etree


def generate_drugbank_file(
        input_file: str = "drugbank_partial.xml",
        output_file: str = "drugbank_partial_and_generated.xml",
        num_new_drugs: int = 19900
):
    """
    Wczytuje plik XML z danymi leków, kopiuje istniejące elementy,
    a następnie generuje nowe elementy <drug> na podstawie losowych pól z już
    istniejących leków. Wynik zapisywany jest do nowego pliku XML.

    :param input_file: Ścieżka do wejściowego pliku XML.
    :param output_file: Ścieżka do wyjściowego pliku XML.
    :param num_new_drugs: Liczba nowych elementów <drug> do wygenerowania.
    """
    namespace = {"ns": "http://www.drugbank.ca"}
    tree = etree.parse(input_file)
    root = tree.getroot()
    drugs = root.findall("ns:drug", namespace)

    if not drugs:
        raise ValueError("Nie znaleziono żadnych elementów <drug> w pliku.")

    fields = {item.tag for drug in drugs for item in drug if "drugbank-id" not in item.tag}

    with open(output_file, "wb") as f:
        f.write(b'<?xml version="1.0" encoding="UTF-8"?>\n')
        f.write(b'<drugbank xmlns="http://www.drugbank.ca">\n')

        # Kopiujemy wszystkie istniejące leki do nowego drzewa
        for drug in drugs:
            f.write(etree.tostring(drug, encoding="utf-8", pretty_print=True))

        # Generujemy nowe leki
        for i in range(num_new_drugs):
            new_drug = etree.Element("{http://www.drugbank.ca}drug")
            id_element = etree.SubElement(new_drug, "{http://www.drugbank.ca}drugbank-id", primary="true")
            id_element.text = str(i)

            # Losowo wybieramy istniejący lek i kopiujemy jego atrybuty do nowego leku
            random_drug = random.choice(drugs)
            new_drug.attrib.update(random_drug.attrib)

            # Dla każdego z ustalonych pól losujemy nowe dane
            for field in fields:
                field_element = etree.SubElement(new_drug, field)
                random_field = random.choice(drugs).find(field)
                if random_field is not None:
                    field_element.text = random_field.text
                    for child in random_field:
                        field_element.append(copy.deepcopy(child))
                    if random_field.attrib:
                        field_element.attrib.update(random_field.attrib)

            # Dodajemy nowo utworzony lek do głównego drzewa
            f.write(etree.tostring(new_drug, encoding="utf-8", pretty_print=True))

        f.write(b'</drugbank>')


if __name__ == "__main__":
    generate_drugbank_file()
