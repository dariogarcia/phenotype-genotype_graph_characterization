"""
This file contains methods to parse the file containing the relations between
genes and drugs. They are meant to operate with methods from ontology_parser.
Methods return the ids of the elements in lists of strings (if elements)
or lists of tuples (string,string) (if links).
"""
from ontology_parser import *


def load_genes_drugs_links(gen_dr_path, genes_ids):
    """
    Read links between genes and drugs from a file.
    Returns only links backed by experimental results.

    Args:
        -gen_dr_path: Complete path to the file linking genes with drugs.
            +Type: str
        -genes_ids: List of relevant genes, previously generated.
            +Type: list[gene]
    Returns:
        -genes_drugs_links: list of links without repetition
            +Type: list[(gene,drug)]
    """
    # Check if the file is the expected one
    import hashlib
    if hashlib.md5(open(gen_dr_path, 'rb').read()).hexdigest() != '647de9dbd724e1f4fb7f792f068af5c7':
        print('WARNING: Hash of file', gen_dr_path, 'does not match expected value.')
        print('WARNING: Parsing may not be correct.')
    # Process the file
    genes_drugs_links = []
    with open(gen_dr_path) as f:
        # Skip the first line of syntax
        next(f)
        for line in f:
            split_line = line.split('\t')
            # Avoid those links for which score is 0
            if float(split_line[16]) > 0:
                continue
            gene = split_line[0]
            drug = split_line[3]
            # Aviod those links where the gene not contained in the list
            if gene in genes_ids:
                genes_drugs_links.append((gene, drug))
    # Return a list without repeated elements
    return list(set(genes_drugs_links))


# Read drug (standard name) from gene-drug links
# Returns list of unique drugs
def load_all_drugs(gen_dr_path):
    """
    Obtains the list of drugs (standard names) appearing in the gene-drug database

    Args:
        -gen_dr_path: Complete path to the file linking genes with drugs.
            +Type: str

    Returns:
        -drug_ids: list of unique drugs appearing in the file
            + Type: list[str]
    """
    # Check if the file is the expected one
    import hashlib
    if hashlib.md5(open(gen_dr_path, 'rb').read()).hexdigest() != '647de9dbd724e1f4fb7f792f068af5c7':
        print('WARNING: Hash of file', gen_dr_path, 'does not match expected value.')
        print('WARNING: Parsing may not be correct.')
    # Process the file
    drugs_ids = set()
    with open(gen_dr_path) as f:
        # Skip the first line of syntax
        next(f)
        for line in f:
            drug = line.split('\t')[3]
            if drug not in drugs_ids:
                drugs_ids.append(drug)
    # Return a list combining both sets, without repeated elements
    return list(drugs_ids)


def load_linked_drugs(gen_dr_path):
    """
    Obtains the list of drugs (standard names) which are linked to at least one gene
    that appears in a previously generated gene list.


    Args:
        -gen_dr_path: Complete path to the file linking genes with drugs.
            +Type: str
        -genes_ids: List of relevant genes, previously generated.
            +Type: list[gene]

    Returns:
        -drug_ids: list of unique drugs linked to at least one gene
            + Type: list[str]
    """
    # Check if the file is the expected one
    import hashlib
    if hashlib.md5(open(gen_dr_path, 'rb').read()).hexdigest() != '647de9dbd724e1f4fb7f792f068af5c7':
        print('WARNING: Hash of file', gen_dr_path, 'does not match expected value.')
        print('WARNING: Parsing may not be correct.')
    # Process the file
    drugs_ids = set()
    with open(gen_dr_path) as f:
        # Skip the first line of syntax
        next(f)
        for line in f:
            # Avoid those links for which score is 0
            if float(line.split('\t')[16]) > 0:
                continue
            gene = line.split('\t')[0]
            drug = line.split('\t')[3]
            if gene in genes_ids and drug not in drugs_ids:
                drugs_ids.append(drug)
    # Return a list combining both sets, without repeated elements
    return list(drugs_ids)
