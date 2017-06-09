from multiprocessing import Pool, cpu_count
import itertools
from graph_methods import *

def get_connected_phenotype_genotype_alternative_paths(phenotypes_ids,\
        genotypes_ids, genes_ids, phenotypes_links, genotypes_links,\
        phenotypes_genes_links, genotypes_genes_links):
    #For each phenotype
    for p_id in phenotypes_ids:
        #Get the list of linked genes
        p_genes = list(set([i[1] for i in phenotypes_genes_links if i[0]==p_id]))
        #And the list of linked genotypes
        p_genotypes = list(set([i[0] for i in genotypes_genes_links if i[1] in p_genes]))
        #For each linked genotype
        pool = Pool(cpu_count())
        list_paths = pool.map(find_phenotype_genotype_alternative_paths,\
                itertools.izip(itertools.repeat(p_id), p_genotypes[:10],\
                itertools.repeat(phenotypes_ids), itertools.repeat(genotypes_ids),\
                itertools.repeat(genes_ids), itertools.repeat(phenotypes_links),\
                itertools.repeat(genotypes_links), itertools.repeat(phenotypes_genes_links),\
                itertools.repeat(genotypes_genes_links)))
        pool.close()
        pool.join()
        return list_paths

def find_phenotype_genotype_alternative_paths(argv):
#def find_paths(p_id, g_id, phenotypes_ids, genotypes_ids, genes_ids, phenotypes_links, genotypes_links, phenotypes_genes_links, genotypes_genes_links):
    """
    Find all paths between a phenotype and a genotype which share a gene.
    Return the number and type of paths without using the phenotype-shared_gene link.

    Args:
        -p_id: Phenotype source of the path
            +Type: str
        -g_id: Genotype target of the path
            +Type: str
        -phenotypes_ids: List of phenotypes to be used as vertices
            +Type: list[str]
        -genotypes_ids: List of genotypes to be used as vertices
            +Type: list[str]
        -genes_ids: List of genes to be used as vertices
            +Type: list[str]
        -phenotypes_links: List of phenotype-phenotype links to be used as edges
            +Type: list[(str,str)]
        -genotypes_links: List of genotype-genotype links to be used as edges
            +Type: list[(str,str)]
        -phenotypes_genes_links: List of phenotype-genes links to be used as edges
            +Type: list[(str,str)]
        -genotypes_genes_links: List of genotype-genes links to be used as edges
            +Type: list[(str,str)]

    Returns:
        -paths_codes: Dictionary containing all path types and their frequency
            +Type: dict{(str,int)}
    """
    #TODO: add directionality of graph as parameter
    
    #TODO: Due to multiprocessing, parameters are passed within a list
    #TODO: and unpacked here. This can be probably fixed.
    p_id = argv[0]
    g_id = argv[1]
    phenotypes_ids = argv[2]
    genotypes_ids = argv[3]
    genes_ids = argv[4]
    phenotypes_links = argv[5]
    genotypes_links = argv[6]
    phenotypes_genes_links = argv[7]
    genotypes_genes_links = argv[8]
    #Find the common genes of this pair
    common_genes = [gene for gene in genes_ids if (p_id,gene) in phenotypes_genes_links \
            and (g_id,gene) in genotypes_genes_links]
    #Remove the links directly linking the source phenotype and the target genotype
    phenotypes_genes_pruned_links = [i for i in phenotypes_genes_links \
            if i[0]!= p_id or i[1] not in common_genes]
    #Create the graph
    graph = build_graph(phenotypes_ids, genotypes_ids, genes_ids, phenotypes_links, genotypes_links, phenotypes_genes_pruned_links, genotypes_genes_links, undirected=True)
    #Find all paths
    paths = find_all_paths(graph, graph.vs.find(p_id).index, graph.vs.find(g_id).index,maxlen=4)
    #Compute the type and frequency of all paths
    paths_codes = {}
    for current_path in paths:
        current_code  = get_phenotype_genotype_path_code(graph,current_path,p_id,common_genes,g_id)
        if current_code in paths_codes:
            paths_codes[current_code]+=1
        else:
            paths_codes[current_code]=1
    return paths_codes 
