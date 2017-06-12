from multiprocessing import Pool, cpu_count
import itertools
from graph_methods import *
import numpy as np
import cPickle as pickle

def persist_partial_results(total_results, list_elems, type_index):
    """
    Persists partial results on disk using numpy

    Args:
        -total_results: paths stored
            +Type: 2D numpy.array. Rows are phenotype-genotype pairs. 
                Columns are path types
        - list_elems: list of computed pairs. Serves as row index.
            +Type: list[(phenotype_id ,genotype_id)]
        -type_index: dictionary of path type and column index
            +Type: dict{path_type, index}

    Returns:
        None. Persists data.
    """
    #TODO: files are currently overwritten. Consider when this may not be desirable
    pickle.dump(total_results, open("../results/total_results.pkl", "wb"))
    pickle.dump(list_elems, open("../results/list_elems.pkl", "wb"))
    pickle.dump(type_index, open("../results/type_index.pkl", "wb"))
    return
    
def merge_alternative_paths(total_results, list_elems, type_index, partial_list):
    """
    Given a set of paths already processed, and a new set, merge both.

    Args:
        -total_results: previous paths stored
            +Type: 2D numpy.array. Rows are phenotype-genotype pairs. 
                Columns are path types
        - list_elems: list of computed pairs. Serves as row index.
            +Type: list[(phenotype_id ,genotype_id)]
        -type_index: dictionary of path type and column index
            +Type: dict{path_type, index}
        -partial_list: new pairs data to be added
            +Type: (phenotype_id, genotype_id, dict{path_type, frequency})

    Returns:
        -total_results: merged paths. added rows (one per pair)
            and maybe some columns (if new path type found)
            +Type: 2D numpy.array. Rows are phenotype-genotype pairs. 
                Columns are path types
        - list_elems: updated list of computed pairs (one new per pair)
            +Type: list[(phenotype_id ,genotype_id)]
        -type_index: dictionary of path type and column index.
            one new entry per new type of path
            +Type: dict{path_type, index}
    """
    #For each pair processed
    for current_paths in partial_list:
        #Get the pair in question
        p_id = current_paths[0]
        g_id = current_paths[1]
        #Append to the list
        list_elems.append((p_id, g_id))
        #Initialize values to zero
        current_np = np.zeros([len(type_index.keys())])
        #For each path
        paths = current_paths[2]
        for k,v in paths.iteritems():
            #If the type of path already exists
            if k in type_index.keys():
                #Set the value
                current_np[type_index[k]] = v
            #New type of path
            else:
                #Assign an index
                type_index[k] = total_results.shape[1]
                #Add a column of zeros to the total 
                empty_col = np.zeros([total_results.shape[0],1])
                total_results = np.hstack((total_results,empty_col))
                #Add the value to current case
                current_np = np.append(current_np,v)
        total_results = np.vstack((total_results,current_np))
    return total_results, list_elems, type_index

def get_connected_phenotype_genotype_alternative_paths(phenotypes_ids,\
        genotypes_ids, genes_ids, phenotypes_links, genotypes_links,\
        phenotypes_genes_links, genotypes_genes_links):
    """
    Given a list of genotypes and phenotypes, which may be linked through genes,
    for every linked genotype-phenotype pair find all the alternative paths
    when removing the linking gene/s.
    WARNING: This method takes a while to compute (i.e., probably weeks).
    For this reason, results are stored periodically on disc, and these are not returned.

    Args:
        -phenotypes_ids: List of phenotypes
            +Type: list[str]
        -genotypes_ids: List of genotypes 
            +Type: list[str]
        -genes_ids: List of genes 
            +Type: list[str]
        -phenotypes_links: List of phenotype-phenotype links 
            +Type: list[(str,str)]
        -genotypes_links: List of genotype-genotype links 
            +Type: list[(str,str)]
        -phenotypes_genes_links: List of phenotype-genes links 
            +Type: list[(str,str)]
        -genotypes_genes_links: List of genotype-genes links 
            +Type: list[(str,str)]

    Returns:
        None (see previous Warning)
    """
    #Initialize structures
    list_elems = []
    type_index = {}
    total_results = np.empty([0,0])
    #For each phenotype
    for p_id in phenotypes_ids:
        #Get the list of linked genes
        p_genes = list(set([i[1] for i in phenotypes_genes_links if i[0]==p_id]))
        #And the list of linked genotypes
        p_genotypes = list(set([i[0] for i in genotypes_genes_links if i[1] in p_genes]))
        #Launch the computation for each linked genotype
        pool = Pool(cpu_count())
        partial_list = pool.map(find_phenotype_genotype_alternative_paths,\
                itertools.izip(itertools.repeat(p_id), p_genotypes,\
                itertools.repeat(phenotypes_ids), itertools.repeat(genotypes_ids),\
                itertools.repeat(genes_ids), itertools.repeat(phenotypes_links),\
                itertools.repeat(genotypes_links), itertools.repeat(phenotypes_genes_links),\
                itertools.repeat(genotypes_genes_links)))
        pool.close()
        pool.join()
        #Merge list with previous results
        total_results, list_elems, type_index = merge_alternative_paths(total_results, list_elems, type_index, partial_list)
        #Persist partial results
        persist_partial_results(total_results, list_elems, type_index)
    return

def find_phenotype_genotype_alternative_paths(argv):
#def find_paths(p_id, g_id, phenotypes_ids, genotypes_ids, genes_ids, phenotypes_links, genotypes_links, phenotypes_genes_links, genotypes_genes_links):
    #TODO: Due to multiprocessing, parameters are passed within a list
    #TODO: and unpacked here. This can be probably fixed.
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
    #TODO: add directionality of graph as parameter
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
    #Return the pair and their alternative paths 
    return p_id,g_id,paths_codes
