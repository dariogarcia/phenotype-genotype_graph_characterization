from multiprocessing import Pool, cpu_count
from extract_paths import *


def get_pair_list_by_connected_prob(phenotypes_ids, genotypes_ids, genes_ids,
                                    phenotypes_links, genotypes_links,
                                    phenotypes_genes_links,
                                    genotypes_genes_links, continuing=False,
                                    total_results_path='../results/total_results.pkl',
                                    list_elems_path='../results/list_elems.pkl',
                                    type_index_path='../results/type_index.pkl'):
    """
    Given a list of genotypes and phenotypes, which may be linked through genes,
    order every genotype-phenotype pair according to the probability that they
    are directly connected by an (undiscovered) gene.
    WARNING: This method is parallelized and will use all available CPUs.
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
        -continuing: Is this process a continuation of a partial execution?
            +Type: bool
        -total_results_path: path to partial and previous total_results
            +Type: str
        -list_elems_path: path to partial and previous list_elems
            +Type: str
        -type_index_path: path to partial and previous type_index
            +Type: str

    Returns:
        None (see previous Warning)
    """
    # Initialize structures
    if not continuing:
        # To be replaced with structures representing the probabilities
        list_elems = []
        type_index = {}
        total_results = np.empty([0, 0])
    # If we are continuing a previous partial execution
    # load the precomputed values
    else:
        list_elems = pickle.load(open(list_elems_path, 'rb'))
        type_index = pickle.load(open(type_index_path, 'rb'))
        total_results = pickle.load(open(total_results_path, 'rb'))
    # For each phenotype
    for p_id in phenotypes_ids:
        # Get the list of linked genes
        p_genes = list(set([i[1] for i in phenotypes_genes_links if i[0] == p_id]))
        # And the list of unlinked genotypes
        # ...but first, find the linked genotypes
        p_linked_genotypes = list(set([i[0] for i in genotypes_genes_links if i[1] in p_genes]))
        # keep the rest
        p_genotypes = list(set(genotypes_ids).difference(p_linked_genotypes))
        # If we are continuing a previous partial execution avoid the already
        # computed pairs
        if continuing:
            # Remove the already computed ones
            p_genotypes = [x for x in p_genotypes if (
                p_id, x) not in list_elems]
        # Launch the computation for each unlinked genotype
        pool = Pool(2)
        prob_connected = pool.map(estimate_connected_probability,
                                  itertools.izip(itertools.repeat(p_id), p_genotypes,
                                                 itertools.repeat(phenotypes_ids),
                                                 itertools.repeat(genotypes_ids),
                                                 itertools.repeat(genes_ids),
                                                 itertools.repeat(phenotypes_links),
                                                 itertools.repeat(genotypes_links),
                                                 itertools.repeat(phenotypes_genes_links),
                                                 itertools.repeat(genotypes_genes_links)))
        pool.close()
        pool.join()
        # Add probability and pair to previous results.
        # N.B.: Results must be ordered according to decreasing connected prob.

        # Persist partial results (?)
        persist_alternative_paths(
            total_results, list_elems, type_index, disconnected=True)
    return


def estimate_connected_probability(argv, statistics_path='../data_files/statistics.txt'):
    """
    Given a phenotype source, a genotype target, and a list of other
    genotypes and phenotypes, which may be linked through genes,
    estimate the probability that a gene connects the source and target.

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
    # Find all paths between phenotype-genotype pair with max length 4.
    path_list = pool.map(find_phenotype_genotype_alternative_paths,
                         itertools.izip(itertools.repeat(p_id), p_genotypes,
                                        itertools.repeat(phenotypes_ids),
                                        itertools.repeat(genotypes_ids),
                                        itertools.repeat(genes_ids),
                                        itertools.repeat(phenotypes_links),
                                        itertools.repeat(genotypes_links),
                                        itertools.repeat(phenotypes_genes_links),
                                        itertools.repeat(genotypes_genes_links)))
    # Estimate probability that the pair is indeed connected.

    with open(ph_gen_path) as f:
        # For each kind of path with mean and std dev.

        # Normalize values and build normal distr.
        # Calculate normal distr. probability for such path.

        # Connected prob = product of all previously calculated prob.
    return
