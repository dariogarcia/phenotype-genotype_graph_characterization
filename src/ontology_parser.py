"""
This file contains methods to parse the phenotype and genotype ontologies,
and the files containing the relations between genes and phenotypes/genotypes.
Methods return the ids of the elements in lists of strings (if elements)
or lists of tuples (string,string) (if links).
"""

def load_phenotypes_ids(po_path):
    """
    Read phenotypes ids from the phenotype ontology.

    Args:
        -po_path: Complete path to the phenotype ontology.
            +Type: str
    Returns:
        -phenotype_ids: list of phenotype ids without repetition
            +Type: list[phenotype]
    """
    #Check if the file is the expected one
    import hashlib
    if hashlib.md5(open(po_path, 'rb').read()).hexdigest() != 'a48f293dfed9eef6aabbe400ac32f032':
        print 'WARNING: Hash of file',po_path,'does not match expected value.'
        print 'WARNING: Parsing may not be correct.'
    #Process the file
    phenotypes_ids = []
    with open(po_path) as f:
        for line in f:
            #Found a new term
            if line=='[Term]\n':
                term_id = next(f).split()[1]
                if 'HP' not in term_id:
                    raise Exception('Wrong phenotype id')
                phenotypes_ids.append(term_id.splitlines()[0])
    #Return a list without repeated elements
    return list(set(phenotypes_ids))


def load_genotypes_ids(go_path):
    """
    Read genotypes ids from the genotype ontology.

    Args:
        -go_path: Complete path to the genotype ontology.
            +Type: str
    Returns:
        -genotype_ids: list of genotype ids without repetition
            +Type: list[genotype]
    """
    #Check if the file is the expected one
    import hashlib
    if hashlib.md5(open(go_path, 'rb').read()).hexdigest() != '569240d38bd319fbf200ff1d09606640':
        print 'WARNING: Hash of file',go_path,'does not match expected value.'
        print 'WARNING: Parsing may not be correct.'
    #Process the file
    genotypes_ids = []
    with open(go_path) as f:
        for line in f:
            #Found a new term
            if line=='[Term]\n':
                genotypes_ids.append(next(f).split()[1].splitlines()[0])
    #Return a list without repeated elements
    return list(set(genotypes_ids))


def load_phenotype_phenotype_links(po_path):
    """
    Read links between phenotypes from the phenotype ontology.

    Args:
        -po_path: Complete path to the phenotype ontology.
            +Type: str
    Returns:
        -phenotype_phenotype_links: list of links without repetition. first element is_a second element
            +Type: list[(phenotype,phenotype)]
    """
    #Check if the file is the expected one
    import hashlib
    if hashlib.md5(open(po_path, 'rb').read()).hexdigest() != 'a48f293dfed9eef6aabbe400ac32f032':
        print 'WARNING: Hash of file',po_path,'does not match expected value.'
        print 'WARNING: Parsing may not be correct.'
    #Process the file
    phenotypes_links = []
    with open(po_path) as f:
        source_phenotype = ''
        target_phenotype = ''
        for line in f:
            #Found a new term
            if line.startswith('id: '):
                source_phenotype = line.split()[1].splitlines()[0]
                if 'HP' not in source_phenotype:
                    raise Exception('Wrong phenotype id ->',line,'<- should start with \'HP\'')
                continue
            #Make sure we have a source at this point
            if source_phenotype == ' ':
                raise Exception('Failed to find a source phenotype id. \'is_a\' found before \'id:\'')
            #Find the target term
            if line.startswith('is_a: '):
                target_phenotype = line.split()[1].splitlines()[0]
                if 'HP' not in target_phenotype:
                    raise Exception('Wrong phenotype id ->',line,'<- should start with \'HP\'')
                phenotypes_links.append((source_phenotype,target_phenotype))
                continue
    #Return a list without repeated elements
    return list(set(phenotypes_links))


def load_genotype_genotype_links(go_path):
    """
    Read links between genotypes from the genotype ontology.

    Args:
        -go_path: Complete path to the genotype ontology.
            +Type: str
    Returns:
        -genotype_genotype_links: list of links without repetition. first element is_a second element
            +Type: list[(genotype,genotype)]
    """
    #Check if the file is the expected one
    import hashlib
    if hashlib.md5(open(go_path, 'rb').read()).hexdigest() != '569240d38bd319fbf200ff1d09606640':
        print 'WARNING: Hash of file',go_path,'does not match expected value.'
        print 'WARNING: Parsing may not be correct.'
    #Process the file
    genotypes_links = []
    with open(go_path) as f:
        source_genotype = ''
        target_genotype = ''
        for line in f:
            #Stop reading when all terms have been read
            if line.startswith('[Typedef]'):
                break
            #Found a new term
            if line.startswith('id: '):
                source_genotype = line.split()[1].splitlines()[0]
                if 'GO' not in source_genotype:
                    raise Exception('Wrong genotype id ->',line,'<- should start with \'GO\'')
                continue
            #Make sure we have a source at this point
            if source_genotype == ' ':
                raise Exception('Failed to find a source genotype id. \'is_a\' found before \'id:\'')
            #Find the target term
            if line.startswith('is_a: '):
                target_genotype = line.split()[1].splitlines()[0]
                if 'GO' not in target_genotype:
                    raise Exception('Wrong genotype id ->',line,'<- should start with \'GO\'')
                genotypes_links.append((source_genotype,target_genotype))
                continue
    #Return a list without repeated elements
    return list(set(genotypes_links))


def load_phenotypes_genes_links(ph_gen_path):
    """
    Read links between phenotypes ids and genes from a file.

    Args:
        -ph_gen_path: Complete path to the file containing the links.
            +Type: str
    Returns:
        -phenotypes_genes_links: list of links without repetition
            +Type: list[(phenotype,gene)]
    """
    #Check if the file is the expected one
    import hashlib
    if hashlib.md5(open(ph_gen_path, 'rb').read()).hexdigest() != 'a5f7d8d8accbd1983cef710ebb586736':
        print 'WARNING: Hash of file',ph_gen_path,'does not match expected value.'
        print 'WARNING: Parsing may not be correct.'
    #Process the file
    phenotypes_genes_links = []
    with open(ph_gen_path) as f:
        #Skip the first line of syntax
        next(f)
        for line in f:
            gene = line.split()[1]
            phenotype = line.split()[-1:][0]
            if 'HP' not in phenotype:
                raise Exception('Wrong phenotype id in link ->',line,'<- should start with \'HP\'')
            phenotypes_genes_links.append((phenotype,gene))
    #Return a list without repeated elements
    return list(set(phenotypes_genes_links))


def load_genotypes_genes_links(go_gen_path):
    """
    Read links between genotypes ids and genes from a file.
    Returns only links backed by experimental results.

    Args:
        -go_gen_path: Complete path to the file containing the links.
            +Type: str
    Returns:
        -genotypes_genes_links: list of links without repetition
            +Type: list[(genotype,gene)]
    """
    #Check if the file is the expected one
    import hashlib
    if hashlib.md5(open(go_gen_path, 'rb').read()).hexdigest() != 'c943dcfbd6a6432761841304ff79d372':
        print 'WARNING: Hash of file',go_gen_path,'does not match expected value.'
        print 'WARNING: Parsing may not be correct.'
    #Process the file
    genotypes_genes_links = []
    with open(go_gen_path) as f:
        for line in f:
            #Skip comment lines
            if line[0] == '!':
                continue
            #In some cases the split needs an offset
            offset = 0 
            if line.split()[6] in ['EXP','IDA','IPI','IMP','IGI','IEP','TAS','ND','NAS','ISS','IKR','IEA','IC','IBA']:
                offset = 1 
            #Avoid those links which are not experimental results
            if line.split()[5+offset] not in ['EXP','IDA','IPI','IMP','IGI','IEP']:
                continue
            gene = line.split()[2+offset]
            genotype = line.split()[3+offset]
            if 'GO:' not in genotype:
                raise Exception('Wrong genotype id in link ->',line,'<- should start with \'GO\'')
            genotypes_genes_links.append((genotype,gene))
    #Return a list without repeated elements
    return list(set(genotypes_genes_links))


#Read genes Ids from phenotypes-gene and genotype-gene links
#Returns list of unique genes Ids
def load_all_genes(ph_gen_path, go_gen_path):
    """
    Obtains the list of genes ids which are linked to at least
    one phenotype or genotype 

    Args:
        - ph_gen_path: path to the file linking phenotypes with genes
            + Type: str
        - go_gen_path: path to the file linking genotypes with genes
            + Type: str
    
    Returns:
        -genes_ids: list of unique genes linked in the ph or go ontologies
            + Type: list[str]
    """
    #Check if the file is the expected one
    import hashlib
    if hashlib.md5(open(go_gen_path, 'rb').read()).hexdigest() != 'c943dcfbd6a6432761841304ff79d372':
        print 'WARNING: Hash of file',go_gen_path,'does not match expected value.'
        print 'WARNING: Parsing may not be correct.'
    if hashlib.md5(open(ph_gen_path, 'rb').read()).hexdigest() != 'a5f7d8d8accbd1983cef710ebb586736':
        print 'WARNING: Hash of file',ph_gen_path,'does not match expected value.'
        print 'WARNING: Parsing may not be correct.'
    #Process the file
    genes_ids = set()
    with open(ph_gen_path) as f:
        #Skip the first line of syntax
        next(f)
        for line in f:
            gene_id = line.split()[1]
            if gene_id not in genes_ids:
                genes_ids.add(gene_id)
    #Read genes from genotype-gene links
    genes_ids2 = set()
    with open(go_gen_path,) as f:
        for line in f:
            #Skip comment lines
            if line[0] == '!':
                continue
            gene_id = line.split()[2]
            genes_ids2.add(gene_id)
    #Return a list combining both sets, without repeated elements
    return list(genes_ids.union(genes_ids2))


def load_shared_genes(ph_gen_path, go_gen_path):
    """
    Obtains the list of genes ids which are shared by both
    the phenotype and the genotype ontology

    Args:
        - po_path: path to the phenotype ontology
            + Type: str
        - go_path: path to the genotype ontology
            + Type: str
        - ph_gen_path: path to the file linking phenotypes with genes
            + Type: str
        - go_gen_path: path to the file linking genotypes with genes
            + Type: str
    
    Returns:
        -genes_ids: list of genes shared by both ph and go ontologies
            + Type: list[str]
    """
    #Get the unique list of genes linked with phenotypes
    genes_phenotypes_links = load_phenotypes_genes_links(ph_gen_path)
    genes_phenotypes = set([i[1] for i in genes_phenotypes_links])
    #Get the unique list of genes linked with genotypes
    genes_genotypes_links = load_genotypes_genes_links(go_gen_path)
    genes_genotypes = set([i[1] for i in genes_genotypes_links])
    #Return the intersection in a list
    return list(genes_genotypes.intersection(genes_phenotypes))


def load_data(po_path, go_path, ph_gen_path, go_gen_path, only_shared_genes=True):
    """
    Loads and returns the phenotypes, go-terms and genes of the PHO and GO ontologies.
    It also returns their links (P-P,P-Gn,GO-GO,GO-Gn)
    Only the genes shared by both ontologies are returned.
    
    Args:
        - po_path: path to the phenotype ontology
            + Type: str
        - go_path: path to the genotype ontology
            + Type: str
        - ph_gen_path: path to the file linking phenotypes with genes
            + Type: str
        - go_gen_path: path to the file linking genotypes with genes
            + Type: str
    Returns:
        - phenotypes_ids: list of phenotype ids from the PHO
            + Type list[str]
        - genotypes_ids: list of genotype ids from the GO
            + Type list[str]
        - genes_ids: list of genes ids shared by the PHO and the GO
            + Type list[str]
        - phenotype_links: list of links between phenotype ids
            + Type list[(source str,target str)]
        - genotype_links: list of links between genotype ids
            + Type list[(source str,target str)]
        - phenotype_gene_links: list of links between phenotypes and genes
            + Type list[(source phenotype str,target gene str)]
        - genotype_gene_links: list of links between genotypes and genes
            + Type list[(source genotype str,target gene str)]
        - 
    """
    #TODO: There is redundancy in this call. load_shared_genes parses() the files
    #TODO: and so do load_phenotypes_genes_links() and load_genotypes_genes_links()
    if only_shared_genes:
        #Get list of genes ids shared by both ontologies
        genes_ids = load_shared_genes(ph_gen_path, go_gen_path)
        #Remove phenotype_gene links which relate to genes not linked by both ontologies
        dirty_phenotype_gene_links = load_phenotypes_genes_links(ph_gen_path)
        phenotype_gene_links = [l for l in dirty_phenotype_gene_links if l[1] in genes_ids]
        #Remove genotype_gene links which relate to genes not linked by both ontologies
        dirty_genotype_gene_links = load_genotypes_genes_links(go_gen_path)
        genotype_gene_links = [l for l in dirty_genotype_gene_links if l[1] in genes_ids]
    else:
        genes_ids = load_all_genes(ph_gen_path, go_gen_path) 
        phenotype_gene_links = load_phenotypes_genes_links(ph_gen_path)
        genotype_gene_links = load_genotypes_genes_links(go_gen_path)
        
    return load_phenotypes_ids(po_path), load_genotypes_ids(go_path), genes_ids, load_phenotype_phenotype_links(po_path), load_genotype_genotype_links(go_path), phenotype_gene_links, genotype_gene_links