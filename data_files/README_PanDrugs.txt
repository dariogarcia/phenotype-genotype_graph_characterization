A brief description of the columns of the document:

(1) gene_symbol
(2) source: database or source where the drug-gene association came from
(3) source_drug_name: name of the drug in the original database
(4) standard_drug_name: standardized name for the drug
(5) show_drug_name: INN name or short name for the drug
(6) family: drug family based on the KEGG target-based classification of drugs and Connectivity Map classification
(7) status: FDA approval status of the drug ("Approved","Withdrawn","Clinical trials","Experimental","","Undefined")
(8) pathology: condition the drug is prescribed for
(9) cancer: cancer type (based on location) for approved drugs in cancer; or 'clinical cancer' label for approved drugs for other conditions but in clinical trials in cancer; or 'cancer' label for drugs in clinical trials in cancer
(10) extra: FDA label for approved drugs in cancer
(11) extra2: type of therapy for approved drugs in cancer
(12) pathways: KEGG pathway ID for the pathways the gene is involved in
(13) target_marker: association type between the drug and the gene:
target: when the gene can be directly targeted by the drug
marker: the gene is a biomarker (when its genetic status is associated with the drug response by clinical or pre-clinical evidences)
(14) resistance: resistance or sensitivity response to the drug based on the alteration in the gene
(15) alteration: type of alteration that causes the sensitivity or resistance response
(16) ind_pathway: genes located upstream the gene_symbol (until a 4 level of distance). This is used to associate treatments in a member pathway manner
(17) score: pre-computed score for the drug-gene association. It gives an idea of the suitability of the treatment when the gene is altered, based on the approval status of the drug, the use in cancer and the drug-gene association

source: http://pandrugs.bioinfo.cnio.es/
