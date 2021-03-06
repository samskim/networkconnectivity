# network connectivity
Kim S, Dai C, Hormozdiari F, van de Geijn B, Gazal S, Park Y, O’Connor L, Amariuta T, Loh PR, Finucane H, Raychaudhuri S, Price AL. Genes with high network connectivity are enriched for disease heritability. Am J Hum Genet. 2019 May 2;104(5):896-913. PMID: 31051114.

Paper link: https://pubmed.ncbi.nlm.nih.gov/31051114/

Annotations available at https://data.broadinstitute.org/alkesgroup/LDSCORE/Kim_pathwaynetwork/

Source codes used for the analysis:

network_calculate_connectivity.py: calculate network connectivity metrics (e.g. closeness centrality)

get_consensus_network.R: make a consensus network by taking intersection of multiple networks

bedToContAnnot.py: make .annot files from .bed files for continuous-valued annotation (weight is not just 0 or 1)

bedToAnnot.py: helping code for bedToBinaryAnnot_folder_window.py

bedToBinaryAnnot_folder_window.py: make .annot files from .bed files for binary annotation

annotToLD_thinannot.py: make .ldscore.gz from .annot files

PartitionHeritability_baselineLD.py: run S-LDSC partition heritability given .annot and .ldscore.gz files

PartitionHeritability_baselineLD_run.sh: run PartitionHeritability_baselineLD.py across multiple summary statistics

PrepareForMetaAnalysis.py: prepare for the meta-analysis by making .sd file (standard deviation of the annotation) for tau* calculation

MetaAnalysis.R: random-effect meta analysis


File formats:

gene networks: gene ID 1, gene ID 2, weight (delimiter = '\t'; we used ENTREZ ID for gene ID; if others are provided, we converted to ENTREZ ID).

.bed: CHR, START, STOP, WEIGHT (optional for continous-valued annotation with delimiter = '\t')

.annot: thin-annot version is used (one column ANNOT)

.ldscore.gz: CHR, SNP, BP, L2	(delimiter = '\t')

For other details including processing gene networks, please read the manuscript.


Update:

We have added 165 reference gene sets analyzed in Kim et al. 2020: 05212020_geneset_list_curated.csv
