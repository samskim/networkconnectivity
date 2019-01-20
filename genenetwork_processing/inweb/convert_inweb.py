import pandas
import numpy

def biomart_search(InWeb_interactor_A, InWeb_interactor_B, gene_set_dict, gene_map):
    interactor_A_gene = numpy.append(numpy.empty(0, int), gene_map.get(InWeb_interactor_A,[]))
    interactor_B_gene = numpy.append(numpy.empty(0, int), gene_map.get(InWeb_interactor_B,[]))
    if (interactor_A_gene.size is not 0) & (interactor_B_gene.size is not 0):
        for gene_a in interactor_A_gene:
            gene_set_dict[gene_a] = gene_set_dict.get(gene_a, set([])).union(interactor_B_gene)
    return None

def convert_inweb(inweb_dataset, biomart_dataset):
	#import InWeb data set, renaming the columns according to documentation
	#see https://www.intomics.com/inbio/map/#downloads under InWeb_InBioMap version 2015_11_02
	inweb = pandas.read_table(inweb_dataset,usecols = [0,1,9,10], header = None)
	inweb.columns = ['interactor A ID','interactor B ID', 'NCBI Taxonomy A', 'NCBI Taxonomy B']

	#confirm that all proteins/genes in the dataset are human related
	assert (''.join(inweb['NCBI Taxonomy A'].unique()) == 'taxid:9606(Homo sapiens)') and (''.join(inweb['NCBI Taxonomy A'].unique()) == 'taxid:9606(Homo sapiens)')

	#filter the InWeb dataset to isolate potential sets of size 10 to 500 (for pathways analysis only)
	filtered = pandas.DataFrame(inweb.loc[:, 'interactor A ID'].str.split(':').str[1])
	filtered['interactor B ID'] = inweb.loc[:, 'interactor B ID'].str.split(':').str[1]
	#filtered = filtered.groupby(filtered.loc[:,'interactor A ID']).filter(lambda x: (len(x) >= 10) and (len(x)<=500))

	#import BioMart data, which will include UniProtKB ID and a map to NCBI GeneID 
	biomart = pandas.read_csv(biomart_dataset, usecols = ['UniProtKB/Swiss-Prot ID','NCBI gene ID'])
	gene_map = biomart.dropna().drop_duplicates()
	gene_map.set_index('UniProtKB/Swiss-Prot ID', inplace = True)
	gene_map = gene_map['NCBI gene ID'].astype(int)

	#initialize dictionary to store the gene sets
	gene_sets_dictionary = {}

	#build the gene sets with numpy vectorize
	vectorized = numpy.vectorize(biomart_search, otypes=[dict], excluded = [2,3])
	vectorized(filtered['interactor A ID'].values, filtered['interactor B ID'].values, gene_sets_dictionary, gene_map)

	#build pandas DataFrame from the resulting gene sets dictionary 
	gene_sets = pandas.Series(gene_sets_dictionary, dtype = object).apply(list)
	taxonomy = pandas.Series([9606]*len(gene_sets.index), index = gene_sets.index)
	inweb_gene_set = pandas.concat([gene_sets, taxonomy], axis =1)
	inweb_gene_set.columns = ['Gene Set', 'Taxonomy ID']
	inweb_gene_set.index.rename('NCBI Gene ID', inplace= True)

	inweb_gene_set.to_csv('inweb_dataset.csv')

	return inweb_gene_set


def main():
	#default inweb dataset file should be named 'core.psimitab', default biomart dataset file name should be 'mart_export.txt'
	convert_inweb('core.psimitab','mart_export.txt')

if __name__ == '__main__':
    main()