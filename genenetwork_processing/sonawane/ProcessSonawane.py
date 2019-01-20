import pandas as pd
import os

dest = '/random/directory/Network_Sonawane/raw_network'
path = '/random/directory/Network_Sonawane'
os.chdir(path)

ref = pd.read_csv("Ensembl_protein_coding_grch37_biocart.csv") #ENSGID, CHR, START, STOP, GENE

df1 = pd.read_csv("Edges.csv")
df1 = df1.drop(df1.columns[0],axis=1)
df2 = pd.read_csv("EdgeWeight_transformed.csv")
df2 = df2.drop(df2.columns[0],axis=1)
df3 = pd.read_csv("EdgeTS.csv")
df3 = df3.drop(df3.columns[0],axis=1)


af = df2 * df3 #only specific edge weight

#for each of 38 tissues
for i in range(38):
	af2 = af[af.columns[i]] #TS edge weight for a given tissue
	tf = pd.concat([df1, af2], axis=1) #TF, Gene, Prior, Edge Weight
	tf = tf.drop('Prior',  axis=1)
	tissue_name = tf.columns[2]
	tf.columns = ['GENE','ENSGID','Weight'] #GENE = TF 
	bf = tf.merge(ref,on='GENE') #merge with reference and convert TF to ENSG ID
	bf = bf[['ENSGID_y','ENSGID_x','Weight']] #TF, GENE, WEIGHT (only protein coding)
	bf = bf[bf.Weight > 0] #remove zero weight and ENSG prefix for centrality
	bf['ENSGID_y'].replace(regex=True, inplace=True, to_replace=r'\ENSG',value=r'')
	bf['ENSGID_x'].replace(regex=True, inplace=True, to_replace=r'\ENSG',value=r'')
	os.chdir(dest)
	bf.to_csv("%s" %tissue_name, sep='\t', index=False, header=False)

	

