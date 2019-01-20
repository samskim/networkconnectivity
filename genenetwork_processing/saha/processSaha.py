import glob
import os, sys
import pandas as pd
import numpy as np
import math

#protein coding only
os.chdir('/random/directory/BIOCART')
rf1 = pd.read_csv("ALLGENES.tsv", sep='\t', names=['CHR','START','STOP','ENST','ENTREZ','Name1'], skiprows=[0]) #CHR, START, STOP, ENST, GENE (GRCH38)
rf2 = pd.read_csv("ALLGENES.tsv", sep='\t', names=['CHR','START','STOP','ENST','ENTREZ','Name2'], skiprows=[0]) #CHR, START, STOP, ENST, GENE (GRCH38)
rf1 = rf1.drop(['ENTREZ'],axis=1)
rf2 = rf2.drop(['ENTREZ'],axis=1)
rf1 = rf1.drop_duplicates(subset='Name1', keep="last")
rf2 = rf2.drop_duplicates(subset='Name2', keep="last")

cf = pd.read_csv("ENSTtoENTREZ_proteincoding.tsv", sep='\t', usecols=[1,5])  #ENSGID, ENTREZ, CHR, START, STOP (GRCH37)
cf = cf.dropna()
cf['ENTREZ'] = cf['ENTREZ'].astype(int)
cf = cf.drop_duplicates(subset='ENTREZ', keep="last")
#cf.to_csv("ENSTtoENTREZ_unique.tsv", sep='\t', index=False)

files = glob.glob('/random/directory/Network_Saha/gtex_portal_networks/twns/*')

for i in range(len(files)):
	#os.chdir(files[i])
	df = pd.read_csv(files[i], sep='\t')
	df['WEIGHT'] = df['Edge weight'].abs()
	df = df.drop(['Edge weight'], axis = 1)
	#df['ENST_ID'] = np.nan
	df1 = df.loc[df.Name1.isin(rf1.Name1)]
	df1 = df1.merge(rf1, on=['Name1'])
	df1 = df1.drop(['Name1'],axis=1)
	
	df2 = df1.loc[df1.Name2.isin(rf2.Name2)]
	df2 = df2.merge(rf2, on=['Name2'])
	df2 = df2.drop(['Name2'],axis=1)
	df = df2
	df['ENST_x'] = df['ENST_x'].str[:15]
	df['ENST_y'] = df['ENST_y'].str[:15]
	df = df.sort_values('ENST_x')
	os.chdir('/random/directory/Network_Saha/raw_network_ENST')
	df.to_csv('%s' %files[i][59:], sep='\t', index=False, header = ['EdgeType','WEIGHT','CHR1','START1', 'STOP1','GENE1','CHR2','START2','STOP2','GENE2'])
	
	#now convert GENE1/GENE2 (ENST ID) to ENTREZ (df.tsv uses ENTREZ)
	cf.columns = ['ENTREZ','ENST_x']
	df1 = df.loc[df.ENST_x.isin(cf.ENST_x)]
	df1 = df1.merge(cf, on=['ENST_x'])
	df1 = df1.drop(['ENST_x'],axis=1)
	cf.columns = ['ENTREZ','ENST_y']
	df2 = df1.loc[df1.ENST_y.isin(cf.ENST_y)]
	df2 = df2.merge(cf, on=['ENST_y'])
	df2 = df2.drop(['ENST_y'],axis=1)		
	
	tf = pd.DataFrame({'GENE1':df2.ENTREZ_x, 'GENE2':df2.ENTREZ_y, 'WEIGHT':df2.WEIGHT})
	#make sure weight is all positive
	tf = tf.sort_values('GENE1')
	os.chdir('/random/directory/Network_Saha/raw_network')
	tf.to_csv('%s' %files[i][59:], sep='\t', index=False, header=False)
	