from sklearn import preprocessing
import pandas as pd
import numpy as np
import os
import glob
import sys

#filter results to protein coding genes
#normalize connectivity to [0, 1]
#input: ENTREZ, unnormalized weight
#output: ENTREZ, normalized_weight
dir = sys.argv[1]
ref = pd.read_csv("/random/directory/proteincoding_genes.csv")
files = glob.glob("%s/*closeness*.txt" %dir)
 
for i in range(len(files)):
	df = pd.read_csv(files[i], names=['ENTREZ','WEIGHT'], sep='\t')
	tissuename = files[i].split("%s/" %dir)[1][:-4]
	tf = df.merge(ref,on=["ENTREZ"])
	tf = tf[["ENTREZ","WEIGHT"]]
	tf.ENTREZ = tf.ENTREZ.astype(int)
	df = tf
	x = df.values #returns a numpy array
	min_max_scaler = preprocessing.MinMaxScaler()
	x_scaled = min_max_scaler.fit_transform(x)
	betweenness = pd.DataFrame(x_scaled, columns=df.columns)
	betweenness.ENTREZ = df.ENTREZ.astype(int)
	betweenness = betweenness.sort_values('WEIGHT', ascending=False)
	betweenness.to_csv("/random/directory/network/%s_normalized.csv" %tissuename, index=False)