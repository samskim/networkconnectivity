from sklearn import preprocessing
import pandas as pd
import numpy as np
import os
import glob
import sys
import math  

#compute the enrichment, enrichment s.e. given the ENTREZ ID, closeness matrix
#input files: dataframe with two columns ["ENTREZ", "WEIGHT"]
#output files: for 10 decile bins, enrichment matrix, enrichment s.e. matrix, summary matrix 
#summary matrix contains for each genes, '1' if found in > 16 gene sets (e.g. high pLI) analyzed
#description of gene sets analyzed is provided in the paper

subfolder = sys.argv[1]
files = glob.glob("/random/directory/network/%s/*normalized*" %subfolder)

def safe_div(x,y):
    if y==0: return 0
    return x/y
    
for i in range(len(files)):
	filename = files[i].split("/random/directory/network/%s/" %subfolder)[1][:-4]
	df = pd.read_csv(files[i])
	ref = pd.read_csv("/random/directory/geneset_ENTREZ.csv")
	for k in range(len(ref.columns)):
		df[ref.columns[k]] = df.ENTREZ.isin(ref[ref.columns[k]]).astype(int)
		df.to_csv("/random/directory/network/constraint_analysis_RESULT/%s_constraintgenes.csv" %filename,index=False)

#10th decile bin
	a =  df[df.WEIGHT > df.WEIGHT.quantile(0.9)]
	af = a.sum()
	af = af.to_frame("top10%")
	af = af.T
	j = 2
	#1-9th decile bin
	for z in (10, 20, 30, 40, 50, 60, 70, 80, 90):
		upper_percentile = (100.0 - z)/100 #0.9
		lower_percentile = upper_percentile - 0.1
		tf = df[(df.WEIGHT < df.WEIGHT.quantile(upper_percentile)) & (df.WEIGHT >= df.WEIGHT.quantile(lower_percentile))]
		af2 = tf.sum()
		af2 = af2.to_frame("top%s0" %j)
		af2 = af2.T
		af = af.append(af2)
		j += 1	
	af = af.drop("ENTREZ",axis=1)
	af = af.drop("WEIGHT",axis=1)
	af.to_csv("/random/directory/network/constraint_analysis_RESULT/%s_constraintgenes_10bins_summary.csv" %filename)
	#calculate enrichment
	#for example:
	#enrichment = (number of ExAC genes captured in bin 1)/(number of genes in bin 1) * (number of all genes)/(number of ExAC genes found in Greene network) = 866/1871*18712/2919 = 2.97
	#enrichment s.e. = sqrt((866/1871)(1-866/1871)/1871)*18712/2919 = 0.0736
	df = af #read bin10 summary data
	enrichment_matrix={}
	enrichmentse_matrix={}
	#number of genes of interest captured in bin 1
	for g in range(0, len(df.columns)):
		geneofinterest = df.columns[g]
		if g == 0:
			geneofinterest = 'Drug_bank' #force
		enrichment_list = [] #for each gene of interest
		enrichmentse_list = []
		for u in range(0,10):
			num1 = df[geneofinterest][u] #number of genes of interest in bin k
			num2 = df["Universe"][u] #number of genes in bin k
			num3 = df.sum()["Universe"] #all genes in the network
			num4 = df.sum()[geneofinterest]
			enrichment = safe_div(num1,num2) * safe_div(num3,num4)
			enrichment_se = math.sqrt(safe_div(num1,num2)*(1-safe_div(num1,num2))/num2)*safe_div(num3,num4)
			enrichment_list.append(enrichment)
			enrichmentse_list.append(enrichment_se)
		enrichment_matrix[geneofinterest]=enrichment_list
		enrichmentse_matrix[geneofinterest]=enrichmentse_list
	enrichment_df = pd.DataFrame(enrichment_matrix) #save enrichment df
	enrichmentse_df = pd.DataFrame(enrichmentse_matrix) #save enrichment s.e. df
	enrichment_df.to_csv("/random/directory/network/constraint_analysis_RESULT/%s_enrichment.csv" %filename)
	enrichmentse_df.to_csv("/random/directory/network/constraint_analysis_RESULT/%s_enrichmentse.csv" %filename)
