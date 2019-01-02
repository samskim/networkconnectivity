#make annot files from .bed for continuous-valued annotation
#bed files should include chr, start, stop, weight for their columns
#weight could be negative as well

import pandas as pd
import numpy as np
import os
from intervaltree import Interval, IntervalTree
import argparse
import sys
import glob

subfolder =sys.argv[1]
basedir = "/random/directory/bed/"
files = glob.glob("/random/directory/bed/%s/*.extend.100kb.bed" %subfolder)
dest = "/random/directory/annot/%s" %subfolder
os.system("mkdir %s" %dest)

def intersect(a, b):
    intersection = max(a[0], b[0]), min(a[1], b[1])
    if intersection[0] > intersection[1]:
        return None
    return intersection

#df1 = reference, df2 = to be filled up
#assume weight is chosen at the first one
def interval_df_intersection(df1, df2):
    #df1 tree created [start, stop, weight] etc
    tree = IntervalTree.from_tuples(zip(
            df1.START.values,
            df1.STOP.values,
            df1.drop(["START", "STOP"], axis=1).values.tolist() #save weight to list
        ))
    intersections = []
    #iterator for df2 to be filled up
    for row in df2.itertuples():
        i1 = Interval(row.START, row.STOP)
        intersections += [list(intersect(i1, i2)) + i2.data for i2 in tree[i1]]
    # Make sure the column names are in the correct order
    data_cols = list(df1.columns)
    data_cols.remove("START")
    data_cols.remove("STOP")
    return pd.DataFrame(intersections, columns=["START", "STOP"] + data_cols)


for j in range(len(files)):
	name = files[j].split("%s%s/" %(basedir,subfolder))[1][:-17]
	bedfile = "./%s.extend.100kb.bed" %name 
	for i in range(1, 23):
		os.chdir(dest)
		if (os.path.isfile("%s.%s.annot.gz" %(name,i))==False):
			os.chdir("%s%s" %(basedir, subfolder)) #bed folder
			my_df = pd.read_csv("%s.extend.100kb.bed.unmerged" %name, sep='\t', names=['CHR','START','STOP','WEIGHT']) #START, STOP, WEIGHT 
			my_df = my_df[my_df['CHR'] == 'chr%s' %i]
			my_df = my_df[['START','STOP','WEIGHT']]
			os.chdir('/random/directory/bed') #BP.1-22; base pair files (tsv format) include 'START','STOP' for list of snps for chr 1 to chr 22
			my_int = pd.read_csv("BP.%s" %i, sep='\t', usecols=[0,1]) #START (BP), STOP (BP+1)
			annotfile_duplicate = ("annot_duplicate.%s" %i)
			annotfile = ("annot_temp.%s" %i)
			annotfile_final = ("%s.%s.annot.gz" %(name,i))
			bfile = ("/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC.%s" %i)
			df_bim = pd.read_csv(bfile + '.bim',delim_whitespace=True,usecols = [0,1,2,3],
                                     names = ['CHR','SNP','CM','BP'])
                	result_df = interval_df_intersection(my_df, my_int)
                	result_df.columns = ['BP', 'STOP', 'ANNOT']
                	result_df = result_df.drop('STOP', 1) #drop BP STOP col
                	annot_max = result_df.groupby('BP').ANNOT.transform(max)
                	result_df = result_df[result_df.ANNOT == annot_max]
                	df_annot = pd.merge(df_bim,result_df,how='left',on='BP') #merge together
                	df_annot.fillna(0,inplace=True) #fill blank with 0
                	df_annot['OTHER'] = 1 - df_annot['ANNOT'] #complement
                	df_annot = df_annot[['CHR','BP','SNP','CM','ANNOT','OTHER']]
                	df_annot.columns = ['CHR', 'BP', 'SNP', 'CM', 'SAM_ANNOT', 'OTHER']
                	df_annot = df_annot.drop('OTHER', 1) #drop BP STOP col
                	df_annot = df_annot.drop_duplicates()
                	os.chdir(dest)
                        df_annot = df_annot.drop(["CHR","BP","SNP","CM"], axis = 1)
                	df_annot.to_csv(annotfile_final,sep='\t',compression='gzip',index=False)
