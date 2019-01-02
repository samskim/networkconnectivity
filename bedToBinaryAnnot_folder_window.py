#from .bed files to .annot files for LDSC
#change ("/random/directory") to working directory as needed

import pandas as pd
import numpy as np
import os
import glob
import sys

subfolder = sys.argv[1] #folder name 
window = sys.argv[2] #variable windows added to the bed file (e.g. input 100 for 100kb)
basedir = "/random/directory"
files = glob.glob("%s/bed/%s/*.sorted.bed" %(basedir,subfolder))
dest = "%s/annot/%s" %(basedir,subfolder)
os.system("mkdir %s" %dest) #make a destination folder

#convert bed to annot
for i in range(len(files)):
	file = files[i]
    os.chdir("%s%s" %(basedir, subfolder))
	index = files[i].split("%s%s/" %(basedir,subfolder))[1][:-11] #remove extension 
    os.chdir(dest)
    #k = chr
    for k in range(1, 23): #if annot is not found in the folder
		if (os.path.isfile("%s.%s.annot.gz" %(index, k)))==False: 
			os.system("python /random/directory/bedToAnnot.py --bedfile-single /random/directory/bed/%s/%s.extend.%skb.bed --bfile /random/directory/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC.%s --annotfile-slim %s.%s.annot.gz" %(subfolder,index,window,k,index,k))
