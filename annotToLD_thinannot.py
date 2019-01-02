#convert .annot to .ldscore.gz for LDSC
#change ("/random/directory") to working directory as needed

import os
import glob
import sys

path = sys.argv[1] #folder name 
files = glob.glob("/random/directory/annot/%s/*.22.annot.gz" %path)

for k in range(len(files)):
    filename = files[k].split('/random/directory/annot/%s' %path)[1][1:-12]
    os.chdir("/random/directory/annot/%s" %path)
    index = filename
    #remove thin annot flag if annot file contains columns other than annotation values (e.g. # cols > 1)
    for j in range(1,23):
        if (os.path.isfile("%s.%s.l2.ldscore.gz" %(index,j)))==False: 
            os.system("python /directory1/ldsc.py --l2\
                     --bfile /directory2/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC.%s\
                     --ld-wind-cm 1 --annot %s.%s.annot.gz --thin-annot  --out %s.%s --print-snps /n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/list.txt" %(j, index, j, index, j))
