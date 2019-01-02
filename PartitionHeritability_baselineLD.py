#Run partition heritability conditioning on baselineLD
#note baselineLD has multiple versions

import os
import sys
import glob

i = sys.argv[1] #sumstats file
subfolder = sys.argv[2]
file = "/random/directory/annot/%s" %subfolder  #annot folder
files = glob.glob("/random/directory/annot/%s/*.22.annot.gz" %subfolder)
os.chdir(file)

for k in range(len(files)):
  j = files[k].split(file)[1][1:-12]
  if (os.path.isfile("%s_%s_550.results" %(i, j)))==False:
    os.system("python /directory1/steven/soft/ldsc/ldsc.py --h2 /directory1/ldsc/sumstats_formatted/%s --ref-ld-chr %s/%s.,/directory1/ldsc/reference_files/1000G_EUR_Phase3/baselineLD_v2.1/baselineLD. --w-ld-chr /directory1/ldsc/reference_files/1000G_EUR_Phase3/weights/weights.hm3_noMHC. --overlap-annot --frqfile-chr /directory1/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC. --out %s_%s_550 --print-delete-vals --print-coefficients" % (i, file, j, i, j))
