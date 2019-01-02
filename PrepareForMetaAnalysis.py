#prepare for meta analysis 
#copy .results and .log files
#generate .sd files of annotation and copy to /random/directory/Meta/

import pandas as pd
import numpy as np
import csv, sys
import os
import glob

folder = '/random/directory/annot'
subfolder = sys.argv[1] #subfolder name
case = sys.argv[2]

#generate .sd file (standard deviation of annotation)
def GenerateAnnotSD():
		os.chdir("%s/%s" %(folder,subfolder))
		if os.path.isfile("%s.sd" %name) == False:
			tf = pd.read_csv("%s.1.annot.gz" %(name), sep='\t') #read 1.annot
			for j in range(2,23):
 				temp = pd.read_csv("%s.%s.annot.gz" %(name, j), sep='\t') #read annot
 				tf = tf.append(temp) #appended df for chr1 to chr22
			tf.to_csv("%s_annot_merged.csv" %name, index=False, header=False) #save merged annot per geneset
			stan_dev = np.std(tf.values) #if sd file does not exist, make one and write to that file
			if os.path.isfile("%s.sd" %name) == False:
				os.system("touch %s.sd" %name) #write name.sd
				file = open("%s.sd" %name, 'w')
				file.write(str(stan_dev))
				file.close()
                                os.system("rm %s_annot_merged.csv" %name)
		
filetypes = glob.glob("/random/directory/annot/%s/*.1.annot.gz" %subfolder)
for i in range(len(filetypes)):
	filetype = filetypes[i].split("/random/directory/annot/%s/" %subfolder)[1][:-11]
        files = glob.glob("%s/%s/*%s*.22.annot.gz" %(folder,subfolder,filetype))
	for k in range(len(files)):
                        name = files[k].split("%s/%s/" %(folder,subfolder))[1][:-12]
			GenerateAnnotSD() #generate name.sd files for different annotations
                        for case in range(int(case), int(case)+1):
                                        dest = '/random/directory/Meta/%s_%s' %(name, case)
                                        os.system("mkdir %s" %dest)
                                        os.chdir("%s/%s" %(folder,subfolder)) #go there and copy needed files for meta analysis
                                        os.system("cp *_%s_%s.results %s" %(name, case, dest))
                                        os.system("cp *_%s_%s.log %s" %(name, case, dest))
                                        os.system("cp %s.sd %s" %(name, dest))
