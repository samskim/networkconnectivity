options(echo=FALSE);
param <- commandArgs(trailingOnly=T)
folder = eval(paste(text=param[1]))
annotation = eval(paste(text=param[2])) #annot name

annot = NULL
annot_temp = NULL
reference_annots = NULL

for (chr in 1:22) {
    temp=read.table(paste("/random/directory/Data/annot/",folder,"/",annotation,".",chr,".annot.gz",sep=""),h=T)
    annot = rbind(annot,temp)
}

#annot = annot[-c(1:4)] uncomment this line if annot is not thin annot 

for (chr in 1:22){
    temp = read.table(paste("/random/directory/reference_files/1000G_EUR_Phase3/baselineLD_v2.0/baselineLD.", chr, ".annot.gz", sep = ""), h = T)
    reference_annots = rbind(reference_annots ,temp[,-c(1:5)])
}


# compute the correlation of my annotations with one another
correlation.matrix = matrix(data= NA, nrow = 1, ncol = 75) 

for (j in 1:75){
	correlation.matrix[1,j] = cor(annot, reference_annots[,j])
}	

write.table(correlation.matrix,
	file=paste("correlation.matrix.all.annotations.",annotation,".csv",sep = ""),
	sep=",", quote = F, col.names = colnames(reference_annots), row.names = colnames(annot))