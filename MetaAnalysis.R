#random-effect meta analysis across traits
library("rmeta")

param = commandArgs(trailingOnly=TRUE)
file = eval(paste(text=param[1]))
filetype = eval(paste(text=param[2])) #example: 500 (extension to the .results file e.g. ~_500.results)
setwd(paste0("/random/directory/Meta/", file, "_", filetype))

#CORRECT VERSION
M       = 5961159

OutPath = paste0("/random/directory/Meta/", file , "_", filetype, "/*.results")
list_pair = Sys.glob(OutPath); #traits with log files
list_pair = tools::file_path_sans_ext(list_pair)
list_pair = gsub(paste0("/random/directory/Meta/", file, "_", filetype, "/"), "", list_pair) #file names (file.log, file.result)

enr        = NULL
enr_sd     = NULL
enrstat    = NULL
enrstat_sd = NULL
tau        = NULL
tau_sd     = NULL

for (trait in list_pair) {
##Get sd_annot by mapping to specific BSID from list_pair instance##
 
#BSID = gsub(".*_", "", trait) #just BSID names
mysd   = read.table(paste0(file, ".sd"))$V1; #read sd_annot: BSID.sd

data         = read.table(paste("",trait,".results",sep=""),h=T)    #modify the .results path
log          = read.table(paste("",trait,".log",sep=""),h=F,fill=T) #modify the .log path
h2g          = as.numeric(as.character(log[which(log$V4=="h2:"),5]))
enr          = cbind(enr   , data$Enrichment)
enr_sd       = cbind(enr_sd, data$Enrichment_std_error)
#
myenrstat    = (h2g/M)*((data$Prop._h2/data$Prop._SNPs)-(1-data$Prop._h2)/(1-data$Prop._SNPs))
myenrstat_z  = qnorm(data$Enrichment_p/2) #step2
myenrstat_sd = myenrstat/myenrstat_z #step3
enrstat      = cbind(enrstat   , myenrstat)
enrstat_sd   = cbind(enrstat_sd, myenrstat_sd)
#
tau          = cbind(tau   , M*mysd*data$Coefficient/h2g)
tau_sd       = cbind(tau_sd, M*mysd*data$Coefficient_std_error/h2g)
} #close for loop

#meta analysis begins here
enr_meta = NULL
tau_meta = NULL
for (i in 1:nrow(enr)){
test1 = meta.summaries(enr[i,],enr_sd[i,],method="random")
if (data$Prop._SNPs[i]==1) {
enr_meta = rbind(enr_meta,c(test1$summary,test1$se.summary,NA)) # case of the base annotation
} else {
test2 = meta.summaries(enrstat[i,],enrstat_sd[i,],method="random")
enr_meta = rbind(enr_meta,c(test1$summary,test1$se.summary,2*pnorm(-abs(test2$summary/test2$se.summary))))
}
test = meta.summaries(tau[i,],tau_sd[i,],method="random")
tau_meta = rbind(tau_meta,c(test$summary,test$se.summary,2*pnorm(-abs(test$summary/test$se.summary))))
}

out = cbind(data[,1:2],enr_meta,tau_meta)

colnames(out)[3:8] = c("Enrichment","Enrichment_std_error","Enrichment_pval","Coefficient","Coefficient_std_error","Coefficient_pval")

out[1,] #print

write.csv(out[1,],"meta.csv")

length(list_pair)
