load("GTEx_PANDA_tissues.RData")
ls()
net_transformed <- log(exp(net)+1)
write.csv(net_transformed, file="EdgeWeight_transformed.csv")