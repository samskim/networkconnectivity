#compute consensus network by taking intersection of given list of gene networks
#gene networks formatted as gene ID 1, gene ID 2, weight
#max, median, mean edge weights could be taken (mean recommended)
#if more than two networks are provided, taking consensus of all networks provided

library('igraph')

#networks = c('./inweb_dedup_edge_list_connected.txt.gz','./thyroid_gland_top.gz')
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
    stop("first arguement is a list of networks separated by commas, second is out file prefix", call.=FALSE)
}

networks = strsplit(args[1],",")[[1]]
out = args[2]

network_file = networks[1]
el=read.table(network_file, header = FALSE, sep = '\t')
el[,1]=as.character(el[,1]) #Because the vertex IDs in this dataset are numbers, we make sure igraph knows these should be treated as characters. Otherwise, it'll create problems (see page on data import)
el[,2]=as.character(el[,2])
el=as.matrix(el)
g=graph.edgelist(el[,1:2], directed = FALSE) #We first greate a network from the first two columns, which has the list of vertices

#We then add the edge weights to this network by assigning an edge attribute called 'weight'.

E(g)$weight=(as.numeric(el[,3]) - min(as.numeric(el[,3])))/(max(as.numeric(el[,3])) - min(as.numeric(el[,3])))# Normalize for Sonawane networks
inter = g

for (network2 in tail(networks, -1)){
    network_file = network2
    el=read.table(network_file, header = FALSE, sep = '\t')
    el[,1]=as.character(el[,1]) #Because the vertex IDs in this dataset are numbers, we make sure igraph knows these should be treated as characters. Otherwise, it'll create problems (see page on data import)
    el[,2]=as.character(el[,2])
    el=as.matrix(el)
    g=graph.edgelist(el[,1:2], directed = FALSE) #We first greate a network from the first two columns, which has the list of vertices

    E(g)$weight_attr=(as.numeric(el[,3]) - min(as.numeric(el[,3])))/(max(as.numeric(el[,3])) - min(as.numeric(el[,3])))
    inter = intersection(inter,g)
}

inter_edge_list = as_long_data_frame(inter)

if (ncol(inter_edge_list) == 6) {
names(inter_edge_list) = c('from','to', 'weight_1','weight_2','gene_1','gene_2')

inter_edge_list$max = apply(inter_edge_list[,c("weight_1","weight_2")],1,max)
inter_edge_list$mean = apply(inter_edge_list[,c("weight_1","weight_2")],1,mean)
inter_edge_list$median = apply(inter_edge_list[,c("weight_1","weight_2")],1,median)
inter_edge_list = apply(inter_edge_list, 2, as.numeric)

inter_edge_list = inter_edge_list[,c('gene_1','gene_2','max','median','mean')]
} else if (ncol(inter_edge_list) == 7) {

names(inter_edge_list) = c('from','to', 'weight_1','weight_2','weight_3','gene_1','gene_2')

inter_edge_list$max = apply(inter_edge_list[,c("weight_1","weight_2","weight_3")],1,max)
inter_edge_list$mean = apply(inter_edge_list[,c("weight_1","weight_2","weight_3")],1,mean)
inter_edge_list$median = apply(inter_edge_list[,c("weight_1","weight_2","weight_3")],1,median)
inter_edge_list = apply(inter_edge_list, 2, as.numeric)

inter_edge_list = inter_edge_list[,c('gene_1','gene_2','max','median','mean')]

} else if (ncol(inter_edge_list) == 8) {

names(inter_edge_list) = c('from','to', 'weight_1','weight_2','weight_3','weight_4','gene_1','gene_2')

inter_edge_list$max = apply(inter_edge_list[,c("weight_1","weight_2","weight_3","weight_4")],1,max)
inter_edge_list$mean = apply(inter_edge_list[,c("weight_1","weight_2","weight_3","weight_4")],1,mean)
inter_edge_list$median = apply(inter_edge_list[,c("weight_1","weight_2","weight_3","weight_4")],1,median)
inter_edge_list = apply(inter_edge_list, 2, as.numeric)

inter_edge_list = inter_edge_list[,c('gene_1','gene_2','max','median','mean')]

}

for (i in c('max','median','mean')){
    write.table(inter_edge_list[,c('gene_1','gene_2',i)], paste(out, '_', i,'.txt', sep = ''), sep = '\t', row.names = FALSE, col.names = FALSE)
}
