import pandas, sys, numpy, multiprocessing

#remove random edges from 10 to 90% 
#input: gene network (gene1, gene2, weight)
#output: gene network with missing edges (same format)
global network_f, percentage, out_folder

network_f = sys.argv[1]
percentage = int(sys.argv[2])
i = int(sys.argv[3])
out_folder = sys.argv[4]

def subset(i, network_file, perc, out):
    network = pandas.read_table(network_f, header = None, sep = '\t').values
    idx = numpy.random.choice(network.shape[0], size = int(numpy.rint(len(network)*percentage/100.0)), replace = False)
    network_out = network[idx, :]

    numpy.savetxt(out_folder + network_f.split('/')[-1].split('.')[0] + '_{0}perc_perm{1}.txt.gz'.format(percentage, i), network_out, fmt = ['%.1f','%.1f','%8f'], delimiter='\t')

subset(i, network_f, percentage, out_folder)
