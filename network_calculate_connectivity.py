#given the input of gzipped gene network (gene ID 1, gene ID 2, weight) with '\t' delimiter, 
#compute different network connectivity metrics
#if list of genes is provided (as .tsv), it searches for neighbors in a given network and uses a subnetwork

import pandas
import numpy
import glob
import scipy
from sklearn import cluster
import graph_tool.all as gt
from matplotlib import pyplot as plt
plt.switch_backend('agg')
import sys

###################### Parsing of Datasets and Preprocessing ######################
def neighbor_graph(GIANT_dataset, core_gene_array = None, core_only = False, peripheral_connections = False, output_type = 'pandas', same_cost = False):
	'''Creates an pandas adjacency matrix of core genes, neighboring peripheral genes, and posterior probability.
		Takes in a tissue network and a set of core genes (OPTIONAL) --> generates a pandas DataFrame identifying
		the set of neighboring genes and posterior probability (-1 if an edge does not exist)

	Args:
		- GIANT_dataset (gzipped, tab-delimited text file): Top Edges file of a tissue from GIANT
		- core_gene_array (list; default = None): an array of core genes;
												if None, all genes in the tissue network are assumed to be relevant
		- core_only: if True, only core-core connections will be evaluated. core_only cannot be combined with peripheral_connections == True.
		- peripheral_connections (default = False): if True, peripheral-peripheral interactions will also be included in the graph
		- return_type (default = 'pandas'): valid entries are 'dict' (returns a nested dictionary) or 'pandas' (DataFrame)
	Returns:
		If core_gene_array == None, returns an adjacency list (pandas Dataframe) with columns:
			['source', 'target', 'weight', 'cost']
		If return_type == 'pandas', returns a pandas Dataframe where columns are the source (core genes) and
			rows are the target (peripheral genes) and values are the posterior proability.
		If return_type == 'dict', returns a nested dictionary where keys are core genes and
			values are a dictionary with peripheral genes (must output as 'dict' if building graph tool graph)
	'''

	#os.chdir('/n/groups/price/sam/Network_Greene/raw_network') #directory that has networks

	#if no core genes are presented, return the entire graph
	if core_gene_array == None:
		print 'No core genes specified; building adjacency list for all genes in the tissue network: ' + str(GIANT_dataset)[:-3]
		tissue_network = pandas.read_table(GIANT_dataset,header = None,
									   names = ['source', 'target', 'weight'])
		if same_cost == True:
			print 'Edge cost is the same as the weights'
			tissue_network['cost'] = tissue_network['weight']
		else:
			tissue_network['cost'] = 1/tissue_network['weight']
		return tissue_network

	else:
		print 'Core genes specified; building adjacency matrix using tissue network: {0}'.format(str(GIANT_dataset)[:-3])
		tissue_network = pandas.read_table(GIANT_dataset,header = None,
										   names = ['gene_id_1', 'gene_id_2', 'posterior_prob'], index_col = [0,1])

		#identify neighboring genes and convert gene_id_2 column and posterior_prob to numpy arrays
		if core_only == True:
			peripheral_connections = False
			core_gene_neighbors = tissue_network.query('gene_id_1 in {0} and gene_id_2 in {0}'.format(list(core_gene_array)))
		else:
			core_gene_neighbors = tissue_network.query('gene_id_1 in {0} or gene_id_2 in {0}'.format(list(core_gene_array)))
		core_peripheral_dict = dict(zip(core_gene_neighbors.index, numpy.ndarray.flatten(core_gene_neighbors.values)))
		if peripheral_connections:
			peripheral_genes = set(core_gene_neighbors.index.get_level_values(1)).union(set(core_gene_neighbors.index.get_level_values(0)))-set(list(core_gene_array))
			peripheral_neighbors = tissue_network.query('gene_id_1 in {0} and gene_id_2 in {0}'.format(list(peripheral_genes)))
			core_gene_neighbors = pandas.concat([core_gene_neighbors,peripheral_neighbors])

		neighbor_series = core_gene_neighbors.groupby(level=0).apply(lambda ind: dict(zip(ind.xs(ind.name).index, \
														numpy.ndarray.flatten(ind.xs(ind.name).values))))
		core_peripheral_dict = dict(zip(neighbor_series.index, numpy.ndarray.flatten(neighbor_series.values)))

		#must output as 'dict' if building graph tool graph
		if output_type == 'dict':
			return core_peripheral_dict

		else:
			core_peripheral = pandas.DataFrame.from_dict(core_peripheral_dict, orient = 'columns')
			core_peripheral.fillna(-1, inplace = True)
			core_peripheral.to_csv(str(GIANT_dataset)[:-3] + '_core_peripheral_graph.csv', sep = '\t')
			return core_peripheral

def build_adj_list(adj_list, core_gene, peripheral_gene, probability):
	'''Helper function for numpy.vectorize for the building of an adjacency list.'''
	adj_list.append((core_gene,peripheral_gene,probability, 1.0/probability))

def convert_to_adj_list(neighbors_adj_matrix):
	'''Converting the output of neighbor graph to an adjacency list matrix that can be accepted by Graph-Tool or NetworkX.
	Weights of the graph are the posterior probability. The cost of an edgeis the inverse of the weight.
	Args:
	   neighbors_adj_matrix: an adjacency matrix with columns as core genes and rows as peripheral genes
	Returns:
		neighbors_adj_list: pandas DataFrame with four columns: ['source', 'target', 'weight', 'cost']
	'''
	print "Converting to adjacency list"
	adj_list = []
	for core_gene in neighbors_adj_matrix.keys():
		vectorized = numpy.vectorize(build_adj_list, excluded = [0,1], otypes=[list])
		vectorized(adj_list, core_gene, neighbors_adj_matrix[core_gene].keys(), neighbors_adj_matrix[core_gene].values())
	return pandas.DataFrame.from_records(adj_list, columns = ['source', 'target', 'weight', 'cost'])


###################### Gene Network Graph ######################

class Gene_interaction_graph(object):
	'''Object representation of a core gene - peripheral gene interaction network'''
	def __init__(self, pandas_adj_list, core_genes = None, name = None):
		'''Initialize the Gene_interaction_graph object.
		Args:
			pandas_adj_list: pandas adjacency list with columns: ['source', 'target', 'weight', 'cost']
			core_genes (OPTIONAL): list of core genes
			name (OPTIONA): name for the graph for outputting
		Returns:
			self.graph: a graph object of the gene interaction network
		'''
		super(Gene_interaction_graph, self).__init__()

		#build graph tool graph
		self.graph = gt.Graph(directed = None)
		self.edge_weights = self.graph.new_edge_property('double')
		self.edge_cost = self.graph.new_edge_property('double')
		self.prop = self.graph.add_edge_list(pandas_adj_list.iloc[:,:4].values, hashed=True,
							   string_vals=True, eprops = [self.edge_weights, self.edge_cost])

		self.adj_list = pandas_adj_list
		self.adj_matrix = gt.adjacency(self.graph, weight = self.edge_weights)

		self.core_genes = core_genes
		self.name = name

	def get_core_peripheral_genes(self):
		'''Return a pandas Series with genes and annotations of core or peripheral'''
		all_genes = self.prop.get_array()
		core_gene_dict = {gene:1.0 for gene in self.core_genes}

		gene_annot = pandas.DataFrame(all_genes)
		gene_annot.columns = ['genes']
		gene_annot['annotation'] = gene_annot['genes'].apply(core_gene_dict.get).fillna(0)
		return gene_annot.set_index('genes')['annotation']

	def betweenness_centrality(self):
		'''Calculates a pandas Series with each gene node and the corresponding normalized betweeness centrality value
		'''

		#get betweeneess centrality for both nodes and edges
		print 'Calculating betweenness centrality'
		vp, ep = gt.betweenness(self.graph, weight = self.edge_cost) #edge cost due to shortest path method

		node_centrality = dict(zip(self.prop.get_array(), vp.get_array()))

		print 'Building betweenness centrality matrix'
		node_centrality_series = pandas.Series(node_centrality).rename('betweenness_centrality')
		if self.name != None:
			node_centrality_series.to_csv(self.name + '_centrality.txt', sep = '\t')
		return node_centrality_series

	def eigenvector_centrality(self, get_top = 1.0):
		'''Calculate the eigenvalue centrality of a graph and identifies the genes with highest eigenvalue centrality
		Args:
			- get_top: porportion of genes with highest eigenvalue centrality to return [0.0 <= get_top <= 1.0]

		Returns:
			- node_eigen_centrality: pandas Series of each entrez gene and its eigenvalue centrality (index, value = Gene ID, eigenvalue centrality)
			- top_nodes: list of Gene ID in decreasing eigenvalue centrality order
		'''
		if not (get_top <= 1.0 and get_top >= 0.0):
			print 'User defined get_top is out of valid range [0.0 <= get_top <= 1.0]'
			return None, None
		else:
			#build graph tool graph
			print 'Calculating eigenvector centrality'
			#calculate eigenvector centrality
			top_eigenvalue, eigenvector = gt.eigenvector(self.graph, weight= self.edge_weights) #edge cost due to represent importance/higher probability

			print 'Building eigenvector centrality matrix'
			node_eigen_centrality = dict(zip(self.prop.get_array(), eigenvector.get_array()))

			#return genes in order of descending eigenvector centrality
			sorted_genes = sorted(node_eigen_centrality, key=node_eigen_centrality.get)[::-1]
			top_nodes = sorted_genes[:int(round(len(sorted_genes)*get_top,0))]

			return pandas.Series(node_eigen_centrality).rename('eigen_centrality'), top_nodes

	def closeness_centrality(self):
		"Calculates a pandas Series mapping each gene vertex to a normalized closeness centrality"

		print 'Calculating closeness centrality'
		closeness = gt.closeness(self.graph, weight = self.edge_cost) #edge cost due to shortest path method
		node_closeness = dict(zip(self.prop.get_array(), closeness.get_array()))

		print 'Building closeness centrality matrix'
		return pandas.Series(node_closeness).rename('closeness_centrality')

	def pagerank(self):
		"Calculates a pandas Series holding the PageRank value for each gene"

		print 'Calculating pagerank'
		pagerank = gt.pagerank(self.graph, weight = self.edge_weights) #edge weight due to importance
		node_pagerank = dict(zip(self.prop.get_array(), pagerank.get_array()))

		print 'Building PageRank matrix'
		return pandas.Series(node_pagerank).rename('PageRank')

	def degree_centrality(self):
		'''Calculate the normalized degree centrality of the weighted graph.
		For a weighted graph, the degree centrality is the sum of the edge weights
		connected to a vertex divided by the total number of possible edges (V-1)

		Returns:
			- degree centrality: pandas Series with the degree centrality of each node
		'''

		print 'Calculating degree centrality'
		degree_centrality = {}
		total_possible_edges = len(self.prop.get_array())-1

		for index in self.graph.vertices():
			degree_centrality[self.prop[index]] = (index.out_degree() + index.in_degree())/float(total_possible_edges)

		return pandas.Series(degree_centrality).rename('Degree Centrality')

	def sum_edge_scores(self):
		'''Calculates the sum of edge weights of vertices in the graph.
		Returns:
			- sum_edge: a pandas Series with the value of the sum of all edge weights for each node
		'''
		print 'Calculating sum of edge scores'
		sum_edge_score = {}

		for index in self.graph.vertices():
			sum_edge_score[self.prop[index]] = index.out_degree(weight = self.edge_weights) + index.in_degree(weight = self.edge_weights)

		return pandas.Series(sum_edge_score).rename('Sum Edge Weight')

        def vertex_degree(self):
                '''Calculates the degree of a vertices in the graph.
                Returns:
                        - vertex_degree: a pandas Series with the value of the sum of all edge weights for each node
                '''
                print 'Calculating degree of vertices'
                degree = {}

                for index in self.graph.vertices():
                        degree[self.prop[index]] = index.out_degree() + index.in_degree()

                return pandas.Series(degree).rename('Degree')

	def mean_edge_scores(self):
		'''Calculates the mean of edge weights of vertices in the graph.
		Returns:
			- mean_edge: a pandas Series with the value of the mean of all edge weights for each node
		'''
		print 'Calculating mean of edge scores'
		mean_edge_score = {}

		for index in self.graph.vertices():
			mean_edge_score[self.prop[index]] = (index.out_degree(weight = self.edge_weights) + index.in_degree(weight = self.edge_weights))/(index.out_degree() + index.in_degree())

		return pandas.Series(mean_edge_score).rename('Mean Edge Weight')

	def max_edge_scores(self):
		'''Calculates the maximum edge weights of all vertices in the graph.
		Returns:
			- max_edge: a pandas Series with the value of the maximum edge weight for each node
		'''
		print 'Calculating maximum of edge scores'
		max_dict = {}
		for index in self.graph.vertices():
			max_dict[self.prop[index]] = max([self.edge_weights[edge] for edge in index.all_edges()])

		return pandas.Series(max_dict).rename('Maximum Edge Weight')

	def local_cluster_coeff(self):
		'''Find the local clustering coefficients for all vertices.
		Returns:
			- cluster_coeff: a pandas Series with the local cluster coefficients for each node
		'''
		print 'Calculating local cluster coefficient'
		local_cluster_coeff = gt.local_clustering(self.graph, undirected = True)
		node_cluster_coeff = dict(zip(self.prop.get_array(), local_cluster_coeff.get_array()))

		return pandas.Series(node_cluster_coeff).rename('Local Clustering')

	def shortest_path_distribution(self, tissue):
		'''Find the distribution of the shortest paths
		Returns:
			- saves a histogram plot as TISSUE+_distance_distribution.png'''
		print 'Finding distribution of shortest paths'
		counts, bins = gt.distance_histogram(self.graph, weight=self.edge_cost)
		print 'Building histogram of shortest paths distribution'
		plt.figure(figsize = (10,5))
		plt.bar((bins[1:] + bins[:-1]) * .5, counts, width=(bins[1] - bins[0]))
		plt.xlim((0,1000))
		plt.title(tissue)
		plt.savefig(tissue+'_distance_distribution.png')
		plt.close()
		return None

	def spectral_clustering(self, k):
		'''Executes spectral clustering on the graph and maps each node to its respective cluster
		Args:
			- k: number of clusters to group the nodes into

		Returns:
			- node_labels: a dictionary mapping each gene to its respective cluster (key, value = Gene ID, cluster)
		'''
		labels = cluster.spectral_clustering(self.adj_matrix, n_clusters= k, eigen_solver='arpack', n_init=10, assign_labels='kmeans')
		return dict(zip(self.prop.get_array(), labels))

	def laplacian_eigendecomposition(self):
		'''Generate the normalized laplacian matrix and calculate its eigenvalues and eigenvectors
		Returns:
			- ranked_w_norm: list of eigenvalues for the normalized laplacian in decreasing order
			- ranked_vr_norm: list of eigenvectors for the normalized laplacian in the same order as their respective eigenvalues
		'''
		cs_graph = scipy.sparse.csgraph.csgraph_from_dense(self.adj_matrix.todense())
		L_norm = scipy.sparse.csgraph.laplacian(cs_graph, normed = True)
		w_norm, vr_norm = scipy.linalg.eig(L_norm.todense(), overwrite_a=True)

		#rank the normalized eigenvalue
		idx = w_norm.argsort()[::-1]
		ranked_w_norm = w_norm[idx]
		ranked_vr_norm = vr_norm[:,idx]

		return ranked_w_norm, ranked_vr_norm

	def hierarchical_clustering_heatmap(self, nodes_to_community_file, labelings_connection_file):
		node_cluster_file = nodes_to_community_file
		node_index_file = labelings_connection_file

		#map Entrez Gene ID with corresponding cluster
		node_index_map = pandas.read_table(node_index_file, sep = ' ', header = None, index_col = 1)
		node_cluster = pandas.read_table(node_cluster_file, sep = ' ', header = None, index_col = 0)
		node_cluster_map = node_index_map.join(node_cluster)
		node_cluster_dict = dict(zip(node_cluster_map.iloc[:,0].values, node_cluster_map.iloc[:,1].values)) #keys = Entrez ID, values = cluster label

		print 'Saving gene-cluster map at: {0}'.format('gene_id_cluster_map.csv')
		node_cluster_map.to_csv('gene_id_cluster_map.csv', index = None, header = None)

		self.heatmap(node_labels = node_cluster_dict)

		return node_cluster_dict

	def heatmap(self, node_labels = None):
		if node_labels == None:
			print 'Saving heatmap at: {0}'.format('heatmap.png')
			plt.imshow(self.adj_matrix.todense(), cmap='hot', interpolation='nearest')
			plt.savefig('heatmap.png')
			plt.close()
			return None

		elif type(node_labels) == dict:
			#map the node indices in graph tool graph with the Gene ID
			property_map = dict(zip(range(len(self.prop.get_array())), self.prop.get_array()))
			property_map = pandas.DataFrame.from_dict(property_map, orient = 'index')
			property_map.columns = ['node_id']

			#sort the nodes in order based on the cluster they're suppose to be in and reset the index for new adjacency matrix
			new_index = pandas.DataFrame.from_dict(node_labels, orient = 'index').sort_values(0).reset_index()
			new_index_dict = dict(zip(new_index.iloc[:,0].values, new_index.index))
			new_property_map = property_map.join(pandas.DataFrame.from_dict(new_index_dict, orient = 'index'), on = ['node_id'])

			#generate new adjacency
			new_adj_index = self.graph.new_vertex_property('int', vals = new_property_map.loc[:,0].values)
			adj_matrix_reordered = gt.adjacency(self.graph, weight = self.edge_weights, index = new_adj_index)

			#generate heatmap
			print 'Saving ordered heatmap at: {0}'.format('heatmap_'+str(new_index.iloc[-1,1]+1)+'_clusters.png')
			plt.imshow(adj_matrix_reordered.todense(), cmap='hot', interpolation='nearest')
			plt.savefig('heatmap_'+str(new_index.iloc[-1,1]+1)+'_clusters.png')
			plt.close()
			return None

		else:
			print 'node_labels arguement must be as a dictionary'
			return None


	def show(self, node_labels = None):
		'''graphs the gene interaction graph
		Args:
			node_labels (OPTIONAL): user defined color for the core gene nodes (other than red), default color is blue
		Returns:
			a image of the graph with core gene nodes as the defined color (or default blue) and the peripheral gene nodes as red
		'''
		print 'Drawing graph, output saved as: gene_interation_graph.png'

		#set color of vertices based on core vs peripheral genes
		if node_labels == None:
			coloring = dict.fromkeys(self.core_genes, 0.85)
			max_value = 0.85
		elif type(node_labels) == dict:
			coloring = node_labels
			max_value = float(max(coloring.values()))
		else:
			print 'node_labels arguement must be as a dictionary'
			return None

		if self.core_genes != None:
			vertex_color = self.graph.new_vertex_property('double')
			vertex_size = self.graph.new_vertex_property('double')

			for index in range(self.graph.num_vertices()):
				if self.prop[index] in self.core_genes:
					vertex_color[index] = coloring.get(self.prop[index], 0.85)/max_value
					vertex_size[index] = 0.1
				else:
					vertex_color[index] = coloring.get(self.prop[index], 0.15)/max_value
					vertex_size[index] = 0.05

			gt.graphviz_draw(self.graph, size = (50,50), layout = 'sfdp', penwidth = 0.1,
						 vcolor = vertex_color, vsize = vertex_size, output="gene_interation_graph.png")
			return None

		else:
			gt.graphviz_draw(self.graph, size = (50,50), layout = 'sfdp', penwidth = 0.1, vsize = 0.025,
							 output="gene_interation_graph.png")
			return None


###################### Past Impementations - For Record Keeping ######################

def max_probability(dictionary, gene_id, posterior_probability):
	'''Helper function. Stores the posterior probability of a neighboring gene in a dictionary.
		If the neighboring gene is already in dictionary, calculate and store the max probability
	Args:
		dictionary: dictionary to store the set of genes
		gene_id: the unique ID of a neighboring gene
		posterior_probability: the probability between a core gene and the neighboring gene
	Returns:
		dictionary containing neighboring genes as keys and the maximum probability as values
	'''
	dictionary[gene_id] = max(dictionary.get(gene_id), posterior_probability)
	return len(dictionary)

def neighbor_search(core_gene_array, GIANT_dataset):
	'''Indentifies the neighboring genes with the highest probability from a core set of genes.
		Takes in a set of core genes and generates a pandas DataFrame identifying
		the set of neighboring genes and the highest posterior probability
	Args:
		core_gene_array (list): an array of core genes
		GIANT_dataset (gzipped, tab-delimited text file): Top Edges file of a tissue from GIANT
	Returns:
		pandas DataFrame with two columns, the first is the set of neighboring genes and the second is the maximum posterior probability
	'''
	#import GIANT .gz files
	#add headers [gene_id_1, gene_id_2, posterior_prob]
	#generate multi-index for the index where layer 0 of the index is gene_id_1 and layer 1 of index is gene_id_2
	os.chdir('/n/groups/price/sam/Network_Greene/raw_network') #directory that has networks
	tissue_network = pandas.read_table(GIANT_dataset,header = None, names = ['gene_id_1', 'gene_id_2', 'posterior_prob'], index_col = [0,1])
	#dictionary to store neighbors and probability: neighbor_gene_id -> posterior probability
	neighbor_genes_dict = {}
	#identify all neighboring genes and store the edge with the maximum probability in neighbor_genes_dict
	for genes in core_gene_array:
		#identify neighboring genes and convert gene_id_2 column and posterior_prob to numpy arrays
		neighbor_ids = tissue_network.loc[genes].reset_index().loc[:,'gene_id_2'].values
		neighbor_prob = tissue_network.loc[genes].reset_index().loc[:,'posterior_prob'].values
		#use numpy vectorize for performance
		numpy.vectorize(max_probability, otypes=[dict], excluded = [0])(neighbor_genes_dict, neighbor_ids, neighbor_prob)
	#pandas DataFrame from dictionary with neighboring genes as first column and probability as second column
	neighbor_genes = pandas.DataFrame.from_dict(neighbor_genes_dict,orient='index').reset_index()
	neighbor_genes.columns = ['gene_id', 'posterior_prob']
	#save as csv, excluding header and index
	neighbor_genes.to_csv(str(GIANT_dataset)[:-3] + '_neighbors.csv', header = False, index = False)
	return neighbor_genes


###################### Main ######################

def read(file):
	data = numpy.array([], dtype=int);
	for line in open(file,'r'):
		data = numpy.append(data, numpy.array([int(ID) for ID in line.split(',')]))
	return list(data)

def parse_core_gene_file(pathway_file):
	return list(set(numpy.ndarray.flatten(pandas.read_table(pathway_file, usecols=['ENTREZ']).values)))

def main():
	if len(sys.argv) == 3:
		tissue = sys.argv[1]
		core_gene_file = sys.argv[2]#.split('.')[0]
		core_genes = parse_core_gene_file(core_gene_file + '.tsv')
		print 'Core Genes for pathway {0} are: {1}'.format(core_gene_file, str(core_genes))
		graph = neighbor_graph(tissue, core_gene_array = core_genes, peripheral_connections = True, output_type = 'dict')

		if len(graph) == 0:
			with open('error.log', 'a') as error:
				error.write('{0} {1} \n'.format(core_gene_file, tissue))
			error.close()
			print "ERROR: {0} not in {1}; see error.log".format(str(core_genes), tissue)
			return
		graph_adj_list = convert_to_adj_list(graph)
		new_graph = Gene_interaction_graph(graph_adj_list, core_genes = core_genes) #To specify core genes, provide a list to the core_genes parameter)
		tissue = tissue.split('/')[-1].split('.')[0]
		core_gene_file = './inweb/18k_pathway_results/' + core_gene_file.split('/')[-1].split('.')[0]

		# graph_genes_annot = new_graph.get_core_peripheral_genes()
		# graph_genes_annot.to_csv(core_gene_file + '_' + tissue + "_gene_annot.txt", sep = '\t')
		graph_eigenvector = new_graph.eigenvector_centrality()[0]
		graph_eigenvector.to_csv(core_gene_file + '_' + tissue + "_pagerank.txt", sep = '\t')
		graph_pagerank = new_graph.pagerank()
		graph_pagerank.to_csv(core_gene_file + '_' + tissue + "_eigenvector_centrality.txt", sep = '\t')
		graph_degree = new_graph.degree_centrality()
		graph_degree.to_csv(core_gene_file + '_' + tissue + "_degree_centrality.txt", sep = '\t')
		graph_sum_edge = new_graph.sum_edge_scores()
		graph_sum_edge.to_csv(core_gene_file + '_' + tissue + "_sum_edge.txt", sep = '\t')
		#graph_mean_edge = new_graph.mean_edge_scores()
		#graph_mean_edge.to_csv(core_gene_file + '_' + tissue + "_mean_edge.txt", sep = '\t')
		graph_max_edge = new_graph.max_edge_scores()
		graph_max_edge.to_csv(core_gene_file + '_' + tissue + "_max_edge.txt", sep = '\t')
		graph_closeness = new_graph.closeness_centrality()
		graph_closeness.to_csv(core_gene_file + '_' + tissue + "_closeness_centrality.txt", sep = '\t')
		graph_betweenness = new_graph.betweenness_centrality()
		graph_betweenness.to_csv(core_gene_file + '_' + tissue + "_betweenness_centrality.txt", sep = '\t')
		#graph_local_cluster_coeff = new_graph.local_cluster_coeff()
		#graph_local_cluster_coeff.to_csv(core_gene_file + '_' + tissue + "_local_cluster_coeff.txt", sep = '\t')
		graph_degree = new_graph.vertex_degree()
		graph_degree.to_csv(core_gene_file + '_' + tissue + '_degree.txt', sep = '\t')


	else:
		tissue = sys.argv[1]
		print 'Inputted tissue file {0} with no core genes'.format(tissue)

		if tissue.split('.')[-1] == 'txt':
			graph = neighbor_graph(tissue + '.gz')
		elif tissue.split('_')[-1] == 'top.gz':
			graph = neighbor_graph(tissue)
		else:
			graph = neighbor_graph(tissue)
		new_graph = Gene_interaction_graph(graph)
		del graph
		tissue = tissue.split('.')[0]
		graph_eigenvector = new_graph.eigenvector_centrality()[0]
		graph_eigenvector.to_csv(tissue + "_pagerank.txt", sep = '\t')
		graph_pagerank = new_graph.pagerank()
		graph_pagerank.to_csv(tissue + "_eigenvector_centrality.txt", sep = '\t')
		graph_degree = new_graph.degree_centrality()
		graph_degree.to_csv(tissue + "_degree_centrality.txt", sep = '\t')
		graph_sum_edge = new_graph.sum_edge_scores()
		graph_sum_edge.to_csv(tissue + "_sum_edge.txt", sep = '\t')
		graph_mean_edge = new_graph.mean_edge_scores()
		graph_mean_edge.to_csv(tissue + "_mean_edge.txt", sep = '\t')
		graph_max_edge = new_graph.max_edge_scores()
		graph_max_edge.to_csv(tissue + "_max_edge.txt", sep = '\t')
		graph_closeness = new_graph.closeness_centrality()
		graph_closeness.to_csv(tissue + "_closeness_centrality.txt", sep = '\t')
		graph_betweenness = new_graph.betweenness_centrality()
		graph_betweenness.to_csv(tissue + "_betweenness_centrality.txt", sep = '\t')
		graph_local_cluster_coeff = new_graph.local_cluster_coeff()
		graph_local_cluster_coeff.to_csv(tissue + "_local_cluster_coeff.txt", sep = '\t')
		graph_degree = new_graph.vertex_degree()
		graph_degree.to_csv(tissue + '_degree.txt', sep = '\t')
		# new_graph.shortest_path_distribution(tissue)
main()
