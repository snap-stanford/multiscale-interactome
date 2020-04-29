from .node_to_node import NodeToNode
from .drug_to_protein import DrugToProtein
from .indication_to_protein import IndicationToProtein
from .protein_to_protein import ProteinToProtein
from .protein_to_functional_pathway import ProteinToFunctionalPathway
from .functional_pathway_to_functional_pathway import FunctionalPathwayToFunctionalPathway

import networkx as nx
import pandas as pd
import os
import scipy
import pickle

DRUG = "drug"
INDICATION = "indication"
PROTEIN = "protein"
FUNCTIONAL_PATHWAY = "functional_pathway"
UP_FUNCTIONAL_PATHWAY = "up_functional_pathway"
DOWN_FUNCTIONAL_PATHWAY = "down_functional_pathway"
WEIGHT = "weight"

DRUG_PROTEIN = "drug-protein"
INDICATION_PROTEIN = "indication-protein"
PROTEIN_PROTEIN = "protein-protein"
PROTEIN_FUNCTIONAL_PATHWAY = "protein-functional_pathway"
FUNCTIONAL_PATHWAY_FUNCTIONAL_PATHWAY = "functional_pathway-functional_pathway"

class MSI():
	def __init__(self, nodes = [DRUG, INDICATION, PROTEIN, FUNCTIONAL_PATHWAY], edges = [DRUG_PROTEIN, INDICATION_PROTEIN, PROTEIN_PROTEIN, PROTEIN_FUNCTIONAL_PATHWAY, FUNCTIONAL_PATHWAY_FUNCTIONAL_PATHWAY], drug2protein_file_path = "data/drug_to_protein.tsv", drug2protein_directed = False, indication2protein_file_path = "data/indication_to_protein.tsv", indication2protein_directed = False, protein2protein_file_path = "data/protein_to_protein.tsv", protein2protein_directed = False, protein2functional_pathway_file_path = "data/protein_to_functional_pathway.tsv", protein2functional_pathway_directed = False, functional_pathway2functional_pathway_file_path = "data/functional_pathway_to_functional_pathway.tsv", functional_pathway2functional_pathway_directed = True):
		# Parameters
		self.nodes = nodes
		self.edges = edges

		# File paths
		self.drug2protein_file_path = drug2protein_file_path
		self.indication2protein_file_path = indication2protein_file_path
		self.protein2protein_file_path = protein2protein_file_path
		self.protein2functional_pathway_file_path = protein2functional_pathway_file_path
		self.functional_pathway2functional_pathway_file_path = functional_pathway2functional_pathway_file_path

		# Directed
		self.drug2protein_directed = drug2protein_directed
		self.indication2protein_directed = indication2protein_directed
		self.protein2protein_directed = protein2protein_directed
		self.protein2functional_pathway_directed = protein2functional_pathway_directed
		self.functional_pathway2functional_pathway_directed = functional_pathway2functional_pathway_directed

	def add_edges(self, edge_list, from_node_type, to_node_type):
		for from_node, to_node in edge_list:
			self.graph.add_edge(from_node, to_node)
			self.graph.node[from_node]["type"] = from_node_type
			self.graph.node[to_node]["type"] = to_node_type

	def merge_one_to_one_dicts(self, dict_list):
		out_dict = dict()
		for dict_ in dict_list:
			for k, v in dict_.items():
				if k in out_dict:
					assert((out_dict[k] == v) or (pd.isnull(out_dict[k]) and pd.isnull(v))) 
				else:
					out_dict[k] = v
		return out_dict

	def load_node2type(self):
		# Merge the node2type of each component (these are one to one dictionaries)
		node2type__list = []
		for node2node__name, node2node__obj in self.components.items():
			node2type__list.append(node2node__obj.node2type)
		node2type = self.merge_one_to_one_dicts(node2type__list)
		self.node2type = node2type

	def load_type2nodes(self):
		type2nodes = {}
		for node, type_ in self.node2type.items():
			if type_ in type2nodes:
				type2nodes[type_].add(node)
			else:
				type2nodes[type_] = set([node])
		self.type2nodes = type2nodes

	def load_node2name(self):
		node2name__list = []
		for node2node__name, node2node__obj in self.components.items():
			node2name__list.append(node2node__obj.node2name)
		node2name = self.merge_one_to_one_dicts(node2name__list)
		self.node2name = node2name

	def load_name2node(self):
		name2node = {v:k for k,v in self.node2name.items()}
		self.name2node = name2node

	def load_graph(self):
		# Initialize graph
		self.graph = nx.Graph()

		# Load components and add edges as appropriate
		self.components = dict()
		if (DRUG in self.nodes) and (DRUG_PROTEIN in self.edges):
			self.components["drug_to_protein"] = DrugToProtein(self.drug2protein_directed, self.drug2protein_file_path)
			self.add_edges(self.components["drug_to_protein"].edge_list, DRUG, PROTEIN)

		if (INDICATION in self.nodes) and (INDICATION_PROTEIN in self.edges):
			self.components["indication_to_protein"] = IndicationToProtein(self.indication2protein_directed, self.indication2protein_file_path)
			self.add_edges(self.components["indication_to_protein"].edge_list, INDICATION, PROTEIN)

		if (PROTEIN in self.nodes) and (PROTEIN_PROTEIN in self.edges):
			self.components["protein_to_protein"] = ProteinToProtein(self.protein2protein_directed, self.protein2protein_file_path)
			self.add_edges(self.components["protein_to_protein"].edge_list, PROTEIN, PROTEIN)

		if (FUNCTIONAL_PATHWAY in self.nodes) and (PROTEIN_FUNCTIONAL_PATHWAY in self.edges):
			self.components["protein_to_functional_pathway"] = ProteinToFunctionalPathway(self.protein2functional_pathway_directed, self.protein2functional_pathway_file_path)
			self.add_edges(self.components["protein_to_functional_pathway"].edge_list, PROTEIN, FUNCTIONAL_PATHWAY)

		if (FUNCTIONAL_PATHWAY in self.nodes) and (FUNCTIONAL_PATHWAY_FUNCTIONAL_PATHWAY in self.edges):
			self.components["functional_pathway_to_functional_pathway"] = FunctionalPathwayToFunctionalPathway(self.functional_pathway2functional_pathway_directed, self.functional_pathway2functional_pathway_file_path)
			self.add_edges(self.components["functional_pathway_to_functional_pathway"].edge_list, FUNCTIONAL_PATHWAY, FUNCTIONAL_PATHWAY)

		# Make graph directional (copy forward and reverse of each edge)
		self.graph = self.graph.to_directed()

	def load_node_idx_mapping_and_nodelist(self):
		# Prepare
		nodes = self.graph.nodes()
		nodelist = []
		node2idx = dict.fromkeys(nodes)
		for idx, node in enumerate(nodes):
			nodelist.append(node)
			node2idx[node] = idx
		idx2node = {v:k for k,v in node2idx.items()}

		# Save
		self.nodelist = nodelist
		self.node2idx = node2idx
		self.idx2node = idx2node

	def load_saved_node_idx_mapping_and_nodelist(self, save_load_file_path):
		# Load node2idx
		node2idx_file = os.path.join(save_load_file_path, "node2idx.pkl")
		assert(os.path.exists(node2idx_file))
		with open(node2idx_file, "rb") as f:
			node2idx = pickle.load(f)
		self.node2idx = node2idx
		
		# Load idx2node
		self.idx2node = {v: k for k, v in self.node2idx.items()}
		
		# Load nodelist
		nodelist = []
		for i in range(0, len(self.idx2node)):
			nodelist.append(self.idx2node[i])
		self.nodelist = nodelist

	def save_node2idx(self, save_load_file_path):
		node2idx_file_path = os.path.join(save_load_file_path, "node2idx.pkl")
		# assert(not(os.path.isfile(node2idx_file_path)))
		with open(node2idx_file_path, "wb") as f:
			pickle.dump(self.node2idx, f)

	def load_drugs_in_graph(self):
		self.drugs_in_graph = list(self.type2nodes[DRUG])

	def load_indications_in_graph(self):
		self.indications_in_graph = list(self.type2nodes[INDICATION])
	
	def load_drug_or_indication2proteins(self):
		drug_or_indication2proteins = dict.fromkeys(self.drugs_in_graph + self.indications_in_graph)
		for drug in self.drugs_in_graph:
			assert(not(nx.is_directed(self.components["drug_to_protein"].graph)))
			drug_or_indication2proteins[drug] = set(self.components["drug_to_protein"].graph.neighbors(drug))
		for indication in self.indications_in_graph:
			assert(not(nx.is_directed(self.components["indication_to_protein"].graph)))
			drug_or_indication2proteins[indication] = set(self.components["indication_to_protein"].graph.neighbors(indication))
		self.drug_or_indication2proteins = drug_or_indication2proteins

	def load(self):
		self.load_graph()
		self.load_node_idx_mapping_and_nodelist()
		self.load_node2type()
		self.load_type2nodes()
		self.load_node2name()
		self.load_name2node()
		self.load_drugs_in_graph()
		self.load_indications_in_graph()
		self.load_drug_or_indication2proteins()

	def save_graph(self, save_load_file_path):
		graph_file_path = os.path.join(save_load_file_path, "graph.pkl")
		# assert(not(os.path.isfile(graph_file_path)))
		with open(graph_file_path, "wb") as f:
			pickle.dump(self.graph, f)

	def add_to_cs_adj_dict(self, node, successor_type, successor):
		if (successor_type in self.cs_adj_dict[node]):
			self.cs_adj_dict[node][successor_type].append(successor)
		else:
			self.cs_adj_dict[node][successor_type] = [successor]

	def create_class_specific_adjacency_dictionary(self):
		self.cs_adj_dict = {node: {} for node in self.graph.nodes()}
		for node in self.graph.nodes():
			node_type = self.node2type[node]
			if node_type == FUNCTIONAL_PATHWAY:
				up_neighbors = list(self.components["functional_pathway_to_functional_pathway"].graph.successors(node))
				down_neighbors = list(self.components["functional_pathway_to_functional_pathway"].graph.predecessors(node))

			successors = self.graph.successors(node)
			for successor in successors:
				successor_type = self.node2type[successor]
				if (node_type == FUNCTIONAL_PATHWAY) and (successor_type == FUNCTIONAL_PATHWAY):
					if successor in up_neighbors:
						self.add_to_cs_adj_dict(node, UP_FUNCTIONAL_PATHWAY, successor)
					elif successor in down_neighbors:
						self.add_to_cs_adj_dict(node, DOWN_FUNCTIONAL_PATHWAY, successor)
					else:
						assert(False)
				else:
					self.add_to_cs_adj_dict(node, successor_type, successor)

	def weight_graph(self, weights):
		self.create_class_specific_adjacency_dictionary()
		for from_node, adj_dict in self.cs_adj_dict.items():
			for node_type, to_nodes in adj_dict.items():
				num_typed_nodes = len(to_nodes)
				for to_node in to_nodes:
					self.graph[from_node][to_node][WEIGHT] = weights[node_type]/float(num_typed_nodes)