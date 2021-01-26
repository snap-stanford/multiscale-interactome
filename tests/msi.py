from msi.msi import MSI
import pickle
import networkx as nx

def test_msi():
	# Compare whether constructed MSI and saved MSI are the same
	with open("data/10_top_msi/di_protein_go_graph.pkl", "rb") as f:
		saved_graph = pickle.load(f)

	msi = MSI()
	msi.load()
	msi_graph = msi.graph

	# Same nodes?
	assert(set(msi_graph.nodes()) == set(saved_graph.nodes()))

	# Same edges?
	assert(set(msi_graph.edges()) == set(saved_graph.edges()))

	# Each node has same type?
	node_mapping = {"go": "biological_function"}
	for node in msi_graph.nodes():
		assert(msi_graph.node[node]["type"] == node_mapping.get(saved_graph.node[node]["type"], saved_graph.node[node]["type"]))
