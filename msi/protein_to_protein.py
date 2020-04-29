from .node_to_node import NodeToNode

class ProteinToProtein(NodeToNode):
	def __init__(self, file_path, sep = "\t"):
		super().__init__(file_path, sep)