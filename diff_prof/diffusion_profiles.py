import os
import pickle
import multiprocessing
import networkx as nx
import math
import time
import copy
import scipy
import numpy as np
WEIGHT = 'weight'

class DiffusionProfiles():
	def __init__(self, alpha, max_iter, tol, weights, num_cores, save_load_file_path):
		self.alpha = alpha
		self.max_iter = max_iter
		self.tol = tol
		self.weights = weights
		self.num_cores = num_cores
		self.save_load_file_path = save_load_file_path

	def get_initial_M(self, msi):
		# This function is adapted from the NetworkX implementation of Personalized PageRank #
		M = nx.to_scipy_sparse_matrix(msi.graph, nodelist = msi.nodelist, weight = WEIGHT, dtype = float)
		N = len(msi.graph)
		if (N == 0):
			assert(False)
		self.initial_M = M

	def convert_M_to_make_all_drugs_indications_sinks_except_selected(self, msi, selected_drugs_and_indications):
		reconstructed_M = copy.deepcopy(self.initial_M)

		# Delete edges INTO selected drugs and indications
		for drug_or_indication in selected_drugs_and_indications:
			proteins_pointing_to_selected_drug_or_indication = msi.drug_or_indication2proteins[drug_or_indication]
			idxs_to_remove_in_edges = [msi.node2idx[protein] for protein in proteins_pointing_to_selected_drug_or_indication]
			reconstructed_M[idxs_to_remove_in_edges, msi.node2idx[drug_or_indication]] = 0

		# Delete edges OUT of nonselected drugs and indications
		idxs_to_remove_out_edges = []
		for drug_or_indication in msi.drugs_in_graph + msi.indications_in_graph:
			if not(drug_or_indication in selected_drugs_and_indications):
				proteins_drug_or_indication_points_to = msi.drug_or_indication2proteins[drug_or_indication]
				for protein in proteins_drug_or_indication_points_to:
					idxs_to_remove_out_edges.append((msi.node2idx[drug_or_indication], msi.node2idx[protein]))
		i, j = zip(*idxs_to_remove_out_edges)
		reconstructed_M[i, j] = 0.0
		return reconstructed_M

	def refine_M_S(self, M):
		# This function is adapted from the NetworkX implementation of Personalized PageRank #
		S = scipy.array(M.sum(axis=1)).flatten()
		S[S != 0] = 1.0 / S[S != 0]
		Q = scipy.sparse.spdiags(S.T, 0, *M.shape, format='csr')
		M = Q * M

		return M, S

	def get_personalization_dictionary(self, nodes_to_start_from, nodelist):
		personalization_dict = dict.fromkeys(nodelist, 0)
		N = len(nodes_to_start_from)
		for node in nodes_to_start_from:
			personalization_dict[node] = 1./N
		return personalization_dict

	def power_iteration(self, M, S, nodelist, per_dict):
		# This function is adapted from the NetworkX implementation of Personalized PageRank #
		# Size of nodelist
		N = len(nodelist)

		# Personalization vector
		missing = set(nodelist) - set(per_dict)
		if missing:
			raise NetworkXError('Personalization dictionary must have a value for every node. Missing nodes %s' % missing)
		p = scipy.array([per_dict[n] for n in nodelist], dtype=float)
		p = p / p.sum()

		# Dangling nodes
		dangling_weights = p
		is_dangling = scipy.where(S == 0)[0]

		# power iteration: make up to max_iter iterations
		x = scipy.repeat(1.0 / N, N) # Initialize; alternatively x = p
		for _ in range(self.max_iter):
			xlast = x
			x = self.alpha * (x * M + sum(x[is_dangling]) * dangling_weights) + (1 - self.alpha) * p
			# check convergence, l1 norm
			err = scipy.absolute(x - xlast).sum()
			if err < N * self.tol:
				return x
		raise NetworkXError('pagerank_scipy: power iteration failed to converge in %d iterations.' % self.max_iter)

	def clean_file_name(self, file_name):
		return "".join([c for c in file_name if c.isalpha() or c.isdigit() or c==' ' or c =="_"]).rstrip()

	def save_diffusion_profile(self, diffusion_profile, selected_drug_or_indication):
		f = os.path.join(self.save_load_file_path, self.clean_file_name(selected_drug_or_indication) + "_p_visit_array.npy")
		np.save(f, diffusion_profile)

	def calculate_diffusion_profile_batch(self, msi, selected_drugs_and_indications_batch):
		for selected_drugs_and_indications in selected_drugs_and_indications_batch:
			self.calculate_diffusion_profile(msi, selected_drugs_and_indications)    		

	def calculate_diffusion_profile(self, msi, selected_drugs_and_indications):
		assert(len(selected_drugs_and_indications) == 1) # Not enabling functionality to run from multiple
		selected_drug_or_indication = selected_drugs_and_indications[0]
		M = self.convert_M_to_make_all_drugs_indications_sinks_except_selected(msi, selected_drugs_and_indications)
		M, S = self.refine_M_S(M)
		per_dict = self.get_personalization_dictionary(selected_drugs_and_indications, msi.nodelist)
		diffusion_profile = self.power_iteration(M, S, msi.nodelist, per_dict)
		self.save_diffusion_profile(diffusion_profile, selected_drug_or_indication)

	def batch_list(self, list_, batch_size = None, num_cores = None):
		if batch_size == float('inf'):
			batched_list = [list_]
		else:
			if batch_size is None:
				batch_size = math.ceil(len(list_) / (float(num_cores)))
			batched_list = []
			for i in range(0, len(list_), batch_size):
				batched_list.append(list_[i:i + batch_size])
		return batched_list

	def calculate_diffusion_profiles(self, msi):
		# Save MSI graph and node2idx
		msi.save_graph(self.save_load_file_path) 
		msi.save_node2idx(self.save_load_file_path)

		# Weight graph
		msi.weight_graph(self.weights)

		# Prepare to run power iteration in parallel
		self.get_initial_M(msi)
		computation_list = [[i] for i in msi.drugs_in_graph + msi.indications_in_graph]

		# Run power iteration in parallel
		computation_list_batches = self.batch_list(computation_list, num_cores = self.num_cores)

		procs = []
		for computation_list_batch in computation_list_batches:
			# Don't launch more processes than the maximum number of cores
			while(len([job for job in procs if job.is_alive()]) == self.num_cores):
				time.sleep(1)

			proc = multiprocessing.Process(target = self.calculate_diffusion_profile_batch, args = (msi, computation_list_batch))
			procs.append(proc)
			proc.start()

		# Wait until all processes done before moving forward
		while(len([job for job in procs if job.is_alive()]) > 0):
			time.sleep(1)

		# Stop all of the jobs and close them
		for job in procs:
			proc.join()
			proc.terminate()

	def load_diffusion_profiles(self, drugs_and_indications):
		assert(not(self.save_load_file_path is None))

		# Load diffusion profiles
		drug_or_indication2diffusion_profile = dict()
		for drug_or_indication in drugs_and_indications:
			file_path = os.path.join(self.save_load_file_path, self.clean_file_name(drug_or_indication) + "_p_visit_array.npy")
			if (os.path.exists(file_path)):
				diffusion_profile = np.load(file_path)
				drug_or_indication2diffusion_profile[drug_or_indication] = diffusion_profile
			else:
				print("Loading failed at " + str(drug_or_indication) + " | " + str(file_path))
		self.drug_or_indication2diffusion_profile = drug_or_indication2diffusion_profile