from diff_prof.diffusion_profiles import DiffusionProfiles
from msi.msi import MSI
import os
import numpy as np

def test_diffusion_profiles(ref_dp_file_path, calculated_dp_file_path):
	# Saved
	dp_saved = DiffusionProfiles(alpha = None, max_iter = None, tol = None, weights = None, num_cores = None, save_load_file_path = ref_dp_file_path)
	msi_saved = MSI()
	msi_saved.load()
	msi_saved.load_saved_node_idx_mapping_and_nodelist(dp_saved.save_load_file_path)
	dp_saved.load_diffusion_profiles(msi_saved.drugs_in_graph + msi_saved.indications_in_graph)

	# Calculated
	dp_calculated = DiffusionProfiles(alpha = None, max_iter = None, tol = None, weights = None, num_cores = None, save_load_file_path = calculated_dp_file_path)
	msi_calculated = MSI()
	msi_calculated.load()
	msi_calculated.load_saved_node_idx_mapping_and_nodelist(dp_calculated.save_load_file_path)
	dp_calculated.load_diffusion_profiles(msi_calculated.drugs_in_graph + msi_calculated.indications_in_graph)

	# Compare
	## Make sure have diffusion profiles for the same drugs and indications
	assert(set(dp_saved.drug_or_indication2diffusion_profile.keys()) == set(dp_calculated.drug_or_indication2diffusion_profile.keys()))
	## Reorder calculated diffusion profile for consistency with saved diffusion profile
	calculated_reorder_idxs = [msi_calculated.node2idx[node] for node in msi_saved.nodelist]

	for drug_or_indication, saved_diffusion_profile in dp_saved.drug_or_indication2diffusion_profile.items():
	    calculated_diffusion_profile = dp_calculated.drug_or_indication2diffusion_profile[drug_or_indication]

	    # Reorder calculated diffusion_profile according to saved
	    calculated_diffusion_profile = calculated_diffusion_profile[calculated_reorder_idxs]

	    # Ensure close enough
	    assert(np.allclose(saved_diffusion_profile, calculated_diffusion_profile))
