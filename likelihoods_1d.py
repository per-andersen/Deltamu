import numpy as np
import matplotlib.pyplot as plt
import CosmoMC_dao as dao


def get_distributions(chain_name, directory, fname):
	cosmomc_n7cpl = dao.CosmoMC_dao(chain_name, directory)

	data, data_labels = cosmomc_n7cpl.read_pickled_chains(fname=fname)
	
	omega_m = data[:,0]
	w0 = data[:,1]
	wa = data[:,2]

	return omega_m, w0, wa

def plot_likelihoods(chain_name, directory, fname, nbins):
	omegam, w0, wa = get_distributions(chain_name, directory, fname)
	plt.figure()
	plt.hist(omegam, nbins)

	plt.figure()
	plt.hist(w0, nbins)

	plt.figure()
	plt.hist(wa, nbins)

def get_mean_std(chain_name, directory, fname):
	omegam, w0, wa = get_distributions(chain_name, directory, fname)
	print np.mean(omegam), np.std(omegam)
	#print np.mean(w0), np.std(w0)
	#print np.mean(wa), np.std(wa)


omegam, w0, wa = get_distributions('n7cpl', '/Users/perandersen/Data/HzSC/Chains/n7CPL/', 'n7cpl_chains.pkl')

#plot_likelihoods('n7cpl', '/Users/perandersen/Data/HzSC/Chains/n7CPL/', 'n7cpl_chains.pkl',70)
#plot_likelihoods('cpl', '/Users/perandersen/Data/HzSC/Chains/CPL/', 'cpl_chains.pkl',70)

get_mean_std('cpl', '/Users/perandersen/Data/HzSC/Chains/CPL/', 'cpl_chains.pkl')
get_mean_std('jbp', '/Users/perandersen/Data/HzSC/Chains/JBP/', 'jbp_chains.pkl')
get_mean_std('n3cpl', '/Users/perandersen/Data/HzSC/Chains/n3CPL/', 'n3cpl_chains.pkl')

cosmomc_lcdm = dao.CosmoMC_dao('lcdm', '/Users/perandersen/Data/HzSC/Chains/LCDM/')
data, data_labels = cosmomc_lcdm.read_pickled_chains(fname='lcdm_chains.pkl')
print np.mean(data), np.std(data)


plt.show()