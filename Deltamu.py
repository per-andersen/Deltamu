import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import pickle as pick
import astropy.cosmology as cosmo
import CosmoMC_dao
import Contour
import FlatJBP
import Flatn3CPL
import Flatn7CPL
import time as ti

def lin_interp_array(X,Y,X_interp):
    #X and Y are arrays to draw interpolated values from
    #X_interp is an array of values we wish to calculate interpoalted Y values for
    if (len(X) != len(Y)): #Sanity check
        print "Mismatched input lengths!", exit()

    Y_interp = np.zeros(len(X_interp)) #Storage of interpolated results
    for j in np.arange(len(X_interp)):
        i = 0
        if ( (X_interp[j] > np.max(X)) or (X_interp[j] < np.min(X)) ):
            print np.min(X), X_interp[j], np.max(X), "Interpolation Impossible!"
        while (X[i] < X_interp[j]):
            i += 1
            if i == 0:
                Y_interp[j] = Y[i] + (Y[i+1]-Y[i]) * (X_interp[j] - X[i]) / (X[i+1] - X[i])
            else:
                Y_interp[j] = Y[i-1] + (Y[i]-Y[i-1]) * (X_interp[j] - X[i-1]) / (X[i] - X[i-1])
    return Y_interp

def chi2(m, redshifts, distmod, redshifts_marg, distmod_marg, sigma_marg):
	return np.sum( (lin_interp_array(redshifts,distmod,redshifts_marg) + m  - distmod_marg )**2 / sigma_marg**2 )

class Deltamu(object):
	"""A contour class for 3d contours
	"""
	def __init__(self, chain_name, cosmology_string, directory='/Users/perandersen/Data/HzSC/',\
		redshifts=np.linspace(0.001,10.,1000), sigmas_marg = None, contour_level = 0.68,\
		tolerance = 0.005, bins_tuple=(20,20,20), do_marg=False, redshift_marg_min=0.1,\
		redshift_marg_max=10., redshift_marg_n=10, redshifts_marg_method='jla',smoothing=0.):

		self.chain_name = chain_name
		self.cosmology_string = cosmology_string
		self.directory = directory
		self.redshifts = redshifts
		self.sigmas_marg = sigmas_marg
		self.contour_level = contour_level
		self.tolerance = tolerance
		self.bins_tuple = bins_tuple
		self.do_marg = do_marg
		self.redshift_marg_min = redshift_marg_min
		self.redshift_marg_max = redshift_marg_max
		self.redshift_marg_n = redshift_marg_n
		self.redshifts_marg_method = redshifts_marg_method
		self.smoothing = smoothing

		if self.do_marg:
			if self.redshifts_marg_method == 'lin':
				self.redshifts_marg = np.linspace(self.redshift_marg_min,self.redshift_marg_max,self.redshift_marg_n)
			elif self.redshifts_marg_method == 'log':
				self.redshifts_marg = np.logspace(self.redshift_marg_min,self.redshift_marg_max,self.redshift_marg_n)
			elif self.redshifts_marg_method == 'jla':
				self.redshifts_marg, self.distmod_marg, self.sigmas_marg = self.read_jla_sample()
				self.redshift_marg_n = len(self.redshifts_marg)
				self.redshift_marg_max = np.max(self.redshifts_marg)
				self.redshift_marg_min = np.min(self.redshifts_marg)

			if type(self.sigmas_marg).__module__ != np.__name__:
				self.sigmas_marg = np.ones(redshift_marg_n) * 0.1			

	def get_marg_file_name(self):
		if self.do_marg:
			if self.chain_name == 'lcdm':
				f_name = "deltamu_lcdm_c" + str(self.contour_level) +\
				"_t" + str(self.tolerance) + "_s" + str(self.smoothing) + "_b" + str(self.bins_tuple) + "_marg_" +\
				self.redshifts_marg_method +"_z" + str(self.redshift_marg_min) +\
				"-" + str(self.redshift_marg_max) + "_n" + str(self.redshift_marg_n) + ".dat"
			else:
				f_name = "deltamu_" + self.chain_name + "_c" + str(self.contour_level) +\
				"_t" + str(self.tolerance) + "_s" + str(self.smoothing) + "_b" + str(self.bins_tuple[0]) + \
				str(self.bins_tuple[1]) + str(self.bins_tuple[2]) + "_marg_" +\
				self.redshifts_marg_method +"_z" + str(self.redshift_marg_min) +\
				"-" + str(self.redshift_marg_max) + "_n" + str(self.redshift_marg_n) + ".dat"
			return f_name
		else:
			print "Not doing marginalisation this run!"
			raise

	def run_prep(self):
		self.lcdm_bestfit = cosmo.FlatLambdaCDM(H0=hubble_const, Om0=omega_m_lcdm_bestfit)
		self.distmod_bestfit_lcdm = np.array(self.lcdm_bestfit.distmod(self.redshifts))
		if self.do_marg:
				self.m_bestfit_lcdm_marg = opt.fmin(chi2,x0=-19.,args=(self.redshifts, self.distmod_bestfit_lcdm, self.redshifts_marg, self.distmod_marg, self.sigmas_marg),xtol=0.005,disp=0)

		self.deltamu = self.get_deltamu()

	def read_jla_sample(self):
		jla_sample = np.genfromtxt(self.directory + 'JLA.txt')
		return jla_sample[:,0], jla_sample[:,1], jla_sample[:,2]

	def write_minmax_deltamuparameters(self):
		'''Writes min/max deltamu and parameters to pickle file
		'''
		self.run_prep()

		if self.chain_name == 'lcdm':
			f_name = "deltamu_lcdm_c" + str(self.contour_level) +\
			"_t" + str(self.tolerance) + "_s" + str(self.smoothing) + "_b" + str(self.bins_tuple) + ".dat"

		else:
			f_name = "deltamu_" + self.chain_name + "_c" + str(self.contour_level) +\
			"_t" + str(self.tolerance) + "_s" + str(self.smoothing) + "_b" + str(self.bins_tuple[0]) + \
			str(self.bins_tuple[1]) + str(self.bins_tuple[2]) + ".dat"

		parameters_min, parameters_max, deltamu_min, deltamu_max = self.get_minmax_deltamuparameters(margi=0.)

		print "Dumping data to file", f_name
		data_dump = self.redshifts, deltamu_min, deltamu_max, parameters_min, parameters_max, self.deltamu, self.marg
		data_file_name = self.directory + "Deltamu/" + f_name
		output = open(data_file_name,'wb')
		pick.dump(data_dump,output)
		output.close()

		if self.do_marg:
			if self.chain_name == 'lcdm':
				f_name = "deltamu_lcdm_c" + str(self.contour_level) +\
				"_t" + str(self.tolerance) + "_s" + str(self.smoothing) + "_b" + str(self.bins_tuple) + "_marg_" +\
				self.redshifts_marg_method +"_z" + str(self.redshift_marg_min) +\
				"-" + str(self.redshift_marg_max) + "_n" + str(self.redshift_marg_n) + ".dat"
			else:
				f_name = "deltamu_" + self.chain_name + "_c" + str(self.contour_level) +\
				"_t" + str(self.tolerance) + "_s" + str(self.smoothing) + "_b" + str(self.bins_tuple[0]) + \
				str(self.bins_tuple[1]) + str(self.bins_tuple[2]) + "_marg_" +\
				self.redshifts_marg_method +"_z" + str(self.redshift_marg_min) +\
				"-" + str(self.redshift_marg_max) + "_n" + str(self.redshift_marg_n) + ".dat"

			parameters_min, parameters_max, deltamu_min, deltamu_max = self.get_minmax_deltamuparameters(margi=self.marg-self.m_bestfit_lcdm_marg)

			print "Dumping data to file", f_name
			data_dump = self.redshifts, deltamu_min, deltamu_max, parameters_min, parameters_max, self.deltamu, self.marg
			data_file_name = self.directory + "Deltamu/"+ f_name
			output = open(data_file_name,'wb')
			pick.dump(data_dump,output)
			output.close()

	def get_deltamu(self):
		'''Returns deltamu values for given cosmology
		'''
		self.parameter_sets = self.get_parameters()
		parameter_sets = self.parameter_sets
		print np.shape(parameter_sets)
		length_loop = np.shape(parameter_sets)[-1]
		if self.do_marg == True:
			self.marg = np.zeros(length_loop)
		print self.chain_name, length_loop
		deltamu = np.zeros((len(self.redshifts),length_loop))
		distmod = np.zeros((len(self.redshifts),length_loop))

		if self.do_marg:
			if self.redshifts_marg_method == 'jla':
				pass
			else:
				self.distmod_marg = np.array(lcdm_bestfit.distmod(self.redshifts_marg))

		print "Calculating deltamu values..."
		for ii in np.arange(length_loop):
			#print ii
			exec(self.cosmology_string)
			distmod[:,ii] = cosmology.distmod(self.redshifts)
			deltamu[:,ii] = distmod[:,ii] - np.array(self.distmod_bestfit_lcdm)
					
			if self.do_marg:
				m_marg = opt.fmin(chi2,x0=-19.,args=(self.redshifts, distmod[:,ii], self.redshifts_marg, self.distmod_marg, self.sigmas_marg),xtol=0.005,disp=0)
				#print m_marg
				self.marg[ii] = m_marg
		return deltamu

	def get_minmax_deltamuparameters(self, margi):
		'''Returns min/max deltamu and associated
		parameter values
		'''
		deltamu_max = np.zeros(len(self.redshifts))
		deltamu_min = np.zeros(len(self.redshifts))
		if len(np.shape(self.parameter_sets)) == 2:
			parameters_max = np.zeros((len(self.redshifts),3))
			parameters_min = np.zeros((len(self.redshifts),3))
			for jj in np.arange(len(self.redshifts)):
				deltamu_marg = self.deltamu[jj,:] + margi
				deltamu_max[jj] = np.max(deltamu_marg)
				deltamu_min[jj] = np.min(deltamu_marg)
				param_max = np.array([ self.parameter_sets[0][np.argmax(deltamu_marg)], self.parameter_sets[1][np.argmax(deltamu_marg)], self.parameter_sets[2][np.argmax(deltamu_marg)] ])
				param_min = np.array([ self.parameter_sets[0][np.argmin(deltamu_marg)], self.parameter_sets[1][np.argmin(deltamu_marg)], self.parameter_sets[2][np.argmin(deltamu_marg)] ])
				parameters_max[jj,:] = param_max
				parameters_min[jj,:] = param_min
		else:
			parameters_max = np.zeros(len(self.redshifts))
			parameters_min = np.zeros(len(self.redshifts))
			for jj in np.arange(len(self.redshifts)):
				deltamu_marg = self.deltamu[jj,:] + margi
				deltamu_max[jj] = np.max(deltamu_marg)
				deltamu_min[jj] = np.min(deltamu_marg)
				param_max = np.array([self.parameter_sets[np.argmax(deltamu_marg)]])
				param_min = np.array([self.parameter_sets[np.argmin(deltamu_marg)]])
				parameters_max[jj] = param_max
				parameters_min[jj] = param_min

		return parameters_min, parameters_max, deltamu_min, deltamu_max

	def get_parameters(self):
		'''Returns parameter sets on contour of given cosmology
		'''
		if self.chain_name == 'lcdm':
			Parameter_contour = Contour.LCDM_Contour(self.chain_name,self.directory,self.contour_level,self.tolerance,self.bins_tuple,self.smoothing)
		else:
			Parameter_contour = Contour.Contour(self.chain_name,self.directory,self.contour_level,self.tolerance,self.bins_tuple,self.smoothing)

		print "Does", self.chain_name ,"contour exist?",
		if Parameter_contour.test_contour_exists():
			print "Yes"
		else:
			print "No"
			Parameter_contour.pickle_contour()

		if self.chain_name == 'lcdm':
			parameter_sets = Parameter_contour.read_pickled_contour()
		else:
			omega_contour, w0_contour, wa_contour = Parameter_contour.read_pickled_contour()
			parameter_sets = omega_contour, w0_contour, wa_contour
		
		return parameter_sets

if __name__ == "__main__":
	#These are comstants set from LCDM CosmoMC marginalised bestfit.
	#Hubble const can be set to any constant value, it is only defined for astropy to work.

	omega_m_lcdm_bestfit, hubble_const = 0.30754277645, 70.

	lcdm_string = 'cosmology=cosmo.FlatLambdaCDM(H0=hubble_const, Om0=parameter_sets[ii])'
	cpl_string = 'cosmology=cosmo.Flatw0waCDM(H0=hubble_const, Om0=parameter_sets[0][ii], w0=parameter_sets[1][ii], wa=parameter_sets[2][ii])'
	jbp_string = 'cosmology=FlatJBP.FlatJBP_CDM(H0=hubble_const, Om0=parameter_sets[0][ii], w0=parameter_sets[1][ii], wa=parameter_sets[2][ii])'
	n3cpl_string = 'cosmology=Flatn3CPL.Flatn3CPL(H0=hubble_const, Om0=parameter_sets[0][ii], w0=parameter_sets[1][ii], wa=parameter_sets[2][ii])'
	n7cpl_string = 'cosmology=Flatn7CPL.Flatn7CPL(H0=hubble_const, Om0=parameter_sets[0][ii], w0=parameter_sets[1][ii], wa=parameter_sets[2][ii])'

	t0 = ti.time()
	
	lcdm_bins = [70, 80, 90]
	for bins in lcdm_bins:
		Deltamu_lcdm = Deltamu('lcdm',lcdm_string,tolerance = 0.01, bins_tuple=bins,do_marg=True)
		Deltamu_lcdm.write_minmax_deltamuparameters()
	

	'''
	cpl_smoothing = [0., 0.6]
	cpl_bins = [(30,30,30),(40,40,40),(50,50,50),(60,60,60)]
	for smoothing in cpl_smoothing:
		for bins in cpl_bins:
			Deltamu_cpl = Deltamu('cpl', cpl_string, bins_tuple=bins, do_marg=True,smoothing=smoothing)
			Deltamu_cpl.write_minmax_deltamuparameters()
	'''

	'''
	jbp_smoothing = [0., 0.6]
	jbp_bins = [(30,30,30),(40,40,40),(50,50,50),(60,60,60)]
	for smoothing in jbp_smoothing:
		for bins in jbp_bins:
			Deltamu_jbp = Deltamu('jbp', jbp_string, bins_tuple=bins, do_marg=True, smoothing=smoothing)
			Deltamu_jbp.write_minmax_deltamuparameters()
	'''

	'''
	n3cpl_smoothing = [0., 0.6]
	n3cpl_bins = [(30,30,30),(40,40,40),(50,50,50),(60,60,60)]
	for smoothing in n3cpl_smoothing:
		for bins in n3cpl_bins:
			Deltamu_n3cpl = Deltamu('n3cpl', n3cpl_string, bins_tuple=bins, do_marg=True, smoothing=smoothing)
			Deltamu_n3cpl.write_minmax_deltamuparameters()
	'''

	'''
	n7cpl_smoothing = [0., 0.2, 0.4]
	n7cpl_bins = [(20,20,20),(30,30,30),(40,40,40)]
	for smoothing in n7cpl_smoothing:
		for bins in n7cpl_bins:
			Deltamu_n7cpl = Deltamu('n7cpl', n7cpl_string, bins_tuple=bins, do_marg=True, smoothing=smoothing)
			Deltamu_n7cpl.write_minmax_deltamuparameters()
	print "Classes done in:", ti.time() - t0, "seconds"
	'''



