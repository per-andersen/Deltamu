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

def chi2(m, redshifts, distmod, redshifts_marg, distmod_lcdm_marg, sigma_marg):
	return np.sum( (lin_interp_array(redshifts,distmod,redshifts_marg) + m  - distmod_lcdm_marg )**2 / sigma_marg**2 )

class Deltamu(object):
	"""A contour class for 3d contours
	"""
	def __init__(self, chain_name, cosmology_string, directory='/Users/perandersen/Data/HzSC/', redshifts=np.linspace(0.001,10.,1000), redshifts_marg = np.linspace(0.1,10.,10), sigmas_marg = None, contour_level = 0.68, tolerance = 0.005, bins_tuple=(50,50,50), do_marg=False):
		self.chain_name = chain_name
		self.cosmology_string = cosmology_string
		self.directory = directory
		self.redshifts = redshifts
		self.redshifts_marg = redshifts_marg
		self.sigmas_marg = sigmas_marg
		self.contour_level = contour_level
		self.tolerance = tolerance
		self.bins_tuple = bins_tuple
		self.do_marg = do_marg
		if self.sigmas_marg == None:
			self.sigmas_marg = np.ones(len(redshifts_marg)) * 0.1


		self.deltamu = self.get_deltamu()
		self.lcdm_bestfit = cosmo.FlatLambdaCDM(H0=hubble_const, Om0=omega_m_lcdm_bestfit)
		self.distmod_bestfit_lcdm = lcdm_bestfit.distmod(self.redshifts)

	def write_minmax_deltamuparameters(self):
		'''Writes min/max deltamu and parameters to pickle file
		'''
		if self.chain_name == 'lcdm':
			f_name = "deltamu_lcdm_c" + str(self.contour_level) +\
			"_t" + str(self.tolerance) + "_b" + str(self.bins_tuple) + ".dat"
		else:
			f_name = "deltamu_" + self.chain_name + "_c" + str(self.contour_level) +\
			"_t" + str(self.tolerance) + "_b" + str(self.bins_tuple[0]) + \
			str(self.bins_tuple[1]) + str(self.bins_tuple[2]) + ".dat"

		parameters_min, parameters_max, deltamu_min, deltamu_max = self.get_minmax_deltamuparameters(margi=0.)

		print "Dumping data to file", f_name
		data_dump = self.redshifts, deltamu_min, deltamu_max, parameters_min, parameters_max
		data_file_name = self.directory + "/Deltamu/" + f_name
		output = open(data_file_name,'wb')
		pick.dump(data_dump,output)
		output.close()

		if self.do_marg:
			if self.chain_name == 'lcdm':
				f_name = "deltamu_lcdm_c" + str(self.contour_level) +\
				"_t" + str(self.tolerance) + "_b" + str(self.bins_tuple) + "_marg.dat"
			else:
				f_name = "deltamu_" + self.chain_name + "_c" + str(self.contour_level) +\
				"_t" + str(self.tolerance) + "_b" + str(self.bins_tuple[0]) + \
				str(self.bins_tuple[1]) + str(self.bins_tuple[2]) + "_marg.dat"

			parameters_min, parameters_max, deltamu_min, deltamu_max = self.get_minmax_deltamuparameters(margi=self.marg)

			print "Dumping data to file", f_name
			data_dump = self.redshifts, deltamu_min, deltamu_max, parameters_min, parameters_max
			data_file_name = self.directory + "/Deltamu/"+ f_name
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
			self.distmod_lcdm_marg = np.array(lcdm_bestfit.distmod(self.redshifts_marg))

		print "Calculating deltamu values..."
		for ii in np.arange(length_loop):
			#print ii
			exec(self.cosmology_string)
			distmod[:,ii] = cosmology.distmod(self.redshifts)
			deltamu[:,ii] = distmod[:,ii] - np.array(distmod_bestfit_lcdm)
					
			if self.do_marg:
				m_marg = opt.fmin(chi2,x0=0.,args=(self.redshifts, distmod[:,ii], self.redshifts_marg, self.distmod_lcdm_marg, self.sigmas_marg),xtol=0.005,disp=0)
				self.marg[ii] = m_marg
		return deltamu

	def get_minmax_deltamuparameters(self, margi=0.):
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
			Parameter_contour = Contour.LCDM_Contour(self.chain_name,self.directory,self.contour_level,self.tolerance,self.bins_tuple)
		else:
			Parameter_contour = Contour.Contour(self.chain_name,self.directory,self.contour_level,self.tolerance,self.bins_tuple)

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

omega_m_lcdm_bestfit, hubble_const = 0.30754277645, 70.
redshifts=np.linspace(0.001,10.,1000)

lcdm_bestfit = cosmo.FlatLambdaCDM(H0=hubble_const, Om0=omega_m_lcdm_bestfit)
distmod_bestfit_lcdm = lcdm_bestfit.distmod(redshifts)

lcdm_string = 'cosmology=cosmo.FlatLambdaCDM(H0=hubble_const, Om0=parameter_sets[ii])'
cpl_string = 'cosmology=cosmo.Flatw0waCDM(H0=hubble_const, Om0=parameter_sets[0][ii], w0=parameter_sets[1][ii], wa=parameter_sets[2][ii])'
jbp_string = 'cosmology=FlatJBP.FlatJBP_CDM(H0=hubble_const, Om0=parameter_sets[0][ii], w0=parameter_sets[1][ii], wa=parameter_sets[2][ii])'

t0 = ti.time()
Deltamu_lcdm = Deltamu('lcdm',lcdm_string,tolerance = 0.01, bins_tuple=100,do_marg=True)
Deltamu_lcdm.write_minmax_deltamuparameters()

#Deltamu_cpl = Deltamu('cpl',cpl_string,do_marg=True)
#Deltamu_cpl.write_minmax_deltamuparameters()

#Deltamu_jbp = Deltamu('jbp',jbp_string, bins_tuple=(20,20,20),do_marg=True)
#Deltamu_jbp.write_minmax_deltamuparameters()
print "Classes done in:", ti.time() - t0, "seconds"



