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

def get_parameters_minmax(parameter_sets, deltamu, marg=0.):
	deltamu_max = np.zeros(len(redshifts))
	deltamu_min = np.zeros(len(redshifts))
	if len(np.shape(parameter_sets)) == 2:
		parameters_max = np.zeros((len(redshifts),3))
		parameters_min = np.zeros((len(redshifts),3))
		for jj in np.arange(len(redshifts)):
			deltamu_marg = deltamu[jj,:] + marg
			deltamu_max[jj] = np.max(deltamu_marg)
			deltamu_min[jj] = np.min(deltamu_marg)
			param_max = np.array([ parameter_sets[0][np.argmax(deltamu_marg)], parameter_sets[1][np.argmax(deltamu_marg)], parameter_sets[2][np.argmax(deltamu_marg)] ])
			param_min = np.array([ parameter_sets[0][np.argmin(deltamu_marg)], parameter_sets[1][np.argmin(deltamu_marg)], parameter_sets[2][np.argmin(deltamu_marg)] ])
			parameters_max[jj,:] = param_max
			parameters_min[jj,:] = param_min
	else:
		parameters_max = np.zeros(len(redshifts))
		parameters_min = np.zeros(len(redshifts))
		for jj in np.arange(len(redshifts)):
			deltamu_marg = deltamu[jj,:] + marg
			deltamu_max[jj] = np.max(deltamu_marg)
			deltamu_min[jj] = np.min(deltamu_marg)
			param_max = np.array([parameter_sets[np.argmax(deltamu_marg)]])
			param_min = np.array([parameter_sets[np.argmin(deltamu_marg)]])
			parameters_max[jj] = param_max
			parameters_min[jj] = param_min

	return parameters_min, parameters_max, deltamu_min, deltamu_max


def write_deltamu(cosmology_string, parameter_sets, length_loop, f_name):
	
	#First we preallocate some memory
	distmod = np.zeros((len(redshifts),length_loop))
	deltamu = np.zeros((len(redshifts),length_loop))
	marg = np.zeros(length_loop)

    #For each parameter set, we define a cosmology and calculate distance moduli.
    #These are then used to calculate deltamu, and if marginalising is chosen
    #the magninalisation constants are calculated and stored.

    #Need to get the max/min deltamu values for each redshift, and parameter set
    #that gave that value
	print "Calculating deltamu values..."
	for ii in np.arange(length_loop):
		#print ii
		exec(cosmology_string)
		distmod[:,ii] = cosmology.distmod(redshifts)
		deltamu[:,ii] = distmod[:,ii] - np.array(distmod_bestfit_lcdm)
		
		if do_marg:
			m_marg = opt.fmin(chi2,x0=0.,args=(redshifts, distmod[:,ii], redshifts_marg, distmod_lcdm_marg, sigma_marg),xtol=0.005,disp=0)
			marg[ii] = m_marg

	print deltamu[:,ii][:10]
    
	print "Getting min/max parameters and deltamus"
	parameters_min, parameters_max, deltamu_min, deltamu_max = get_parameters_minmax(parameter_sets, deltamu)
    
	print "Dumping data to file"
	data_dump = redshifts, deltamu_min, deltamu_max, parameters_min, parameters_max
	data_file_name = f_name
	output = open(data_file_name,'wb')
	pick.dump(data_dump,output)
	output.close()

	if do_marg:
		parameters_min, parameters_max, deltamu_min, deltamu_max = get_parameters_minmax(parameter_sets, deltamu, marg)
		data_dump = redshifts, deltamu_min, deltamu_max, parameters_min, parameters_max
		data_file_name = f_name + '_marg'
		output = open(data_file_name,'wb')
		pick.dump(data_dump,output)
		output.close()

	print np.shape(parameters_max)
	print np.shape(parameters_min)
		
def define_redshifts_marg(zmin=0.1, zmax=10., nz=10, sigma=0.1):
	redshifts_marg = np.linspace(zmin, zmax, nz)
	sigma_marg = np.ones(len(redshifts_marg)) * sigma
	return redshifts_marg, sigma_marg

def define_contours(lcdm_contour_level = 0.68, lcdm_tolerance = 0.01, lcdm_bins_tuple=100,cpl_contour_level = 0.68, cpl_tolerance = 0.005, cpl_bins_tuple=(50,50,50),jbp_contour_level = 0.68, jbp_tolerance = 0.005, jbp_bins_tuple=(50,50,50)):
	#Set contours we need and check that they exist
	LCDM_Contour = Contour.LCDM_Contour(chain_name='lcdm',directory='/Users/perandersen/Data/HzSC/',contour_level=lcdm_contour_level,tolerance=lcdm_tolerance,bins_tuple=lcdm_bins_tuple)
	print "Does LCDM Contour exist?",
	if LCDM_Contour.test_contour_exists():
		print "Yes"
	else:
		print "No"
		LCDM_Contour.pickle_contour()

	CPL_Contour = Contour.Contour(chain_name='cpl',directory='/Users/perandersen/Data/HzSC/',contour_level=cpl_contour_level,tolerance=cpl_tolerance,bins_tuple=cpl_bins_tuple)
	print "Does CPL Contour exist?",
	if CPL_Contour.test_contour_exists():
		print "Yes"
	else:
		print "No"
		CPL_Contour.pickle_contour()
	JBP_Contour = Contour.Contour(chain_name='jbp',directory='/Users/perandersen/Data/HzSC/',contour_level=jbp_contour_level,tolerance=jbp_tolerance,bins_tuple=jbp_bins_tuple)
	print "Does JBP Contour exist?",
	if JBP_Contour.test_contour_exists():
		print "Yes"
	else:
		print "No"
		JBP_Contour.pickle_contour()

	#Read in contours
	print "Reading in contours...",
	lcdm_parameter_sets, lcdm_contour_level, lcdm_tolerance, lcdm_bins_tuple = LCDM_Contour.read_pickled_contour()

	cpl_omega_contour, cpl_w0_contour, cpl_wa_contour, cpl_contour_level, cpl_tolerance, cpl_bins_tuple = CPL_Contour.read_pickled_contour()
	cpl_parameter_sets = cpl_omega_contour, cpl_w0_contour, cpl_wa_contour

	jbp_omega_contour, jbp_w0_contour, jbp_wa_contour, jbp_contour_level, jbp_tolerance, jbp_bins_tuple = JBP_Contour.read_pickled_contour()
	jbp_parameter_sets = jbp_omega_contour, jbp_w0_contour, jbp_wa_contour
	print "Done!"
	return lcdm_parameter_sets, lcdm_contour_level, lcdm_tolerance, lcdm_bins_tuple, cpl_parameter_sets, cpl_contour_level, cpl_tolerance, cpl_bins_tuple, jbp_parameter_sets, jbp_contour_level, jbp_tolerance, jbp_bins_tuple


lcdm_parameter_sets,lcdm_contour_level,lcdm_tolerance,lcdm_bins_tuple,cpl_parameter_sets,cpl_contour_level,cpl_tolerance,cpl_bins_tuple,jbp_parameter_sets,jbp_contour_level,jbp_tolerance,jbp_bins_tuple = define_contours(jbp_bins_tuple=(20,20,20))



#Define the best fit LCDM parameter
omega_m_lcdm_bestfit, hubble_const = 0.30754277645, 70.

#Decide if we want to do marginalization run too (time consuming)
do_marg = True

#Define redshift range
#(Needs to be dense enough to get fine detail at lower redshifts)
redshifts = np.linspace(0.001, 10, 1000)

#Create redshift range to marginalise with, need to experiment with this
if do_marg:
	redshifts_marg, sigma_marg = define_redshifts_marg()


#Defining cosmology for best fitting parameters and derive distance moduli
#for the redshift range we look at
lcdm_bestfit = cosmo.FlatLambdaCDM(H0=hubble_const, Om0=omega_m_lcdm_bestfit)
distmod_bestfit_lcdm = lcdm_bestfit.distmod(redshifts)

if do_marg:    
    #Distance moduli for marginalising over
    distmod_lcdm_marg = np.array(lcdm_bestfit.distmod(redshifts_marg))

print "Function begin"
t0 = ti.time()
lcdm_cosmo_string = 'cosmology=cosmo.FlatLambdaCDM(H0=hubble_const, Om0=lcdm_parameter_sets[ii])'
cpl_cosmo_string = 'cosmology=cosmo.Flatw0waCDM(H0=hubble_const, Om0=cpl_parameter_sets[0][ii], w0=cpl_parameter_sets[1][ii], wa=cpl_parameter_sets[2][ii])'
jbp_cosmo_string = 'cosmology=FlatJBP.FlatJBP_CDM(H0=hubble_const, Om0=jbp_parameter_sets[0][ii], w0=jbp_parameter_sets[1][ii], wa=jbp_parameter_sets[2][ii])'

root_dir = '/Users/perandersen/Data/HzSC/Deltamu/'
lcdm_fname = "deltamu_lcdm_c" + str(lcdm_contour_level) +\
"_t" + str(lcdm_tolerance) + "_b" + str(lcdm_bins_tuple) + ".dat"

cpl_fname = "deltamu_cpl_c" + str(cpl_contour_level) +\
"_t" + str(cpl_tolerance) + "_b" + str(cpl_bins_tuple[0]) + \
str(cpl_bins_tuple[1]) + str(cpl_bins_tuple[2]) + ".dat"

jbp_fname = "deltamu_jbp_c" + str(jbp_contour_level) +\
"_t" + str(jbp_tolerance) + "_b" + str(jbp_bins_tuple[0]) + \
str(jbp_bins_tuple[1]) + str(jbp_bins_tuple[2]) + ".dat"

write_deltamu(lcdm_cosmo_string, lcdm_parameter_sets, len(lcdm_parameter_sets), root_dir + lcdm_fname)
write_deltamu(cpl_cosmo_string, cpl_parameter_sets, np.shape(cpl_parameter_sets)[1], root_dir + cpl_fname)
write_deltamu(jbp_cosmo_string, jbp_parameter_sets, np.shape(jbp_parameter_sets)[1], root_dir + jbp_fname)

print "Function end in:", ti.time() - t0, "seconds"

#----Get the max/min deltamu values and the parameter sets, as function of redshift

#----Save values, redshifts, and parameter sets to file


