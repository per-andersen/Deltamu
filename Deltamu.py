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

def write_deltamu(cosmology_string, parameter_sets, length_loop, f_name):
	
	#First we preallocate some memory
	distmod = np.zeros((len(redshifts),length_loop))
	deltamu = np.zeros((len(redshifts),length_loop))
	marg = np.zeros(length_loop)
	deltamu_max = np.zeros(len(redshifts))
	deltamu_min = np.zeros(len(redshifts))

	if do_marg:
		sigma_marg = np.ones(len(redshifts_marg)) * 0.1

    #For each parameter set, we define a cosmology and calculate distance moduli.
    #These are then used to calculate deltamu, and if marginalising is chosen
    #the magninalisation constants are calculated and stored.

    #Need to get the max/min deltamu values for each redshift, and parameter set
    #that gave that value

	for ii in np.arange(length_loop):
		exec(cosmology_string)
		distmod[:,ii] = cosmology.distmod(redshifts)
		deltamu[:,ii] = distmod[:,ii] - np.array(distmod_bestfit_lcdm)
		
		if do_marg:
			m_marg = opt.fmin(chi2,x0=0.,args=(redshifts, distmod[:,ii], redshifts_marg, distmod_lcdm_marg, sigma_marg),xtol=0.005,disp=0)
			marg[ii] = m_marg
			#plt.figure()
			#plt.plot(redshifts_marg, distmod_lcdm_marg / distmod_lcdm_marg,'rx')
			#plt.plot(redshifts,(distmod[:,ii] + m_marg) / distmod_bestfit_lcdm,'b')
			#plt.plot(redshifts,(distmod[:,ii]) / distmod_bestfit_lcdm,'g')
			#plt.show()
    
	if len(np.shape(parameter_sets)) == 2:
		parameters_max = np.zeros((len(redshifts),3))
		parameters_min = np.zeros((len(redshifts),3))
		for jj in np.arange(len(redshifts)):
			deltamu_max[jj] = np.max(deltamu[jj,:])
			deltamu_min[jj] = np.min(deltamu[jj,:])
			param_max = np.array([ parameter_sets[0][np.argmax(deltamu[jj,:])], parameter_sets[1][np.argmax(deltamu[jj,:])], parameter_sets[2][np.argmax(deltamu[jj,:])] ])
			param_min = np.array([ parameter_sets[0][np.argmin(deltamu[jj,:])], parameter_sets[1][np.argmin(deltamu[jj,:])], parameter_sets[2][np.argmin(deltamu[jj,:])] ])
			parameters_max[jj,:] = param_max
			parameters_min[jj,:] = param_min
	else:
		parameters_max = np.zeros(len(redshifts))
		parameters_min = np.zeros(len(redshifts))
		for jj in np.arange(len(redshifts)):
			deltamu_max[jj] = np.max(deltamu[jj,:])
			deltamu_min[jj] = np.min(deltamu[jj,:])
			param_max = np.array([parameter_sets[np.argmax(deltamu[jj,:])]])
			param_min = np.array([parameter_sets[np.argmin(deltamu[jj,:])]])
		parameters_max[jj] = param_max
		parameters_min[jj] = param_min
    
	print np.shape(parameters_max)
	print np.shape(parameters_min)
	
	data_dump = (redshifts, deltamu_min, deltamu_max, parameters_min, parameters_max)
	data_file_name = f_name
	output = open(data_file_name,'wb')
	pick.dump(data_dump,output)
	output.close()

	#ToDo:
	#----Write the above redshifts, min/max deltamu values, min/max parameters to pickle file
	#----Do the same process for marginalised deltamu values, and store those in a pickle file as well
		

#Set contours we need and check that they exist
LCDM_Contour = Contour.LCDM_Contour(chain_name='lcdm',directory='/Users/perandersen/Data/HzSC/')
print "Does LCDM Contour exist?",
if LCDM_Contour.test_contour_exists():
	print "Yes"
else:
	print "No"
	LCDM_Contour.pickle_contour()

CPL_Contour = Contour.Contour(chain_name='cpl',directory='/Users/perandersen/Data/HzSC/')
print "Does CPL Contour exist?",
if CPL_Contour.test_contour_exists():
	print "Yes"
else:
	print "No"
	CPL_Contour.pickle_contour()
JBP_Contour = Contour.Contour(chain_name='jbp',directory='/Users/perandersen/Data/HzSC/', bins_tuple=(20,20,20))
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

#Define the best fit LCDM parameter
omega_m_lcdm_bestfit, hubble_const = 0.31, 70.

#Decide if we want to do marginalization run too (time consuming)
do_marg = False

#Define redshift range
#(Needs to be dense enough to get fine detail at lower redshifts)
redshifts = np.linspace(0.001, 10, 1000)

#Create redshift range to marginalise with, need to experiment with this
if do_marg:
	redshifts_marg = np.linspace(0.1, 10, 11)
	sigma_marg = np.ones(len(redshifts_marg)) * 0.1


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
write_deltamu(lcdm_cosmo_string, lcdm_parameter_sets, len(lcdm_parameter_sets), root_dir + 'lcdm.pkl')
#write_deltamu(cpl_cosmo_string, cpl_parameter_sets, np.shape(cpl_parameter_sets)[1])
#write_deltamu(jbp_cosmo_string, jbp_parameter_sets, np.shape(jbp_parameter_sets)[1])
print ""
write_deltamu(cpl_cosmo_string, cpl_parameter_sets, 20, root_dir + 'cpl.pkl')
print ""
write_deltamu(jbp_cosmo_string, jbp_parameter_sets, 20, root_dir + 'jbp.pkl')

print "Function end in:", ti.time() - t0, "seconds"

#----Get the max/min deltamu values and the parameter sets, as function of redshift

#----Save values, redshifts, and parameter sets to file


