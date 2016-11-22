import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import astropy.cosmology as cosmo
import CosmoMC_dao
import Contour
import FlatJBP
import Flatn3CPL
import Flatn7CPL

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
lcdm_parameter_sets = LCDM_Contour.read_pickled_contour()
cpl_parameter_sets = CPL_Contour.read_pickled_contour()
jbp_parameter_sets = JBP_Contour.read_pickled_contour()
print "Done!"

#Define the best fit LCDM parameter
omega_m_lcdm_bestfit = 0.31
hubble_const = 70.
do_marg = True

#Define redshift range
redshifts = np.linspace(0.001, 10, 1000)

#Calculate best fit LCDM distance modulus values
lcdm_bestfit = cosmo.FlatLambdaCDM(H0=hubble_const, Om0=omega_m_lcdm_bestfit)
distance_mod_bestfit_lcdm = lcdm_bestfit.distmod(redshifts)

#For all parameter sets on contours, calculate distance modulus
#and unmarginalised distance moduli
distance_mod_lcdm = np.zeros((len(redshifts),len(lcdm_parameter_sets)))
distance_mod_cpl = np.zeros((len(redshifts),np.shape(cpl_parameter_sets)[1]))
distance_mod_jbp = np.zeros((len(redshifts),np.shape(jbp_parameter_sets)[1]))

deltamu_lcdm = np.zeros((len(redshifts),len(lcdm_parameter_sets)))
deltamu_cpl = np.zeros((len(redshifts),np.shape(cpl_parameter_sets)[1]))
deltamu_jbp = np.zeros((len(redshifts),np.shape(jbp_parameter_sets)[1]))

if do_marg:
    deltamu_lcdm_marg = np.zeros((len(redshifts),len(lcdm_parameter_sets)))
    deltamu_cpl_marg = np.zeros((len(redshifts),np.shape(cpl_parameter_sets)[1]))
    deltamu_jbp_marg = np.zeros((len(redshifts),np.shape(jbp_parameter_sets)[1]))
    
    marg_lcdm = np.zeros(len(lcdm_parameter_sets))
    marg_cpl = np.zeros(np.shape(cpl_parameter_sets)[1])
    marg_jbp = np.zeros(np.shape(jbp_parameter_sets)[1])
    
    redshifts_marg = np.linspace(0.1, 10, 11)
    sigma_marg = np.ones(len(redshifts_marg)) * 0.1
    distmod_lcdm_marg = np.array(lcdm_bestfit.distmod(redshifts_marg))

for ii in np.arange(len(lcdm_parameter_sets)):
	lcdm = cosmo.FlatLambdaCDM(H0=hubble_const, Om0=lcdm_parameter_sets[ii])
	distance_mod_lcdm[:,ii] = lcdm.distmod(redshifts)
	deltamu_lcdm[:,ii] = distance_mod_lcdm[:,ii] - np.array(distance_mod_bestfit_lcdm)
	if do_marg:
		chi2 = lambda m: np.sum( (lin_interp_array(redshifts,distance_mod_lcdm[:,ii],redshifts_marg) + m  - distmod_lcdm_marg )**2 / sigma_marg**2 )
		#print distmod_lcdm_marg
		#print lin_interp_array(redshifts,distance_mod_lcdm[:,ii],redshifts_marg)
		m_marg = opt.fmin(chi2,0.,xtol=0.005,disp=0)
		print m_marg, lcdm_parameter_sets[ii]
		marg_lcdm[ii] = m_marg
		#plt.figure()
		#plt.plot(redshifts_marg, distmod_lcdm_marg + m_marg,'rx')
		#plt.plot(redshifts,distance_mod_lcdm[:,ii],'b')
		#plt.show()


for ii in np.arange(np.shape(cpl_parameter_sets)[1]):
	cpl = cosmo.Flatw0waCDM(H0=hubble_const, Om0=cpl_parameter_sets[0][ii], w0=cpl_parameter_sets[1][ii], wa=cpl_parameter_sets[2][ii])
	distance_mod_cpl[:,ii] = cpl.distmod(redshifts)
	deltamu_cpl[:,ii] = distance_mod_cpl[:,ii] - np.array(distance_mod_bestfit_lcdm)
	if do_marg:
		distmod_cpl_marg = np.array(cpl.distmod(redshifts_marg))
		chi2 = lambda m: np.sum( (lin_interp_array(redshifts,distance_mod_cpl[:,ii],redshifts_marg) + m  - distmod_lcdm_marg )**2 / sigma_marg**2 )
		m_marg = opt.fmin(chi2,0.,xtol=0.005,disp=0)
		print m_marg, cpl_parameter_sets[0][ii],cpl_parameter_sets[1][ii],cpl_parameter_sets[2][ii]
		marg_cpl[ii] = m_marg
		#plt.figure()
		#plt.plot(redshifts_marg, distmod_lcdm_marg / distmod_lcdm_marg,'rx')
		#plt.plot(redshifts,(distance_mod_cpl[:,ii] + m_marg) / distance_mod_bestfit_lcdm,'b')
		#plt.show()

print np.shape(jbp_parameter_sets)[1]
for ii in np.arange(np.shape(jbp_parameter_sets)[1]):
	jbp = FlatJBP.FlatJBP_CDM(H0=hubble_const, Om0=jbp_parameter_sets[0][ii], w0=jbp_parameter_sets[1][ii], wa=jbp_parameter_sets[2][ii])
	distance_mod_jbp[:,ii] = jbp.distmod(redshifts)
	deltamu_jbp[:,ii] = distance_mod_jbp[:,ii] - np.array(distance_mod_bestfit_lcdm)
	if do_marg:
		distmod_jbp_marg = np.array(jbp.distmod(redshifts_marg))
		chi2 = lambda m: np.sum( (lin_interp_array(redshifts,distance_mod_jbp[:,ii],redshifts_marg) + m  - distmod_lcdm_marg )**2 / sigma_marg**2 )
		m_marg = opt.fmin(chi2,0.,xtol=0.005,disp=0)
		print m_marg, jbp_parameter_sets[0][ii],jbp_parameter_sets[1][ii],jbp_parameter_sets[2][ii]
		marg_jbp[ii] = m_marg
		#plt.figure()
		#plt.plot(redshifts_marg, distmod_lcdm_marg / distmod_lcdm_marg,'rx')
		#plt.plot(redshifts,(distance_mod_jbp[:,ii] + m_marg) / distance_mod_bestfit_lcdm,'b')
		#plt.plot(redshifts,(distance_mod_jbp[:,ii]) / distance_mod_bestfit_lcdm,'g')
		#plt.show()


#Calculate marginalising M offset for all parameter sets

#For both marginalised and unmarginalised distance moduli do the following

#----For all parameter sets on contours, calculate deltamu

#----Get the max/min deltamu values and the parameter sets, as function of redshift

#----Save values, redshifts, and parameter sets to file


