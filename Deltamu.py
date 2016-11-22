import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
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

def write_deltamu(cosmology_string, parameter_sets, length_loop):
	distmod = np.zeros((len(redshifts),length_loop))
	deltamu = np.zeros((len(redshifts),length_loop))
	marg = np.zeros(length_loop)
	if do_marg:
		sigma_marg = np.ones(len(redshifts_marg)) * 0.1
	for ii in np.arange(length_loop):
		exec(cosmology_string)
		distmod[:,ii] = cosmology.distmod(redshifts)
		deltamu[:,ii] = distmod[:,ii] - np.array(distmod_bestfit_lcdm)
		
		if do_marg:
			#chi2 = lambda m: np.sum( (lin_interp_array(redshifts,distmod[:,ii],redshifts_marg) + m  - distmod_lcdm_marg )**2 / sigma_marg**2 )
			#print distmod_lcdm_marg
			#print lin_interp_array(redshifts,distmod[:,ii],redshifts_marg)
			m_marg = opt.fmin(chi2,x0=0.,args=(redshifts, distmod[:,ii], redshifts_marg, distmod_lcdm_marg, sigma_marg),xtol=0.005,disp=0)
			#print m_marg
			marg[ii] = m_marg
			#plt.figure()
			#plt.plot(redshifts_marg, distmod_lcdm_marg / distmod_lcdm_marg,'rx')
			#plt.plot(redshifts,(distmod[:,ii] + m_marg) / distmod_bestfit_lcdm,'b')
			#plt.plot(redshifts,(distmod[:,ii]) / distmod_bestfit_lcdm,'g')
			#plt.show()
		

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

if do_marg:
	redshifts_marg = np.linspace(0.1, 10, 11)
	sigma_marg = np.ones(len(redshifts_marg)) * 0.1

#Calculate best fit LCDM distance modulus values
lcdm_bestfit = cosmo.FlatLambdaCDM(H0=hubble_const, Om0=omega_m_lcdm_bestfit)
distmod_bestfit_lcdm = lcdm_bestfit.distmod(redshifts)

#For all parameter sets on contours, calculate distance modulus
#and unmarginalised distance moduli
distmod_lcdm = np.zeros((len(redshifts),len(lcdm_parameter_sets)))
distmod_cpl = np.zeros((len(redshifts),np.shape(cpl_parameter_sets)[1]))
distmod_jbp = np.zeros((len(redshifts),np.shape(jbp_parameter_sets)[1]))

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
    
    distmod_lcdm_marg = np.array(lcdm_bestfit.distmod(redshifts_marg))

print "Function begin"
t0 = ti.time()
write_deltamu('cosmology=cosmo.FlatLambdaCDM(H0=hubble_const, Om0=lcdm_parameter_sets[ii])', lcdm_parameter_sets, len(lcdm_parameter_sets))
test_len = 500
#write_deltamu('cosmology=cosmo.Flatw0waCDM(H0=hubble_const, Om0=cpl_parameter_sets[0][ii], w0=cpl_parameter_sets[1][ii], wa=cpl_parameter_sets[2][ii])', cpl_parameter_sets, np.shape(cpl_parameter_sets)[1])
#write_deltamu('cosmology=FlatJBP.FlatJBP_CDM(H0=hubble_const, Om0=jbp_parameter_sets[0][ii], w0=jbp_parameter_sets[1][ii], wa=jbp_parameter_sets[2][ii])', jbp_parameter_sets, np.shape(jbp_parameter_sets)[1])
#write_deltamu('cosmology=cosmo.Flatw0waCDM(H0=hubble_const, Om0=cpl_parameter_sets[0][ii], w0=cpl_parameter_sets[1][ii], wa=cpl_parameter_sets[2][ii])', cpl_parameter_sets, test_len)
#write_deltamu('cosmology=FlatJBP.FlatJBP_CDM(H0=hubble_const, Om0=jbp_parameter_sets[0][ii], w0=jbp_parameter_sets[1][ii], wa=jbp_parameter_sets[2][ii])', jbp_parameter_sets, test_len)


print "Function end in:", ti.time() - t0, "seconds"
print ""
print "Old begin"
t1 = ti.time()
for ii in np.arange(len(lcdm_parameter_sets)):
	lcdm = cosmo.FlatLambdaCDM(H0=hubble_const, Om0=lcdm_parameter_sets[ii])
	distmod_lcdm[:,ii] = lcdm.distmod(redshifts)
	deltamu_lcdm[:,ii] = distmod_lcdm[:,ii] - np.array(distmod_bestfit_lcdm)
	
	if do_marg:
		chi2 = lambda m: np.sum( (lin_interp_array(redshifts,distmod_lcdm[:,ii],redshifts_marg) + m  - distmod_lcdm_marg )**2 / sigma_marg**2 )
		#print distmod_lcdm_marg
		#print lin_interp_array(redshifts,distance_mod_lcdm[:,ii],redshifts_marg)
		m_marg = opt.fmin(chi2,0.,xtol=0.005,disp=0)
		#print m_marg
		marg_lcdm[ii] = m_marg
		#plt.figure()
		#plt.plot(redshifts_marg, distmod_lcdm_marg + m_marg,'rx')
		#plt.plot(redshifts,distance_mod_lcdm[:,ii],'b')
		#plt.show()

#for ii in np.arange(test_len):
for ii in np.arange(np.shape(cpl_parameter_sets)[1]):
	cpl = cosmo.Flatw0waCDM(H0=hubble_const, Om0=cpl_parameter_sets[0][ii], w0=cpl_parameter_sets[1][ii], wa=cpl_parameter_sets[2][ii])
	distmod_cpl[:,ii] = cpl.distmod(redshifts)
	deltamu_cpl[:,ii] = distmod_cpl[:,ii] - np.array(distmod_bestfit_lcdm)
	if do_marg:
		distmod_cpl_marg = np.array(cpl.distmod(redshifts_marg))
		chi2 = lambda m: np.sum( (lin_interp_array(redshifts,distmod_cpl[:,ii],redshifts_marg) + m  - distmod_lcdm_marg )**2 / sigma_marg**2 )
		m_marg = opt.fmin(chi2,0.,xtol=0.005,disp=0)
		#print m_marg
		marg_cpl[ii] = m_marg
		#plt.figure()
		#plt.plot(redshifts_marg, distmod_lcdm_marg / distmod_lcdm_marg,'rx')
		#plt.plot(redshifts,(distmod_cpl[:,ii] + m_marg) / distmod_bestfit_lcdm,'b')
		#plt.show()

#print np.shape(jbp_parameter_sets)[1]
#for ii in np.arange(test_len):
for ii in np.arange(np.shape(jbp_parameter_sets)[1]):
	jbp = FlatJBP.FlatJBP_CDM(H0=hubble_const, Om0=jbp_parameter_sets[0][ii], w0=jbp_parameter_sets[1][ii], wa=jbp_parameter_sets[2][ii])
	distmod_jbp[:,ii] = jbp.distmod(redshifts)
	deltamu_jbp[:,ii] = distmod_jbp[:,ii] - np.array(distmod_bestfit_lcdm)
	if do_marg:
		distmod_jbp_marg = np.array(jbp.distmod(redshifts_marg))
		chi2 = lambda m: np.sum( (lin_interp_array(redshifts,distmod_jbp[:,ii],redshifts_marg) + m  - distmod_lcdm_marg )**2 / sigma_marg**2 )
		m_marg = opt.fmin(chi2,0.,xtol=0.005,disp=0)
		#print m_marg, jbp_parameter_sets[0][ii],jbp_parameter_sets[1][ii],jbp_parameter_sets[2][ii]
		marg_jbp[ii] = m_marg
		#plt.figure()
		#plt.plot(redshifts_marg, distmod_lcdm_marg / distmod_lcdm_marg,'rx')
		#plt.plot(redshifts,(distance_mod_jbp[:,ii] + m_marg) / distance_mod_bestfit_lcdm,'b')
		#plt.plot(redshifts,(distance_mod_jbp[:,ii]) / distance_mod_bestfit_lcdm,'g')
		#plt.show()

print "Old end in:", ti.time() - t1, "seconds"
#Calculate marginalising M offset for all parameter sets

#For both marginalised and unmarginalised distance moduli do the following

#----For all parameter sets on contours, calculate deltamu

#----Get the max/min deltamu values and the parameter sets, as function of redshift

#----Save values, redshifts, and parameter sets to file


