import numpy as np
import scipy.integrate as integrate

def dist_mod_marg(redshift, omega_matter, w0, wz):
    """Calculates the marginalised distance modulus, ignoring constant offsets
    """
    return 5.*np.log10( (1.+redshift)*metric_dist_unitless(redshift, omega_matter, w0, wz) )

def metric_dist_unitless(redshift, omega_matter, w0, wz):
    """The unitless metric distance
    """
    omega_lambda = 1. - omega_matter
    E_inv = lambda zz: ( omega_matter*(1.+zz)**3 + omega_lambda*f_w(zz,w0,wz) )**(-0.5)
    return integrate.quad(E_inv,0,redshift)[0]

def f_w(z,w0,wz):
   #For CPL cosmology 

	return np.exp(-(3.*wz*z)/(1.+z)) * (1.+z)**(3.*(1.+w0+wz))

#def f_w(z,w0,wz):
#	#For JBP cosmology
#	return np.exp(1.5*wz*z**2 / (1.+z)**2) * (1.+z)**(3.*(1+w0))

#def f_w(z,w0,wz):
#	#For n=3 nCPL cosmology
#	aa = 1. / (1. + z)
#	return np.exp(-0.5*wz*(1.-aa)*(2*aa**2 - 7*aa + 11)) * (1.+z)**(3.*(1+w0+wa))

#def f_w(z,w0,wz):
#	#For n=7 nCPL cosmology
#	aa = 1. / (1. + z)
#	return np.exp(-wz*(1.-aa)*(60*aa**6 - 430*aa**5 + 1334*aa**4 - 2341*aa**3 + 2559*aa**2 - 1851*aa + 1089)/140.) * (1.+z)**(3.*(1+w0+wa))

if __name__ == "__main__":

	'''Here we test the CPL cosmology module against astropy values
	!!!make sure to set the correct f_w function above!!!'''
	'''
	import astropy.cosmology as cosmo
	w0, wa = -0.9, -0.8
	hubble_const = 70.
	rz = 0.1
	o_m = 0.1
	o_x = 1. - o_m
	cpl = cosmo.w0waCDM(H0=hubble_const, Om0=o_m, Ode0=o_x, w0=w0, wa=wa)
	print "The below numbers should agree fairly well"
	print 2.9979e5*metric_dist_unitless(rz, o_m, w0, wa) * (1.+rz) / hubble_const, cpl.luminosity_distance(rz)

	o_m = 0.2
	o_x = 1. - o_m
	cpl = cosmo.w0waCDM(H0=hubble_const, Om0=o_m, Ode0=o_x, w0=w0, wa=wa)
	print 2.9979e5*metric_dist_unitless(rz, o_m, w0, wa) * (1.+rz) / hubble_const, cpl.luminosity_distance(rz)

	o_m = 0.9
	o_x = 1. - o_m
	cpl = cosmo.w0waCDM(H0=hubble_const, Om0=o_m, Ode0=o_x, w0=w0, wa=wa)
	print 2.9979e5*metric_dist_unitless(rz, o_m, w0, wa) * (1.+rz) / hubble_const, cpl.luminosity_distance(rz)

	rz = 5.5
	cpl = cosmo.w0waCDM(H0=hubble_const, Om0=o_m, Ode0=o_x, w0=w0, wa=wa)
	print 2.9979e5*metric_dist_unitless(rz, o_m, w0, wa) * (1.+rz) / hubble_const, cpl.luminosity_distance(rz)

	o_m = 0.2
	o_x = 1. - o_m
	cpl = cosmo.w0waCDM(H0=hubble_const, Om0=o_m, Ode0=o_x, w0=w0, wa=wa)
	print 2.9979e5*metric_dist_unitless(rz, o_m, w0, wa) * (1.+rz) / hubble_const, cpl.luminosity_distance(rz)

	o_m = 0.9
	o_x = 1. - o_m
	cpl = cosmo.w0waCDM(H0=hubble_const, Om0=o_m, Ode0=o_x, w0=w0, wa=wa)
	print 2.9979e5*metric_dist_unitless(rz, o_m, w0, wa) * (1.+rz) / hubble_const, cpl.luminosity_distance(rz)
	'''

	'''Here we test the JBP cosmology module against astropy values
	!!!make sure to set the correct f_w function above!!!'''
	
	'''
	import FlatJBP

	w0, wa = -0.5, -1.8
	hubble_const = 70.
	rz = 2.5
	o_m = 0.1
	jbp = FlatJBP.FlatJBP_CDM(H0=hubble_const, Om0=o_m, w0=w0, wa=wa)
	print "The below numbers should agree fairly well"
	print 2.9979e5*metric_dist_unitless(rz, o_m, w0, wa) * (1.+rz) / hubble_const, jbp.luminosity_distance(rz)
	print dist_mod_marg(rz, o_m, w0, wa) - dist_mod_marg(rz-1., o_m, w0, wa), jbp.distmod(rz) - jbp.distmod(rz-1)
	w0, wa = -0.9, -0.8
	rz = 5.5
	o_m = 0.1
	jbp = FlatJBP.FlatJBP_CDM(H0=hubble_const, Om0=o_m, w0=w0, wa=wa)
	print 2.9979e5*metric_dist_unitless(rz, o_m, w0, wa) * (1.+rz) / hubble_const, jbp.luminosity_distance(rz)
	print dist_mod_marg(rz, o_m, w0, wa) - dist_mod_marg(rz-1., o_m, w0, wa), jbp.distmod(rz) - jbp.distmod(rz-1)
	'''
	
	'''Here we test the n=3 nCDM cosmology module against astropy values
	!!!make sure to set the correct f_w function above!!!'''
	
	'''
	import Flatn3CPL

	w0, wa = -0.5, -1.8
	hubble_const = 70.
	rz = 2.5
	o_m = 0.1
	n3cpl = Flatn3CPL.Flatn3CPL(H0=hubble_const, Om0=o_m, w0=w0, wa=wa)
	print "The below numbers should agree fairly well"
	print 2.9979e5*metric_dist_unitless(rz, o_m, w0, wa) * (1.+rz) / hubble_const, n3cpl.luminosity_distance(rz)
	print dist_mod_marg(rz, o_m, w0, wa) - dist_mod_marg(rz-1., o_m, w0, wa), n3cpl.distmod(rz) - n3cpl.distmod(rz-1)
	
	w0, wa = -0.9, -0.8
	rz = 5.5
	o_m = 0.1
	n3cpl = Flatn3CPL.Flatn3CPL(H0=hubble_const, Om0=o_m, w0=w0, wa=wa)
	print 2.9979e5*metric_dist_unitless(rz, o_m, w0, wa) * (1.+rz) / hubble_const, n3cpl.luminosity_distance(rz)
	print dist_mod_marg(rz, o_m, w0, wa) - dist_mod_marg(rz-1., o_m, w0, wa), n3cpl.distmod(rz) - n3cpl.distmod(rz-1)
	'''

	'''Here we test the n=7 nCDM cosmology module against astropy values
	!!!make sure to set the correct f_w function above!!!'''
	
	'''
	import Flatn7CPL

	w0, wa = -0.5, -1.8
	hubble_const = 70.
	rz = 2.5
	o_m = 0.1
	n7cpl = Flatn7CPL.Flatn7CPL(H0=hubble_const, Om0=o_m, w0=w0, wa=wa)
	print "The below numbers should agree fairly well"
	print 2.9979e5*metric_dist_unitless(rz, o_m, w0, wa) * (1.+rz) / hubble_const, n7cpl.luminosity_distance(rz)
	print dist_mod_marg(rz, o_m, w0, wa) - dist_mod_marg(rz-1., o_m, w0, wa), n7cpl.distmod(rz) - n7cpl.distmod(rz-1)
	
	w0, wa = -0.9, -0.8
	rz = 5.5
	o_m = 0.1
	n7cpl = Flatn7CPL.Flatn7CPL(H0=hubble_const, Om0=o_m, w0=w0, wa=wa)
	print 2.9979e5*metric_dist_unitless(rz, o_m, w0, wa) * (1.+rz) / hubble_const, n7cpl.luminosity_distance(rz)
	print dist_mod_marg(rz, o_m, w0, wa) - dist_mod_marg(rz-1., o_m, w0, wa), n7cpl.distmod(rz) - n7cpl.distmod(rz-1)
	print np.array(dist_mod_marg(rz, o_m, w0, wa) - dist_mod_marg(rz-1., o_m, w0, wa)) - np.array(n7cpl.distmod(rz) - n7cpl.distmod(rz-1))
	'''


	'''Here we test the performance of functions against astropy for CPL'''
	
	import astropy.cosmology as cosmo
	import time as ti
	redshift_range = np.linspace(0.001, 10, 50)
	dist_mods_functions = np.zeros(len(redshift_range))

	o_m_array = np.linspace(0.2,0.4,1000)
	w0, wa = -0.5, 1.
	hubble_const = 70.
	
	print "Beginning"
	t0 = ti.time()
	for jj in np.arange(len(o_m_array)):
		o_m = o_m_array[jj]
		o_x = 1. - o_m
		for ii in np.arange(len(redshift_range)):
			dist_mods_functions[ii] = dist_mod_marg(redshift_range[ii], o_m, w0, wa)
	print "Functions done in", ti.time() - t0, "seconds"
	print dist_mods_functions[-10:] - dist_mods_functions[-11:-1]

	t0 = ti.time()
	for ii in np.arange(len(o_m_array)):
	    cpl = cosmo.w0waCDM(H0=hubble_const, Om0=o_m_array[ii], Ode0=o_x, w0=w0, wa=wa)
	    dist_mods_astropy = cpl.distmod(redshift_range)
	print "Astropy done in", ti.time() - t0, "seconds"

	print dist_mods_astropy[-10:] - dist_mods_astropy[-11:-1]

	print (dist_mods_functions[:10] - dist_mods_functions[1:11]) - np.array(dist_mods_astropy[:10] - dist_mods_astropy[1:11])
	

	'''Here we test the performance of functions against astropy for n=3 nCPL'''
	'''
	import Flatn3CPL
	import time as ti
	redshift_range = np.linspace(0.001, 10, 30)
	dist_mods_functions = np.zeros(len(redshift_range))

	o_m_array = np.linspace(0.2,0.4,1000)
	w0, wa = -0.5, 1.
	hubble_const = 70.
	
	print "Beginning n=3"
	t0 = ti.time()
	for jj in np.arange(len(o_m_array)):
		o_m = o_m_array[jj]
		o_x = 1. - o_m
		for ii in np.arange(len(redshift_range)):
			dist_mods_functions[ii] = dist_mod_marg(redshift_range[ii], o_m, w0, wa)
	print "Functions done in", ti.time() - t0, "seconds"
	print dist_mods_functions[-10:] - dist_mods_functions[-11:-1]

	t0 = ti.time()
	for ii in np.arange(len(o_m_array)):
	    cpl = Flatn3CPL.Flatn3CPL(H0=hubble_const, Om0=o_m_array[ii], w0=w0, wa=wa)
	    dist_mods_astropy = cpl.distmod(redshift_range)
	print "Astropy done in", ti.time() - t0, "seconds"

	print dist_mods_astropy[-10:] - dist_mods_astropy[-11:-1]

	print (dist_mods_functions[:10] - dist_mods_functions[1:11]) - np.array(dist_mods_astropy[:10] - dist_mods_astropy[1:11])
	'''
	
	'''Here we test the performance of functions against astropy for n=7 nCPL'''
	'''
	import Flatn7CPL
	import time as ti
	redshift_range = np.linspace(0.001, 10, 20)
	dist_mods_functions = np.zeros(len(redshift_range))

	o_m_array = np.linspace(0.2,0.4,500)
	w0, wa = -0.5, 1.
	hubble_const = 70.
	
	print "Beginning n=7"
	t0 = ti.time()
	for jj in np.arange(len(o_m_array)):
		o_m = o_m_array[jj]
		o_x = 1. - o_m
		for ii in np.arange(len(redshift_range)):
			dist_mods_functions[ii] = dist_mod_marg(redshift_range[ii], o_m, w0, wa)
	print "Functions done in", ti.time() - t0, "seconds"
	print dist_mods_functions[-10:] - dist_mods_functions[-11:-1]

	t0 = ti.time()
	for ii in np.arange(len(o_m_array)):
	    cpl = Flatn7CPL.Flatn7CPL(H0=hubble_const, Om0=o_m_array[ii], w0=w0, wa=wa)
	    dist_mods_astropy = cpl.distmod(redshift_range)
	print "Astropy done in", ti.time() - t0, "seconds"

	print dist_mods_astropy[-10:] - dist_mods_astropy[-11:-1]

	print (dist_mods_functions[:10] - dist_mods_functions[1:11]) - np.array(dist_mods_astropy[:10] - dist_mods_astropy[1:11])
	'''

	'''Here we test the performance of functions against astropy for JBP'''
	'''
	import FlatJBP
	import time as ti
	redshift_range = np.linspace(0.001, 10, 15)
	dist_mods_functions = np.zeros(len(redshift_range))

	o_m_array = np.linspace(0.2,0.4,500)
	w0, wa = -0.5, 1.
	hubble_const = 70.
	
	print "Beginning JBP"
	t0 = ti.time()
	for jj in np.arange(len(o_m_array)):
		o_m = o_m_array[jj]
		o_x = 1. - o_m
		for ii in np.arange(len(redshift_range)):
			dist_mods_functions[ii] = dist_mod_marg(redshift_range[ii], o_m, w0, wa)
	print "Functions done in", ti.time() - t0, "seconds"
	print dist_mods_functions[-10:] - dist_mods_functions[-11:-1]

	t0 = ti.time()
	for ii in np.arange(len(o_m_array)):
	    jbp = FlatJBP.FlatJBP_CDM(H0=hubble_const, Om0=o_m_array[ii], w0=w0, wa=wa)
	    dist_mods_astropy = jbp.distmod(redshift_range)
	print "Astropy done in", ti.time() - t0, "seconds"

	print dist_mods_astropy[-10:] - dist_mods_astropy[-11:-1]

	print (dist_mods_functions[:10] - dist_mods_functions[1:11]) - np.array(dist_mods_astropy[:10] - dist_mods_astropy[1:11])
    '''



