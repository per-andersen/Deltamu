import numpy as np
import matplotlib.pyplot as plt
import pickle as pick
import Contour

'''
Program to analyse and plot output from the Deltamu functions. 
ToDo:
----Read in pickled data
----Determine how many unique sets of parameters are needed to get min/max delta over redshift space
----Plot the min/max deltamu as a function of redshift
'''

def read_pickled_deltamu(fname):
    pkl_data_file = open(fname,'rb')
    data = pick.load(pkl_data_file)
    pkl_data_file.close()
    return data

def unique_par_sets_1d(parameters):
    unique_par_sets = np.array([])
    for pars in parameters:
        if pars in unique_par_sets:
            pass
        else:
            unique_par_sets = np.append(unique_par_sets,pars)
    return unique_par_sets

def unique_par_sets_3d(parameters):
    unique_par_sets = []
    for pars in parameters:
        keep = True
        for accepted_sets in unique_par_sets:
            if np.abs(np.sum(pars - accepted_sets)) < 0.00001:
                keep = False
        if keep:
            unique_par_sets.append(pars)
    return unique_par_sets

def get_unique_parameter_sets(fname):
    redshifts, deltamu_min, deltamu_max, parameters_min, parameters_max = read_pickled_deltamu(fname)
    if len(np.shape(parameters_min)) == 1:
        unique_par_sets_min = unique_par_sets_1d(parameters_min)
        unique_par_sets_max = unique_par_sets_1d(parameters_max)
    elif len(np.shape(parameters_min)) == 2:
        unique_par_sets_min = unique_par_sets_3d(parameters_min)
        unique_par_sets_max = unique_par_sets_3d(parameters_max)




LCDM_Contour = Contour.LCDM_Contour(chain_name='lcdm',directory='/Users/perandersen/Data/HzSC/')
CPL_Contour = Contour.Contour(chain_name='cpl',directory='/Users/perandersen/Data/HzSC/')
JBP_Contour = Contour.Contour(chain_name='jbp',directory='/Users/perandersen/Data/HzSC/', bins_tuple=(20,20,20))

lcdm_parameter_sets, lcdm_contour_level, lcdm_tolerance, lcdm_bins_tuple = LCDM_Contour.read_pickled_contour()
cpl_omega_contour, cpl_w0_contour, cpl_wa_contour, cpl_contour_level, cpl_tolerance, cpl_bins_tuple = CPL_Contour.read_pickled_contour()
jbp_omega_contour, jbp_w0_contour, jbp_wa_contour, jbp_contour_level, jbp_tolerance, jbp_bins_tuple = JBP_Contour.read_pickled_contour()

root_dir = '/Users/perandersen/Data/HzSC/Deltamu/'
lcdm_fname = "deltamu_lcdm_c" + str(lcdm_contour_level) +\
"_t" + str(lcdm_tolerance) + "_b" + str(lcdm_bins_tuple) + ".dat"

cpl_fname = "deltamu_cpl_c" + str(cpl_contour_level) +\
"_t" + str(cpl_tolerance) + "_b" + str(cpl_bins_tuple[0]) + \
str(cpl_bins_tuple[1]) + str(cpl_bins_tuple[2]) + ".dat"

jbp_fname = "deltamu_jbp_c" + str(jbp_contour_level) +\
"_t" + str(jbp_tolerance) + "_b" + str(jbp_bins_tuple[0]) + \
str(jbp_bins_tuple[1]) + str(jbp_bins_tuple[2]) + ".dat"

get_unique_parameter_sets(root_dir + lcdm_fname)
get_unique_parameter_sets(root_dir + cpl_fname)