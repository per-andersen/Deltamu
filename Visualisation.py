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
            if np.abs(np.sum(pars - accepted_sets)) < 0.0001:
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
    return unique_par_sets_min, unique_par_sets_max

def get_file_name(chain_name, bins_tuple, contour_level=0.68, tolerance=0.005, do_marg=False,\
    redshifts_marg_method='lin', redshift_marg_min=0.1,redshift_marg_max=10., redshift_marg_n=10):
    if chain_name == 'lcdm':
        f_name = "deltamu_lcdm_c" + str(contour_level) +\
        "_t" + str(tolerance) + "_b" + str(bins_tuple) + ".dat"
    else:
        f_name = "deltamu_" + chain_name + "_c" + str(contour_level) +\
        "_t" + str(tolerance) + "_b" + str(bins_tuple[0]) + \
        str(bins_tuple[1]) + str(bins_tuple[2]) + ".dat"

    if do_marg:
        if chain_name == 'lcdm':
            f_name = "deltamu_lcdm_c" + str(contour_level) +\
            "_t" + str(tolerance) + "_b" + str(bins_tuple) + "_marg_" +\
            redshifts_marg_method +"_z" + str(redshift_marg_min) +\
            "-" + str(redshift_marg_max) + "_n" + str(redshift_marg_n) + ".dat"
        else:
            f_name = "deltamu_" + chain_name + "_c" + str(contour_level) +\
            "_t" + str(tolerance) + "_b" + str(bins_tuple[0]) + \
            str(bins_tuple[1]) + str(bins_tuple[2]) + "_marg_" +\
            redshifts_marg_method +"_z" + str(redshift_marg_min) +\
            "-" + str(redshift_marg_max) + "_n" + str(redshift_marg_n) + ".dat"
    return f_name

def plot_minmax_deltamu(fname, title, legend_string=None):
    redshifts, deltamu_min, deltamu_max, parameters_min, parameters_max = read_pickled_deltamu(fname)
    #plt.figure()
    plt.plot(redshifts, deltamu_min,label=legend_string)
    plt.plot(redshifts, deltamu_max)
    plt.title(title)

def plot_min_deltamu(fname, title, legend_string=None):
    redshifts, deltamu_min, deltamu_max, parameters_min, parameters_max = read_pickled_deltamu(fname)
    #plt.figure()
    plt.plot(redshifts, deltamu_min,label=legend_string)
    plt.title(title)

def plot_max_deltamu(fname, title, legend_string=None):
    redshifts, deltamu_min, deltamu_max, parameters_min, parameters_max = read_pickled_deltamu(fname)
    #plt.figure()
    plt.plot(redshifts, deltamu_max,label=legend_string)
    plt.title(title)

root_dir = '/Users/perandersen/Data/HzSC/Deltamu/'

lcdm_fname = get_file_name(chain_name='lcdm',bins_tuple=100,tolerance=0.01)
lcdm_marg_fname = get_file_name(chain_name='lcdm',bins_tuple=100,tolerance=0.01,do_marg=True)
jbp_fname = get_file_name(chain_name='jbp',bins_tuple=(20,20,20))
jbp_marg_fname = get_file_name(chain_name='jbp',bins_tuple=(20,20,20),do_marg=True)
cpl_fname = get_file_name(chain_name='cpl',bins_tuple=(50,50,50))
cpl_marg_fname = get_file_name(chain_name='cpl',bins_tuple=(50,50,50),do_marg=True)

#get_unique_parameter_sets(root_dir + lcdm_fname)
#get_unique_parameter_sets(root_dir + cpl_fname)
#plot_minmax_deltamu(root_dir + lcdm_fname, 'lcdm')
#plot_minmax_deltamu(root_dir + cpl_fname, 'cpl')
#plot_minmax_deltamu(root_dir + jbp_fname, 'jbp')

#plot_minmax_deltamu(root_dir + lcdm_fname, 'lcdm')
#plot_minmax_deltamu(root_dir + lcdm_marg_fname, 'lcdm')
#plot_minmax_deltamu(root_dir + jbp_fname, 'jbp')
#plot_minmax_deltamu(root_dir + jbp_marg_fname, 'jbp')
#plot_minmax_deltamu(root_dir + cpl_fname, 'cpl')
#plot_minmax_deltamu(root_dir + cpl_marg_fname, 'cpl')

'''
cpl_bins_tuples = [(20,20,20),(30,30,30),(40,40,40),(50,50,50),(60,60,60),(70,70,70)]
plt.figure()
for cpl_bins in cpl_bins_tuples:
    cpl_fname = get_file_name(chain_name='cpl',bins_tuple=cpl_bins)
    plot_min_deltamu(root_dir + cpl_fname, 'cpl', legend_string=str(cpl_bins[0]))
plt.legend()

plt.figure()
for cpl_bins in cpl_bins_tuples:
    cpl_fname = get_file_name(chain_name='cpl',bins_tuple=cpl_bins)
    plot_max_deltamu(root_dir + cpl_fname, 'cpl', legend_string=str(cpl_bins[0]))
plt.legend()
'''
lcdm_bins_tuples = [50,70,90,110,120,130]
for lcdm_bins in lcdm_bins_tuples:
    lcdm_fname = get_file_name(chain_name='lcdm',bins_tuple=lcdm_bins)
    plot_min_deltamu(root_dir + lcdm_fname, 'lcdm', legend_string=str(lcdm_bins))
plt.legend()
plt.show()

