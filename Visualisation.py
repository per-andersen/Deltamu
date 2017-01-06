import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import pickle as pick
import Contour
import Deltamu

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
    for ii in np.arange(len(unique_par_sets_min)):
        print unique_par_sets_min[ii]
    return unique_par_sets_min, unique_par_sets_max

def get_file_name(chain_name, bins_tuple, contour_level=0.68, tolerance=0.005, do_marg=False,\
    redshifts_marg_method='jla', redshift_marg_min=0.1,redshift_marg_max=10., redshift_marg_n=10):
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
    #redshifts, deltamu_min, deltamu_max, parameters_min, parameters_max = read_pickled_deltamu(fname)
    redshifts, deltamu_min, deltamu_max, parameters_min, parameters_max, deltamu, marg = read_pickled_deltamu(fname)
    #plt.figure()
    plt.plot(redshifts, deltamu_min,label=legend_string)
    plt.plot(redshifts, deltamu_max)
    plt.title(title)
    plt.show()

def plot_min_deltamu(fname, title, legend_string=None):
    redshifts, deltamu_min, deltamu_max, parameters_min, parameters_max, deltamu, marg = read_pickled_deltamu(fname)
    #plt.figure()
    plt.plot(redshifts, deltamu_min,label=legend_string)
    plt.title(title)

def plot_max_deltamu(fname, title, legend_string=None):
    redshifts, deltamu_min, deltamu_max, parameters_min, parameters_max, deltamu, marg = read_pickled_deltamu(fname)
    #plt.figure()
    plt.plot(redshifts, deltamu_max,label=legend_string)
    plt.title(title)

def plot_deltamu(fname, title, legend_string=None):
    redshifts, deltamu_min, deltamu_max, parameters_min, parameters_max, deltamu, marg = read_pickled_deltamu(fname)
    #plt.figure()
    print np.shape(deltamu)
    print np.shape(deltamu)[1]
    for ii in np.arange(np.shape(deltamu)[1]):
        plt.plot(redshifts, deltamu[:,ii],label=legend_string)
    plt.title(title)
    plt.show()

def oplot_deltamus(chain_name, bins,smoothings):
    plt.figure()
    for ii in np.arange(len(bins)):
        for jj in np.arange(len(smoothings)):
            deltamus = Deltamu.Deltamu(chain_name,'',do_marg=True,bins_tuple=bins[ii],smoothing=smoothings[jj])
            fname = root_dir + deltamus.get_marg_file_name()
            redshifts, deltamu_min, deltamu_max, parameters_min, parameters_max, deltamu, marg = read_pickled_deltamu(fname)
            for kk in np.arange(np.shape(deltamu)[1]):
                plt.plot(redshifts, deltamu[:,kk])
    plt.title(chain_name)


def oplot_3d_contours():
    '''
    This function tests if the contours produced with different binning
    and smoothing settings agree visually. It is really messy, and could
    use some cleaning up, if it is to be used more than just the one time.
    '''
    CPL_Contour = Contour.Contour(chain_name='cpl', directory='/Users/perandersen/Data/HzSC/',bins_tuple=(30,30,30),tolerance = 0.001, smoothing=0.6)
    #CPL_Contour_2 = Contour.Contour(chain_name='cpl', directory='/Users/perandersen/Data/HzSC/',bins_tuple=(40,40,40),tolerance = 0.001, smoothing=0.6)
    CPL_Contour_2 = Contour.Contour(chain_name='n3cpl', directory='/Users/perandersen/Data/HzSC/',bins_tuple=(60,60,60),tolerance = 0.001, smoothing=0.6)
    CPL_Contour_3 = Contour.Contour(chain_name='n3cpl', directory='/Users/perandersen/Data/HzSC/',bins_tuple=(50,50,50),tolerance = 0.001, smoothing=0.6)
    x_contour, y_contour, z_contour = CPL_Contour.read_pickled_contour()
    x_contour_2, y_contour_2, z_contour_2 = CPL_Contour_2.read_pickled_contour()
    x_contour_3, y_contour_3, z_contour_3 = CPL_Contour_3.read_pickled_contour()
    print len(x_contour)
    print len(x_contour_2)
    print len(x_contour_3)
    fig_scatter = plt.figure()
    ax_scatter = fig_scatter.add_subplot(111, projection='3d')
    #ax_scatter.scatter(x_contour, y_contour, z_contour, color='g')
    ax_scatter.scatter(x_contour_2, y_contour_2, z_contour_2)
    ax_scatter.scatter(x_contour_3, y_contour_3, z_contour_3, color='r')

root_dir = '/Users/perandersen/Data/HzSC/Deltamu/'

deltamu_cpl = Deltamu.Deltamu('cpl','',do_marg=True,bins_tuple=(60,60,60),smoothing=0.6)
cpl_marg_fname = deltamu_cpl.get_marg_file_name()
#plot_deltamu(root_dir + cpl_marg_fname,'cpl')
oplot_deltamus('cpl', [(60,60,60), (50,50,50), (40,40,40), (30,30,30)],[0.6])
oplot_deltamus('cpl', [(60,60,60), (50,50,50), (40,40,40), (30,30,30)],[0.2])
oplot_deltamus('cpl', [(60,60,60), (50,50,50), (40,40,40), (30,30,30)],[0.0])
plt.show()

