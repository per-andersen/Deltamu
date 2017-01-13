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

def unique_par_sets_1d(parameters, redshifts):
    unique_par_sets = np.array([])
    unique_par_redshifts = np.array([])
    for ii in np.arange(len(parameters)):
        pars = parameters[ii]
        redshift = redshifts[ii]
        if pars in unique_par_sets:
            pass
        else:
            unique_par_sets = np.append(unique_par_sets,pars)
            unique_par_redshifts = np.append(unique_par_redshifts, redshift)

    return unique_par_sets, unique_par_redshifts

def unique_par_sets_3d(parameters, redshifts):
    unique_par_sets = []
    unique_par_redshifts = np.array([])
    for ii in np.arange(len(parameters)):
        pars = parameters[ii]
        redshift = redshifts[ii]
        keep = True
        for accepted_sets in unique_par_sets:
            if np.abs(np.sum(pars - accepted_sets)) < 0.0001:
                keep = False
        if keep:
            unique_par_sets.append(pars)
            unique_par_redshifts = np.append(unique_par_redshifts,redshift)
    return unique_par_sets, unique_par_redshifts

def get_unique_parameter_sets(fname):
    redshifts, deltamu_min, deltamu_max, parameters_min, parameters_max, deltamu, marg = read_pickled_deltamu(fname)
    if len(np.shape(parameters_min)) == 1:
        unique_par_sets_min, unique_par_redshifts_min = unique_par_sets_1d(parameters_min, redshifts)
        unique_par_sets_max, unique_par_redshifts_max = unique_par_sets_1d(parameters_max, redshifts)
    elif len(np.shape(parameters_min)) == 2:
        unique_par_sets_min, unique_par_redshifts_min = unique_par_sets_3d(parameters_min, redshifts)
        unique_par_sets_max, unique_par_redshifts_max = unique_par_sets_3d(parameters_max, redshifts)
    #for ii in np.arange(len(unique_par_sets_min)):
    #    print unique_par_sets_min[ii], unique_par_redshifts_min[ii]
    #for ii in np.arange(len(unique_par_sets_max)):
    #    print unique_par_sets_max[ii], unique_par_redshifts_max[ii]
    return unique_par_sets_min, unique_par_sets_max, unique_par_redshifts_min, unique_par_redshifts_max

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
    for ii in np.arange(np.shape(deltamu)[1]):
        plt.plot(redshifts, deltamu[:,ii],label=legend_string)
    plt.title(title)
    plt.show()

def oplot_deltamus(chain_name, bins, smoothings, tolerance=0.005):
    plt.figure()
    plt.ylim((-0.1,0.1))
    plt.xlabel('Redshift')
    plt.ylabel('Deltamu')
    for ii in np.arange(len(bins)):
        for jj in np.arange(len(smoothings)):
            deltamus = Deltamu.Deltamu(chain_name,'',do_marg=True,bins_tuple=bins[ii],smoothing=smoothings[jj],tolerance=tolerance)
            fname = root_dir + deltamus.get_marg_file_name()
            redshifts, deltamu_min, deltamu_max, parameters_min, parameters_max, deltamu, marg = read_pickled_deltamu(fname)
            for kk in np.arange(np.shape(deltamu)[1]):
                #plt.plot(redshifts, deltamu[:,kk])
                plt.plot(redshifts, deltamu[:,kk] + marg[kk] + 19.14820557)
    plt.title(chain_name)

def plot_3d_contours(chain_name, bins, smoothing, tolerance=0.005,labels = ['x','y','z']):
    redshift_cut = 0.4

    fig_scatter = plt.figure()
    ax_scatter = fig_scatter.add_subplot(111, projection='3d')
    ax_scatter.set_xlabel(labels[0])
    ax_scatter.set_ylabel(labels[1])
    ax_scatter.set_zlabel(labels[2])

    contour = Contour.Contour(chain_name=chain_name,directory='/Users/perandersen/Data/HzSC/',bins_tuple=bins[0],tolerance=tolerance,smoothing=smoothing)
    try:
        x_contour, y_contour, z_contour = contour.read_pickled_contour()
    except:
        contour.pickle_contour()
        x_contour, y_contour, z_contour = contour.read_pickled_contour()

    ax_scatter.scatter(x_contour, y_contour, z_contour,color='b',s=1.,depthshade=False)
    for ii in np.arange(len(bins)):
        deltamus = Deltamu.Deltamu(chain_name,'',do_marg=True,bins_tuple=bins[ii],smoothing=smoothing,tolerance=tolerance)
        fname = root_dir + deltamus.get_marg_file_name()
        unique_par_sets_min, unique_par_sets_max, unique_par_redshifts_min, unique_par_redshifts_max = get_unique_parameter_sets(fname)
        unique_par_sets_min = np.array(unique_par_sets_min)
        unique_par_sets_max = np.array(unique_par_sets_max)

        unique_par_sets_min = unique_par_sets_min[unique_par_redshifts_min > redshift_cut]
        unique_par_sets_max = unique_par_sets_max[unique_par_redshifts_max > redshift_cut]

        unique_par_redshifts_min = unique_par_redshifts_min[unique_par_redshifts_min > redshift_cut]
        unique_par_redshifts_max = unique_par_redshifts_max[unique_par_redshifts_max > redshift_cut]

        print unique_par_sets_min
        print unique_par_redshifts_min

        om_min = unique_par_sets_min[:,0]
        w0_min = unique_par_sets_min[:,1]
        wa_min = unique_par_sets_min[:,2]
        om_max = unique_par_sets_max[:,0]
        w0_max = unique_par_sets_max[:,1]
        wa_max = unique_par_sets_max[:,2]

        color_min = np.zeros((len(unique_par_redshifts_min),3))
        for jj in np.arange(len(color_min)):
            color_min[jj,0] = 1.
            color_min[jj,1] = float(jj) / len(color_min)
            color_min[jj,2] = float(jj) / len(color_min)

        color_min = color_min[::-1]

        color_max = np.zeros((len(unique_par_redshifts_max),3))
        for jj in np.arange(len(color_max)):
            color_max[jj,0] = 1.
            color_max[jj,1] = float(jj) / len(color_max)
            color_max[jj,2] = float(jj) / len(color_max)
            #print jj, color_max[jj,0], color_max[jj,1], color_max[jj,2]

        color_max = color_max[::-1]

        
        ax_scatter.scatter(om_min[:], w0_min[:], wa_min[:], color=color_min,s=100.,depthshade=False, edgecolor='k')
        ax_scatter.scatter(om_max[:], w0_max[:], wa_max[:], color=color_max,s=100.,depthshade=False, marker='^', edgecolor='k')

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

deltamu_cpl = Deltamu.Deltamu('cpl','',do_marg=True,bins_tuple=(50,50,50),smoothing=0.6)
cpl_marg_fname = deltamu_cpl.get_marg_file_name()

plot_3d_contours('jbp', [(50,50,50)], 0.6)


#oplot_deltamus('n3cpl', [(30,30,30),(40,40,40),(50,50,50)],[0.6])
#oplot_deltamus('n3cpl', [(30,30,30),(40,40,40),(50,50,50)],[0.0])
#oplot_deltamus('jbp', [(30,30,30),(40,40,40),(50,50,50)],[0.6])
#oplot_deltamus('cpl', [(30,30,30),(40,40,40),(50,50,50)],[0.6])
#oplot_deltamus('lcdm', [70,80,90,100],[0.6],tolerance=0.01)
plt.show()

