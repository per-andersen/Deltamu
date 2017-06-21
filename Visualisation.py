import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.patheffects as patheffects
from matplotlib.font_manager import FontProperties
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
    redshifts, deltamu_min, deltamu_max, parameters_min, parameters_max, deltamu, marg, m_bestfit_lcdm_marg = read_pickled_deltamu(fname)
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

def eos_from_chain_name(chain_name, w0, wa, redshifts):
    scale_factor = 1. / (1. + redshifts)
    if chain_name == 'cpl':
        #print 'CPL'
        return w0 + wa * (1. - scale_factor)
    elif chain_name == 'jbp':
        #print 'JBP'
        return w0 + wa * (1. - scale_factor) * scale_factor
    elif chain_name == 'n3cpl':
        #print 'n3CPL'
        return w0 + wa * (1. - scale_factor)**3
    elif chain_name == 'n7cpl':
        #print 'n7CPL'
        return w0 + wa * (1. - scale_factor)**7
    else:
        print "STRANGER DANGER!"
        exit()

def is_phantom(chain_name, w0, wa, redshifts):
    eos = eos_from_chain_name(chain_name, w0, wa, redshifts)
    
    if np.min(eos) < -1:
        return True
    else:
        #print np.min(eos)
        return False

def get_color(shade):
    if shade == 'green':
        colors = ['g', 'forestgreen','lime']
    elif shade == 'red':
        colors = ['r', 'sienna','tomato']
    return np.random.choice(colors,size=1)[0]

    return 
def plot_minmax_deltamu(fname, title, legend_string=None):
    #redshifts, deltamu_min, deltamu_max, parameters_min, parameters_max = read_pickled_deltamu(fname)
    redshifts, deltamu_min, deltamu_max, parameters_min, parameters_max, deltamu, marg, m_bestfit_lcdm_marg = read_pickled_deltamu(fname)
    #plt.figure()
    plt.plot(redshifts, deltamu_min,label=legend_string)
    plt.plot(redshifts, deltamu_max)
    plt.title(title)
    plt.show()

def plot_min_deltamu(fname, title, legend_string=None):
    redshifts, deltamu_min, deltamu_max, parameters_min, parameters_max, deltamu, marg, m_bestfit_lcdm_marg = read_pickled_deltamu(fname)
    #plt.figure()
    plt.plot(redshifts, deltamu_min,label=legend_string)
    plt.title(title)

def plot_max_deltamu(fname, title, legend_string=None):
    redshifts, deltamu_min, deltamu_max, parameters_min, parameters_max, deltamu, marg, m_bestfit_lcdm_marg = read_pickled_deltamu(fname)
    #plt.figure()
    plt.plot(redshifts, deltamu_max,label=legend_string)
    plt.title(title)

def plot_deltamu(fname, title, legend_string=None):
    redshifts, deltamu_min, deltamu_max, parameters_min, parameters_max, deltamu, marg, m_bestfit_lcdm_marg = read_pickled_deltamu(fname)
    #plt.figure()
    for ii in np.arange(np.shape(deltamu)[1]):
        plt.plot(redshifts, deltamu[:,ii],label=legend_string)
    plt.title(title)
    plt.show()

def oplot_deltamus(chain_name, bins, smoothings, tolerance=0.005, label='CPL', thinning=10, ignore_k=False):
    if individual_plots:
        plt.figure()
        plt.xlabel('Redshift',size='x-large')
        plt.ylabel(r'$\Delta \mu$',size='x-large')

    plt.ylim((-0.1,0.1))
    ii_counter = 0

    deltamu_max_global = np.zeros(1000)
    deltamu_min_global = np.zeros(1000)

    deltamu_nonphantom = np.zeros((1,1000))

    for ii in np.arange(len(bins)):
        for jj in np.arange(len(smoothings)):
            deltamus = Deltamu.Deltamu(chain_name,'',do_marg=True,bins_tuple=bins[ii],smoothing=smoothings[jj],tolerance=tolerance)
            fname = root_dir + deltamus.get_marg_file_name()
            redshifts, deltamu_min, deltamu_max, parameters_min, parameters_max, deltamu, marg, m_bestfit_lcdm_marg = read_pickled_deltamu(fname)
            parameters = deltamus.get_parameters(verbose=False)
            
            for kk in np.arange(len(deltamu_min_global)):
                if deltamu_max[kk] > deltamu_max_global[kk]:
                    deltamu_max_global[kk] = deltamu_max[kk]
                if deltamu_min[kk] < deltamu_min_global[kk]:
                    deltamu_min_global[kk] = deltamu_min[kk]

            for kk in np.arange(np.shape(deltamu)[1]):
                ii_counter += 1
                w0 = parameters[1][kk]
                wa = parameters[2][kk]
                
                if is_phantom(chain_name, w0, wa, redshifts) == False:
                    if ignore_k:
                        deltamu_nonphantom = np.concatenate((deltamu_nonphantom,np.array([deltamu[:,kk]])))
                    else:
                        deltamu_nonphantom = np.concatenate((deltamu_nonphantom,np.array([deltamu[:,kk] + marg[kk] - m_bestfit_lcdm_marg])))
                
                if ii_counter == thinning:
                    if ignore_k:
                        if is_phantom(chain_name, w0, wa, redshifts):
                            plt.plot(redshifts, deltamu[:,kk], c=get_color(shade='green'))
                            #pass
                        else:
                            plt.plot(redshifts, deltamu[:,kk], c=get_color(shade='red'))
                            #pass
                    else:
                        if is_phantom(chain_name, w0, wa, redshifts):
                            plt.plot(redshifts, deltamu[:,kk] + marg[kk] - m_bestfit_lcdm_marg,c=get_color(shade='green'))
                            #pass
                        else:
                            plt.plot(redshifts, deltamu[:,kk] + marg[kk] - m_bestfit_lcdm_marg,c=get_color(shade='red'))
                            #pass
                    ii_counter = 0
    
    
    deltamu_nonphantom_max = np.zeros(1000)
    deltamu_nonphantom_min = np.zeros(1000)
    for ii in np.arange(len(deltamu_nonphantom_max)):
        deltamu_nonphantom_max[ii] = np.max(deltamu_nonphantom[:,ii])
        deltamu_nonphantom_min[ii] = np.min(deltamu_nonphantom[:,ii])
    

    if ignore_k==False:
        dashes = [20,10]
        lmax,=plt.plot(redshifts, deltamu_max_global,'k',ls='--',lw=3)
        lmin,=plt.plot(redshifts, deltamu_min_global,'k',ls='--',lw=3)
        lmax.set_dashes(dashes)
        lmin.set_dashes(dashes)

        dashes_nonphantom = [5,5]
        llmax,=plt.plot(redshifts, deltamu_nonphantom_max,'r',ls='--',lw=3)
        llmin,=plt.plot(redshifts, deltamu_nonphantom_min,'r',ls='--',lw=3)
        llmax.set_dashes(dashes_nonphantom)
        llmin.set_dashes(dashes_nonphantom)

    plt.text(7.5, 0.08, label,size='x-large')

def oplot_deltamu_test(chain_name,bins, smoothings, tolerance=0.005, label='CPL'):
    if individual_plots:
        fig = plt.figure()
        plt.xlabel('Redshift',size='x-large')
        plt.ylabel(r'$\Delta \mu$',size='x-large')
    plt.ylim((-0.1,0.1))

    deltamu_max_global = np.zeros(1000)
    deltamu_min_global = np.zeros(1000)
    deltamu_max_global_test = np.zeros(1000)
    deltamu_min_global_test = np.zeros(1000)

    for ii in np.arange(len(bins)):
        for jj in np.arange(len(smoothings)):
            deltamus = Deltamu.Deltamu(chain_name,'',do_marg=True,bins_tuple=bins[ii],smoothing=smoothings[jj],tolerance=tolerance)
            deltamus_test = Deltamu.Deltamu(chain_name,'',do_marg=True,bins_tuple=bins[ii],smoothing=smoothings[jj],tolerance=tolerance,testcase=True)
            fname = root_dir + deltamus.get_marg_file_name()
            fname_test = root_dir + deltamus_test.get_marg_file_name()
            redshifts, deltamu_min, deltamu_max, parameters_min, parameters_max, deltamu, marg, m_bestfit_lcdm_marg = read_pickled_deltamu(fname)
            redshifts_test, deltamu_min_test, deltamu_max_test, parameters_min_test, parameters_max_test, deltamu_test, marg_test, m_bestfit_lcdm_marg_test = read_pickled_deltamu(fname_test)
            for kk in np.arange(len(deltamu_min_global)):
                if deltamu_max[kk] > deltamu_max_global[kk]:
                    deltamu_max_global[kk] = deltamu_max[kk]
                if deltamu_min[kk] < deltamu_min_global[kk]:
                    deltamu_min_global[kk] = deltamu_min[kk]

                if deltamu_max_test[kk] > deltamu_max_global_test[kk]:
                    deltamu_max_global_test[kk] = deltamu_max_test[kk]
                if deltamu_min_test[kk] < deltamu_min_global_test[kk]:
                    deltamu_min_global_test[kk] = deltamu_min_test[kk]
    plt.fill_between(redshifts, deltamu_max_global,deltamu_min_global,color='darkgrey',label=label)
    plt.fill_between(redshifts_test, deltamu_max_global_test,deltamu_min_global_test,color='lightgrey',hatch='X',label=label + r", $w_0=-1, w_a=0$",edgecolor='darkgrey')
    plt.legend(frameon=False)
    if individual_plots:
        plt.savefig('Figures/test.pdf',format='pdf',dpi=fig.dpi)
    #plt.plot(redshifts, deltamu_max_global,c='b')
    #plt.plot(redshifts, deltamu_min_global,c='b')
    #plt.plot(redshifts_test, deltamu_max_global_test,c='g')
    #plt.plot(redshifts_test, deltamu_min_global_test,c='g')

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

    ax_scatter.scatter(x_contour, y_contour, z_contour,color='b',s=1.,depthshade=True)
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

'''
def plot_equation_of_state(wa_1,wa_2):
    redshifts = np.linspace(0.0001,10.,1000)
    scale_factor = 1. / (1. + redshifts)
    linestyles = ['-','--','-.',':']
    fig = plt.figure()
    plt.xlabel('Redshift',size='x-large')
    plt.ylabel(r'$w(z)$',size='xx-large')
    for ii in np.arange(len(wa_1)):
        eos_cpl = wa_1[ii][0] + wa_1[ii][1]*(1.-scale_factor)
        plt.plot(redshifts, eos_cpl,c='b',ls=linestyles[ii],lw=3, label=r'$w_0$ : ' + str(wa_1[ii][0]) + r', $w_a$ : ' + str(wa_1[ii][1]))
    for ii in np.arange(len(wa_2)):
        eos_cpl = wa_2[ii][0] + wa_2[ii][1]*((1.-scale_factor)**7)
        plt.plot(redshifts, eos_cpl,c='g',ls=linestyles[ii],lw=3, label=r'$w_0$ : ' + str(wa_2[ii][0]) + r', $w_a$ : ' + str(wa_2[ii][1]))
    plt.ylim((-1,-0.5))
    plt.legend(frameon=False, loc=2, fontsize=17)
    plt.xticks(size='x-large')
    plt.yticks(size='x-large')
    #plt.text(7.5,-0.8,'Thawing\n    CPL',size='xx-large',color='b')
    #plt.text(6,-0.6,'Freezing\n  n7CPL',size='xx-large',color='g')
    plt.text(7.5,-0.8,'Thawing',size='xx-large',color='b')
    plt.text(6,-0.6,'Freezing',size='xx-large',color='g')
    fig.set_tight_layout('True')
    plt.savefig('Figures/equationofstate.pdf',format='pdf')
'''

def plot_equation_of_state(wa_1,wa_2):
    redshifts = np.linspace(0.0001,10.,1000)
    scale_factor = 1. / (1. + redshifts)
    linestyles = ['-','--','-.',':']
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.xlabel('Redshift',size='x-large')
    plt.ylabel(r'$w(z)$',size='xx-large')
    
    for ii in np.arange(len(wa_1)):
        eos_cpl = wa_1[ii][0] + wa_1[ii][1]*(1.-scale_factor)
        plt.plot(redshifts, eos_cpl,c='b',ls=linestyles[ii],lw=3, label=r'$w_0$ : ' + str(wa_1[ii][0]) + r', $w_a$ : ' + str(wa_1[ii][1]))

    for ii in np.arange(len(wa_2)):
        eos_cpl = wa_2[ii][0] + wa_2[ii][1]*((1.-scale_factor)**7)
        plt.plot(redshifts, eos_cpl,c='g',ls=linestyles[ii],lw=3, label=r'$w_0$ : ' + str(wa_2[ii][0]) + r', $w_a$ : ' + str(wa_2[ii][1]))

    plt.ylim((-1.5,-0.5))
    plt.xlim((0,4))
    plt.legend(frameon=False, loc=9, fontsize=17,handlelength=2.3)
    plt.xticks(size='x-large')
    plt.yticks(size='x-large')
    plt.text(.5,-0.82,'Convex',size='x-large',color='b',rotation=-14)
    plt.text(.5,-1.18,'Concave',size='x-large',color='b',rotation=14)
    plt.text(1.6,-0.95,'Convex',size='x-large',color='g',rotation=10)
    plt.text(1.6,-1.06,'Concave',size='x-large',color='g',rotation=-8)
    ax.add_patch(patches.Rectangle((0,-1.5),4.,0.5,color='grey',alpha=0.2))
    txt = plt.text(1.2, -1.3,'Phantom regime',size='xx-large',color='darkgrey')
    txt.set_path_effects([patheffects.withStroke(linewidth=0.5,foreground='k')])
    fig.set_tight_layout('True')
    plt.savefig('Figures/equationofstate.pdf',format='pdf')

def combined_plot():
    f, (ax1, ax2, ax3) = plt.subplots(3,2,figsize=(8,10))

    log_x_axis = False
    
    plt.sca(ax1[0])
    oplot_deltamus('cpl', [(30,30,30),(40,40,40),(50,50,50)],[0.6],label='CPL',ignore_k=True,thinning=80)
    ax1[0].set_ylabel(r'$\mathbf{\Delta \mu}$',size='x-large')
    ax1[0].set_yticks([-0.08, -0.04, 0., 0.04, 0.08])
    ax1[0].set_xticks([0.])
    ax1[0].set_xticklabels([''])
    ax1[0].text(0.3,0.08,'(a)',size='x-large')
    if log_x_axis:
        ax1[0].set_xscale("log", nonposx='clip')
    
    plt.sca(ax1[1])
    oplot_deltamu_test('cpl', [(30,30,30),(40,40,40),(50,50,50)],[0.6],label='CPL')
    ax1[1].set_yticks([0.])
    ax1[1].set_yticklabels([''])
    ax1[1].set_xticks([0.])
    ax1[1].set_xticklabels([''])
    ax1[1].text(0.3,0.08,'(b)',size='x-large')
    if log_x_axis:
        ax1[1].set_xscale("log", nonposx='clip')

    plt.sca(ax2[0])
    oplot_deltamus('cpl', [(30,30,30),(40,40,40),(50,50,50)],[0.6],label='CPL',ignore_k=False,thinning=100)
    ax2[0].set_ylabel(r'$\mathbf{\Delta \mu}$',size='x-large')
    ax2[0].set_yticks([-0.08, -0.04, 0., 0.04, 0.08])
    ax2[0].set_xticks([0.])
    ax2[0].set_xticklabels([''])
    ax2[0].text(0.3,0.08,'(c)',size='x-large')
    if log_x_axis:
        ax2[0].set_xscale("log", nonposx='clip')

    plt.sca(ax2[1])
    oplot_deltamus('jbp', [(30,30,30),(40,40,40),(50,50,50)],[0.6],label='JBP',ignore_k=False,thinning=100)
    ax2[1].set_yticks([0.])
    ax2[1].set_yticklabels([''])
    ax2[1].set_xticks([0.])
    ax2[1].set_xticklabels([''])
    ax2[1].text(0.3,0.08,'(d)',size='x-large')
    if log_x_axis:
        ax2[1].set_xscale("log", nonposx='clip')

    plt.sca(ax3[0])
    oplot_deltamus('n3cpl', [(30,30,30),(40,40,40),(50,50,50)],[0.6],label='n3CPL',ignore_k=False,thinning=100)
    ax3[0].set_ylabel(r'$\mathbf{\Delta \mu}$',size='x-large')
    ax3[0].set_yticks([-0.08, -0.04, 0., 0.04, 0.08])
    ax3[0].set_xlabel('Redshift',size='x-large')
    ax3[0].text(0.3,0.08,'(e)',size='x-large')
    
    ax3[0].add_patch(patches.Rectangle((1.,-0.068),0.5,0.015,color='g'))
    ax3[0].text(1.6,-0.066,'Phantom',size='large')

    ax3[0].add_patch(patches.Rectangle((1.,-0.093),0.5,0.015,color='r'))
    ax3[0].text(1.6,-0.091,'Non-phantom',size='large')

    if log_x_axis:
        ax3[0].set_xscale("log", nonposx='clip')


    plt.sca(ax3[1])
    oplot_deltamus('n7cpl', [(30,30,30),(40,40,40)],[0.4],label='n7CPL',ignore_k=False,thinning=90)
    ax3[1].set_xlabel('Redshift',size='x-large')
    ax3[1].set_yticks([0.])
    ax3[1].set_yticklabels([''])
    ax3[1].set_xticks([2,4,6,8,10])
    ax3[1].text(0.3,0.08,'(f)',size='x-large')
    if log_x_axis:
        ax3[1].set_xscale("log", nonposx='clip')


    plt.subplots_adjust(left=0.11,bottom=0.05,right=0.98,top=0.98,wspace=0., hspace=0.)
    plt.savefig('Figures/combinedplot.pdf',format='pdf')

def additional_plots():
    f, (ax1, ax2, ax3) = plt.subplots(3,2,figsize=(8,10))
    
    plt.sca(ax1[0])
    oplot_deltamu_test('jbp', [(30,30,30),(40,40,40),(50,50,50)],[0.6],label='JBP')
    ax1[0].set_ylabel(r'$\mathbf{\Delta \mu}$',size='x-large')
    ax1[0].set_yticks([-0.08, -0.04, 0., 0.04, 0.08])
    ax1[0].set_xticks([0.])
    ax1[0].set_xticklabels([''])
    ax1[0].text(0.3,0.08,'(a)',size='x-large')
    
    plt.sca(ax1[1])
    oplot_deltamu_test('n3cpl', [(30,30,30),(40,40,40),(50,50,50)],[0.6],label='n3CPL')
    ax1[1].set_yticks([0.])
    ax1[1].set_yticklabels([''])
    ax1[1].set_xticks([0.])
    ax1[1].set_xticklabels([''])
    ax1[1].text(0.3,0.08,'(b)',size='x-large')

    plt.sca(ax2[0])
    oplot_deltamu_test('n7cpl', [(30,30,30),(40,40,40)],[0.4],label='n7CPL')
    ax2[0].set_ylabel(r'$\mathbf{\Delta \mu}$',size='x-large')
    ax2[0].set_yticks([-0.08, -0.04, 0., 0.04, 0.08])
    ax2[0].set_xticks([0.])
    ax2[0].set_xticklabels([''])
    ax2[0].text(0.3,0.08,'(c)',size='x-large')

    plt.sca(ax2[1])
    oplot_deltamus('jbp', [(30,30,30),(40,40,40),(50,50,50)],[0.6],label='JBP',ignore_k=True,thinning=80)
    ax2[1].set_yticks([0.])
    ax2[1].set_yticklabels([''])
    ax2[1].set_xticks([0.])
    ax2[1].set_xticklabels([''])
    ax2[1].text(0.3,0.08,'(d)',size='x-large')

    plt.sca(ax3[0])
    oplot_deltamus('n3cpl', [(30,30,30),(40,40,40),(50,50,50)],[0.6],label='n3CPL',ignore_k=True,thinning=80)
    ax3[0].set_ylabel(r'$\mathbf{\Delta \mu}$',size='x-large')
    ax3[0].set_yticks([-0.08, -0.04, 0., 0.04, 0.08])
    ax3[0].set_xlabel('Redshift',size='x-large')
    ax3[0].text(0.3,0.08,'(e)',size='x-large')


    plt.sca(ax3[1])
    oplot_deltamus('n7cpl', [(30,30,30),(40,40,40)],[0.4],label='n7CPL',ignore_k=True,thinning=80)
    ax3[1].set_xlabel('Redshift',size='x-large')
    ax3[1].set_yticks([0.])
    ax3[1].set_yticklabels([''])
    ax3[1].set_xticks([2,4,6,8,10])
    ax3[1].text(0.3,0.08,'(f)',size='x-large')


    plt.subplots_adjust(left=0.11,bottom=0.05,right=0.98,top=0.98,wspace=0., hspace=0.)
    plt.savefig('Figures/additionalplots.pdf',format='pdf')

def oplot_deltamu_extrema(chain_names, bins_list, smoothings_list, labels, tolerance=0.005):
    if individual_plots:
        fig = plt.figure()
        ax = plt.subplot(1,1,1)
        plt.xlabel('Redshift',size='x-large')
        plt.ylabel(r'$\Delta \mu$',size='x-large')

    deltamu_max_global = np.zeros((len(chain_names),1000))
    deltamu_min_global = np.zeros((len(chain_names),1000))

    deltamu_nonphantom = np.zeros((1,1000))

    colors = ['g', 'lawngreen', 'limegreen','b']
    colors_nonphantom = ['orange','orange','orange','r']
    linestyles = ['-','--',':','-']
    for ll in np.arange(len(chain_names)):
        chain_name = chain_names[ll]
        bins = bins_list[ll]
        smoothings = smoothings_list[ll]
        label = labels[ll]
        color = colors[ll]
        linestyle = linestyles[ll]
        for ii in np.arange(len(bins)):
            for jj in np.arange(len(smoothings)):
                deltamus = Deltamu.Deltamu(chain_name,'',do_marg=True,bins_tuple=bins[ii],smoothing=smoothings[jj],tolerance=tolerance)
                fname = root_dir + deltamus.get_marg_file_name()
                redshifts, deltamu_min, deltamu_max, parameters_min, parameters_max, deltamu, marg, m_bestfit_lcdm_marg = read_pickled_deltamu(fname)
                parameters = deltamus.get_parameters(verbose=False)

                for kk in np.arange(len(deltamu_min_global[ll])):
                    w0 = parameters[1][kk]
                    wa = parameters[2][kk]

                    if is_phantom(chain_name, w0, wa, redshifts) == False:
                        deltamu_nonphantom = np.concatenate((deltamu_nonphantom,np.array([deltamu[:,kk] + marg[kk] - m_bestfit_lcdm_marg])))

                    if deltamu_max[kk] > deltamu_max_global[ll][kk]:
                        deltamu_max_global[ll][kk] = deltamu_max[kk]
                    if deltamu_min[kk] < deltamu_min_global[ll][kk]:
                        deltamu_min_global[ll][kk] = deltamu_min[kk]
        
        deltamu_max_nonphantom = np.zeros(1000)
        deltamu_min_nonphantom = np.zeros(1000)
        for ii in np.arange(len(deltamu_max_nonphantom)):
            deltamu_max_nonphantom[ii] = np.max(deltamu_nonphantom[:,ii])
            deltamu_min_nonphantom[ii] = np.min(deltamu_nonphantom[:,ii])

        plt.plot(redshifts, deltamu_max_global[ll],lw=3,color=color,ls=linestyle,label=label)
        plt.plot(redshifts, deltamu_min_global[ll],lw=3,color=color,ls=linestyle)
        if chain_name == 'n7cpl' or chain_name == 'cpl':
            plt.plot(redshifts, deltamu_max_nonphantom,lw=3,color=colors_nonphantom[ll],ls=linestyle,label=label)
            plt.plot(redshifts, deltamu_min_nonphantom,lw=3,color=colors_nonphantom[ll],ls=linestyle)

    handles, label_names = ax.get_legend_handles_labels()
    handles.insert(6,plt.Line2D(redshifts,deltamu_min_nonphantom,linestyle='none',marker='None'))
    label_names.insert(6,'')
    
    
    handles_plot, labels_plot = ['','','','','','',''], ['']*7
    handles_plot[0] = handles[0]
    handles_plot[1] = handles[2]
    handles_plot[2] = handles[3]
    handles_plot[3] = handles[4]
    handles_plot[4] = handles[1]
    handles_plot[5] = handles[5]
    handles_plot[6] = handles[6]
    labels_plot[0] = label_names[0]
    labels_plot[1] = label_names[2]
    labels_plot[2] = label_names[3]
    labels_plot[3] = label_names[4]
    labels_plot[4] = label_names[1]
    labels_plot[5] = label_names[5]
    labels_plot[6] = label_names[6]
    
    font0 = FontProperties()
    font0.set_size('xx-large')
    font0.set_weight('bold')
    plt.text(2.3,0.007,'Phantom',fontproperties=font0)
    plt.text(6.3,0.007,'Non-phantom',fontproperties=font0)
    ax.legend(handles_plot,labels_plot,frameon=False, fontsize=20, ncol=2, bbox_to_anchor=(0.95,0.55), columnspacing=4.)
    #ax.legend(handles,label_names,frameon=False, loc=5, fontsize=20, ncol=2)
    plt.ylim((-0.06,0.06))
    plt.yticks([-0.06, -0.03, 0, 0.03, 0.06],size='x-large')
    plt.xticks(size='x-large')

    if individual_plots:
        fig.set_tight_layout('True')
        plt.savefig('Figures/deltamus_extrema.pdf',format='pdf')

root_dir = '/Users/perandersen/Data/HzSC/Deltamu/'

individual_plots = False
combined_plot()
additional_plots()


#oplot_deltamu_extrema(['cpl', 'jbp', 'n3cpl','n7cpl'],\
#[[(30,30,30),(40,40,40),(50,50,50)], [(30,30,30),(40,40,40),(50,50,50)],[(30,30,30),(40,40,40),(50,50,50)],[(30,30,30),(40,40,40)]],\
#[[0.6],[0.6],[0.6],[0.4]], ['CPL','JBP','n3CPL','n7CPL'])

#plot_3d_contours('n7cpl', [(40,40,40)], 0.4)

#oplot_deltamus('n7cpl', [(30,30,30),(40,40,40)],[0.4],label='n7CPL',ignore_k=True,thinning=10)
#oplot_deltamus('n3cpl', [(30,30,30),(40,40,40),(50,50,50)],[0.6],label='n3CPL')
#oplot_deltamus('jbp', [(30,30,30),(40,40,40),(50,50,50)],[0.6],label='JBP')
#oplot_deltamus('cpl', [(30,30,30),(40,40,40),(50,50,50)],[0.6],label='CPL')
#oplot_deltamus('lcdm', [70,80,90,100],[0.6],tolerance=0.01)

#plt.sca(ax1[1])
#oplot_deltamu_test('n7cpl', [(30,30,30),(40,40,40)],[0.4],label='n7CPL')
#oplot_deltamu_test('n3cpl', [(30,30,30),(40,40,40),(50,50,50)],[0.3],label='n3CPL')
#oplot_deltamu_test('jbp', [(30,30,30),(40,40,40),(50,50,50)],[0.6],label='JBP')
#oplot_deltamu_test('cpl', [(30,30,30),(40,40,40),(50,50,50)],[0.6],label='CPL')

#plot_equation_of_state([(-1.,0.1), (-1.,0.2), (-1.,0.3)],[(-1.,0.7), (-1.,0.8), (-1.,.9)])
#plot_equation_of_state([(-0.6,-0.4), (-1.4,0.4)],[(-1.,0.4), (-1.,-0.4)])
#plt.show()


