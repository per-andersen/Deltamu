from __future__ import division
import numpy as np
import pickle as pick
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from scipy.ndimage.filters import gaussian_filter

class Contour(object):
    """A contour class for 3d contours
    """

    def __init__(self, chain_name, directory, contour_level = 0.6827, tolerance = 0.001, bins_tuple=(50,50,50),smoothing=0):
        self.chain_name = chain_name
        self.directory = directory
        self.contour_level = contour_level
        self.tolerance = tolerance
        self.bins_tuple = bins_tuple
        self.smoothing = smoothing

    #These first methods deal with handling data from CosmoMC_dao, creating
    #histogram and contours and finally writing to .pkl file
    def histogram_points(self):
        points, data_labels = self.read_pickled_chains()
        hist, edges = np.histogramdd(points, bins = self.bins_tuple)

        if self.smoothing != 0:
            hist = gaussian_filter(hist, self.smoothing)

        xbins = edges[0]
        ybins = edges[1]
        zbins = edges[2]
        
        shape = hist.shape
        
        hist = hist.ravel() # Returns a flattened array
        
        i_sort = np.argsort(hist)[::-1] # Returns indices sorting hist from largest to smallest
        i_unsort = np.argsort(i_sort)# Returns indices undoing sort above
        
        hist_cumsum = hist[i_sort].cumsum() # Gets the cumulative sum, adding most filled indices first
        
        hist_cumsum /= hist_cumsum[-1]
        
        xbins = 0.5 * (xbins[1:] + xbins[:-1])
        ybins = 0.5 * (ybins[1:] + ybins[:-1])
        zbins = 0.5 * (zbins[1:] + zbins[:-1])
        
        sigma = hist_cumsum[i_unsort].reshape(shape)
        #plt.figure()
        #plt.imshow(sigma[:,:,int(self.bins_tuple[0] / 2)].T, extent=(0,1,0,1), origin='lower', interpolation='none')
        #plt.colorbar()
        #plt.figure()
        #plt.contour(xbins, ybins, sigma[:,:,int(self.bins_tuple[0] / 2)].T,levels=[0.683])
        #plt.title('Original')

        return xbins, ybins, zbins, sigma

    def histogram_to_contours(self):
        xbins, ybins, zbins, sigma = self.histogram_points()
        likelihoods = np.zeros(len(xbins)**3)
        x_contour = np.zeros(len(xbins)**3)
        y_contour = np.zeros(len(xbins)**3)
        z_contour = np.zeros(len(xbins)**3)

        for ii in np.arange(len(xbins)):
            for jj in np.arange(len(xbins)):
                for kk in np.arange(len(xbins)):
                    x_contour[kk + jj*len(xbins) + ii*len(xbins)**2] = xbins[ii]
                    y_contour[kk + jj*len(xbins) + ii*len(xbins)**2] = ybins[jj]
                    z_contour[kk + jj*len(xbins) + ii*len(xbins)**2] = zbins[kk]
                    likelihoods[kk + jj*len(xbins) + ii*len(xbins)**2] = sigma[ii,jj,kk]
        
        idx = np.where((likelihoods - self.contour_level) < self.tolerance)
        x_contour = x_contour[idx]
        y_contour = y_contour[idx]
        z_contour = z_contour[idx]

        print x_contour
        return x_contour, y_contour, z_contour

    def histogram_to_contours_paths(self, all_directions=True):
        '''
        This is an alternative to the first version of 'histogram_to_contours'
        where we use the get_paths() method to get the contours
        '''
        xbins, ybins, zbins, sigma = self.histogram_points()
        
        x_save = np.array([])
        y_save = np.array([])
        z_save = np.array([])

        if all_directions:
            
            fig1, ax1 = plt.subplots()
            for ii in np.arange(len(zbins)):
                cs = ax1.contour(xbins, ybins, sigma[:,:,ii].T,levels=[self.contour_level])
                go = 1
                path_number = 0
                while (go == 1):
                    try:
                        p = cs.collections[0].get_paths()[path_number]
                        v = p.vertices
                        xx = v[:,0]
                        yy = v[:,1]
                        #print zbins[ii]
                        #plt.plot(xx,yy,'bx')
                        #plt.show()
                        path_number += 1
                        for jj in np.arange(len(xx)):
                            #print xx[jj], yy[jj], zbins[ii]
                            x_save = np.append(x_save,xx[jj])
                            y_save = np.append(y_save,yy[jj])
                            z_save = np.append(z_save,zbins[ii])
                    except:
                        go = 0
                #print "i:", ii, " paths:", path_number, " z:", zbins[ii]
            plt.clf()
            plt.close()

            fig3, ax3 = plt.subplots()
            for ii in np.arange(len(ybins)):
                cst = ax3.contour(xbins, zbins, sigma[:,ii,:].T,levels=[self.contour_level])
                go = 1
                path_number = 0
                while (go == 1):
                    try:
                        p = cst.collections[0].get_paths()[path_number]
                        v = p.vertices
                        xx = v[:,0]
                        zz = v[:,1]
                        #print ybins[ii]
                        #plt.plot(xx,zz,'bx')
                        #plt.show()
                        path_number += 1
                        for jj in np.arange(len(xx)):
                            #print xx[jj], yy[jj], zbins[ii]
                            x_save = np.append(x_save,xx[jj])
                            y_save = np.append(y_save,ybins[ii])
                            z_save = np.append(z_save,zz[jj])
                    except:
                        go = 0
                #print "i:", ii, " paths:", path_number, " y:", ybins[ii]
            plt.clf()
            plt.close()
            
            
            fig2, ax2 = plt.subplots()
            for ii in np.arange(len(xbins)):
                cst = ax2.contour(ybins, zbins, sigma[ii,:,:].T,levels=[self.contour_level])
                go = 1
                path_number = 0
                while (go == 1):
                    try:
                        p = cst.collections[0].get_paths()[path_number]
                        v = p.vertices
                        yy = v[:,0]
                        zz = v[:,1]
                        #plt.plot(yy,zz,'bx')
                        #plt.show()
                        path_number += 1
                        for jj in np.arange(len(yy)):
                            #print xx[jj], yy[jj], zbins[ii]
                            x_save = np.append(x_save,xbins[ii])
                            y_save = np.append(y_save,yy[jj])
                            z_save = np.append(z_save,zz[jj])
                    except:
                        go = 0
                #print "i:", ii, " paths:", path_number, " y:", ybins[ii]
            plt.clf()
            plt.close()
            
            
            #print len(x_save), len(y_save), len(z_save)
            print np.min(x_save), np.max(x_save)
            print np.min(y_save), np.max(y_save)
            print np.min(z_save), np.max(z_save)
            return x_save, y_save, z_save

        else:
            fig3, ax3 = plt.subplots()
            for ii in np.arange(len(xbins)):
                cs = ax3.contour(xbins, ybins, sigma[:,:,ii].T,levels=[self.contour_level])
                go = 1
                path_number = 0
                while (go == 1):
                    try:
                        p = cs.collections[0].get_paths()[path_number]
                        v = p.vertices
                        xx = v[:,0]
                        yy = v[:,1]
                        #plt.plot(xx,yy,'bx')
                        path_number += 1
                        for jj in np.arange(len(xx)):
                            #print xx[jj], yy[jj], zbins[ii]
                            x_save = np.append(x_save,xx[jj])
                            y_save = np.append(y_save,yy[jj])
                            z_save = np.append(z_save,zbins[ii])
                    except:
                        go = 0
                #print "i:", ii, " paths:", path_number, " z:", zbins[ii]
                #plt.show()
            plt.clf()
            plt.close()
            #print len(x_save), len(y_save), len(z_save)
            return x_save, y_save, z_save 

    def pickle_contour(self):
        x_contour, y_contour, z_contour = self.histogram_to_contours_paths()
        
        output = x_contour, y_contour, z_contour
        
        fname = "contour_" + self.chain_name + "_c" + str(self.contour_level) +\
        "_t" + str(self.tolerance) + "_s" + str(self.smoothing) +\
        "_b" + str(self.bins_tuple[0]) + str(self.bins_tuple[1]) +\
        str(self.bins_tuple[2]) + ".dat"

        data_file_name = self.directory + "Contours/" + fname
        try:
            output_file = open(data_file_name,'wb')
            pick.dump(output,output_file)
            output_file.close()
        except:
            print "Can not write pickled file to disk!"

    def read_pickled_contour(self):
        fname = "contour_" + self.chain_name + "_c" + str(self.contour_level) +\
        "_t" + str(self.tolerance) + "_s" + str(self.smoothing) +\
        "_b" + str(self.bins_tuple[0]) + str(self.bins_tuple[1]) +\
        str(self.bins_tuple[2]) + ".dat"

        pkl_data_file = open(self.directory + "Contours/" + fname,'rb')
        data = pick.load(pkl_data_file)
        pkl_data_file.close()
        return data

    def read_pickled_chains(self):
        pkl_data_file = open(self.directory + "Chains/Pickle/" + self.chain_name + "_chains.pkl",'rb')
        data = pick.load(pkl_data_file)
        pkl_data_file.close()
        return data
    
    #These methods below use the pickled contours and do stuff with them
    def plot_contour(self, labels = ['omegam','w0','wa']):
        try:
            x_contour, y_contour, z_contour = self.read_pickled_contour()
        except:
            self.pickle_contour()
            x_contour, y_contour, z_contour = self.read_pickled_contour()
        
        plt.figure()
        plt.plot(x_contour, y_contour,'.b')
        plt.xlabel(labels[0])
        plt.ylabel(labels[1])
        plt.figure()
        plt.plot(x_contour, z_contour,'.b')
        plt.xlabel(labels[0])
        plt.ylabel(labels[2])
        plt.figure()
        plt.plot(z_contour, y_contour,'.b')
        plt.xlabel(labels[2])
        plt.ylabel(labels[1])
        #plt.show()

    def plot_contour_slice(self):
        try:
            x_contour, y_contour, z_contour = self.read_pickled_contour()
        except:
            self.pickle_contour()
            x_contour, y_contour, z_contour = self.read_pickled_contour()

        x_median = x_contour[z_contour == np.median(z_contour)]
        y_median = y_contour[z_contour == np.median(z_contour)]
        plt.figure()
        plt.title('median')
        plt.plot(x_median,y_median,'.b')

    def plot_contour_3d(self, labels = ['x','y','z']):
        try:
            x_contour, y_contour, z_contour = self.read_pickled_contour()
        except:
            self.pickle_contour()
            x_contour, y_contour, z_contour = self.read_pickled_contour()
        print len(x_contour), len(y_contour), len(z_contour)
        print np.min(x_contour), np.max(x_contour)
        print np.min(y_contour), np.max(y_contour)
        print np.min(z_contour), np.max(z_contour)
        fig_scatter = plt.figure()
        ax_scatter = fig_scatter.add_subplot(111, projection='3d')
        ax_scatter.scatter(x_contour, y_contour, z_contour)
        ax_scatter.set_xlabel(labels[0])
        ax_scatter.set_ylabel(labels[1])
        ax_scatter.set_zlabel(labels[2])
        #plt.show()


    def test_contour_exists(self):
        try:
            self.read_pickled_contour()
            return True
        except:
            return False

class LCDM_Contour(Contour):
    """A contour class for 1d contours, specifically LCDM contours
    This class is really hacky, as it tries to do contour stuff on
    a 1d array. Also, this produces the unmarginalised distribution
    of omega_m values, which doesn't take into account covariances
    between omega_m, w0, wa, and H0. It also looks at the entire
    distribution, not just the boundary (which would be two points).
    This way we will not have to show that the extreme value theorem
    holds true for this particular case.
    """

    def __init__(self, chain_name, directory, contour_level = 0.6827, tolerance = 0.01, bins_tuple=100,smoothing=0.5):
        self.chain_name = chain_name
        self.directory = directory
        self.contour_level = contour_level
        self.tolerance = tolerance
        self.bins_tuple = bins_tuple
        self.smoothing = smoothing

    def histogram_points(self):
        points, data_labels = self.read_pickled_chains()
        H, edges = np.histogram(points, bins = self.bins_tuple)
        
        xbins = edges
        
        shape = H.shape
        
        if self.smoothing != 0:
            H = gaussian_filter(H, self.smoothing)

        H = H.ravel() # Returns a flattened array
        
        i_sort = np.argsort(H)[::-1] # Returns indices sorting H from largest to smallest
        i_unsort = np.argsort(i_sort)# Returns indices undoing sort above
        
        H_cumsum = H[i_sort].cumsum() # Gets the cumulative sum, adding most filled indices first
        
        H_cumsum = H_cumsum / H_cumsum[-1]
        
        
        xbins = 0.5 * (xbins[1:] + xbins[:-1])
        
        sigma = H_cumsum[i_unsort].reshape(shape)
        #print xbins
        #print sigma
        return xbins, sigma

    def histogram_to_contours(self):
        xbins, sigma = self.histogram_points()
        
        
        idx = np.where(np.abs(sigma) < self.contour_level)
        x_contour = xbins[idx]

        return x_contour

    def pickle_contour(self):
        x_contour = self.histogram_to_contours()
        
        output = x_contour
        
        fname = "contour_" + self.chain_name + "_c" + str(self.contour_level) +\
        "_t" + str(self.tolerance) + "_s" + str(self.smoothing) + "_b" + str(self.bins_tuple) + ".dat"

        data_file_name = self.directory + "Contours/" + fname
        try:
            output_file = open(data_file_name,'wb')
            pick.dump(output,output_file)
            output_file.close()
        except:
            print "Can not write pickled file to disk!"

    def read_pickled_contour(self):
        fname = "contour_" + self.chain_name + "_c" + str(self.contour_level) +\
        "_t" + str(self.tolerance) + "_s" + str(self.smoothing) + "_b" + str(self.bins_tuple) + ".dat"

        pkl_data_file = open(self.directory + "Contours/" + fname,'rb')
        data = pick.load(pkl_data_file)
        pkl_data_file.close()
        return data

    #These methods below use the pickled contours and do stuff with them
    def plot_contour(self, labels = ['x']):
        xbins, sigma = self.histogram_points()
        x_contour = self.histogram_to_contours()
        print x_contour

        plt.figure()
        plt.plot(xbins, -sigma)
        plt.xlabel(labels[0])
        plt.ylabel('-sigma')
        plt.show()


if __name__ == "__main__":
    #CPL_Contour = Contour(chain_name='cpl', directory='/Users/perandersen/Data/HzSC/')

    #CPL_Contour.test_contour_exists()
    #CPL_Contour.pickle_contour()
    #CPL_Contour.plot_contour(labels=['omega_m','w0','wa'])
    #CPL_Contour.plot_contour_3d(labels=['omega_m','w0','wa'])

    #JBP_Contour = Contour(chain_name='jbp', directory='/Users/perandersen/Data/HzSC/', bins_tuple=(20,20,20))
    #JBP_Contour.test_contour_exists()
    #JBP_Contour.pickle_contour()
    #JBP_Contour.plot_contour(labels=['omega_m','w0','wa'])
    #JBP_Contour.plot_contour_3d(labels=['omega_m','w0','wa'])

    #LCDM_Contour = LCDM_Contour(chain_name='lcdm', directory='/Users/perandersen/Data/HzSC/')
    #LCDM_Contour.pickle_contour()
    #print LCDM_Contour.test_contour_exists()
    #LCDM_Contour.plot_contour()
    '''
    cpl_bins_tuples = [(20,20,20),(30,30,30),(40,40,40),(50,50,50),(60,60,60),(70,70,70),(80,80,80)]
    for cpl_bins in cpl_bins_tuples:
        print cpl_bins
        CPL_Contour = Contour(chain_name='cpl', directory='/Users/perandersen/Data/HzSC/',bins_tuple=cpl_bins,tolerance=0.001)
        CPL_Contour.plot_contour()
    '''

    #CPL_Contour = Contour(chain_name='cpl', directory='/Users/perandersen/Data/HzSC/',bins_tuple=(30,30,30),tolerance = 0.001)
    #CPL_Contour.pickle_contour()
    #CPL_Contour.plot_contour_3d()
    #CPL_Contour.plot_contour_slice()
    
    #CPL_Contour = Contour(chain_name='cpl', directory='/Users/perandersen/Data/HzSC/',bins_tuple=(40,40,40),tolerance = 0.001, smoothing=0.6)
    #CPL_Contour.pickle_contour()
    #CPL_Contour.plot_contour_3d()
    #CPL_Contour.plot_contour_slice()

    #JBP_Contour = Contour(chain_name='jbp', directory='/Users/perandersen/Data/HzSC/',bins_tuple=(30,30,30),tolerance = 0.001)
    #JBP_Contour.pickle_contour()
    #JBP_Contour.plot_contour_3d()
    #JBP_Contour.plot_contour_slice()
    
    #JBP_Contour = Contour(chain_name='jbp', directory='/Users/perandersen/Data/HzSC/',bins_tuple=(30,30,30),tolerance = 0.001, smoothing=0.4)
    #JBP_Contour.pickle_contour()
    #JBP_Contour.plot_contour_3d()
    #JBP_Contour.plot_contour_slice()

    #n3CPL_Contour = Contour(chain_name='n3cpl', directory='/Users/perandersen/Data/HzSC/',bins_tuple=(60,60,60),tolerance = 0.001, smoothing = 0.6)
    #n3CPL_Contour.pickle_contour()
    #n3CPL_Contour.plot_contour_3d()
    #n3CPL_Contour.plot_contour_slice()

    #n3CPL_Contour = Contour(chain_name='n3cpl', directory='/Users/perandersen/Data/HzSC/',bins_tuple=(30,30,30),tolerance = 0.001, smoothing = 0.6)
    #n3CPL_Contour.pickle_contour()
    #n3CPL_Contour.plot_contour_3d()
    #n3CPL_Contour.plot_contour_slice()

    #n3CPL_Contour = Contour(chain_name='n3cpl', directory='/Users/perandersen/Data/HzSC/',bins_tuple=(30,30,30),tolerance = 0.001, smoothing=0.6)
    #n3CPL_Contour.pickle_contour()
    #n3CPL_Contour.plot_contour_3d(labels=['omegam','w0','wa'])
    #n3CPL_Contour.plot_contour_slice()

    #n7CPL_Contour = Contour(chain_name='n7cpl', directory='/Users/perandersen/Data/HzSC/',bins_tuple=(40,40,40),tolerance = 0.001)
    #n7CPL_Contour.pickle_contour()
    #n7CPL_Contour.plot_contour_3d()
    #n7CPL_Contour.plot_contour_slice()

    n7CPL_Contour = Contour(chain_name='n7cpl', directory='/Users/perandersen/Data/HzSC/',bins_tuple=(10,10,30),tolerance = 0.001, smoothing=0.4)
    #n7CPL_Contour.pickle_contour()
    #n7CPL_Contour.plot_contour_3d(labels=['omegam','w0','wa'])
    #n7CPL_Contour.plot_contour_slice()

    plt.show()