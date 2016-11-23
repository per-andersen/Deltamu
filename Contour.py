from __future__ import division
import numpy as np
import pickle as pick
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

class Contour(object):
    """A contour class for 3d contours
    """

    def __init__(self, chain_name, directory, contour_level = 0.68, tolerance = 0.005, bins_tuple=(50,50,50)):
        self.chain_name = chain_name
        self.directory = directory
        self.contour_level = contour_level
        self.tolerance = tolerance
        self.bins_tuple = bins_tuple

    #These first methods deal with handling data from CosmoMC_dao, creating
    #histogram and contours and finally writing to .pkl file
    def histogram_points(self):
        points, data_labels = self.read_pickled_chains()
        H, edges = np.histogramdd(points, bins = self.bins_tuple)
        xbins = edges[0]
        ybins = edges[1]
        zbins = edges[2]
        
        shape = H.shape
        
        H = H.ravel() # Returns a flattened array
        
        i_sort = np.argsort(H)[::-1] # Returns indices sorting H from largest to smallest
        i_unsort = np.argsort(i_sort)# Returns indices undoing sort above
        
        H_cumsum = H[i_sort].cumsum() # Gets the cumulative sum, adding most filled indices first
        
        H_cumsum /= H_cumsum[-1]
        
        xbins = 0.5 * (xbins[1:] + xbins[:-1])
        ybins = 0.5 * (ybins[1:] + ybins[:-1])
        zbins = 0.5 * (zbins[1:] + zbins[:-1])
        
        sigma = H_cumsum[i_unsort].reshape(shape)
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
        
        idx = np.where((np.abs(likelihoods) - self.contour_level) < self.tolerance)
        x_contour = x_contour[idx]
        y_contour = y_contour[idx]
        z_contour = z_contour[idx]

        return x_contour, y_contour, z_contour

    def pickle_contour(self):
        x_contour, y_contour, z_contour = self.histogram_to_contours()
        
        output = x_contour, y_contour, z_contour, self.contour_level, self.tolerance, self.bins_tuple
        
        fname = "contour_" + self.chain_name + "_c" + str(self.contour_level) +\
        "_t" + str(self.tolerance) + "_b" + str(self.bins_tuple[0]) + \
        str(self.bins_tuple[1]) + str(self.bins_tuple[2]) + ".dat"

        data_file_name = self.directory + "Contours/" + fname
        try:
            output_file = open(data_file_name,'wb')
            pick.dump(output,output_file)
            output_file.close()
        except:
            print "Can not write pickled file to disk!"

    def read_pickled_contour(self):
        fname = "contour_" + self.chain_name + "_c" + str(self.contour_level) +\
        "_t" + str(self.tolerance) + "_b" + str(self.bins_tuple[0]) + \
        str(self.bins_tuple[1]) + str(self.bins_tuple[2]) + ".dat"

        pkl_data_file = open(self.directory + "Contours/" + fname,'rb')
        data = pick.load(pkl_data_file)
        pkl_data_file.close()
        return data

    def read_pickled_chains(self):
        pkl_data_file = open(self.directory + "Chains/Pickle/" + self.chain_name + "_contours.pkl",'rb')
        data = pick.load(pkl_data_file)
        pkl_data_file.close()
        return data
    
    #These methods below use the pickled contours and do stuff with them
    def plot_contour(self, labels = ['x','y','z']):
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
        plt.show()

    def plot_contour_3d(self, labels = ['x','y','z']):
        try:
            x_contour, y_contour, z_contour = self.read_pickled_contour()
        except:
            self.pickle_contour()
            x_contour, y_contour, z_contour = self.read_pickled_contour()
        
        fig_scatter = plt.figure()
        ax_scatter = fig_scatter.add_subplot(111, projection='3d')
        ax_scatter.scatter(x_contour, y_contour, z_contour)
        plt.show()

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

    def __init__(self, chain_name, directory, contour_level = 0.68, tolerance = 0.01, bins_tuple=100):
        self.chain_name = chain_name
        self.directory = directory
        self.contour_level = contour_level
        self.tolerance = tolerance
        self.bins_tuple = bins_tuple

    def histogram_points(self):
        points, data_labels = self.read_pickled_chains()
        H, edges = np.histogram(points, bins = self.bins_tuple)
        
        xbins = edges
        
        shape = H.shape
        
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
        
        output = x_contour, self.contour_level, self.tolerance, self.bins_tuple
        
        fname = "contour_" + self.chain_name + "_c" + str(self.contour_level) +\
        "_t" + str(self.tolerance) + "_b" + str(self.bins_tuple) + ".dat"

        data_file_name = self.directory + "Contours/" + fname
        try:
            output_file = open(data_file_name,'wb')
            pick.dump(output,output_file)
            output_file.close()
        except:
            print "Can not write pickled file to disk!"

    def read_pickled_contour(self):
        fname = "contour_" + self.chain_name + "_c" + str(self.contour_level) +\
        "_t" + str(self.tolerance) + "_b" + str(self.bins_tuple) + ".dat"

        pkl_data_file = open(self.directory + "Contours/" + fname,'rb')
        data = pick.load(pkl_data_file)
        pkl_data_file.close()
        return data

    #These methods below use the pickled contours and do stuff with them
    def plot_contour(self, labels = ['x']):
        xbins, sigma = self.histogram_points()

        plt.figure()
        plt.plot(xbins, -sigma)
        plt.xlabel(labels[0])
        plt.ylabel('-sigma')
        plt.show()


if __name__ == "__main__":
    CPL_Contour = Contour(chain_name='cpl', directory='/Users/perandersen/Data/HzSC/')
    #CPL_Contour.test_contour_exists()
    CPL_Contour.pickle_contour()
    #CPL_Contour.plot_contour(labels=['omega_m','w0','wa'])
    #CPL_Contour.plot_contour_3d(labels=['omega_m','w0','wa'])

    JBP_Contour = Contour(chain_name='jbp', directory='/Users/perandersen/Data/HzSC/', bins_tuple=(20,20,20))
    #JBP_Contour.test_contour_exists()
    JBP_Contour.pickle_contour()
    #JBP_Contour.plot_contour(labels=['omega_m','w0','wa'])
    #JBP_Contour.plot_contour_3d(labels=['omega_m','w0','wa'])

    LCDM_Contour = LCDM_Contour(chain_name='lcdm', directory='/Users/perandersen/Data/HzSC/')
    LCDM_Contour.pickle_contour()
    #print LCDM_Contour.test_contour_exists()
    #LCDM_Contour.plot_contour()
