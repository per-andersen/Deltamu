import os
import numpy as np
import pickle as pick

class CosmoMC_dao(object):
    """A data access object for CosmoMC output
    """

    def __init__(self, chain_name, directory):
        self.chain_name = chain_name
        self.directory = directory


    def print_parameter_names(self, find_names=[]):
        files =  os.listdir(self.directory)
        param_names_file = [s for s in files if "paramnames" in s]
        
        if len(param_names_file) != 1:
            print "Too many parameter name files!"
            raisej
        
        param_names_file = param_names_file[0]

        index = 2
        print "\nIndices to pull parameter from CosmoMC chain .txt file"
        print "(First two (0 and 1) are weight and likelihood)\n"
        if (len(find_names) == 0):
            with open(self.directory + param_names_file) as f:
                for line in f:
                    print index, ":", line
                    index += 1
        else:
            for param_name in find_names:
                with open(self.directory + param_names_file) as f:
                    for line in f:
                        if (param_name in line):
                            print index, ":", line
                        index += 1
                    index = 2

    def get_parameter_names(self, indices=[]):
        '''Returns list of strings with parameter names
        NOTE: For consistency we use index from raw cosmomc chain files,
        which include a weight and likelihood parameter. Therefore we
        subtract two from the index.
        '''
        files =  os.listdir(self.directory)
        param_names_file = [s for s in files if "paramnames" in s]
        
        if len(param_names_file) != 1:
            print "Too many parameter name files!"
            raise
        
        param_names_file = param_names_file[0]

        param_names = []
        with open(self.directory + param_names_file) as f:
                for line in f:
                    param_names.append(line)

        param_indices = []
        for ii in np.arange(len(indices)):
            param_indices.append(param_names[indices[ii]-2]) #Subtract two from index for consistency with CosmoMC format
        return param_indices


    def get_parameter_chains(self,burn,indices=[]):
        files =  os.listdir(self.directory) #Get all files in dir
        files = [s for s in files if self.chain_name in s] #Find ones associated with chain
        files = [s for s in files if 'txt' in s]#Find only data files

        if (len(files) == 0):
            print "No files with chain name in directory!"
            raise
        
        #We read the first file to get the number of parameters in chain
        first_data_file = np.genfromtxt(self.directory + files[0])
        number_of_params = np.shape(first_data_file)[1]
        data_all = np.zeros((0,number_of_params))

        #We read full data from all files into one array
        print "Reading files..."
        for ff in files:
            data_from_file = np.genfromtxt(self.directory + ff)
            print ff, np.shape(data_from_file)
            data_all = np.concatenate((data_all, data_from_file[burn:,:]))
        
        #Data from
        weights = data_all[:,0]
        number_iterations = int(np.sum(weights))

        data_indices = np.zeros((number_iterations,len(indices)))
        index_counter = 0
        #We loop over all parameters we wish to grab
        for ii in np.arange(len(indices)):
            #This is the length of the chain output from CosmoMC
            for jj in np.arange(len(weights)):
                #To account for the weights we loop over each value at a time
                for kk in np.arange(weights[jj]):
                    data_indices[index_counter,ii] = data_all[jj,indices[ii]]
                    index_counter+=1
            index_counter = 0
        return data_indices

    def pickle_sliced_chains(self,indices=[],fname="", burn=2000):
        data = self.get_parameter_chains(indices=indices,burn=burn)
        
        data_labels = self.get_parameter_names(indices=indices)
        output = [data, data_labels]
        data_file_name = self.directory + "../Pickle/" + fname
        try:
            output_file = open(data_file_name,'wb')
            pick.dump(output,output_file)
            output_file.close()
        except:
            print "Can not write pickled file to disk!"


    def read_pickled_chains(self, fname):
        pkl_data_file = open(self.directory + "../Pickle/" + fname,'rb')
        data = pick.load(pkl_data_file)
        pkl_data_file.close()
        return data



if __name__ == "__main__":

    #The code below tests some basic functionality
    cosmomc_lcdm = CosmoMC_dao('lcdm', '/Users/perandersen/Data/HzSC/Chains/LCDM/')
    cosmomc_cpl = CosmoMC_dao('cpl', '/Users/perandersen/Data/HzSC/Chains/CPL/')
    cosmomc_jbp = CosmoMC_dao('jbp', '/Users/perandersen/Data/HzSC/Chains/JBP/')
    cosmomc_jbp_ns = CosmoMC_dao('jbp-ns', '/Users/perandersen/Data/HzSC/Chains/JBP-ns/')
    cosmomc_n3cpl = CosmoMC_dao('n3cpl', '/Users/perandersen/Data/HzSC/Chains/n3CPL/')
    cosmomc_n7cpl = CosmoMC_dao('n7cpl', '/Users/perandersen/Data/HzSC/Chains/n7CPL/')

    '''Test printing of parameter names. Should be run one at a time'''
    #cosmomc_lcdm.print_parameter_names(["omega", "chi"])
    #cosmomc_lcdm.print_parameter_names(["chi", "omega"])
    #cosmomc_lcdm.print_parameter_names(["chi", "omega"])
    #cosmomc_lcdm.print_parameter_names(["omega"])
    #cosmomc_cpl.print_parameter_names(["omega", "w"])
    #cosmomc_jbp.print_parameter_names(["omega", "w"])
    #cosmomc_jbp_ns.print_parameter_names(["omega", "w"])
    #cosmomc_n3cpl.print_parameter_names(["omega", "w"])
    #cosmomc_n7cpl.print_parameter_names(["omega", "w"])

    '''Below tests getting data, gets and plots it'''
    
    '''
    params = cosmomc_cpl.get_parameter_chains(indices=[31, 6, 7, 30])
    params_lcdm = cosmomc_lcdm.get_parameter_chains(indices=[29])
    import matplotlib.pyplot as plt
    print np.shape(params)
    print np.mean(params[:,3]), np.min(params[:,3]), np.max(params[:,3])
    print np.mean(params[:,0]), np.min(params[:,0]), np.max(params[:,0])
    print np.mean(params[:,1]), np.min(params[:,1]), np.max(params[:,1])
    print np.mean(params[:,2]), np.min(params[:,2]), np.max(params[:,2])
    plt.figure()
    plt.hist(params[:,0], bins = 200)
    plt.figure()
    plt.hist(params[:,1], bins = 200)
    plt.figure()
    plt.hist(params[:,2], bins = 200)
    plt.figure()
    plt.hist(params[:,3], bins = 200)
    plt.figure()
    plt.hist(params_lcdm[:,0], bins = 200)
    plt.show()
    '''
    

    '''Below tests writing and reading data'''
    #cosmomc_lcdm.pickle_sliced_chains(indices=[29], fname='lcdm_chains.pkl')
    #cosmomc_cpl.pickle_sliced_chains(indices=[31, 6, 7], fname='cpl_chains.pkl')
    #cosmomc_jbp.pickle_sliced_chains(indices=[31, 6, 7], fname='jbp_chains.pkl')
    cosmomc_jbp_ns.pickle_sliced_chains(indices=[31, 6, 7], fname='jbp-ns_chains.pkl')
    cosmomc_n3cpl.pickle_sliced_chains(indices=[31, 6, 7], fname='n3cpl_chains.pkl')
    cosmomc_n7cpl.pickle_sliced_chains(indices=[31, 6, 7], fname='n7cpl_chains.pkl')
    
    #data, data_labels = cosmomc_cpl.read_pickled_chains(fname='cpl_contours.pkl')
    #print cosmomc_cpl.get_parameter_names(indices=[31, 6, 7])
    #print np.shape(data)
    #print data_labels




