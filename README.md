# [Deltamu](https://github.com/per-andersen/Deltamu)
The code to perform analysis of CosmoMC output chains. Final results is plot of deltamu values versus redshift for the cosmology given as input to CosmoMC. See below for the steps necessary to produce the final results. This code is released under the GNU general public license, version 3. In order to run the code the numpy, scipy, matplotlib, pickle, and astropy libraries must be installed. The easiest way is to use pip, e.g. run "pip install numpy" for all the needed libraries. 

### CosmoMC_DAO
First run CosmoMC_DAO.py. This takes the raw chains from CosmoMC, extracts the values we need and save it in Pickle format for faster reading.

### Contour
This class produces the omega_, w_0, w_a likelihood contours needed to calculate deltamu. Running Contour.py can be skipped, as deltamu.py will automatically check if the contour it needs exist, and it not produce it.

### Deltamu
Run Deltamu.py. This produces the deltamu values and the marginalisation constant K and saves them to file.

### Visualisation
Contains all plotting function to produce plots from the output of Deltamu.py.