# solar-c-
Updated solar neutrino code

This is the code for the simulation of SNO and Borexino data.

The codes are inside the codes folder. There is one dependency, it is 
gsl. gsl has been used for interpolation. 

To compile the code you need to run make. 

To run the code you need to run ./solar

You can edit the input parameters via a default card file : "input_card.txt"

To run into your machine you need to set the path of the repository

All the output will be saved to a folder called data_files inside the main folder 

Sample input_card.txt file below

Presently the code is only giving the event spectra. Developing the 
chi2 code now...
