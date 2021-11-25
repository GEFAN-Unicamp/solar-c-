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

#############################################################################
#					Se5t level of verbosity
############################################################################

@verbosity : 2

#############################################################################
#					Information about path and output
############################################################################
@file_path : /home/dipyaman/Documents/Github/solar-c-
@out_file : test_out.dat

#############################################################################
#					Initial oscillation parameters
############################################################################
@th12 : 33.56
@th13 : 8.5
@th23 :0
@dcp : 0
@dm21 : 7.5e-5
@dm31 : 0

@tau : 1e-6
@delta : 0.9

#############################################################################
#					Which Experiment
############################################################################
@name : Borexino

#############################################################################
#					Borexino Parameters
############################################################################


#############################################################################
#					SNO Parameters
############################################################################
