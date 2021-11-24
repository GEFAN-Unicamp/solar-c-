//
//  main.cpp
//  solar_neutrino
//
//  Created by Dipyaman Pramanik on 23/08/21.
//

#include <iostream>
#include "wrapper.hpp"


using namespace solar;




int main(int argc, const char * argv[]) {

	Wrapper_class _wrap;
	
	char option;
   
	std::cout<<"Please put the correct option ...\n";
	std::cout<<"A : Probability, B : Event Spectra, C: 1d chi2, D: 2d chi2\n";
	std::cin>>option;

	
	switch(option)
	{
		case 'A' :
				std::cout<<"Calculating probability...\n";
				_wrap.Calculate_probability();
			break;
			
		case 'B' :
				std::cout<<"Calculating Event spectra...\n";
				_wrap.Calculate_Event();
			break;
			
		case 'C' :
				std::cout<<"Calculating 1d chi2 projection...\n";
				_wrap.Calculate_1d_chi2();
			break;
			
		case 'D' :
				std::cout<<"Calculating 2d chi2 projection...\n";
				_wrap.Calculate_2d_chi2();
			break;
			
		default:
			std::cout<<"Invalid option...\n";
			abort();
		
	}
	


    
    return 0;
}


