#ifndef wrapper_hpp
#define wrapper_hpp


#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <vector>
#include <functional>
#include "read_files.hpp"
#include "numerical.hpp"
#include "probability.hpp"
#include "event.hpp"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

namespace solar {


class Wrapper_class
{
private:
    double th12,th13,th23,dcp,dm21,dm31,tau,delta;
    
    
    std::string file_path;
    std::string out_file;
    std::string Experiment;
    
    std::string _Input;
    int verb;
    vec params;


public:
	Wrapper_class();
	Wrapper_class(std::string);
	
	int Calculate_probability();
	int Calculate_Event();
	int Calculate_1d_chi2();
	int Calculate_2d_chi2();
	
private:
	int Init_wrapper();
	
};

}

#endif /* wrapper_hpp */
