//
//  numerical.hpp
//  solar_neutrino
//
//  Created by Dipyaman Pramanik on 23/08/21.
//

#ifndef numerical_hpp
#define numerical_hpp

#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>
#include <functional>
#include <cmath>
#include <gsl/gsl_spline.h>

typedef std::vector<double> vec;

double square(double);

class Cspline_interpolator
{
private:
    gsl_interp_accel *acc;
    gsl_spline *spline;
    
public:
    int free_interpolator();
    int set_interpolator(vec,vec,bool);
    double interpolate(double);
};

#endif /* numerical_hpp */
