//
//  numerical.cpp
//  solar_neutrino
//
//  Created by Dipyaman Pramanik on 23/08/21.
//

#include "numerical.hpp"
#include <gsl/gsl_spline.h>

int Cspline_interpolator::set_interpolator(vec x,vec y, bool logarithmic)
{
    int size_d = int(x.size())-1;
    
    
    double *x_d,*y_d;
    
    x_d = new double[size_d];
    y_d = new double[size_d];
    
    if(logarithmic!=true)
    {
        for(int i=0;i<size_d;i++)
        {
            x_d[i] = x[i];
            y_d[i] = y[i];
        }
        
        
    }
    else
    {
        for(int i=0;i<size_d;i++)
        {
            x_d[i] = x[i];
            y_d[i] = log10(y[i]);
        }
        
        
    }
    
    acc = gsl_interp_accel_alloc();
    spline =gsl_spline_alloc(gsl_interp_cspline,size_d);
    
    gsl_spline_init(spline, x_d, y_d,size_d);
    
    
    return 0;
}

int Cspline_interpolator::free_interpolator()
{
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    
    return 0;
}

double Cspline_interpolator::interpolate(double x)
{
    double inter=0;
    
    inter = gsl_spline_eval(spline,x,acc);
    
    return inter;
}



double square(double x)
{
    return x*x;
    
}
