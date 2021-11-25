//
//  probability.hpp
//  solar_neutrino
//
//  Created by Dipyaman Pramanik on 24/08/21.
//

#ifndef probability_hpp
#define probability_hpp

#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <vector>
#include <functional>
#include "read_files.hpp"
#include "numerical.hpp"



typedef std::vector<double> vec;

namespace solar{

#define SOL_YES 1
#define SOL_NO 0

class Probability
{
private:
    int regeneration_earth();
    int interpolate_data();

    
public:
    int wrap_oscparams();
    int Init_probability_engine();
    double Calculate_probability(double,int);
    int free_data();
    int prob_inside_sun();

    
public:
    vec osc_params;
    double Energy;
    std::string prob_path;

    double *pday,*pmed;
    double Th12,Th13,Dm21,Dm31;
    double pe1[8],pe2[8];
private:
    vec *ff;
    double *rhosol;
    int fin_flav;
    int prod_data_n_row;
    
public:
    double L;
    double P_ij;

};

class weighted_differential
{
public:
    
    weighted_differential()
    {
        
    }
    
    weighted_differential(const char* _which)
        : which_type(_which)
    {
        
    }
    double weighted_rate(double,double,double,int);
public:
    std::string which_type;
};


class dec_probability: public Probability, public weighted_differential
{
public:
    dec_probability()
    {
        
        
    }
    
    dec_probability(int& _source,const char* which)
        : source(_source),_which(which),regeneration(SOL_NO)
    {
        which_type = which;
        L =1.0;
        
    }
    int init_interpolate_flux();

    double Calculate_decayed_flux(double,int,int,int);
private:
    int interpolate_data();
    double integrate_flux();
    double Propagation(double);
    double survival(double);
    double integrate();
    double integrand(double);
private:
    std::string _which;
    double x_max,E_max;
    double life;
    int particle;
    int _channel;
    int source;
    int regeneration;
    Cspline_interpolator flux_interpolator;
    file_reader flux;
    double flux_min,flux_max;
    
};



}

#endif /* probability_hpp */
