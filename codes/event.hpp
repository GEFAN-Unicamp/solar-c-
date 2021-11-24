//
//  event.hpp
//  solar_neutrino
//
//  Created by Dipyaman Pramanik on 24/09/21.
//

#ifndef event_hpp
#define event_hpp


#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <vector>
#include <functional>
#include "read_files.hpp"
#include "numerical.hpp"
#include "probability.hpp"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

namespace solar {



#define SOL_YES 1
#define SOL_NO 0

double poiss_likelihood(double,double);
double gauss_likelihood(double,double);


class SNO
{
public:
    
    int InitEventGenerator();
    int FreeEventGenerator();
    
    int Set_fluxfile(std::string,int);
    int Set_crossfile(std::string,int);
    int Set_datafile(std::string);
    int Set_probability_engine(Probability);
    int decayed_events();
    double ReconsEvents();
    double Chi_nosys(vec);
    double Chi_sys(vec);
    SNO()
    {
        smearing_matrix = SOL_NO;
        Man_bins = SOL_YES;
    }
    
    
    
private:
    int InitFluxInterpolator();
    int InitCrossInterpolator();
    int SetProbability();
    int create_bins();
    double resFunc(double,double);
    

    double Propagation(double);
    double Survival(double);
    double Kernel_recons(double);
    
    
public:
    int smearing_matrix,Man_bins;
    Probability proba;
    double resolution[3];
    int n_bins;
    double bin_0,bin_w,bin_n,norm;
    double e_min,e_max;
    vec bin_i,bin_f,bin_center,event_vector;
    vec manual_bins;
    double Normalization;
    double off_set;

private:
    int cross_file_stat;
    Cspline_interpolator flux_interpolator;
    Cspline_interpolator cross_interpolator;
    file_reader back,exp_data;
    file_reader flux,cross;
    double life;
    double flux_min,flux_max;
    double cross_min,cross_max;
    double **smear_mat;
    int cross_ind;
    int flux_ind;
    std::string flux_file;
    std::string cross_file;
    int n_samples;
    vec smp_points;
    double sam_width;
    weighted_differential w;
    
    
};

class Borexino
{
public:
    Borexino()
        :n_bins(140),bin_0(200),bin_n(1600),_standard(SOL_YES),norm(0.52e47),n_samplings(1000),Sigma(50)
    {
        bin_w = (bin_n-bin_0)/n_bins;
        
        for(int i=0;i<n_bins;i++)
        {
            bin_i.push_back(bin_0+i*bin_w);
            bin_f.push_back(bin_0+(i+1)*bin_w);
            bin_center.push_back(bin_0+(i+0.5)*bin_w);
        }

//
        back.read_file("/home/dipyaman/Documents/solar/solar_nu/exp_data/Borexino_be.dat");
        
    }
    
//    ~Borexino()
//    {
//        delete smear_mat;
//        samplings.clear();
//    }
    
    Probability proba;
    double standard_events(double);
    double decayed_events(double,double);
    int InitEventGenerator();
    int FreeEventGenerator();
    int Set_probability(Probability);
    int Reconstructed_events();
    int Reconstructed_events(double);

    double Po_210_background(double);
    double beta(double,double);
    double Chi_nosys(vec);
    double Chi_sys(vec,double);
    
    
private:
    double resFunc(double, double);
    double elastic_cross(double,double);
    double integrand(double,double);
    double integrand(double,double,double);
    double integrate_cross(double,double);
    file_reader back;


    double Propagation(double);
    double Survival(double);
    double Kernel_recons(double);   


    
    
public:
    int n_bins;
    double bin_0,bin_w,bin_n,norm;
    vec bin_i,bin_f,bin_center,event_vector;
    
private:
    int _standard;
    int n_samplings;
    double** smear_mat;
    vec samplings;
    double Sigma;
    double life;
    double* inv_kern;
    double** vis_kern;
    weighted_differential w;
    std::vector<double> smp_points;
    double sam_width;

};


}


#endif /* event_hpp */
