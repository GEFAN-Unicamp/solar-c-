//
//  probability.cpp
//  solar_neutrino
//
//  Created by Dipyaman Pramanik on 24/08/21.
//

#include "probability.hpp"

namespace solar {



int Probability::wrap_oscparams()
{
    osc_params.clear();
    
    osc_params.push_back(Th12);
    osc_params.push_back(Th13);
    osc_params.push_back(0.0);
    osc_params.push_back(0.0);
    osc_params.push_back(Dm21);
    osc_params.push_back(Dm31);
    
    
    return 0;
}

int Probability::Init_probability_engine()
{
	
    interpolate_data();
    
    pday = new double[8];
    pmed = new double[8];
    
    
    return 0;
    
}

int Probability::interpolate_data()
{

    file_reader distri;
    std::string distri_path = prob_path + "/prod_dist/bs05opflux2.dat";
    distri.read_file(distri_path);

    
    vec *prod_data;
    int distri_n_col = distri.n_col;
    prod_data = new vec[distri_n_col];
    
    for(int i=0;i<distri.n_row-1;i++)
    {
        for(int j=0;j<distri.n_col;j++)
        {
            prod_data[j].push_back(distri.data[j][i]);

        }
    }
    
    
    
    distri.clean_data();
    
    
    Cspline_interpolator *inter_prod;
    inter_prod = new Cspline_interpolator[distri_n_col];
    
    for(int i=1;i<distri_n_col;i++)
    {
        inter_prod[i].set_interpolator(prod_data[0], prod_data[i], false);
    }
    
    file_reader rho_sol;
    std::string rho_path = prob_path + "/solar_density/bs2005op.dat";
    rho_sol.read_file(rho_path);
    
    vec rho_sol_x,rho_sol_y;
    
    rho_sol_x.push_back(0.0);
    rho_sol_y.push_back(rho_sol.data[1][0]);
    
    
    for(int i=0;i<rho_sol.data[0].size()-1;i++)
    {
        rho_sol_x.push_back(rho_sol.data[0][i]);
        rho_sol_y.push_back(rho_sol.data[1][i]);
        
    }
    
    Cspline_interpolator sol_density;
    sol_density.set_interpolator(rho_sol_x, rho_sol_y, false);
    
    prod_data_n_row = int(prod_data[0].size());
    rhosol = new double[prod_data_n_row];
    
    ff = new vec[8];
    
    for(int i=1;i<prod_data_n_row-1;i++)
    {
        rhosol[i] = sol_density.interpolate(prod_data[0][i]);//*prod_data[0][i]);
        
        for(int j=0;j<8;j++)
        {
            double X = prod_data[0][i];
            ff[j].push_back(inter_prod[j+1].interpolate(fabs(X)));

        }
        
    }
    
    rho_sol.clean_data();
    
    delete [] inter_prod;
    delete [] prod_data;
    
    
    return 0;
}

double Probability::Calculate_probability(double E, int flav)
{
    Energy = E;
    fin_flav = flav;
    prob_inside_sun();
    
    return  pday[4];
}


int Probability::prob_inside_sun()
{
    
 //   wrap_oscparams();

    
    double c12 = cos(osc_params[0]);
    double s12 = sin(osc_params[0]);
    
    double c13 = cos(osc_params[1]);
 //   double s13 = sin(osc_params[1]);
    
  //  double s2_th13 = s13*s13;
    
    double c_2th12 = c12*c12-s12*s12;
    double s_2th12 = 1-c_2th12;
    
    
    
    double *fmed;
    fmed = new double[8];
    pmed = new double[8];
    double pe1med[8];
    double pe2med[8];
    double pmed1[8];
    
    
    for(int i=0;i<8;i++)
    {
        fmed[i] = 0.0;
        pmed[i] = 0.0;
        pe1med[i]=0.0;
        pe2med[i]=0.0;
        pmed1[i]=0.0;
    }
    
    
    double E = Energy;
    
    for(int i=0;i<prod_data_n_row-1;i++)
    {
        double temp_rhosol = square(c13)*pow(10,rhosol[i]);
        
        double aa = osc_params[4]*c_2th12-4.0*E*1e6*temp_rhosol;
        double bb = osc_params[4]*s_2th12;
        
        double c_2tm = aa/(sqrt(aa*aa+bb*bb));
        
        double s_tm = 0.5*(1-c_2tm);
        double c_tm = 0.5*(1+c_2tm);
//
        double pee = 0.5*(1+c_2tm*c_2th12);
        
        
        for(int j=0;j<8;j++)
        {
            fmed[j] = fmed[j] + ff[j][i];
            pe1med[j] = pe1med[j] + ff[j][i]*c_tm;
            pe2med[j] = pe2med[j] + ff[j][i]*s_tm;
            pmed[j] = pmed[1] + ff[j][i]*pee;
            

        }
        
    }

    
    
    for(int i=0;i<8;i++)
    {
        pday[i] = pmed[i]/fmed[i];
        pe1[i] = pe1med[i]/fmed[i];
        pe2[i] = pe2med[i]/fmed[i];
    }
    
    delete [] fmed;
    delete [] pmed; 
    return 0;
}


int Probability::free_data()
{
    delete [] ff;
    delete [] rhosol;
    delete [] pday;
    delete [] pmed;
    
    
    return 0;
}


int dec_probability::init_interpolate_flux()
{
    if(source==4)
    {
		std::string flux_path = prob_path + "/flux/b8spectrum.txt";
        flux.read_file(flux_path);
    }
    
    int size_v = int(flux.data[0].size()-2);
    flux_min = flux.data[0][0];
    flux_max = flux.data[0][size_v];
    
    flux_interpolator.set_interpolator(flux.data[0],flux.data[1], false);
    
    return 0;
}


double dec_probability::integrate()
{
    x_max = flux_max;
    E_max = Energy/square(osc_params[7]);
    
    if(E_max>x_max)
    {
        E_max = x_max;
    }
    
    double h = (E_max-Energy)/4.0;
    
    
    double integral = (integrand(Energy)+integrand(E_max));
    
    for(int i=0;i<3;i++)
    {
        double term = 2.0*integrand(Energy+(i+1)*h);
        integral = integral + term;
    }
    
    return h/2.0*integral;
}

double dec_probability::integrand(double E)
{
    pe2[source] = 0;
    pe1[source] = 0;
    prob_inside_sun();
    
    double E_beta = Energy;
    double E_alpha = E;
    
    
    /*ONLY BORON IS BEING CONSIDERED FOR SIMPLICITY. LATER WE CAN ADD ALL SOURCES*/
    
    double res=0.0;
    
//    double w;

  /*  if(particle==1)
    {
        w = E_beta/square(E_alpha);
    }
    else
    {
        w = 1/E_alpha*(1-E_beta/E_alpha);
    }
*/
    if(E<x_max)
    {
        if(_channel==1)
        {
            res = flux_interpolator.interpolate(E_alpha)*pe2[source]*Propagation(E_alpha)*weighted_rate(osc_params[7],E_alpha,E_beta,particle);
        }
        else if(_channel==2)
        {
            res = flux_interpolator.interpolate(E_alpha)*Propagation(E_alpha)*weighted_rate(osc_params[7],E_alpha,E_beta,particle);
        }
        else if(_channel==3)
        {
        
     //       std::cout<<"Energy  =====  "<<E<<"\n";

            res = flux_interpolator.interpolate(E_alpha)*Propagation(E_alpha)*weighted_rate(osc_params[7],E_alpha,E_beta,particle);
        }
        else if (_channel==4)
        {
            res = flux_interpolator.interpolate(E_alpha)*pe2[source]*Propagation(E_alpha)*weighted_rate(osc_params[7],E_alpha,E_beta,particle);


        }
        else if(_channel==5)
        {

            res = flux_interpolator.interpolate(E_alpha)*pe1[source]*Propagation(E_alpha)*weighted_rate(osc_params[7],E_alpha,E_beta,particle);
        }
       // std::cout<<w<<"\t"<<W_rate.weighted_rate(osc_params[3],E_alpha,E_beta,particle)<<"\n";
        
     //   std::cout<<"Res = "<<res<<std::endl;
    }
    else
    {
        res = 0.0;
    }
    //~ std::cout<<osc_params[3]<<"\n";
    
    
    
    return res;
}


double dec_probability::Calculate_decayed_flux(double E, int fin_flav,int parti,int channel)
{
    _channel = channel;
    
    Energy = E;
    particle = parti;
    life = osc_params[6];

    
    prob_inside_sun();
    
    double c13 = cos(osc_params[1]);
    double s13 = sin(osc_params[1]);
    double s12 = sin(osc_params[0]);
    double c12 = cos(osc_params[0]);
    
    double s13_4 = s13*s13*s13*s13;
    double c13_4 = c13*c13*c13*c13;
    
    //double s13_2 = s13*s13;
    double c13_2 = c13*c13;
    
    std::cout<<c12<<std::endl;
    
    double s12_2 = s12*s12;
    double c12_2 = c12*c12;
    
    double c13_2_s13_2 = c13*c13*s13*s13;
    
    double Flux_dec=0;
    
    if(channel==1) /* nu2->nu1 */
    {
        if(fin_flav==0)
        {
            if(parti==1)
            {
            
                double P_inv_prime = pe1[source]*c12_2+pe2[source]*s12_2*survival(Energy);
            
                double flux_0 = flux_interpolator.interpolate(Energy);
                double flux_inv_prime = flux_0*P_inv_prime;
            
                double flux_inv = c13_4*flux_inv_prime + s13_4*flux_0;
                
                double flux_vis=0;
                
                if(regeneration!=SOL_YES)
                {
                    flux_vis = c13_4*c12_2*integrate();
                    
                }
                
                Flux_dec = flux_inv + flux_vis;
                
            
            }
            if(parti==-1)
            {
                double flux_vis=0;
                
                if(regeneration!=SOL_YES)
                {
                    flux_vis = c13_4*c12_2*integrate();
                }
                
                Flux_dec = flux_vis;
            }
            
            
        }
    }
    else if(channel==2) /* nu3->nu2 */
    {
        if(fin_flav==0)
        {
            if(parti==1)
            {
                double P_inv_prime = pe1[source]*c12_2+pe2[source]*s12_2;
                    
                double flux_0 = flux_interpolator.interpolate(Energy);
                    
                double flux_inv_prime = flux_0*P_inv_prime;
                    
                double flux_inv = c13_4*flux_inv_prime + s13_4*survival(Energy);
                    
                double flux_vis=0;
                    
                if(regeneration!=SOL_YES)
                {
                    flux_vis = c13_2_s13_2*s12*integrate();
                }
                    
                Flux_dec = flux_inv + flux_vis;
            }
            if(parti==-1)
            {
                double flux_vis=0;
            
                if(regeneration!=SOL_YES)
                {
                    flux_vis = c13_2_s13_2*s12_2*integrate();
                }
                    
                Flux_dec = flux_vis;
                    
            }
        }
    }
    else if(channel==3)  /* nu3->nu1 */
    {
        if(fin_flav==0)
        {
            if(parti==1)
            {
                double P_inv_prime = pe1[source]*c12_2 + pe2[source]*s12_2;
                
                double flux_0 = flux_interpolator.interpolate(Energy);
                    
                double flux_inv_prime = flux_0*P_inv_prime;
                    
                double flux_inv = c13_4*flux_inv_prime+s13_4*survival(Energy)*flux_0;
                    
                double flux_vis=0;
                    
                if(regeneration!=SOL_YES)
                {
                    flux_vis = c13_2_s13_2*c12_2*integrate();
                }
                
                Flux_dec = flux_inv + flux_vis;
            }
            if(parti==-1)
            {
                double flux_vis=0;
                    
                if(regeneration!=SOL_YES)
                {
                    flux_vis = c13_2_s13_2*c12_2*integrate();
                }
                    
                Flux_dec = flux_vis;
            }
        }
    }
    else if(channel==4)
    {
        if(fin_flav==0)
        {
            if(parti==1)
            {
                double P_inv_prime = pe1[source]*c12_2 + pe2[source]*s12_2*survival(Energy);
                
                double flux_0  = flux_interpolator.interpolate(Energy);
                    
                double flux_inv_prime = flux_0*P_inv_prime;
                    
                double flux_inv = c13_4*flux_inv_prime + s13_4*flux_0;
                    
                double flux_vis=0;
                    
                if(regeneration!=SOL_YES)
                {
                    flux_vis = c13_2_s13_2*integrate();
                }
                    
                Flux_dec = flux_inv + flux_vis;
            }
            if(parti==-1)
            {
                if(regeneration!=SOL_YES)
                {
                    Flux_dec = c13_2_s13_2*integrate();
                }
            }
        }
    }
    else if(channel==5)
    {
        if(fin_flav==0)
        {
            if(parti==1)
            {
                double P_inv_prime = pe1[source]*c12_2*survival(Energy)+pe2[source]*s12_2;
                
                double flux_0 = flux_interpolator.interpolate(Energy);
                
                double flux_inv_prime = flux_0*P_inv_prime;
                
                double flux_inv = c13_4*flux_inv_prime + s13_4*flux_0;
                
                double flux_vis=0;
                
                if(regeneration!=SOL_YES)
                {
                    flux_vis = c13_2_s13_2*integrate();
                }
                
                Flux_dec = flux_inv + flux_vis;
            }
            if(parti==-1)
            {
                if(regeneration!=SOL_YES)
                {
                    Flux_dec = c13_2_s13_2*integrate();
                }
            }
        }
    }
    else
    {
        std::cerr<<"ERROR!! Invalid Channel..."<<std::endl;
    }
    
    
    
    return Flux_dec;
    
    
}


double dec_probability::survival(double E)
{
    double Gamma_i;
    
    if(life>1e-16)
    {
        Gamma_i = 1/(E*life);
    }
    else
    {
        Gamma_i = 0;
    }
    
    double p = exp(-Gamma_i*L*4.98e-4);
    
    return p;
}

double dec_probability::Propagation(double E)
{
    double Gamma_i;
    
    if(life>1e-16)
    {
        Gamma_i = 1/(E*life);
        
        
    }
    else
    {
        Gamma_i = 0.0;
    }
    
    double p = (1.0-exp(-Gamma_i*L*4.98e-4));
    
    P_ij = p;
    
    return p;
}


double weighted_differential::weighted_rate(double delta,double E_alpha,double E_beta,int parti)
{
    double w=0;


    
    if(parti==1)
    {
        double factor2 = 1.0/E_alpha;
        double factor3 = (1.0+square(delta)-E_beta/E_alpha-square(delta)*E_alpha/E_beta);
        
        double term1 = 1.0/((1-square(delta))*E_alpha);
        
        if(which_type=="Scalar")
        {
            double factor1 = 1.0/((1.0-square(delta))*square(1.0+delta));
            
//            std::cout<<term1<<"\t"<<factor<<"\n";

            w = term1 - factor1*factor2*factor3;
            
//            std::cout<<E_beta/sqr(E_alpha)<<"\t"<<(1/E_alpha - factor1*factor2*factor3
  //                  )<<"\n";
        }
        else if(which_type=="Pseudo")
        {
            double factor1 = 1.0/((1-square(delta))*square(1-delta));
            
            w = term1 - factor1*factor2*factor3;
            
        }
        else if(which_type=="Mixed")
        {
            double factor1 = 1.0/((1-square(delta*delta)));
            
            w = term1 - factor1*factor2*factor3;
        }
        else
        {
            std::cerr<<"ERROR!! Invalid option\n";
            abort();
        }
        
    }
    else if(parti==-1)
    {

        double factor2 = 1.0/E_alpha;
        double factor3 = (1.0+square(delta)-E_beta/E_alpha-square(delta)*E_alpha/E_beta);
        
    
        
        if(which_type=="Scalar")
        {
            double factor1 = 1.0/((1.0-square(delta))*square(1.0+delta));
            
            w = factor1*factor2*factor3;
            

      //      std::cout<<E_alpha/E_beta<<"\t"<<E_alpha<<"\t"<<E_beta<<"\t"<<factor3<<"\n";
        }
        else if(which_type=="Pseudo")
        {
            double factor1 = 1.0/((1-square(delta))*square(1-delta));
            
            w = factor1*factor2*factor3;
            
        }
        else if(which_type=="Mixed")
        {
            double factor1 = 1.0/((1-square(delta*delta)));
            
            w = factor1*factor2*factor3;
        }
        else
        {
            std::cerr<<"ERROR!! Invalid option\n";
            abort();
        }
        

    }
    
    
    return w;
}




}
