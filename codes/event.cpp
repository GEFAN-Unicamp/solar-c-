//
//  event.cpp
//  solar_neutrino
//
//  Created by Dipyaman Pramanik on 24/09/21.
//

#include "event.hpp"


namespace solar {





double elastic_cross(double E,double* par)
{
    double g1 = 0.73;
    double g2 = 0.23;
    
    double me = 0.511;
    
    double sigma_0 = 25*88.06e-46;
    
    double T_max = 2.0*square(E)/(me*me+2.0*E);
    
    double term1 = (square(g1)+square(g2))*(T_max-5.5);
    double term2 = (square(g2)+g1*g2*me/(2.0*E))*(square(T_max)-5.5*5.5)/E;
    double term3 = 1.0/3.0 *square(g2)*(T_max*T_max*T_max-5.5*5.5*5.5)/(E*E);
    
    
    
    double sigma = sigma_0/me*(term1-term2+term3);
    
 //   std::cout<<sigma<<"\n";
    
    return sigma;
}


double SNO::resFunc(double E, double E_p)
{
    double Sigma=0.0;
    
    if(smearing_matrix!=SOL_YES)
    {
        double alpha = resolution[0];
        double beta = resolution[1];
        double gamma = resolution[2];
        
        
        Sigma = alpha+beta*sqrt(E)+gamma*E;
    }
    if(fabs(Sigma)<1e-25)
    {
        std::cout<<"Solar error!! Resolution not set!!!\n";
        abort();
    }
    
    double res = 1.0/(sqrt(2.0*M_PI)*Sigma)*exp(-0.5*square((E-E_p)/Sigma));
    return res;
    
}

double SNO::Survival(double E)
{
    double Gamma_i;

    if(life>1e-16)
    {
        Gamma_i = 1/(E*life);
    }
    else
    {
        Gamma_i =0.0;

    }

  //  std::cout<<life<<std::endl;
    
    double p = exp(-Gamma_i*4.98e-4);


    return p;

}

double SNO::Propagation(double E)
{
        double Gamma_i;

        if(life>1e-16)
        {
            Gamma_i = 1/(life*E);

        }
        else
        {
            Gamma_i = 0.0;
        }

    
    
        double p = (1.0-exp(-Gamma_i*4.98e-4));
    
    
        return p;


}

int SNO::InitEventGenerator()
{
    
    InitFluxInterpolator();
    InitCrossInterpolator();
    
    
    double a = 5.3;
    double b = 15.0;
    
    n_samples = 30;
   
    create_bins();
    

    
    smear_mat = new double*[n_samples];
    for(int i=0;i<n_samples;i++)
    {
        smear_mat[i] = new double[n_bins];

    }

    sam_width = (b-a)/n_samples;
    
    
    for(int i=0;i<n_samples;i++)
    {
        smp_points.push_back(a+i*sam_width);
    }

    for(int i=0;i<n_samples;i++)
    {
        for(int j=0;j<n_bins;j++)
        {
            double _events = 0.0;

            _events = 0.5*(bin_f[j]-bin_i[j])*(resFunc(smp_points[i],bin_i[j])+resFunc(smp_points[i],bin_f[j]));
            smear_mat[i][j] = _events;
            
        }
        
    }

	

    
    
    return 0;
}

int SNO::FreeEventGenerator()
{
    for(int i=0;i<n_samples;i++)
    {
        delete smear_mat[i];
    }
    delete [] smear_mat;
    smp_points.clear();
    return 0;
}

int SNO::decayed_events()
{
    w.which_type = "Mixed";

    
    double c12_2 = square(cos(proba.osc_params[0]));
    double s12_2 = 1 - c12_2;
    
    life = proba.osc_params[6];
    double _delta = proba.osc_params[7];
    //~ life = 1e-4;

    
    double c13 = cos(proba.osc_params[1]);
    double c13_4 = c13*c13*c13*c13;
    
    
    
    
    for(int k=0;k<n_bins;k++)
    {
    
        double sum_inv=0.0;
        for(int i=0;i<n_samples-1;i++)
        {
            proba.Calculate_probability(smp_points[i]+off_set,4);
        
            double prob_inv1 = c13_4*(proba.pe1[4]*c12_2+proba.pe2[4]*s12_2*Survival(smp_points[i]+off_set));

            proba.Calculate_probability(smp_points[i+1]+off_set,4);
        
            double prob_inv2 = c13_4*(proba.pe1[4]*c12_2+proba.pe2[4]*s12_2*Survival(smp_points[i+1]+off_set));

        
            double event_inv1 = 0.0;
            double event_inv2 = 0.0;
            if((smp_points[i]+off_set)>=flux_min && (smp_points[i]+off_set)<=flux_max)
            {
                if((smp_points[i]+off_set)>=cross_min && (smp_points[i]+off_set)<=cross_max)
                {
                    event_inv1 = flux_interpolator.interpolate(smp_points[i]+off_set)*prob_inv1*cross_interpolator.interpolate(smp_points[i]+off_set)*smear_mat[i][k];
                    
//                    std::cout<<smear_mat[i][k]<<std::endl;
                    
                }

            }
            
            if((smp_points[i+1]+off_set)>=flux_min && (smp_points[i+1]+off_set)<=flux_max)
            {
                if((smp_points[i+1]+off_set)>=cross_min && (smp_points[i+1]+off_set)<=cross_max)
                {
            
                    event_inv2 = flux_interpolator.interpolate(smp_points[i+1]+off_set)*prob_inv2*cross_interpolator.interpolate(smp_points[i+1]+off_set)*smear_mat[i+1][k];
                }
            }
            
            
            double event_inv = 0.5*sam_width*(event_inv1+event_inv2);
            
            
            
            sum_inv = sum_inv+event_inv;
        
        }
   
        double sum_vis = 0.0;
        for(int i=0;i<n_samples;i++)
        {
            proba.Calculate_probability(smp_points[i],4);

            double prob_vis = 0.0;
            
            if((smp_points[i])>=flux_min && (smp_points[i])<flux_max)
            {
                prob_vis = flux_interpolator.interpolate(smp_points[i])*c13_4*proba.pe2[4]*c12_2*Propagation(smp_points[i]);
            }
            
            double ini_be = (smp_points[i])*square(_delta)-off_set;
            double fin_be = (smp_points[i])-off_set;
        
            double wid_bet = (fin_be-ini_be)/n_samples;
        
            vec smp_bet;
        
            for(int j=0;j<n_samples;j++)
            {
                smp_bet.push_back(ini_be+j*wid_bet);
            }
        
            double kern=0.0;
        
            for(int j=0;j<n_samples;j++)
            {
                double term = 0.0;
                
                if((smp_bet[j])>=cross_min && (smp_bet[j])<=cross_max)
                {
                    //std::cout<<(smp_bet[j]+off_set)<<"\t"<<cross_max<<std::endl;

                    term = w.weighted_rate(_delta,smp_points[i],smp_bet[j],1)*cross_interpolator.interpolate(smp_bet[j])*smear_mat[j][k];
                    
                }
                
                kern = kern + wid_bet*term;
                
            
            }
            
            smp_bet.clear();
            
           // std::cout<<prob_vis<<std::endl;
            
            double term = sam_width*prob_vis*kern;
            
        sum_vis = sum_vis + term;
        }
       
        double events = sum_inv+sum_vis;
        
        event_vector.push_back(Normalization*events);
    }
    
    

    return 0;
}

int SNO::InitFluxInterpolator()
{
    flux.read_file(flux_file);
    int size_v = int(flux.data[0].size()-2);
    flux_min = flux.data[0][0];
    flux_max = flux.data[0][size_v];
    
    vec data_x,data_y;
    flux_interpolator.set_interpolator(flux.data[0],flux.data[1], false);
    return 0;

}

int SNO::Set_datafile(std::string file_name)
{
    exp_data.read_file(file_name);
    
    return 0;
}

int SNO::InitCrossInterpolator()
{
    cross.read_file(cross_file);
    
    cross.read_file(cross_file);
    int size_v = int(cross.data[0].size()-2);
    cross_min = cross.data[0][0];
    cross_max = cross.data[0][size_v];
    
    
    cross_interpolator.set_interpolator(cross.data[0],cross.data[cross_ind],false);
    
    return 0;
    
}

int SNO::Set_fluxfile(std::string _file,int col_ind)
{
    flux_ind = col_ind;
    flux_file = _file;
    
    return 0;
}

int SNO::Set_crossfile(std::string _file, int col_ind)
{
    cross_ind = col_ind;
    cross_file = _file;
    
    cross_file_stat = SOL_YES;
    
    return cross_file_stat;
    
    
}


double SNO::Chi_nosys(vec params)
{
    proba.osc_params = params;
    
    
    decayed_events();
    
    double res = 0.0;
    
    for(int i=0;i<n_bins;i++)
    {
            res = res + poiss_likelihood(exp_data.data[1][i],event_vector[i]);
    }
    
    event_vector.clear();
    
    return res;
}


double SNO::Chi_sys(vec params)
{
    double res = 0.0;
    
    return res;
}


int SNO::create_bins()
{
    if(Man_bins!=SOL_YES)
    {
        bin_w = (e_max-e_min)/n_bins;
        
        for(int i=0;i<n_bins;i++)
        {
            bin_i.push_back(e_min+i*bin_w);
            bin_f.push_back(e_min+(i+1)*bin_w);
            
            bin_center.push_back(e_min+(0.5+i)*bin_w);
        }
    }
    if(Man_bins==SOL_YES)
    {
        for(int i=0;i<manual_bins.size()-1;i++)
        {
            bin_i.push_back(manual_bins[i]);
            bin_f.push_back(manual_bins[i+1]);
            bin_center.push_back((manual_bins[i]+manual_bins[i+1])/2.0);
            

        }

        n_bins = int(manual_bins.size()-1);

        e_min = manual_bins[0];
        e_max = manual_bins[manual_bins.size()-1];
    }

    
    return 0;

    
    
}

int SNO::Set_probability_engine(Probability _proba)
{
    proba = _proba;
    return 0;
}




double Borexino::Po_210_background(double E)
{
    return 80.0*exp(-square(E-430.0)/square(40.0));
}

double Borexino::beta(double E, double Q)
{
    if(E<Q)
    {
        return sqrt(E*(E+2.0*511.0))*square(E-Q);
    }
    
    
        return 0.0;
    
}


double Borexino::integrand(double T , double T_eff)
{
    double res = elastic_cross(T,862)*resFunc(T,T_eff);
    return res;
}

double Borexino::integrand(double E,double T , double T_eff)
{
    double res = elastic_cross(T,E)*resFunc(T,T_eff);
    
//    if(res<0)
//    {
//        std::cout<<"res = "<<res<<"; E_nu = "<<E<<"; T = "<<T<<"; T_eff = "<<T_eff<<std::endl;
//    }
    return res;
}
double Borexino::elastic_cross(double E, double E_nu)
{
    double g1 = -0.27;
    double g2 = 0.23;
    
    double me = 511;
    
    double sigma_0 = 88.06e-46;
    
    double T_max = 2.0*square(E_nu)/(me+2.0*E_nu);
    
    double term1 = square(g1);
    double term2 = square(g2)*square(1-E/E_nu);
    double term3 = g1*g2*me*E/square(E_nu);
    
    
    
    double sigma = sigma_0/me*(term1+term2-term3);
    
    if(E>T_max || E<200.0)
    {
         sigma = 0.0;
     }
    

    return sigma;
}

double Borexino::resFunc(double E, double E_p)
{
    
    
    double res = 1.0/(sqrt(2.0*M_PI)*Sigma)*exp(-0.5*square((E-E_p)/Sigma));
    return res;
}



double Borexino::standard_events(double T)
{
    double a = 10;
    double b = 1600.0;//2.0*square(862)/(me+2.0*862);

    int n_samples = 50;
   

    std::vector<double> smp_points;
;
    
    double sam_width = (b-a)/n_samples;
    
    for(int i=0;i<n_samples;i++)
    {
        smp_points.push_back(a+i*sam_width);
    }
    
    double sum = 0.0;
    
    for(int i=0;i<n_samples;i++)
    {
        double term=0.0;

            term = sam_width/2.0*(integrand(smp_points[i],T)+integrand(smp_points[i]+sam_width,T));
 
        sum = sum + term;
    }
    

    smp_points.clear();

    return sum;
    
}


double Borexino::decayed_events(double E, double T)
{
    double a = 200;
    double me = 511.;
    double b = 2.0*square(E)/(me+2.0*E);

    int n_samples = 50;
   

    std::vector<double> Smp_points;
    

    double Sam_width = (b-a)/n_samples;
    
    for(int i=0;i<n_samples;i++)
    {
        Smp_points.push_back(a+i*Sam_width);
    }
    
    double sum = 0.0;
    
    for(int i=0;i<n_samples;i++)
    {
        double term=0.0;
        if(b>a)
        {
            term = Sam_width/2.0*(integrand(E,Smp_points[i],T)+integrand(E,Smp_points[i]+Sam_width,T));
        
        }
        else
        {
            term = 0.0;
        }
        
        sum = sum + term;
    }
    
    Smp_points.clear();

    return sum;
    
}

int Borexino::InitEventGenerator()
{

    double a = 10;
    double b = 1600.0;
    
    int n_samples = 1000;
    n_samplings = 1000;
   
	vis_kern = new double*[n_samples];
	for(int i=0;i<n_samples;i++)
	{
		vis_kern[i] = new double[n_bins];

	}

	inv_kern = new double[n_bins];

    
    
    for(int j=0;j<n_bins;j++)
    {
        double _events = 0.0;

        _events = bin_w/2.0*(decayed_events(862.0,bin_i[j])+decayed_events(862.0,bin_f[j]));
        inv_kern[j] = _events;
    }

    
    sam_width = (b-a)/n_samples;
    
    for(int i=0;i<n_samples;i++)
    {
        smp_points.push_back(a+i*sam_width);
    }

	for(int i=0;i<n_samplings;i++)
	{
		for(int j=0;j<n_bins;j++)
		{
			double _events = 0.0;

			_events = bin_w/2.0*(decayed_events(smp_points[i],bin_i[j])+decayed_events(smp_points[i],bin_f[j]));
            
			vis_kern[i][j] = _events;
            
		}	
		
    
		
	}

    back.read_file(back_file);



	return 0;

}

int Borexino::FreeEventGenerator()
{
	delete inv_kern;
	for(int i=0;i<50;i++)
	{
		delete vis_kern[i];
	}

	delete vis_kern;
	return 0;

}



int Borexino::Reconstructed_events()
{
    double probability = proba.Calculate_probability(0.862,3);
    
    
    double c12_2 = square(cos(proba.osc_params[0]));
    double s12_2 = 1 - c12_2;
    
    life = proba.osc_params[6];
    double _delta = proba.osc_params[7];

    
    double c13 = cos(proba.osc_params[1]);
    double c13_4 = c13*c13*c13*c13;
    
    double prob_inv = c13_4*(proba.pe1[3]*c12_2+proba.pe2[3]*s12_2*Survival(0.862));
    
    
    double prob_vis = c13_4*proba.pe2[3]*c12_2*Propagation(0.862);
    
    w.which_type = "Mixed";
    


    
        for(int i=0;i<n_bins;i++)
        {
            double _events=0.0;
        
            
            _events = prob_inv*inv_kern[i];
        
            double vis_events = 0.0;
            
            for(int j=0;j<n_samplings-1;j++)
            {
                double E_beta = smp_points[j];
                
                double  term = 0.0;
                
                
               if(E_beta<862.0 && E_beta>=862.0*square(_delta))
                {

                    
                    vis_events = vis_events + prob_vis*sam_width/2.0*(w.weighted_rate(_delta,862,E_beta,1)*vis_kern[j][i]+w.weighted_rate(_delta,862,(E_beta+sam_width),1)*vis_kern[j+1][i]);
                    
              
                
                }
                else
                {
                    vis_events+=0.0;
                }
             //
            }
            


            _events=_events+vis_events;
            norm = 0.55e47;
            event_vector.push_back(norm*_events);//+Po_210_background(bin_center[i])+0.08e-8*beta(bin_center[i],1160.0));
            
            
        }


    return 0;
}

int Borexino::Reconstructed_events(double _po_s)
{
    double probability = proba.Calculate_probability(0.862,0);
    
    
    life = proba.osc_params[6];
    
    double c12_2 = square(cos(proba.osc_params[0]));
    double s12_2 = 1 - c12_2;
    
    double c13 = cos(proba.osc_params[1]);
    double c13_4 = c13*c13*c13*c13;
    
    double delta = proba.osc_params[7];
    
    double prob_inv = c13_4*(proba.pe1[3]*c12_2+proba.pe2[3]*s12_2*Survival(0.862));

    double prob_vis = c13_4*proba.pe2[3]*s12_2*Propagation(0.862);

	w.which_type = "Scalar";
    

	for(int i=0;i<n_bins;i++)
	{
		double _events =0.0;
		
		_events = prob_inv*bin_w/2.0*(integrate_cross(862.0,bin_i[i])+integrate_cross(862.0,bin_f[i]));
		
		double vis1=0.;
		double vis2=0.0;
		
        
        double a = 862.0*delta*delta;
		double b = 862.0;
		
        
		int n_samp = 50;
		vec samp;
		
		double width = (b-a)/n_samp;
		
		for(int j=0;j<n_samp+1;j++)
		{
			samp.push_back(a+j*width);
		}
		
		for(int j=0;j<n_samp;j++)
		{
            
                double term1=0.0;
                vis1  += width/2.0*(w.weighted_rate(delta,862.0,samp[j],1)*integrate_cross(samp[j],bin_i[i])+w.weighted_rate(delta,862.0,samp[j+1],1)*integrate_cross(samp[j+1],bin_i[i]));
                double term2 = 0.0;
                
                vis2 +=width/2.0*(w.weighted_rate(delta,862.0,samp[j],1)*integrate_cross(samp[j],bin_f[i])+w.weighted_rate(delta,862.0,samp[j+1],1)*integrate_cross(samp[j+1],bin_f[i]));
                

  
		}
		
        _events = _events + bin_w/2.0*(vis1+vis2);
		
        norm = 0.55e47;

		event_vector.push_back(norm*_events);
		
		samp.clear();
	}	
	
		



    return 0;
}

double Borexino::integrate_cross(double E_nu,double T)
{
	double me = 511.0;
	double a = 200.0;
    double b = 2.0*square(E_nu)/(me+2.0*E_nu);

	vec points;
	
	
	
	int n_points = 50;
	
	double _width = (b-a)/n_points;
	
	for(int i=0;i<n_points+1;i++)
	{
		points.push_back(a+i*_width);
	}
	
	double sum=0.0;
	for(int  i=0;i<n_points;i++)
	{
        
		double term = 0.0;
        
        if(b>a)
        {
            term = _width/2.0*(integrand(E_nu,points[i],T)+integrand(E_nu,points[i+1],T));
        }
		sum = sum +term;
	}


	points.clear();

	return sum;
}



int Borexino::Set_probability(Probability _proba)
{
    proba = _proba;
    return 0;
}

double Borexino::Chi_nosys(vec params)
{
    proba.osc_params = params;
    
    
    Reconstructed_events();
    
    double res = 0.0;
    

    
    
    for(int i=0;i<n_bins;i++)
    {
        if(bin_center[i]>350)
        {
            res = res + poiss_likelihood(570.4*back.data[1][i],570.4*event_vector[i]);
        }
    }
    
    event_vector.clear();
    
    return res;
}

double Borexino::Chi_sys(vec params,double Po_s)
{
    proba.osc_params = params;
    Reconstructed_events(Po_s);
    
    
    double res = 0.0;
    

    
    
    for(int i=0;i<n_bins;i++)
    {
        res = res + poiss_likelihood(back.data[1][i],event_vector[i]);

    }
    
    event_vector.clear();
    
    return res;
}




double Borexino::Survival(double E)
{
	double Gamma_i;

	if(life>1e-16)
	{
		Gamma_i = 1/(E*life);
	}
	else
	{
		Gamma_i =0.0;

	}

    
	double p = exp(-Gamma_i*4.98e-4);


	return p;

}

double Borexino::Propagation(double E)
{
		double Gamma_i;

		if(life>1e-16)
		{
			Gamma_i = 1/(life*E);
		}
		else
		{
			Gamma_i = 0.0;
		}

    
    
		double p = (1.0-exp(-Gamma_i*4.98e-4));
    
    
		return p;


}



double poiss_likelihood(double n_true, double n_test)
{
    double res;
    
    if(n_test<0)
    {
        n_test=0;
    }
    
    if(n_true>0)
    {
        res = 2*(n_test-n_true+n_true*log(n_true/n_test));
    }
    else
    {
        res = 2*n_test;
    }
    return res;
    
}

double gauss_likelihood(double n_true, double n_test)
{
    double res;
    
    if(n_test<0)
    {
        n_test=0;
    }
    
    if(n_true>0)
    {
        res =square(n_true-n_test)/n_true;
    }
    else
    {
        res = 2*n_test;
    }
    return res;
    
}

}

