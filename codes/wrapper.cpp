
#include "wrapper.hpp"


namespace solar {

Wrapper_class::Wrapper_class()
{
	_Input = "input_card.txt";
	
	Init_wrapper();
	
}

Wrapper_class::Wrapper_class(std::string _input_file)
{
	_Input = _input_file;
	
	Init_wrapper();
	
	
}

int Wrapper_class::Init_wrapper()
{
	if(_Input.size()==0)
	{
		std::cout<<"Give the name of the input card ...\n";
		std::cin>>_Input;
		
	}
    std::ifstream ifl;
    ifl.open(_Input);
    if(!ifl)
    {
		std::cerr<<"Error!! File does not exists\n";
		return 0;
	}

    
  
	if(verb==2)
	{
		std::cout<<"Reading the input card..."<<std::endl;
	}
	while(ifl)
	{
		std::string s1,t1;
		
		getline(ifl,s1,'@');
		ifl>>t1;
		
		std::string s,t;
		getline(ifl,s,':');
		ifl>>t;
		
		if(t1=="verbosity")
		{
			verb = std::stoi(t);
			
		}
		
		if(t1=="name")
		{
			Experiment = t;
			if(verb==2)
			{
				std::cout<<"Current Experiment is : "<<Experiment<<std::endl;
			}
		}
		
		if(t1=="file_path")
		{
			file_path = t;
			if(verb==2)
			{
				std::cout<<"The path is : "<<file_path<<std::endl;
			}
			
		}
		
		if(t1=="out_file")
		{
			out_file = t;
			if(verb==2)
			{
				std::cout<<"The output will be saved to : "<<out_file<<std::endl;
			}
			
		}

		

		
		if(t1=="th12")
		{
			th12 = std::stod(t);
			if(verb==2)
			{
				std::cout<<"th12 = "<<th12<<std::endl;
			}
			
		}
		
		if(t1=="th13")
		{
			th13 = std::stod(t);
			if(verb==2)
			{
				std::cout<<"th13 = "<<th13<<std::endl;
			}
			
		}
		
		if(t1=="th23")
		{
			th23 = std::stod(t);
			if(verb==2)
			{
				std::cout<<"th23 = "<<th23<<std::endl;
			}
			
		}
		
		if(t1=="dcp")
		{
			dcp = std::stod(t);
			if(verb==2)
			{
				std::cout<<"delta_cp = "<<dcp<<std::endl;
			}
			
		}
		
		
		if(t1=="dm21")
		{
			dm21 = std::stod(t);
			if(verb==2)
			{
				std::cout<<"dm21 = "<<dm21<<std::endl;
			}
			
		}
		
		if(t1=="dm31")
		{
			dm31 = std::stod(t);
			if(verb==2)
			{
				std::cout<<"dm31 = "<<dm31<<std::endl;
			}
			
		}
		
		if(t1=="tau")
		{
			tau = std::stod(t);
			if(verb==2)
			{
				std::cout<<"tau = "<<tau<<std::endl;
			}
			
		}
		
		if(t1=="delta")
		{
			delta = std::stod(t);
			if(verb==2)
			{
				std::cout<<"delta = "<<delta<<std::endl;
			}
			
		}		
	}
	
	ifl.close();
	if(verb==2)
	{
		std::cout<<"The input file is closed\n";
	}
    
    
    
    params.push_back(th12*(M_PI/180.0));
    params.push_back(th13*(M_PI/180.0));
    params.push_back(th23*(M_PI/180.0));
    params.push_back(dcp*(M_PI/180.0));
    params.push_back(dm21);
    params.push_back(dm31);
    params.push_back(tau);
    params.push_back(delta);
    
	
	return 0;
}

int Wrapper_class::Calculate_probability()
{
	return 0;
}

int Wrapper_class::Calculate_Event()
{
	Probability proba;
	proba.prob_path = file_path;
	proba.Init_probability_engine();
	proba.osc_params = params;
    
	
	if(Experiment=="Borexino")
	{
		
		Borexino _borex;
		_borex.Set_probability(proba);
		_borex.proba.osc_params = params;
		_borex.back_file = file_path + "/exp_data/Borexino_be.dat";
		_borex.InitEventGenerator();

		std::ofstream ofl;
		
		std::string output = file_path + "/data_files/event_Borexino_" + out_file;
		ofl.open(output);

   
		_borex.Reconstructed_events();

		for(int i=0;i<_borex.n_bins;i++)
		{
			if(verb>1)
			{
				std::cout<<_borex.bin_center[i]<<"\t"<<_borex.event_vector[i]<<std::endl;
			}
			ofl<<_borex.bin_center[i]<<"\t"<<_borex.event_vector[i]<<std::endl;
			
		}

		ofl.close();
        
        std::cout<<"The output is written into ../data_files/event_Borexino_"<<out_file<<std::endl;
        
	}
	else if(Experiment=="SNO")
	{
	    std::string exp_file = file_path + "/exp_data/SNO_CC_data.dat";
    
    
		std::string flux_file_name = file_path + "/flux/b8spectrum.txt";
		std::string cross_file_name = file_path + "/cross_section/SNO_NueDcs.dat";
    
		double _bins[17] = {5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0,12.5,13.0,13.5};

		file_reader exp_data;
		exp_data.read_file(exp_file);

    
		SNO _sno;
    
		for(int i=0;i<17;i++)
		{
			_sno.manual_bins.push_back(_bins[i]);
		}
    
		_sno.resolution[0] = -0.4620;
		_sno.resolution[1] = 0.5470;
		_sno.resolution[2] = 0.008722;
		_sno.off_set = 1.44+1.0;//+2.0;//-1.1;
		_sno.Set_probability_engine(proba);
		_sno.Set_fluxfile(flux_file_name,1);
		_sno.Set_crossfile(cross_file_name,3);
    
		_sno.Normalization = 7.8e45;

		_sno.InitEventGenerator();

	    _sno.decayed_events();
    
		std::ofstream ofl;
		
		std::string output = file_path + "/data_files/event_SNO_" + out_file;
		ofl.open(output);
		//std::string _data = file_path + "exp_data/SNO_CC_data.dat";
    
		//_sno.Set_datafile(_data);
    
		
	   for(int i=0;i<_sno.n_bins;i++)
 	   {
		  std::cout<<_sno.bin_center[i]<<"\t"<<_sno.event_vector[i]<<std::endl;
		  ofl<<_sno.bin_center[i]<<"\t"<<_sno.event_vector[i]<<std::endl;

	   }
		
	   std::cout<<"The output is written into ../data_files/event_SNO_"<<out_file<<std::endl;
       
		
	}
	else
	{
		std::cerr<<"Error!! Invalid experiment\n";
		abort();
	}
	
	return 0;
}

int Wrapper_class::Calculate_1d_chi2()
{
	return 0;
}

int Wrapper_class::Calculate_2d_chi2()
{
	return 0;
}

}
