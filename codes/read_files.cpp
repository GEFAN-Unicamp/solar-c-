//
//  read_files.cpp
//  solar_neutrino
//
//  Created by Dipyaman Pramanik on 23/08/21.
//

#include "read_files.hpp"


file_reader::file_reader()
{
    
}

file_reader::~file_reader()
{
    
    
}

int file_reader::clean_data()
{
    delete[] data;
    return 0;
}

int file_reader::count_column(std::string file_name)
{
    std::ifstream ifile;
    ifile.open(file_name);
    std::string line;
    
    std::getline(ifile,line);
    std::stringstream ss;
    ss<<line;
    
    int count = 0;
    double value;
    while(ss>>value) count++;
    
    col_count = count;
    if(count==0)
    {
        std::cout<<"ERROR!! Number of column is zero.\n";
        abort();
        return -1;
    }
    return 0;
    
}


int file_reader::read_file(std::string file_name)
{
    count_column(file_name);
    n_col = col_count;
    if(col_count==0)
    {
        std::cout<<"ERROR!! Invalid  file...\n";
        exit(-1);
        return  -1;
        
    }
    
    data = new vec[col_count];
    
    double *col;
    col = new double[col_count];
    
    std::ifstream ifile;
    ifile.open(file_name);
    
    while(ifile)
    {
        for(int i=0;i<col_count;i++)
        {
            ifile>>col[i];
            data[i].push_back(col[i]);
        }
    }
    
    n_row = int(data[0].size());
    
    delete [] col;
    
    ifile.close();
    
    for(int i=0;i<col_count;i++)
    {
        data[i][data[0].size()-1]=0;
    }
    
    if(data==NULL)
    {
        std::cout<<"ERROR!! File reading failed\n";
    }
    
    
    return 0;
}
