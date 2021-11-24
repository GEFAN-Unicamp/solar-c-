//
//  read_files.hpp
//  solar_neutrino
//
//  Created by Dipyaman Pramanik on 23/08/21.
//

#ifndef read_files_hpp
#define read_files_hpp

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cstring>
#include <sstream>

typedef std::vector<double>  vec;

class file_reader {
private:
    int count_column(std::string);
    
    
public:
    int n_row;
    int n_col;
    vec *data;
    int col_count;
    file_reader();
    ~file_reader();
    int clean_data();
    int read_file(std::string);
};

#endif /* read_files_hpp */
