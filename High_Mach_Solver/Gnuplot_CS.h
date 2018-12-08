#ifndef GNUPLOT_CS_H
#define GNUPLOT_CS_H

#include <iostream>
#include <fstream>
#include <cstdarg>
#include <vector>

template <typename data_type>
void gnuplot_CS(std::string file_name, std::vector<std::vector<data_type>> data, std::string first_line){

  std::ofstream gnu_input_file;
  gnu_input_file.open (file_name + ".dat");
  gnu_input_file << first_line << std::endl;
  for(size_t i = 0; i < data[0].size(); ++i) {
    for(size_t j = 0; j < data.size();++j) {
        gnu_input_file << data[j][i] << " ";
      }
      gnu_input_file << std::endl;
  }
  gnu_input_file.close();
}

#endif
