//
// Author : Toru Shiozaki
// Date   : May 2009
//

#ifndef __src_global_h
#define __src_global_h

#include <iostream>
#include <fstream>
#include <cassert>
#include <stdexcept>
#include <boost/regex.hpp>

void print_header() {
  std::cout << std::endl; 
  std::cout << "  ===============================================================" << std::endl;
  std::cout << "    Code name: White Rabit (since 2009)                          " << std::endl;
  std::cout << "        \"Oh dear! Oh dear! I shall be too late!\"               " << std::endl;
  std::cout << "  ===============================================================" << std::endl;
  std::cout << std::endl;
}


void print_footer() {
  std::cout << std::endl;
  std::cout << "  " << std::endl;
  std::cout << "  ===============================================================" << std::endl;
  std::cout << std::endl;
}


const int count_string(const std::string inputfile, const std::string keyword) {
  try {
    std::ifstream ifs;
    ifs.open(inputfile.c_str());
    if (!ifs.is_open())
      throw std::runtime_error("input file cannot be opened.");

    boost::smatch what;
    boost::regex reg(keyword);
    int out = 0;
    while(true) {
      std::string sline;
      if (!getline(ifs, sline)) break;
      std::string::const_iterator start = sline.begin();
      std::string::const_iterator end   = sline.end();
      if (boost::regex_search(start, end, what, reg)) ++out;
    }

    ifs.close();
    return out;
  } catch (...) {
    throw;
  }
}

#endif

