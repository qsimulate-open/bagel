//
// Author : Toru Shiozaki
// Date   : August 2009
//

#include <src/util/filename.h>
#include <sstream>

using namespace std;

static unsigned int counter_ = 0;

Filename::Filename() {

}


Filename::~Filename() {

}


const string Filename::filename_next() const {
  stringstream ss; 
  ss << "temp_file_" << counter_ << ".data";
  ++counter_;

  return ss.str();
} 

