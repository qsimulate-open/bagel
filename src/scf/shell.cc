//
// Author: Toru Shiozaki
// Date  : April 2009
// 

#include <src/scf/shell.h>
#include <iostream>
#include <sstream>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;

Shell::Shell(const bool sph, vector<double> _position, int _ang, vector<double> _expo, 
                       vector<vector<double> > _contr,  vector<pair<int, int> > _range)
 : spherical_(sph), position_(_position), angular_number_(_ang),
   exponents_(_expo), contractions_(_contr), contraction_ranges_(_range) {

  for (vector<pair<int, int> >::iterator piter = _range.begin(); piter != _range.end(); ++piter) {
    contraction_lower_.push_back(piter->first);  
    contraction_upper_.push_back(piter->second);  
  }


  if (spherical_)
    nbasis_ = (angular_number_*2+1) * num_contracted();
  else
    nbasis_ = (angular_number_+1) * (angular_number_+2) / 2 * num_contracted();

}


Shell::~Shell() {

}


boost::shared_ptr<Shell> Shell::move_atom(const vector<double>& displacement) {
  boost::shared_ptr<Shell> out(new Shell(*this));
  out->position_[0] += displacement[0]; 
  out->position_[1] += displacement[1]; 
  out->position_[2] += displacement[2]; 
  return out;
}


boost::shared_ptr<Shell> Shell::move_atom(const double* displacement) {
  boost::shared_ptr<Shell> out(new Shell(*this));
  out->position_[0] += displacement[0]; 
  out->position_[1] += displacement[1]; 
  out->position_[2] += displacement[2]; 
  return out;
}


const string Shell::show() const {
  stringstream ss;
  ss << "position: ";
  ss << position_[0] << " " << position_[1] << " "  << position_[2] << endl;
  ss << "angular: "  << angular_number_ << endl;
  ss << "exponents: ";
  for (int i = 0; i != exponents_.size(); ++i) {
    ss << " " << exponents_[i];
  }
  ss << endl;
  for (int i = 0; i != contractions_.size(); ++i) {
    ss << " (" << contraction_ranges_[i].first << "," << contraction_ranges_[i].second << ") ";
    for (int j = contraction_ranges_[i].first; j != contraction_ranges_[i].second; ++j)
      ss << contractions_[i][j] << " ";
  }

  return ss.str();
}

