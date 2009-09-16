//
// Author: Toru Shiozaki
// Date  : April 2009
// 

#include <src/scf/shell.h>
#include <iostream>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;

Shell::Shell(const bool sph, vector<double> _position, int _ang, vector<double> _expo, 
                       vector<vector<double> > _contr,  vector<pair<int, int> > _range)
 : spherical_(sph), position_(_position), angular_number_(_ang), exponents_(_expo), contractions_(_contr), contraction_ranges_(_range) {

  for (vector<pair<int, int> >::iterator piter = _range.begin(); piter != _range.end(); ++piter) {
    contraction_lower_.push_back(piter->first);  
    contraction_upper_.push_back(piter->second);  
  }


  if (spherical_)
    nbasis_ = (angular_number_ * 2 + 1) * num_contracted();
  else
    nbasis_ = (angular_number_ + 1) * (angular_number_ + 2) / 2 * num_contracted();

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
  string out;
  out += "position: ";
  out += lexical_cast<string>(position_[0]) + " " + lexical_cast<string>(position_[1]) + " " +lexical_cast<string>(position_[2]) + "\n";
  out += "angular: " + lexical_cast<string>(angular_number_) + "\n"; 
  out += "exponents: ";
  for (int i = 0; i != exponents_.size(); ++i) {
    out += " " + lexical_cast<string>(exponents_[i]);
  }
  out += "\n"; 
  for (int i = 0; i != contractions_.size(); ++i) {
    out += " (" + lexical_cast<string>(contraction_ranges_[i].first) + "," + lexical_cast<string>(contraction_ranges_[i].second) + ") "; 
    for (int j = contraction_ranges_[i].first; j != contraction_ranges_[i].second; ++j)
      out += lexical_cast<string>(contractions_[i][j]) + " "; 
  } 
  return out;
}

