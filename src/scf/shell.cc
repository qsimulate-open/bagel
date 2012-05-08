//
// Newint - Parallel electron correlation program.
// Filename: shell.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include <src/scf/shell.h>
#include <iostream>
#include <sstream>

using namespace std;

Shell::Shell(const bool sph, vector<double> _position, int _ang, vector<double> _expo, 
                       vector<vector<double> > _contr,  vector<pair<int, int> > _range)
 : spherical_(sph), position_(_position), angular_number_(_ang),
   exponents_(_expo), contractions_(_contr), contraction_ranges_(_range), dummy_(false) {

  for (auto piter = _range.begin(); piter != _range.end(); ++piter) {
    contraction_lower_.push_back(piter->first);  
    contraction_upper_.push_back(piter->second);  
  }

  if (spherical_)
    nbasis_ = (angular_number_*2+1) * num_contracted();
  else
    nbasis_ = (angular_number_+1) * (angular_number_+2) / 2 * num_contracted();

}


Shell::Shell(const bool sph) : spherical_(sph), position_(3,0.0), angular_number_(0), exponents_(1,0.0), contraction_ranges_(1,make_pair(0,1)),
        dummy_(true) {
  contractions_.push_back(vector<double>(1,1.0));
  contraction_lower_.push_back(0);
  contraction_upper_.push_back(1);
}


Shell::~Shell() {

}


std::shared_ptr<Shell> Shell::move_atom(const vector<double>& displacement) {
  std::shared_ptr<Shell> out(new Shell(*this));
  out->position_[0] += displacement[0]; 
  out->position_[1] += displacement[1]; 
  out->position_[2] += displacement[2]; 
  return out;
}


std::shared_ptr<Shell> Shell::move_atom(const double* displacement) {
  std::shared_ptr<Shell> out(new Shell(*this));
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

