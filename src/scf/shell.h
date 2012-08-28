//
// Newint - Parallel electron correlation program.
// Filename: shell.h
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


#ifndef __newint_scf_shell_h
#define __newint_scf_shell_h

#include <cassert>
#include <array>
#include <vector>
#include <string>
#include <memory>

namespace bagel {

class Shell {

  protected:
    bool spherical_;

    std::array<double,3> position_;
    int angular_number_;
    std::vector<double> exponents_;     // length of primitive basis function
    std::vector<std::vector<double> > contractions_;  // length of contracted basis function
    std::vector<std::pair<int, int> > contraction_ranges_;

    const bool dummy_;
    std::vector<int> contraction_upper_;
    std::vector<int> contraction_lower_;

    int nbasis_;

  public:
    Shell(const bool, std::array<double,3>, int, std::vector<double>,
        std::vector<std::vector<double> >, std::vector<std::pair<int, int> >);
    // default constructor for adding null basis
    Shell(const bool sph);

    ~Shell();

    bool dummy() const { return dummy_; };
    bool spherical() const { return spherical_; };
    int num_primitive() const { return exponents_.size(); };
    int num_contracted() const { return contractions_.size(); };

    double position(const int i) const { return position_[i]; };
    const std::array<double,3> position() const { return position_; };
    int angular_number() const { return angular_number_; };
    double exponents(const int i) const { return exponents_[i]; }; 
    const std::vector<double>& exponents() const { return exponents_; };
    const double* exponents_pointer() const { return &(exponents_[0]); };
    const std::vector<double> contractions(const int i) const { return contractions_[i]; };
    const std::vector<std::vector<double> >& contractions() const { return contractions_; };
    const std::pair<int, int>& contraction_ranges(const int i) const { return contraction_ranges_[i]; }; 
    const std::vector<std::pair<int, int> >& contraction_ranges() const { return contraction_ranges_; }; 

    const std::vector<int>& contraction_upper() const { return contraction_upper_; }; 
    const std::vector<int>& contraction_lower() const { return contraction_lower_; }; 

    const std::string show() const;
    int nbasis() const { return nbasis_; }; 

    std::shared_ptr<const Shell> move_atom(const std::array<double,3>&) const;
    std::shared_ptr<const Shell> move_atom(const double*) const;

    bool operator==(const Shell& o) const;
};

}

#endif

