//
// Newint - Parallel electron correlation program.
// Filename: atom.h
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


#ifndef __scf_atom_h
#define __scf_atom_h

#include <vector>
#include <string>
#include <src/scf/shell.h>
#include <memory>

class Atom {
  protected:

    bool spherical_;

    std::string name_;
    std::vector<double> position_;
    std::vector<std::shared_ptr<Shell> > shells_; 
    int atom_number_;
    int nbasis_;
    int lmax_;

    // This function sets shell_ and lmax_
    // in : a vector of an angular label, exponents, and coefficients. 
    void construct_shells(std::vector<std::tuple<std::string, std::vector<double>, std::vector<std::vector<double> > > > in);
    void construct_shells(std::vector<std::tuple<std::string, std::vector<double>, std::vector<double> > > in);

    void common_init();
     
  public:
    Atom(const bool spherical, const std::string name, const std::vector<double>& position, const std::string basisfile);
    Atom(const std::string name, const std::vector<std::shared_ptr<Shell> > shell);
    Atom(const Atom&, const std::vector<double>&);
    Atom(const Atom&, const double*);
    ~Atom() {};

    const std::string name() const { return name_; };
    const int atom_number() const { return atom_number_;};
    const std::vector<double> position() const { return position_; };
    const double position(const unsigned int i) const { return position_[i]; };
    const std::vector<std::shared_ptr<Shell> >& shells() const { return shells_; };
    const int nshell() const { return shells_.size(); };

    const int nbasis() const { return nbasis_; };
    const int lmax() const { return lmax_; };
    const bool spherical() const { return spherical_; };

    void print_basis() const;
    void print() const;
};

#endif

