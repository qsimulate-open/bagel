//
// Newint - Parallel electron correlation program.
// Filename: rysint.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
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


// Base class for the Rys-type integral evaluator

#ifndef __src_rysint_rysint_h
#define __src_rysint_rysint_h

#include <src/rysint/vrrlist.h>
#include <src/rysint/hrrlist.h>
#include <src/rysint/sortlist.h>
#include <src/scf/shell.h>
#include <vector>
#include <memory>

class RysInt {
  protected:
    VRRList vrr_; 
    HRRList hrr_; 
    SortList sort_;

    bool spherical_;

    std::vector<std::shared_ptr<Shell> > basisinfo_;
    double *data_;
    double *data2_;
    unsigned int size_final_;

    /// info for Rys quadruture
    double *roots_; 
    double *weights_;
    int rank_;

    // for screening
    int* screening_;
    int screening_size_;

  public:
    RysInt(const std::vector<std::shared_ptr<Shell> >);
    ~RysInt();

    virtual void compute() {};

    /// retrieve a batch of integrals
    const double* data() const { return data_; };
    const double* data2() const { return data2_; };
    const bool data2_exists() const { return data2_ != NULL; };
    const unsigned int data_size() const { return size_final_; };

};

#endif

