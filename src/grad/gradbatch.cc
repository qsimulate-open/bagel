//
// Newint - Parallel electron correlation program.
// Filename: gradbatch.cc
// Copyright (C) 2012 Toru Shiozaki
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

#define PITWOHALF 17.493418327624862
#define PIMHALF 0.564189583547756
#define SQRTPI2 0.886226925452758013649083741671

#include <cmath>
#include <cassert>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <stdexcept>
#include <src/grad/gradbatch.h>
#include <src/grad/gvrrlist.h>
#include <src/rysint/inline.h>
#include <src/rysint/macros.h>
#include <src/stackmem.h>

using namespace std;

typedef std::shared_ptr<Shell> RefShell;

// This object lives in main.cc
extern StackMem* stack;


GradBatch::GradBatch(const vector<RefShell> shells, const double max_density, const double dummy, const bool dum) :  RysInt(shells) { 
  centers_ = 4;  
  for (auto i = shells.begin(); i != shells.end(); ++i) if ((*i)->dummy()) --centers_;
  vrr_ = shared_ptr<VRRListBase>(dynamic_cast<VRRListBase*>(new GVRRList()));

  // a member variable in RysInt <- ERIBatch <- GradBatch.
  deriv_rank_ = 1;

  const double integral_thresh = 0.0;

  // determins if we want to swap shells
  set_swap_info(true);

  // stores AB and CD
  set_ab_cd();

  // set primsize_ and contsize_, as well as relevant members
  set_prim_contsizes();

  // sets angular info
  int asize_final, csize_final, asize_final_sph, csize_final_sph;
  tie(asize_final, csize_final, asize_final_sph, csize_final_sph) = set_angular_info();

  allocate_data(asize_final, csize_final, asize_final_sph, csize_final_sph);
  
  allocate_arrays(primsize_);

  compute_ssss(integral_thresh);

  root_weight(primsize_);

  set_exponents();
}


void GradBatch::set_exponents() {
  exponents_ = unique_ptr<double[]>(new double[primsize_*4]);
  double* tmp = exponents_.get();
  for (auto i0 = basisinfo_[0]->exponents().begin(); i0 != basisinfo_[0]->exponents().end(); ++i0) {
    for (auto i1 = basisinfo_[1]->exponents().begin(); i1 != basisinfo_[1]->exponents().end(); ++i1) {
      for (auto i2 = basisinfo_[2]->exponents().begin(); i2 != basisinfo_[2]->exponents().end(); ++i2) {
        for (auto i3 = basisinfo_[3]->exponents().begin(); i3 != basisinfo_[3]->exponents().end(); ++i3, tmp += 4) {
          tmp[0] = *i0;
          tmp[1] = *i1;
          tmp[2] = *i2;
          tmp[3] = *i3;
        }
      }
    }
  }
};


GradBatch::~GradBatch() {
}


