//
// Newint - Parallel electron correlation program.
// Filename: naibatch.cc
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

#include <cmath>
#include <cassert>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <src/rysint/inline.h>
#include <src/rysint/naibatch.h>
#include <src/util/f77.h>
#include <src/rysint/f77.h>
#include <src/rysint/macros.h>
#include <src/stackmem.h>
#define PI 3.1415926535897932
#define SQRTPI2 0.886226925452758013649083741671

using namespace std;

typedef std::shared_ptr<Geometry> RefGeometry;
typedef std::shared_ptr<Atom> RefAtom;
typedef std::shared_ptr<Shell> RefShell;

extern StackMem* stack;

NAIBatch::NAIBatch(const vector<RefShell> _info, const RefGeometry gm, const int L, const double A)
 :  RysInt(_info), geom_(gm), L_(L), A_(A) {

  const double integral_thresh = PRIM_SCREEN_THRESH; 

  if (_info.size() != 2) throw logic_error("NAIBatch should be called with shell pairs");
  shared_ptr<Shell> dum(new Shell(_info[0]->spherical()));
  basisinfo_.push_back(dum); basisinfo_.push_back(dum);

  // natom_ 
  natom_ = geom_->natom() * (2 * L + 1);

  // cartesian to spherical tranformation has not yet been implemented.
  const bool spherical = false;
  assert(_info[0]->spherical());

  // swap 01 indices when needed: Larger angular momentum function comes first
  set_swap_info();

  AB_[0] = basisinfo_[0]->position(0) - basisinfo_[1]->position(0);
  AB_[1] = basisinfo_[0]->position(1) - basisinfo_[1]->position(1);
  AB_[2] = basisinfo_[0]->position(2) - basisinfo_[1]->position(2);

  set_prim_contsizes();

  int asize_intermediate, asize_final, dum0, dum1;
  tie(asize_intermediate, dum0, asize_final, dum1) = set_angular_info();

  const unsigned int size_start = asize_ * primsize_; 
  size_final_ = asize_final * contsize_;
  const unsigned int size_intermediate = asize_intermediate * primsize_;
  size_alloc_ = max(max(size_start, size_intermediate), size_final_);
  data_ = stack->get(size_alloc_);
  fill(data_, data_ + size_alloc_, 0.0);
  assert(size_final_ <= size_alloc_);

  buff_ = stack->get((rank_ * 2 + 6) * primsize_ * natom_);
  double* pointer = buff_;
  p_ = pointer;     pointer += primsize_ * natom_ * 3;
  xp_ = pointer;    pointer += primsize_ * natom_;
  coeff_ = pointer; pointer += primsize_ * natom_;
  T_ = pointer;     pointer += primsize_ * natom_;

  screening_ = (int*)(stack->get(primsize_ * natom_));
  screening_size_ = 0;

  vector<double>::const_iterator expi0, expi1, expi2, expi3;
  const vector<double> exp0 = basisinfo_[0]->exponents();
  const vector<double> exp1 = basisinfo_[1]->exponents();

  int index = 0;
  vector<RefAtom> atoms = geom_->atoms();

  const double onepi2 = 1.0 / (PI * PI);
  const double sqrtpi = ::sqrt(PI);
  for (expi0 = exp0.begin(); expi0 != exp0.end(); ++expi0) { 
    for (expi1 = exp1.begin(); expi1 != exp1.end(); ++expi1) { 
      for (auto aiter = atoms.begin(); aiter != atoms.end(); ++aiter, ++index) {
        double Z = static_cast<double>((*aiter)->atom_number()); 
        const double cxp = *expi0 + *expi1;
        xp_[index] = cxp;
        const double ab = *expi0 * *expi1; 
        const double cxp_inv = 1.0 / cxp;
        p_[index * 3    ] = (basisinfo_[0]->position(0) * *expi0 + basisinfo_[1]->position(0) * *expi1) * cxp_inv;
        p_[index * 3 + 1] = (basisinfo_[0]->position(1) * *expi0 + basisinfo_[1]->position(1) * *expi1) * cxp_inv;
        p_[index * 3 + 2] = (basisinfo_[0]->position(2) * *expi0 + basisinfo_[1]->position(2) * *expi1) * cxp_inv;
        const double Eab = ::exp(-(AB_[0] * AB_[0] + AB_[1] * AB_[1] + AB_[2] * AB_[2]) * (ab * cxp_inv) );
        coeff_[index] = - 2 * Z * PI * cxp_inv * Eab; 
        const double PCx = p_[index * 3    ] - (*aiter)->position(0);
        const double PCy = p_[index * 3 + 1] - (*aiter)->position(1);
        const double PCz = p_[index * 3 + 2] - (*aiter)->position(2);
        const double T = cxp * (PCx * PCx + PCy * PCy + PCz * PCz); 
        const double sqrtt = ::sqrt(T);
        const double ss = - coeff_[index] * ::pow(4.0 * ab * onepi2, 0.75) * (T > 1.0e-15 ? sqrtpi * ::erf(sqrtt) / sqrtt * 0.5 : 1.0); 
        if (ss > integral_thresh) {
          T_[index] = T;
          screening_[screening_size_] = index;
          ++screening_size_;
        } else {
          T_[index] = -1.0;
          coeff_[index] = 0.0;
        }
      }
    }
  }

  roots_ = pointer; pointer += rank_ * primsize_ * natom_; 
  weights_ = pointer;

  root_weight(primsize_ * natom_);

}


NAIBatch::~NAIBatch() {
  stack->release(size_alloc_ + (rank_ * 2 + 6) * primsize_ * natom_+primsize_ * natom_);
}


