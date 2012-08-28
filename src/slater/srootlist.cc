//
// BAGEL - Parallel electron correlation program.
// Filename: srootlist.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include <src/slater/srootlist.h>

extern "C" {
 // slater integrals
 void root1_(const double*, const double*, double*, double*, const int*);
 void root2_(const double*, const double*, double*, double*, const int*);
 void root3_(const double*, const double*, double*, double*, const int*);
 void root4_(const double*, const double*, double*, double*, const int*);
 void root5_(const double*, const double*, double*, double*, const int*);
 void root6_(const double*, const double*, double*, double*, const int*);
 void root7_(const double*, const double*, double*, double*, const int*);
 void root8_(const double*, const double*, double*, double*, const int*);
 void root9_(const double*, const double*, double*, double*, const int*);
 void root10_(const double*, const double*, double*, double*, const int*);
 void root11_(const double*, const double*, double*, double*, const int*);
 void root12_(const double*, const double*, double*, double*, const int*);
 void root13_(const double*, const double*, double*, double*, const int*);
}

SRootList::SRootList() {

  srfunc[1] = &root1_;
  srfunc[2] = &root2_;
  srfunc[3] = &root3_;
  srfunc[4] = &root4_;
  srfunc[5] = &root5_;
  srfunc[6] = &root6_;
  srfunc[7] = &root7_;
  srfunc[8] = &root8_;
  srfunc[9] = &root9_;
  srfunc[10] = &root10_;
  srfunc[11] = &root11_;
  srfunc[12] = &root12_;
  srfunc[13] = &root13_;
}


SRootList::~SRootList() {

}
