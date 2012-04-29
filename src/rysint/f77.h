//
// Newint - Parallel electron correlation program.
// Filename: f77.h
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


#ifndef __src_rysint_f77_h
#define __src_rysint_f77_h

extern "C" {
 void rysroot_(const double*, double*, double*, const int*, const int*);
 void eriroot1_(const double*, double*, double*, const int*);
 void eriroot2_(const double*, double*, double*, const int*);
 void eriroot3_(const double*, double*, double*, const int*);
 void eriroot4_(const double*, double*, double*, const int*);
 void eriroot5_(const double*, double*, double*, const int*);
 void eriroot6_(const double*, double*, double*, const int*);
 void eriroot7_(const double*, double*, double*, const int*);
 void eriroot8_(const double*, double*, double*, const int*);
 void eriroot9_(const double*, double*, double*, const int*);
 void eriroot10_(const double*, double*, double*, const int*);
 void eriroot11_(const double*, double*, double*, const int*);
 void eriroot12_(const double*, double*, double*, const int*);
 void eriroot13_(const double*, double*, double*, const int*);
};

#endif
