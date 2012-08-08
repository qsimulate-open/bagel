//
// Newint - Parallel electron correlation program.
// Filename: int2d.cc
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


#include <src/rysint/int2d.h>
#include <src/util/f77.h>
#include <cstring>

using namespace std;

Int2D::Int2D(const array<double,11>& dparam, const double* roots, const int rank, const int datasize, double* data_pointer,
   void (*vrrfunc)(double*, const double*, const double*, const double*, const double*, const double*))
 : rank_(rank), data_(data_pointer), datasize_(datasize) {

  const double P = dparam[0];
  const double Q = dparam[1];
  const double A = dparam[2];
  const double B = dparam[3];
  const double C = dparam[4];
  const double D = dparam[5];
  const double xp = dparam[6];
  const double xq = dparam[7];

  const double one_2p = dparam[8];
  const double one_2q = dparam[9];

  const double one_pq = dparam[10];
  const double xqopq = xq * one_pq;
  const double xpopq = xp * one_pq;

  const double c00i0 = P - A;
  const double c00i1 = (P - Q) * xqopq; 
  const double d00i0 = Q - C;
  const double d00i1 = (P - Q) * xpopq;
  const double b00i0 = 0.5 * one_pq;
  const double b10i0 = xqopq * one_2p; 
  const double b01i0 = xpopq * one_2q; 

  for (int i = 0; i != rank_; ++i) {
    const double tsq = roots[i];
    C00_[i] = c00i0 - c00i1 * tsq; 
    D00_[i] = d00i0 + d00i1 * tsq; 
    B00_[i] = b00i0 * tsq; 
    B10_[i] = one_2p - b10i0 * tsq; 
    B01_[i] = one_2q - b01i0 * tsq; 
  }

  vrrfunc(data_, C00_, D00_, B00_, B01_, B10_);

}


Int2D::~Int2D() {

}

