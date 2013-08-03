//
// BAGEL - Parallel electron correlation program.
// Filename: dimer_jop.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
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

#include <cmath>
#include <src/dimer/dimer_jop.h>

using namespace std;
using namespace bagel;

DimerJop::DimerJop(const shared_ptr<const Reference> ref, const int nstart, const int nfenceA, const int nfenceB,
  const shared_ptr<const Coeff> coeff)
: Jop(ref, nstart, nfenceB, coeff, string("HZ")) {

  const int norbA = nfenceA - nstart;
  const int norbB = nfenceB - nfenceA;
  nact_ = make_pair(norbA, norbB);
  const int norb = norbA + norbB;

  /************************************************************
  * Repackage mo1e integrals for monomers                     *
  ************************************************************/

  unique_ptr<double[]> mo1eA(new double[norbA*norbA]);
  unique_ptr<double[]> mo1eB(new double[norbB*norbB]);

  {
    double* modata = mo1eA.get();
    for (int i = 0; i < norbA; ++i) {
      for (int j = 0; j < norbA; ++j, ++modata) {
        *modata = ( i < j ? mo1e(i,j) : mo1e(j,i) );
      }
    }
  }

  {
    double* modata = mo1eB.get();
    for (int i = 0; i < norbB; ++i) {
      for (int j = 0; j < norbB; ++j, ++modata) {
        *modata = ( i < j ? mo1e(i+norbA,j+norbA) : mo1e(j+norbA,i+norbA) );
      }
    }
  }

  monomer_mo1es_ = make_pair( move(mo1eA), move(mo1eB) );

  /************************************************************
  * Repackage mo2e integrals for monomers                     *
  ************************************************************/

  unique_ptr<double[]> mo2eA(new double[norbA * norbA * norbA * norbA]);
  unique_ptr<double[]> mo2eB(new double[norbB * norbB * norbB * norbB]);

  {
    double* Adata = mo2eA.get();
    double* modata = mo2e_ptr();

    for (int i = 0; i < norbA; ++i) {
      for (int j = 0; j < norbA; ++j) {
        for (int k = 0; k < norbA; ++k) {
          Adata = copy(modata, modata + norbA, Adata);
          modata += norb;
        }
        modata += norb*norbB;
      }
      modata += norb*norb*norbB;
    }
  }

  {
    double* Bdata = mo2eB.get();
    double* modata = mo2e_ptr() + norb*norb*norb*norbA;

    for (int i = 0; i < norbB; ++i) {
      modata += norb*norb*norbA;
      for (int j = 0; j < norbB; ++j) {
        modata += norb*norbA;
        for (int k = 0; k < norbB; ++k) {
          modata += norbA;
          Bdata = copy(modata, modata + norbB, Bdata);
          modata += norbB;
        }
      }
    }
  }

  monomer_mo2es_ = make_pair( move(mo2eA), move(mo2eB) );

  /************************************************************
  * Package cross_mo1e integrals into a matrix                *
  ************************************************************/

  auto cross_mo1e = make_shared<Matrix>(norbA, norbB);

  {
    double* modata = cross_mo1e->data();
    for (int i = 0; i < norbB; ++i) {
      for (int j = 0; j < norbA; ++j, ++modata) {
        *modata = mo1e(j,i+norbA);
      }
    }
  }

  cross_mo1e_ = cross_mo1e;
}
