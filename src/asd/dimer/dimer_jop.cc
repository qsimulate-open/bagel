//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dimer_jop.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <cmath>
#include <src/asd/dimer/dimer_jop.h>

using namespace std;
using namespace bagel;

DimerJop::DimerJop(const shared_ptr<const Reference> ref, const int nstart, const int nfenceA, const int nfenceB,
  const shared_ptr<const Coeff> coeff)
: Jop(ref, nstart, nfenceB, coeff, false, string("HZ")) {

  const int norbA = nfenceA - nstart;
  const int norbB = nfenceB - nfenceA;
  common_init(norbA, norbB);
}

DimerJop::DimerJop(const int norbA, const int norbB, shared_ptr<CSymMatrix> mo1e, shared_ptr<Matrix> mo2e) : Jop(mo1e, mo2e) {
  common_init(norbA, norbB);
}

void DimerJop::common_init(const int norbA, const int norbB) {
  nact_ = {norbA, norbB};
  const int norb = norbA + norbB;

  /************************************************************
  * Repackage mo1e integrals for monomers                     *
  ************************************************************/

  auto mo1eA = make_shared<Matrix>(norbA, norbA);
  auto mo1eB = make_shared<Matrix>(norbB, norbB);

  {
    double* modata = mo1eA->data();
    for (int i = 0; i < norbA; ++i) {
      for (int j = 0; j < norbA; ++j, ++modata) {
        *modata = ( i < j ? mo1e(i,j) : mo1e(j,i) );
      }
    }
  }

  {
    double* modata = mo1eB->data();
    for (int i = 0; i < norbB; ++i) {
      for (int j = 0; j < norbB; ++j, ++modata) {
        *modata = ( i < j ? mo1e(i+norbA,j+norbA) : mo1e(j+norbA,i+norbA) );
      }
    }
  }

  /************************************************************
  * Repackage mo2e integrals for monomers                     *
  ************************************************************/

  auto mo2eA = make_shared<Matrix>(norbA*norbA, norbA*norbA);
  auto mo2eB = make_shared<Matrix>(norbB*norbB, norbB*norbB);

  {
    double* Adata = mo2eA->data();
    double* modata = mo2e_->data();

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
    double* Bdata = mo2eB->data();
    double* modata = mo2e_->data() + norb*norb*norb*norbA;

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

  auto c1eA = make_shared<CSymMatrix>(mo1eA);
  auto c1eB = make_shared<CSymMatrix>(mo1eB);

  jops_ = {make_shared<Jop>(c1eA, mo2eA), make_shared<Jop>(c1eB, mo2eB)};

  /************************************************************
  * Package cross_mo1e integrals into a matrix                *
  ************************************************************/

  auto cross_mo1e = make_shared<Matrix>(norbA, norbB, true);

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

