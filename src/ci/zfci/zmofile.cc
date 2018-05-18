//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: zmofile.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <src/util/prim_op.h>
#include <src/ci/zfci/zmofile.h>

using namespace std;
using namespace bagel;

ZMOFile::ZMOFile(const shared_ptr<const Geometry> geom, shared_ptr<const ZCoeff_Block> co)
 : geom_(geom), coeff_(co) {
  // density fitting is assumed
  assert(geom_->df());
}


void ZMOFile::init(const int nstart, const int nfence, const bool store_c, const bool store_g) {
  // first compute all the AO integrals in core
  nbasis_ = geom_->nbasis();
  nocc_ = (nfence - nstart)/2;
  assert((nfence - nstart) % 2 == 0);
  assert(geom_->dfs());

  // calculates the core fock matrix
  shared_ptr<const ZMatrix> hcore = compute_hcore();

  if (nstart != 0) {
    shared_ptr<const ZMatrix> den = coeff_->distmatrix()->form_density_rhf(nstart)->matrix();
    core_fock_ = compute_fock(hcore, nstart, store_c, store_g);
    const complex<double> prod = (*den * (*hcore+*core_fock_)).trace();
    if (fabs(prod.imag()) > 1.0e-12) {
      stringstream ss; ss << "imaginary part of energy is nonzero!! Perhaps Fock is not Hermite for some reasons " << setprecision(10) << prod.imag();
      cout << ss.str() << endl;
    }
    core_energy_ = 0.5*prod.real();
  } else {
    core_fock_ = hcore;
    core_energy_ = 0.0;
  }

  kramers_coeff_ = coeff_->kramers_active();

  // calculate 1-e MO integrals
  shared_ptr<Kramers<2,ZMatrix>> buf1e = compute_mo1e(kramers_coeff_);

  // calculate 2-e MO integrals
  shared_ptr<Kramers<4,ZMatrix>> buf2e = compute_mo2e(kramers_coeff_);

  // compress and set mo1e_ and mo2e_
  compress_and_set(buf1e, buf2e);
}


void ZMOFile::compress_and_set(shared_ptr<Kramers<2,ZMatrix>> buf1e, shared_ptr<Kramers<4,ZMatrix>> buf2e) {
  mo1e_ = buf1e;
  mo2e_ = make_shared<Kramers<4,ZMatrix>>();

  // Harrison requires <ij|kl> = (ik|jl)
  for (auto& mat : *buf2e) {
    shared_ptr<ZMatrix> tmp = mat.second->clone();
    sort_indices<0,2,1,3,0,1,1,1>(mat.second->data(), tmp->data(), nocc_, nocc_, nocc_, nocc_);
    bitset<4> s = mat.first.tag();
    s[2] = mat.first.tag()[1];
    s[1] = mat.first.tag()[2];
    mo2e_->emplace(s, tmp);
  }
}
