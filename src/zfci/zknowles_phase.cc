//
// BAGEL - Parallel electron correlation program.
// Filename: zknowles_phase.cc
// Copyright (C) 2013 Michael Caldwell
//
// Author: Michael Caldwell <caldwell@u.northwestern.edu>>
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

#include <cmath>
#include <random>
#include <src/math/algo.h>
#include <src/zfci/zknowles.h>

using namespace std;
using namespace bagel;

void ZKnowlesHandy::mult_phase_factor() {
  //jop_->print(3);

  cout << "         "  << "Applying random phase factor to integrals" << endl << endl;
  const size_t norb2 = norb_*norb_;
  const size_t norb3 = norb2*norb_;
  const size_t norb4 = norb3*norb_;

  ZMatrix phase(norb_, norb_);
  for (int i = 0; i != norb_; ++i) {
    const double ran = rand();
   // const complex<double> fac(cos(ran), sin(ran));
    const complex<double> fac(1.0, 0.0);
    phase(i,i) = fac;
  }
//phase(norb_-1, norb_-1) = complex<double>(0.0, 1.0);

  //transform 1e integrals.
  auto mo1e = make_shared<ZMatrix>(norb_, norb_);

  copy_n(jop_->mo1e_ptr(), norb2, mo1e->data());
#if 0
 //reverting unphased h'kl = hkl - sum k (jk|ki)
  int ij = 0;
  for (int i = 0; i != norb_; ++i) {
    for (int j = 0; j != norb_; ++j, ++ij) {
      for (int k = 0; k != norb_; ++k)
        mo1e->data(ij) += 0.5 * jop_->mo2e_kh(j, k, k, i);
    }
  }
#endif
  *mo1e = phase % *mo1e * phase;

  //transforming 4 index mo2e
  unique_ptr<complex<double>[]> mo2e(new complex<double>[norb4]);
  unique_ptr<complex<double>[]> tmp(new complex<double>[norb4]);
  unique_ptr<complex<double>[]> trans(new complex<double>[norb4]);

  // 1) make (i*, j, k, l) (left multiply by phase*)
  zgemm3m_("c","n", norb_, norb3, norb_, 1.0, phase.data(), norb_, jop_->mo2e_ptr(), norb_, 0.0, trans.get(), norb_);

  // 2) make (i*, j, k, l') (right multiply by phase)
  zgemm3m_("n","n", norb3, norb_, norb_, 1.0, trans.get(), norb3, phase.data(), norb_, 0.0, tmp.get(), norb3);

  // 3) transpose to make (k, l', i*, j) now (kl|ij)
  mytranspose_complex_(tmp.get(), norb2, norb2, trans.get());

  // 4) make (k*, l', i*, j) (left multiply by phase*)
  zgemm3m_("c","n", norb_, norb3, norb_, 1.0, phase.data(), norb_, trans.get(), norb_, 0.0, tmp.get(), norb_);

  // 5) make (k', l', i*, j') (right multiply by phase)
  zgemm3m_("n","n", norb3, norb_, norb_, 1.0, tmp.get(), norb3, phase.data(), norb_, 0.0, trans.get(), norb3);

  // 6) transpose to return indices to (ij|kl)
  mytranspose_complex_(trans.get(), norb2, norb2, tmp.get());


  // applying phased h'kl = hkl - sum k (jk|ki)
#if 0
  ij = 0;
  for (int i = 0; i != norb_; ++i) {
    for (int j = 0; j != norb_; ++j, ++ij) {
      for (int k = 0; k != norb_; ++k)
        mo1e->data(ij) -= 0.5 * tmp[j+k*norb_+k*norb2+i*norb3];
    }
  }
#endif
  //pass complex conjugate of 1 and 2 body integrals to account for replacement of 8fold symmetry with 4fold symmetry in complex KH
  mytranspose_complex_conjg_(tmp.get(), norb2, norb2, mo2e.get());
//apply transformation to hamiltonian
#if 1
  jop_->set_moints(mo1e->transpose_conjg(), mo2e);
#else
  cout << "one body" << endl;
  complex<double>* a = jop_->mo1e_ptr();
  auto pt = mo1e->transpose_conjg();
  complex<double>* b = pt->data();
  cout << setprecision(10) << zdotc_(norb2, a, 1, a, 1) - 2*zdotc_(norb2, a, 1, b, 1).real() + zdotc_(norb2, b, 1, b, 1) << endl;
  cout << "two body" << endl;
  a = jop_->mo2e_ptr();
  b = mo2e.get();
  cout << setprecision(10) << zdotc_(norb4, a, 1, a, 1) - 2*zdotc_(norb4, a, 1, b, 1).real() + zdotc_(norb4, b, 1, b, 1) << endl;
  assert(false);

  //jop_->print(3);
#endif
}
