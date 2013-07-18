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

  cout << "         "  << "Applying random phase factor to integrals" << endl << endl;
  const size_t norb2 = norb_*norb_;
  const size_t norb3 = norb2*norb_;
  const size_t norb4 = norb3*norb_;

  ZMatrix phase(norb_, norb_);
  for (int i = 0; i != norb_; ++i) {
    const double ran = rand();
    const complex<double> fac(cos(ran), sin(ran));
    phase(i,i) = fac;
  }

  ZMatrix test(norb_, norb_);
  for (int i = 0; i!=norb_; ++i) test(i,i) = 1.0;

  //transform 1e integrals.
  auto mo1e = make_shared<ZMatrix>(norb_, norb_);

  copy_n(jop_->mo1e_nodelta_ptr(), norb2, mo1e->data());
  *mo1e = phase % (*mo1e * phase);

  //transforming 4 index mo2e
#if 0
  unique_ptr<complex<double>[]> mo2e(new complex<double>[norb4]);
  unique_ptr<complex<double>[]> tmp(new complex<double>[norb4]);
  unique_ptr<complex<double>[]> trans(new complex<double>[norb4]);
  unique_ptr<complex<double>[]> out(new complex<double>[norb4]);

  // 1) make -> out(n, n, n, n')
  zgemm3m_("n","n", norb3, norb_, norb_, 1.0, jop_->mo2e_ptr(), norb3, phase.data(), norb_, 0.0, tmp.get(), norb3);

  // 2) make -> out(n, n, n'(dagger), n')
  for (int l = 0; l != norb_; ++l)
    zgemm3m_("n","n", norb2, norb_, norb_, 1.0, tmp.get()+l*norb3, norb2, phase.data(), norb_, 0.0, trans.get()+l*norb3, norb2);
  copy_n(trans.get(), norb4, mo2e.get());

  // 3) transpose to make -> tmp(n', n'(dagger), n, n)
  mytranspose_complex_(jop_->mo2e_ptr(), norb2 , norb2, trans.get());

  // 4) repeat step 1 to make -> out(n', n'(dagger), n, n')
  zgemm3m_("n","c", norb3, norb_, norb_, 1.0, trans.get(), norb3, phase.data(), norb_, 0.0, tmp.get(), norb3);

  // 5) repeat step 2 to make -> mo2e(n', n', n'(dagger), n')
  for (int l = 0; l != norb_; ++l)
    zgemm3m_("n","c", norb2, norb_, norb_, 1.0, tmp.get()+l*norb3, norb2, phase.data(), norb_, 0.0, trans.get()+l*norb3, norb2);

  // 6) transpose again to return to original
  mytranspose_complex_(trans.get(), norb2, norb2, mo2e.get());
#endif

//apply transformation to hamiltonian
  jop_ =  make_shared<ZHtilde>(ref_, ncore_, ncore_+norb_, mo1e, move(mo2e));
}
