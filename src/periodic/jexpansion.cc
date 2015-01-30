//
// BAGEL - Parallel electron correlation program.
// Filename: jexpansion.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#include <src/periodic/jexpansion.h>
#include <src/util/constants.h>

using namespace std;
using namespace bagel;

static const double pisq__ = pi__ * pi__;

JExpansion::JExpansion(const array<shared_ptr<const Shell>, 4>& sh, const int lmax, const int ws)
 : basisinfo_(sh), lmax_(lmax), ws_(ws) {

  init();
}


bool JExpansion::is_well_separated() {

  bool out = false;
  const double r = sqrt(r12_[0] * r12_[0] + r12_[1] * r12_[1] + r12_[2] * r12_[2]);
  cout << " r = " << setprecision(9) << r << " min sepapration = " << (1.0 + ws_) * (extent0_ + extent1_) << endl;
  if (r > (1.0 + ws_) * (extent0_ + extent1_)) out = true;

  return out;
}


void JExpansion::init() {

  centre0_ = distribution_centre(array<shared_ptr<const Shell>, 2>{{basisinfo_[0], basisinfo_[1]}});
  centre1_ = distribution_centre(array<shared_ptr<const Shell>, 2>{{basisinfo_[2], basisinfo_[3]}});

  extent0_ = distribution_extent(array<shared_ptr<const Shell>, 2>{{basisinfo_[0], basisinfo_[1]}});
  extent1_ = distribution_extent(array<shared_ptr<const Shell>, 2>{{basisinfo_[2], basisinfo_[3]}});

  r12_[0] = centre0_[0] - centre1_[0];
  r12_[1] = centre0_[1] - centre1_[1];
  r12_[2] = centre0_[2] - centre1_[2];

  num_multipoles_ = (lmax_ + 1) * (lmax_ + 1);
}


array<double, 3> JExpansion::distribution_centre(array<shared_ptr<const Shell>, 2> shells) {

  const vector<double> exp0 = shells[0]->exponents();
  const vector<double> exp1 = shells[1]->exponents();
  array<double, 3> out = {{0.0, 0.0, 0.0}};
  for (auto& expi0 : exp0) {
    for (auto& expi1 : exp1) {
      const double cxp_inv = 1.0 / (expi0 + expi1);
      out[0] += (shells[0]->position(0) * expi0 + shells[1]->position(0) * expi1) * cxp_inv;
      out[1] += (shells[0]->position(1) * expi0 + shells[1]->position(1) * expi1) * cxp_inv;
      out[2] += (shells[0]->position(2) * expi0 + shells[1]->position(2) * expi1) * cxp_inv;
    }
  }
  const int denom = shells[0]->exponents().size() * shells[1]->exponents().size();
  out[0] /= denom;
  out[1] /= denom;
  out[2] /= denom;

  return out;
}


double JExpansion::distribution_extent(array<shared_ptr<const Shell>, 2> shells, const double thresh) {

  const array<double, 3> centre = distribution_centre(shells);

  const vector<double> exp0 = shells[0]->exponents();
  const vector<double> exp1 = shells[1]->exponents();
  array<double, 3> AB;
  AB[0] = shells[0]->position(0) - shells[1]->position(0);
  AB[1] = shells[0]->position(1) - shells[1]->position(1);
  AB[2] = shells[0]->position(2) - shells[1]->position(2);
  const double rsq = AB[0] * AB[0] + AB[1] * AB[1] + AB[2] * AB[2];
  const double lnthresh = log(thresh);

  double out = 0.0;
  for (auto& expi0 : exp0) {
    for (auto& expi1 : exp1) {
      const double cxp_inv = 1.0 / (expi0 + expi1);
      const double expi01 = expi0 * expi1;
      const double lda_kl = sqrt((- lnthresh - expi01 * rsq * cxp_inv + 0.75 * log(4.0 * expi01 / pisq__)) * cxp_inv);

      array<double, 3> tmp;
      tmp[0] = (shells[0]->position(0) * expi0 + shells[1]->position(0) * expi1) * cxp_inv - centre[0];
      tmp[1] = (shells[0]->position(1) * expi0 + shells[1]->position(1) * expi1) * cxp_inv - centre[1];
      tmp[2] = (shells[0]->position(2) * expi0 + shells[1]->position(2) * expi1) * cxp_inv - centre[2];

      const double extent = sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2]) + lda_kl;
      if (extent > out) out = extent;
    }
  }

  return out;
}


shared_ptr<const ZMatrix> JExpansion::compute(shared_ptr<const Matrix> density) {

  assert(is_well_separated());
  vector<shared_ptr<const ZMatrix>> multipoles0(num_multipoles_);
  vector<shared_ptr<const ZMatrix>> multipoles1(num_multipoles_);

  const int dimb1 = basisinfo_[0]->nbasis();
  const int dimb0 = basisinfo_[1]->nbasis();
  {
    MultipoleBatch mpole0(array<shared_ptr<const Shell>, 2>{{basisinfo_[0], basisinfo_[1]}}, centre0_, lmax_);
    mpole0.compute();

    for (int i = 0; i != num_multipoles_; ++i) {
      ZMatrix tmp0(dimb1, dimb0);
      tmp0.copy_block(0, 0, dimb1, dimb0, mpole0.data(i));
      multipoles0[i] = make_shared<const ZMatrix>(tmp0);
    }
  }

  const int dimb3 = basisinfo_[2]->nbasis();
  const int dimb2 = basisinfo_[3]->nbasis();
  assert(density->ndim() == dimb3 && density->mdim() == dimb2);
  {
    MultipoleBatch mpole1(array<shared_ptr<const Shell>, 2>{{basisinfo_[2], basisinfo_[3]}}, centre1_, lmax_);
    mpole1.compute();

    for (int i = 0; i != num_multipoles_; ++i) {
      ZMatrix tmp1(dimb3, dimb2);
      tmp1.copy_block(0, 0, dimb3, dimb2, mpole1.data(i));
      multipoles1[i] = make_shared<const ZMatrix>(tmp1);
    }
  }

  LocalExpansion local(r12_, multipoles1, lmax_);
  vector<shared_ptr<const ZMatrix>> lmoments = local.local_moments();

  ZMatrix out(dimb1, dimb0);
  for (int i = 0; i != num_multipoles_; ++i) {
    complex<double> contract = 0.0;
    for (int j = 0; j != dimb2; ++j)
      for (int k = 0; k != dimb3; ++k)
        contract += lmoments[i]->element(k, j) * density->element(k, j);

    out += contract * *multipoles0[i];
  }
  out.print("Multipole expansion approximation");

  return make_shared<const ZMatrix>(out);
}


void JExpansion::map_lm_index() {

  int cnt = 0;
  for (int l = 0; l <= lmax_; ++l)
    for (int m = 0; m <= 2 * l; ++m, ++cnt)
      lm_map_.push_back(make_pair(l, m-l));
}
