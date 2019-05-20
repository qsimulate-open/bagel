//
// BAGEL - Parallel electron correlation program.
// Filename: shellpair.cc
// Copyright (C) 2016 Toru Shiozaki
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

#include <src/molecule/shellpair.h>
#include <src/integral/os/multipolebatch.h>
#include <src/integral/rys/eribatch.h>
#include <src/integral/libint/libint.h>

using namespace std;
using namespace bagel;

const static double pisq__ = pi__ * pi__;

ShellPair::ShellPair(const array<shared_ptr<const Shell>, 2>& sh, const array<int, 2>& ofs, const pair<int, int>& ind, const string ext, const double thr)
 : shells_(sh), offset_(ofs), shell_ind_(ind), extent_type_(ext), thresh_(thr) {
  init();
}


vector<shared_ptr<const ZMatrix>> ShellPair::multipoles(const int lmax, const array<double, 3>& Q) const {

  const int nmult =  (lmax + 1) * (lmax + 1);
  vector<shared_ptr<const ZMatrix>> mult(nmult);

  const int isize = shells_[0]->nbasis();
  const int jsize = shells_[1]->nbasis();

  {
    MultipoleBatch mpole(shells_, Q, lmax);
    mpole.compute();
    ZMatrix tmp(jsize, isize);
    for (int im = 0; im != nmult; ++im) {
      tmp.copy_block(0, 0, jsize, isize, mpole.data(im));
      mult[im] = make_shared<const ZMatrix>(tmp);
    }
  }

  return mult;
}


bool ShellPair::is_neighbour(shared_ptr<const ShellPair> sp, const double ws) const {

  array<double, 3> v;
  v[0] = centre_[0] - sp->centre(0);
  v[1] = centre_[1] - sp->centre(1);
  v[2] = centre_[2] - sp->centre(2);
  const double r = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  return (r <= (1.0 + ws) * (extent_ + sp->extent()));
}


void ShellPair::init() {

  shared_ptr<const Shell> b0 = shells_[0];
  shared_ptr<const Shell> b1 = shells_[1];
  nbasis0_ = b0->nbasis();
  nbasis1_ = b1->nbasis();
  assert(b0->angular_number() < 7 && b1->angular_number() < 7);

  array<double, 3> AB;
  AB[0] = b0->position(0) - b1->position(0);
  AB[1] = b0->position(1) - b1->position(1);
  AB[2] = b0->position(2) - b1->position(2);
  const double rsq = AB[0] * AB[0] + AB[1] * AB[1] + AB[2] * AB[2];

  const double tol = 20.0 / (log10(exp(1))*rsq);
  // centre
  centre_ = {{0.0, 0.0, 0.0}};
  for (auto& expi0 : b0->exponents()) {
    for (auto& expi1 : b1->exponents()) {
      const double cxp_inv = 1.0 / (expi0 + expi1);
      centre_[0] += (b0->position(0) * expi0 + b1->position(0) * expi1) * cxp_inv;
      centre_[1] += (b0->position(1) * expi0 + b1->position(1) * expi1) * cxp_inv;
      centre_[2] += (b0->position(2) * expi0 + b1->position(2) * expi1) * cxp_inv;
    }
  }
  const int nexp = b0->exponents().size()*b1->exponents().size();
  centre_[0] /= nexp;
  centre_[1] /= nexp;
  centre_[2] /= nexp;

  // extent
//  array<double, 7> scale = {{1.0, 1.1781, 1.3333, 1.4726, 1.7181, 1.8286}};
  const double lnthresh = log(thresh_);
  if (extent_type_ == "sierka" || extent_type_ == "sierka-alhrichs" || extent_type_ == "alhrichs") {
    double extentA = 0.0;
    double extentB = 0.0;
    for (auto& expi0 : b0->exponents()) {
      for (auto& expi1 : b1->exponents()) {
        const double cxp_inv = 1.0 / (expi0 + expi1);
        const double expi01 = expi0 * expi1;
        if (expi01*cxp_inv > tol) continue;
        const double s01 = pow(4.0 * expi01 * cxp_inv * cxp_inv, 0.75) * exp(-expi01 * cxp_inv * rsq);
        const double r01sq = (-lnthresh + log(s01) + 0.5 * log(expi0 + expi1)) * cxp_inv;

        array<double, 3> P;
        P[0] = (b0->position(0) * expi0 + b1->position(0) * expi1) * cxp_inv;
        P[1] = (b0->position(1) * expi0 + b1->position(1) * expi1) * cxp_inv;
        P[2] = (b0->position(2) * expi0 + b1->position(2) * expi1) * cxp_inv;

        array<double, 3> r0;
        r0[0] = P[0] - b0->position(0);
        r0[1] = P[1] - b0->position(1);
        r0[2] = P[2] - b0->position(2);
        const double d0 = sqrt(r0[0]*r0[0] + r0[1]*r0[1] + r0[2]*r0[2]);

        array<double, 3> r1;
        r1[0] = P[0] - b1->position(0);
        r1[1] = P[1] - b1->position(1);
        r1[2] = P[2] - b1->position(2);
        const double d1 = sqrt(r1[0]*r1[0] + r1[1]*r1[1] + r1[2]*r1[2]);

        if (r01sq > 0) {
          const double r01 = sqrt(r01sq);
          extentA = max(r01 + d0, extentA);
          extentB = max(r01 + d1, extentB);
        } else {
          extentA = max(1.0e-2, extentA);
          extentB = max(1.0e-2, extentB);
        }
      }
    }
    extent_ = min(extentA, extentB);
  } else if (extent_type_ == "perez-jorda" || extent_type_ == "pj-yang" || extent_type_ == "yang") {
    extent_ = 0.0;
    const double lntau = log(1.0e-6);
    //const double lntau = log(sqrt(thresh_));
    for (auto& expi0 : b0->exponents()) {
      for (auto& expi1 : b1->exponents()) {
        const double cxp_inv = 1.0 / (expi0 + expi1);
        const double expi01 = expi0 * expi1;
        if (expi01*cxp_inv > tol) continue;
        const double tmp = (-lntau- expi01*cxp_inv*rsq + 0.75*log(4.0*expi01/pisq__))*cxp_inv;
        const double r01sq = max(0.0, tmp);

        array<double, 3> ld;
        ld[0] = (b0->position(0) * expi0 + b1->position(0) * expi1) * cxp_inv - centre_[0];
        ld[1] = (b0->position(1) * expi0 + b1->position(1) * expi1) * cxp_inv - centre_[1];
        ld[2] = (b0->position(2) * expi0 + b1->position(2) * expi1) * cxp_inv - centre_[2];
        const double d = sqrt(ld[0]*ld[0] + ld[1]*ld[1] + ld[2]*ld[2]);
        extent_ = max(extent_, sqrt(r01sq)+d);
      }
    }
  }
//  extent_ *= scale[b0->angular_number()] * scale[b1->angular_number()];

  // schwarz
  array<shared_ptr<const Shell>,4> input = {{b1, b0, b1, b0}};
#ifdef LIBINT_INTERFACE
  Libint eribatch(input);
#else
  ERIBatch eribatch(input, 0.0);
#endif
  eribatch.compute();
  const double* eridata = eribatch.data();
  const int datasize = eribatch.data_size();
  schwarz_ = 0.0;
  for (int xi = 0; xi != datasize; ++xi, ++eridata) {
    const double absed = sqrt(abs(*eridata));
    if (absed > schwarz_) schwarz_ = absed;
  }
}
