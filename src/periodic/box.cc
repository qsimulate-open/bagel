//
// BAGEL - Parallel electron correlation program.
// Filename: box.cc
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


#include <src/util/f77.h>
#include <src/periodic/box.h>
#include <src/periodic/multipolebatch.h>
#include <src/integral/os/overlapbatch.h>
#include <src/periodic/localexpansion.h>
#include <src/integral/rys/eribatch.h>

using namespace bagel;
using namespace std;

static const double pisq__ = pi__ * pi__;
const static Legendre plm;

void Box::init() {

  centre_ = {{0, 0, 0}};
  ndim_ = 0;
  nbasis0_ = 0;   nbasis1_ = 0;
  for (auto& i : sp_) {
    centre_[0] += i->centre(0);
    centre_[1] += i->centre(1);
    centre_[2] += i->centre(2);
    nbasis0_ += i->nbasis0();
    nbasis1_ += i->nbasis1();
    ndim_ += i->nbasis0() * i->nbasis1();
  }
  centre_[0] /= nsp();
  centre_[1] /= nsp();
  centre_[2] /= nsp();

  extent_ = 0;
  for (auto& i : sp_) {
    if (i->schwarz() < 1e-20) continue;
    double tmp = 0;
    for (int j = 0; j != 3; ++j)
      tmp += pow(i->centre(j)-centre_[j], 2.0);
    const double ei = sqrt(tmp) + i->extent();
    assert(ei > 0);
    if (extent_ < ei) extent_ = ei;
  }

#if 0
  ndim_ = 0;
  if (nchild() == 0) {
    for (auto& i : sp_)
      ndim_ += i->nbasis0() * i->nbasis1();
  } else {
    for (int n = 0; n != nchild(); ++n) {
      shared_ptr<const Box> c = child(n);
      ndim_ += c->ndim();
    }
  }
#endif

  nmult_ = (lmax_ + 1) * (lmax_ + 1);
  multipole_.resize(nmult_);
  offset0_.resize(nsp());
  fill_n(offset0_.begin(), nsp(), 0);
  offset1_.resize(nsp());
  fill_n(offset1_.begin(), nsp(), 0);
}


void Box::insert_sp(vector<shared_ptr<const ShellPair>> sp) {

  const int nsp = sp_.size();
  sp_.resize(nsp + sp.size());
  for (int i = 0; i != sp.size(); ++i)
    sp_[i + nsp] = sp[i];

}


void Box::insert_child(shared_ptr<const Box> child) {

  const int nchild = child_.size();
  if (child) {
    child_.resize(nchild + 1);
    child_[nchild] = child;
  }
}


void Box::get_neigh(vector<shared_ptr<Box>> box, const int ws) {

  neigh_.resize(box.size());
  int nn = 0;
  for (auto& b : box) {
    if (b->rank() == rank_ && is_neigh(b, ws)) {
      neigh_[nn] = b;
      ++nn;
    }
  }
  neigh_.resize(nn);
}


bool Box::is_neigh(shared_ptr<const Box> box, const int ws) const {

  double rr = 0;
  for (int i = 0; i != 3; ++i)
    rr += pow(centre_[i] - box->centre(i), 2);

  const bool out = (sqrt(rr) < (1+ws)*(extent_ + box->extent()));
  return out;
}


void Box::get_inter(vector<shared_ptr<Box>> box, const int ws) {

  inter_.resize(box.size());
  int ni = 0;
  for (auto& b : box) {
    if (b->rank() == rank_ && !is_neigh(b, ws) && (!parent_ || parent_->is_neigh(b->parent(), ws))) {
      inter_[ni] = b;
      ++ni;
    }
  }
  inter_.resize(ni);
}


void Box::compute_multipoles() {

  if (nchild() == 0) { // leaf = shift sp's multipoles
    int isp = 1;
    int offset = 0;
    for (auto& v : sp_) {
      vector<shared_ptr<const ZMatrix>> vmult = v->multipoles(lmax_);
      array<double, 3> r12;
      r12[0] = centre_[0] - v->centre(0);
      r12[1] = centre_[1] - v->centre(1);
      r12[2] = centre_[2] - v->centre(2);
      LocalExpansion shift(r12, vmult, lmax_);
      vector<shared_ptr<const ZMatrix>> smoment = shift.compute_shifted_multipoles();

      for (int i = 0; i != nmult_; ++i) {
        multipole_[i].resize(ndim_);

        for (int j = 0; j != v->nbasis0(); ++j)
          for (int k = 0; k != v->nbasis1(); ++k) {
            const int pos = offset + j * v->nbasis1() + k;
            multipole_[i][pos] = smoment[i]->element(j, k);
          }
      }

      offset0_[isp] = offset0_[isp-1] + v->nbasis0();
      offset1_[isp] = offset1_[isp-1] + v->nbasis1();
      offset += v->nbasis0() * v->nbasis1();
      ++isp;
    }
    assert(multipole_[0].size() == ndim_);
  } else { // shift children's multipoles
    int ich = 1;
    for (int n = 0; n != nchild(); ++n) {
      shared_ptr<const Box> c = child(n);
      array<double, 3> r12;
      r12[0] = centre_[0] - c->centre(0);
      r12[1] = centre_[1] - c->centre(1);
      r12[2] = centre_[2] - c->centre(2);
      vector<vector<complex<double>>> smoment = shift_multipoles(c->multipole(), r12);

      for (int i = 0; i != nmult_; ++i)
        multipole_[i].insert(multipole_[i].end(), smoment[i].begin(), smoment[i].end());

      offset0_[ich] = offset0_[ich-1] + c->nbasis0();
      offset1_[ich] = offset1_[ich-1] + c->nbasis1();
      ++ich;
    }
    assert(multipole_[0].size() == ndim_);
  }
}


#if 0
void Box::compute_local_expansions(shared_ptr<const Matrix> density) {

  for (auto& it : inter_) {
    array<double, 3> r12;
    r12[0] = centre_[0] - it->centre(0);
    r12[1] = centre_[1] - it->centre(1);
    r12[2] = centre_[2] - it->centre(2);
    vector<double> taylor = compute_local_expansions(r12, subden);
  }
}
#endif

shared_ptr<const ZMatrix> Box::compute_node_energy(shared_ptr<const Matrix> density, vector<double> max_den, const double schwarz_thresh) const {

  assert(nchild() == 0);
  auto out = make_shared<ZMatrix>(density->ndim(), density->ndim());
  out->zero();


  // FF: add local expansions

  int ninteg = 0;
  // NF: 4c integrals
  const int shift = sizeof(int) * 4;
  const double* density_data = density->data();
  const int nsh = sqrt(max_den.size());
  assert (nsh*nsh == max_den.size());

  for (auto& v01 : sp_) {
    shared_ptr<const Shell> b0 = v01->shell(1);
    const int i0 = v01->shell_ind(1);
    const int b0offset = v01->offset(1);
    const int b0size = b0->nbasis();

    shared_ptr<const Shell> b1 = v01->shell(0);
    const int i1 = v01->shell_ind(0);
    if (i1 < i0) continue;
    const int b1offset = v01->offset(0);
    const int b1size = b1->nbasis();

    const int i01 = i0 * nsh + i1;
    const double density_01 = max_den[i01] * 4.0;

    for (auto& neigh : neigh_) {
      for (auto& v23 : neigh->sp()) {
        shared_ptr<const Shell> b2 = v23->shell(1);
        const int i2 = v23->shell_ind(1);
        if (i2 < i0) continue;
        const int b2offset = v23->offset(1);
        const int b2size = b2->nbasis();

        shared_ptr<const Shell> b3 = v23->shell(0);
        const int i3 = v23->shell_ind(0);
        if (i3 < i2) continue;
        const int b3offset = v23->offset(0);
        const int b3size = b3->nbasis();

        const int i23 = i2 * nsh + i3;
        if (i23 < i01) continue;

        const double density_23 = max_den[i23] * 4.0;
        const double density_02 = max_den[i0 * nsh + i2];
        const double density_03 = max_den[i0 * nsh + i3];
        const double density_12 = max_den[i1 * nsh + i2];
        const double density_13 = max_den[i1 * nsh + i3];

        const double mulfactor = max(max(max(density_01, density_02),
                                         max(density_03, density_12)),
                                         max(density_13, density_23));
        const double integral_bound = mulfactor * v01->schwarz() * v23->schwarz();
        const bool skip_schwarz = integral_bound < schwarz_thresh;
        if (skip_schwarz) continue;

        ++ninteg;
        array<shared_ptr<const Shell>,4> input = {{b3, b2, b1, b0}};
        ERIBatch eribatch(input, mulfactor);
        eribatch.compute();
        const double* eridata = eribatch.data();

        for (int j0 = b0offset; j0 != b0offset + b0size; ++j0) {
          const int j0n = j0 * density->ndim();
          for (int j1 = b1offset; j1 != b1offset + b1size; ++j1) {
            if (j1 < j0) {
              eridata += b2size * b3size;
              continue;
            }
            const int j1n = j1 * density->ndim();
            const unsigned int nj01 = (j0 << shift) + j1;
            const double scale01 = (j0 == j1) ? 0.5 : 1.0;
            for (int j2 = b2offset; j2 != b2offset + b2size; ++j2) {
              const int j2n = j2 * density->ndim();
              for (int j3 = b3offset; j3 != b3offset + b3size; ++j3, ++eridata) {
                if (j3 < j2) continue;
                const int j3n = j3 * density->ndim();
                const unsigned int nj23 = (j2 << shift) + j3;
                if (nj23 < nj01 && i01 == i23) continue;
                const double scale23 = (j2 == j3) ? 0.5 : 1.0;
                const double scale = (nj01 == nj23) ? 0.25 : 0.5;

                const double eri = *eridata;
                const double intval = eri * scale * scale01 * scale23;
                const double intval4 = 4.0 * intval;
                out->element(j1, j0) += density_data[j2n + j3] * intval4;
                out->element(j3, j2) += density_data[j0n + j1] * intval4;
                out->element(max(j2, j0), min(j2,j0)) -= density_data[j1n + j3] * intval;
                out->element(j3, j0) -= density_data[j1n + j2] * intval;
                out->element(max(j1,j2), min(j1,j2)) -= density_data[j0n + j3] * intval;
                out->element(max(j1,j3), min(j1,j3)) -= density_data[j0n + j2] * intval;
              }
            }
          }
        }

      }
    }
  }
  //cout << "ninteg = " << ninteg << endl;

  return out;
}


void Box::print_box() const {
  cout << "Box " << boxid_ << " Rank = " << rank_ << " *** nchild = " << nchild() << " *** nsp = " << nsp()
       << " *** nneigh = " << nneigh() << " *** ninter = " << ninter() << " *** extent = " << extent_
       << " *** centre = " << setprecision(3) << centre_[0] << "  " << centre_[1] << "  " << centre_[2] << endl;

//  cout << " *** Shell pairs at ***" << endl;
//  for (int i = 0; i != nsp(); ++i)
//    cout << setprecision(5) << sp(i)->centre(0) << "  " << sp(i)->centre(1) << "  " << sp(i)->centre(2) << endl;
}


vector<complex<double>> Box::get_Mlm(array<double, 3> r12, vector<double> den) const {

  const double r = sqrt(r12[0]*r12[0] + r12[1]*r12[1] + r12[2]*r12[2]);
  const double ctheta = (r > numerical_zero__) ? r12[2]/r : 0.0;
  const double phi = atan2(r12[1], r12[0]);

  vector<complex<double>> out(nmult_);

  int i1 = 0;
  for (int l = 0; l <= lmax_; ++l) {
    for (int m = 0; m <= 2 * l; ++m, ++i1) {

      vector<complex<double>> local(ndim_);
      int i2 = 0;
      for (int j = 0; j <= lmax_; ++j) {
        for (int k = 0; k <= 2 * j; ++k, ++i2) {

          const int a = l + j;
          const int b = m - l + k - j;

          double prefactor = plm.compute(a, abs(b), ctheta) / pow(r, a + 1);
          double ft = 1.0;
          for (int i = 1; i <= a - abs(b); ++i) {
            prefactor *= ft;
            ++ft;
          }
          const double real = (b >= 0) ? (prefactor * cos(abs(b) * phi)) : (-1.0 * prefactor * cos(abs(b) * phi));
          const double imag = prefactor * sin(abs(b) * phi);
          const complex<double> coeff(real, imag);

          transform(multipole_[i2].begin(), multipole_[i2].end(), local.begin(), bind1st(multiplies<complex<double>>(), coeff));
        }
      }
      assert(local.size() == den.size());
      //out[i1] = inner_product(local.begin(), local.end(), den.begin(), 0);
      out[i1] = 0;
      for (int i = 0; i != ndim_; ++i)
        out[i1] += local[i]*den[i];
    }
  }

  return out;
}


// given O(a) and R(b-a) get O(b)
vector<vector<complex<double>>> Box::shift_multipoles(vector<vector<complex<double>>> oa, array<double, 3> rab) const {

  const double r = sqrt(rab[0]*rab[0] + rab[1]*rab[1] + rab[2]*rab[2]);
  const double ctheta = (r > numerical_zero__) ? rab[2]/r : 0.0;
  const double phi = atan2(rab[1], rab[0]);

  vector<vector<complex<double>>> ob(nmult_);

  int i1 = 0;
  for (int l = 0; l <= lmax_; ++l) {
    for (int m = 0; m <= 2 * l; ++m, ++i1) {

      vector<complex<double>> ob_lm(oa[0].size());
      int i2 = 0;
      for (int j = 0; j <= lmax_; ++j) {
        for (int k = 0; k <= 2 * j; ++k, ++i2) {

          const int a = l - j;
          const int b = m - l - k + j;
          if (abs(b) <= a && a >= 0) {
            double prefactor = pow(r, a) * plm.compute(a, abs(b), ctheta);
            double ft = 1.0;
            for (int i = 1; i <= a + abs(b); ++i) {
              prefactor /= ft;
              ++ft;
            }
            const double real = (b >= 0) ? (prefactor * cos(abs(b) * phi)) : (-1.0 * prefactor * cos(abs(b) * phi));
            const double imag = prefactor * sin(abs(b) * phi);
            const complex<double> coeff(real, imag);

            if (abs(coeff) > numerical_zero__)
              for (int ib = 0; ib != ob_lm.size(); ++ib)
                ob_lm[ib] = coeff * oa[i2][ib];
          }
        }
      }
      ob[i1] = ob_lm;
    }
  }

  return ob;
}
