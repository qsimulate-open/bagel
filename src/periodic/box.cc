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
#include <src/integral/os/multipolebatch.h>
#include <src/integral/os/overlapbatch.h>
#include <src/periodic/localexpansion.h>
#include <src/integral/rys/eribatch.h>
#include <src/util/taskqueue.h>
#include <mutex>
#include <src/util/timer.h>

using namespace bagel;
using namespace std;

static const double pisq__ = pi__ * pi__;
const static Legendre plm;

void Box::init() {

  centre_ = {{0, 0, 0}};
  nbasis0_ = 0;   nbasis1_ = 0;
  for (auto& i : sp_) {
    centre_[0] += i->centre(0);
    centre_[1] += i->centre(1);
    centre_[2] += i->centre(2);
    nbasis0_ += i->nbasis0();
    nbasis1_ += i->nbasis1();
  }
  centre_[0] /= nsp();
  centre_[1] /= nsp();
  centre_[2] /= nsp();

  extent_ = 0;
  for (auto& i : sp_) {
    if (i->schwarz() < 1e-15) continue;
    double tmp = 0;
    for (int j = 0; j != 3; ++j)
      tmp += pow(i->centre(j)-centre_[j], 2.0);
    const double ei = sqrt(tmp) + i->extent();
    assert(ei > 0);
    if (extent_ < ei) extent_ = ei;
  }

  nmult_ = (lmax_ + 1) * (lmax_ + 1);
  multipole_.resize(nmult_); fill_n(multipole_.begin(), nmult_, 0.0);
  localJ_.resize(nmult_);    fill_n(localJ_.begin(), nmult_, 0.0);
}


void Box::insert_sp(const vector<shared_ptr<const ShellPair>>& sp) {

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


void Box::insert_parent(shared_ptr<const Box> parent) {

  assert(!parent_);
  parent_ = parent;
}


void Box::get_neigh(const vector<shared_ptr<Box>>& box, const int ws) {

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

#if 1
  double rr = 0;
  for (int i = 0; i != 3; ++i)
    rr += pow(centre_[i] - box->centre(i), 2);

  const bool out = (sqrt(rr) <= (1+ws)*(extent_ + box->extent()));
#endif

#if 0
  bool out = false;
  for (auto& v : sp_)
    for (auto& v1 : box->sp())
      if (v1->is_neighbour(v, ws)) {
        out = true;
        return out;
      }
#endif

  return out;
}


void Box::get_inter(const vector<shared_ptr<Box>>& box, const int ws) {

  inter_.resize(box.size());
  int ni = 0;
  for (auto& b : box) {
    if ((b->rank() == rank_) && (!is_neigh(b, ws)) && (!parent_ || parent_->is_neigh(b->parent(), ws))) {
      inter_[ni] = b;
      ++ni;
    }
  }
  inter_.resize(ni);
}


void Box::sort_sp() {

  if (nchild() == 0) return;

  assert(nchild() > 0);
  vector<shared_ptr<const ShellPair>> newsp;
  for (int n = 0; n != nchild(); ++n) {
    shared_ptr<const Box> c = child(n);
    const vector<shared_ptr<const ShellPair>> csp = c->sp();
    newsp.insert(newsp.end(), csp.begin(), csp.end());
  }
  assert(newsp.size() == nsp());
  sp_ = newsp;
}


void Box::compute_M2M(shared_ptr<const Matrix> density) {

  fill(multipole_.begin(), multipole_.end(), 0.0);
  fill(localJ_.begin(), localJ_.end(), 0.0);

  if (nchild() == 0) { // leaf
    TaskQueue<function<void(void)>> tasks(nsp());
    mutex jmutex;
    for (auto& v : sp_) {
      if (v->schwarz() < 1e-15) continue;

      tasks.emplace_back(
        [this, &v, &density, &jmutex]() {
          vector<complex<double>> olm(nmult_);
          MultipoleBatch mpole(v->shells(), centre_, lmax_);
          mpole.compute();
          const int dimb0 = v->shell(0)->nbasis();
          const int dimb1 = v->shell(1)->nbasis();
          vector<const complex<double>*> dat(nmult_);
          for (int i = 0; i != nmult_; ++i)
            dat[i] = mpole.data() + mpole.size_block()*i;

          lock_guard<mutex> lock(jmutex);
          for (int i = v->offset(1); i != dimb1 + v->offset(1); ++i)
            for (int j = v->offset(0); j != dimb0 + v->offset(0); ++j)
              for (int k = 0; k != mpole.num_blocks(); ++k)
                olm[k] += conj(*dat[k]++ * density->element(i,j));

          transform(multipole_.begin(), multipole_.end(), olm.begin(), multipole_.begin(), std::plus<complex<double>>());
        }
      );
    }
    tasks.compute();
  } else { // shift children's multipoles
    for (int n = 0; n != nchild(); ++n) {
      shared_ptr<const Box> c = child(n);
      array<double, 3> r12;
      r12[0] = c->centre(0) - centre_[0];
      r12[1] = c->centre(1) - centre_[1];
      r12[2] = c->centre(2) - centre_[2];
      vector<complex<double>> smoment = shift_multipoles(c->multipole(), r12);
      blas::ax_plus_y_n(1.0, smoment.data(), nmult_, multipole_.data());
    }
  }
}


void Box::compute_M2M_X(const VectorB& ocoeff_i, const VectorB& ocoeff_j) {

  mlm_ji_.resize(nmult_); fill(mlm_ji_.begin(), mlm_ji_.end(), complex<double>(0.0, 0.0));
  olm_ji_.resize(nmult_); fill(olm_ji_.begin(), olm_ji_.end(), complex<double>(0.0, 0.0));

  if (nchild() == 0) { // leaf
    VectorB ocoeff_si(nsize_);
    int start = 0;
    for (auto& ci : coffsets_s_) {
      const int ndim = ci.second;
      const VectorB c_si = ocoeff_i.slice(ci.first, ci.first + ndim);
      for (int idx = 0; idx != ndim; ++idx)
        ocoeff_si[idx+start] = c_si[idx];
      start += ndim;
    }
    assert(start == nsize_);
    VectorB ocoeff_uj(msize_);
    start = 0;
    for (auto& cj : coffsets_u_) {
      const int mdim = cj.second;
      const VectorB c_uj = ocoeff_j.slice(cj.first, cj.first + mdim);
      for (int idx = 0; idx != mdim; ++idx)
        ocoeff_uj[idx+start] = c_uj[idx];
      start += mdim;
    }
    assert(start == msize_);
    for (int k = 0; k != nmult_; ++k) {
      const VectorB olm_sj_real = *box_olm_[k]->get_real_part() * ocoeff_uj;
      const VectorB olm_sj_imag = *box_olm_[k]->get_imag_part() * ocoeff_uj;
      const double tmp1 = olm_sj_real % ocoeff_si;
      const double tmp2 = olm_sj_imag % ocoeff_si;
      olm_ji_[k] = complex<double>(tmp1, tmp2);
    }
  } else { // shift children's multipoles
    for (int n = 0; n != nchild(); ++n) {
      shared_ptr<const Box> c = child(n);
      array<double, 3> r12;
      r12[0] = c->centre(0) - centre_[0];
      r12[1] = c->centre(1) - centre_[1];
      r12[2] = c->centre(2) - centre_[2];

      vector<complex<double>> smoment = shift_multipoles(c->olm_ji(), r12);
      blas::ax_plus_y_n(1.0, smoment.data(), nmult_, olm_ji_.data());
    }
  }
}


void Box::compute_multipolesX() {

  assert(nchild() == 0);
  get_offsets();
  box_olm_.resize(nmult_);
  for (int k = 0; k != nmult_; ++k)
    box_olm_[k] = make_shared<ZMatrix>(nsize_, msize_);

  TaskQueue<function<void(void)>> tasks(nsp());
  int isp = 0;
  for (auto& v : sp_) {
    if (v->schwarz() < 1e-15) continue;

    tasks.emplace_back(
      [this, &v, isp]() {
        MultipoleBatch mpole(v->shells(), centre_, lmax_);
        mpole.compute();
        for (int k = 0; k != nmult_; ++k) {
          const int nstart = offsets_[isp][0];
          const int mstart = offsets_[isp][1];
          box_olm_[k]->add_block(1.0, nstart, mstart, v->shell(0)->nbasis(), v->shell(1)->nbasis(), mpole.data(k));
        }
      }
    );
    ++isp;
  }
  tasks.compute();
}


void Box::compute_L2L_X() {

  // from parent
  if (parent_) {
    array<double, 3> r12;
    r12[0] = centre_[0] - parent_->centre(0);
    r12[1] = centre_[1] - parent_->centre(1);
    r12[2] = centre_[2] - parent_->centre(2);
    vector<complex<double>> slocal = shift_localL(parent_->mlm_ji(), r12);
    blas::ax_plus_y_n(1.0, slocal.data(), nmult_, mlm_ji_.data());
  }
}


void Box::compute_L2L() {

  // from parent
  if (parent_) {
    array<double, 3> r12;
    r12[0] = centre_[0] - parent_->centre(0);
    r12[1] = centre_[1] - parent_->centre(1);
    r12[2] = centre_[2] - parent_->centre(2);
    vector<complex<double>> slocal = shift_localL(parent_->localJ(), r12);
    blas::ax_plus_y_n(1.0, slocal.data(), nmult_, localJ_.data());
  }
}


void Box::compute_M2L_X() {

  // from interaction list
  for (auto& it : inter_) {
    array<double, 3> r12;
    r12[0] = centre_[0] - it->centre(0);
    r12[1] = centre_[1] - it->centre(1);
    r12[2] = centre_[2] - it->centre(2);
    vector<complex<double>> slocal = shift_localM(it->olm_ji(), r12);
    blas::ax_plus_y_n(1.0, slocal.data(), nmult_, mlm_ji_.data());
  }
}


void Box::compute_M2L() {

  // from interaction list
  for (auto& it : inter_) {
    array<double, 3> r12;
    r12[0] = centre_[0] - it->centre(0);
    r12[1] = centre_[1] - it->centre(1);
    r12[2] = centre_[2] - it->centre(2);
    vector<complex<double>> slocal = shift_localM(it->multipole(), r12);
    blas::ax_plus_y_n(1.0, slocal.data(), nmult_, localJ_.data());
  }
}


shared_ptr<const ZMatrix> Box::compute_Fock_nf(shared_ptr<const Matrix> density, vector<double>& max_den, const double schwarz_thresh) const {

  assert(nchild() == 0);
  auto out = make_shared<ZMatrix>(density->ndim(), density->ndim());
  out->zero();

  // NF: 4c integrals
  const int shift = sizeof(int) * 4;
  const int nsh = sqrt(max_den.size());
  assert (nsh*nsh == max_den.size());

  int ntask = 0;
  for (auto& v01 : sp_) {
    if (v01->schwarz() < 1e-15) continue;
    const int i0 = v01->shell_ind(1);

    const int i1 = v01->shell_ind(0);
    if (i1 < i0) continue;

    const int i01 = i0 * nsh + i1;
    const double density_01 = max_den[i01] * 4.0;

    for (auto& neigh : neigh_) {
      for (auto& v23 : neigh->sp()) {
        if (v23->schwarz() < 1e-15) continue;
        const int i2 = v23->shell_ind(1);
        if (i2 < i0) continue;
        const int i3 = v23->shell_ind(0);
        if (i3 < i2) continue;
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
        ++ntask;
      }
    }
  }

  TaskQueue<function<void(void)>> tasks(ntask);
  mutex jmutex;

  for (auto& v01 : sp_) {
    if (v01->schwarz() < 1e-15) continue;
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
        if (v23->schwarz() < 1e-15) continue;
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

        array<shared_ptr<const Shell>,4> input = {{b3, b2, b1, b0}};

        tasks.emplace_back(
          [this, &out, &density, input, b0offset, i01, i23, b0size, b1offset, b1size, b2offset, b2size, b3offset, b3size, mulfactor, &jmutex]() {

            ERIBatch eribatch(input, mulfactor);
            eribatch.compute();
            const double* eridata = eribatch.data();

            const double* density_data = density->data();
            lock_guard<mutex> lock(jmutex);
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
        );

      }
    }
  }

  tasks.compute();

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


// M2L: given O(a) and R(a-b) get M(b)
vector<complex<double>> Box::shift_localM(const vector<complex<double>>& olm, array<double, 3> r12) const {

  const double r = sqrt(r12[0]*r12[0] + r12[1]*r12[1] + r12[2]*r12[2]);
  const double ctheta = (r > numerical_zero__) ? r12[2]/r : 0.0;
  const double phi = atan2(r12[1], r12[0]);

  vector<complex<double>> mb(nmult_, 0);

  for (int l = 0; l <= lmax_; ++l) {
    for (int m = 0; m <= 2 * l; ++m) {

      for (int j = 0; j <= lmax_; ++j) {
        for (int k = 0; k <= 2 * j; ++k) {

          const int a = l + j;
          const int b = m - l + k - j;

          double prefactor = plm.compute(a, abs(b), ctheta) / pow(r, a + 1);
          double ft = 1.0;
          for (int i = 1; i <= a - abs(b); ++i) {
            prefactor *= ft;
            ++ft;
          }
          const double real = (b >= 0) ? (prefactor * cos(abs(b) * phi)) : (pow(-1.0, b) * prefactor * cos(abs(b) * phi));
          const double imag = (b >= 0) ? (prefactor * sin(abs(b) * phi)) : (pow(-1.0, b+1) * prefactor * sin(abs(b) * phi));
          const complex<double> Mab(real, imag);
          mb[l*l+m] += pow(-1.0, l) * Mab * olm[j*j+k];
        }
      }
    }
  }

  return mb;
}


// M2M: given O(a) and R(b-a) get O(b)
vector<complex<double>> Box::shift_multipoles(const vector<complex<double>>& oa, array<double, 3> rab) const {

  const double r = sqrt(rab[0]*rab[0] + rab[1]*rab[1] + rab[2]*rab[2]);
  const double ctheta = (r > numerical_zero__) ? rab[2]/r : 0.0;
  const double phi = atan2(rab[1], rab[0]);

  vector<complex<double>> ob(nmult_, 0);

  for (int l = 0; l <= lmax_; ++l) {
    for (int m = 0; m <= 2 * l; ++m) {

      for (int j = 0; j <= lmax_; ++j) {
        for (int k = 0; k <= 2 * j; ++k) {

          const int a = l - j;
          const int b = m - l - k + j;
          if (abs(b) <= a && a >= 0) {
            double prefactor = pow(r, a) * plm.compute(a, abs(b), ctheta);
            double ft = 1.0;
            for (int i = 1; i <= a + abs(b); ++i) {
              prefactor /= ft;
              ++ft;
            }
            const double real = (b >= 0) ? (prefactor * cos(abs(b) * phi)) : (pow(-1.0, b) * prefactor * cos(abs(b) * phi));
            const double imag = (b >= 0) ? (prefactor * sin(abs(b) * phi)) : (pow(-1.0, b+1) * prefactor * sin(abs(b) * phi));
            const complex<double> Oab(real, -imag);

            ob[l*l+m] += Oab * oa[j*j+k];
          }

        }
      }
    }
  }

  return ob;
}


// L2L: given M(r) and R(b) get M(r-b)
vector<complex<double>> Box::shift_localL(const vector<complex<double>>& mr, array<double, 3> rb) const {

  const double r = sqrt(rb[0]*rb[0] + rb[1]*rb[1] + rb[2]*rb[2]);
  const double ctheta = (r > numerical_zero__) ? rb[2]/r : 0.0;
  const double phi = atan2(rb[1], rb[0]);

  vector<complex<double>> mrb(nmult_, 0);

  for (int l = 0; l <= lmax_; ++l) {
    for (int m = 0; m <= 2 * l; ++m) {

      for (int j = 0; j <= lmax_; ++j) {
        for (int k = 0; k <= 2 * j; ++k) {

          const int a = j - l;
          const int b = k - j - m + l;
          if (abs(b) <= a && a >= 0) {
            double prefactor = pow(r, a) * plm.compute(a, abs(b), ctheta);
            double ft = 1.0;
            for (int i = 1; i <= a + abs(b); ++i) {
              prefactor /= ft;
              ++ft;
            }
            const double real = (b >= 0) ? (prefactor * cos(abs(b) * phi)) : (pow(-1.0, b) * prefactor * cos(abs(b) * phi));
            const double imag = (b >= 0) ? (prefactor * sin(abs(b) * phi)) : (pow(-1.0, b+1) * prefactor * sin(abs(b) * phi));
            const complex<double> Oab(real, -imag);

            mrb[l*l+m] += Oab * mr[j*j+k];
          }

        }
      }
    }
  }

  return mrb;
}


shared_ptr<const ZMatrix> Box::compute_Fock_ff(shared_ptr<const Matrix> density) const {

  assert(nchild() == 0);
  auto out = make_shared<ZMatrix>(density->ndim(), density->mdim());

  // compute multipoles again - not efficient - for now
  TaskQueue<function<void(void)>> tasks(nsp());
  mutex jmutex;
  for (auto& v : sp_) {
    if (v->schwarz() < 1e-15) continue;

    tasks.emplace_back(
      [this, &out, &v, &density, &jmutex]() {
        MultipoleBatch mpole(v->shells(), centre_, lmax_);
        mpole.compute();
        const int dimb0 = v->shell(0)->nbasis();
        const int dimb1 = v->shell(1)->nbasis();
        vector<const complex<double>*> dat(nmult_);
        for (int i = 0; i != nmult_; ++i)
          dat[i] = mpole.data() + mpole.size_block()*i;

        lock_guard<mutex> lock(jmutex);
        for (int i = v->offset(1); i != dimb1 + v->offset(1); ++i)
          for (int j = v->offset(0); j != dimb0 + v->offset(0); ++j)
            for (int k = 0; k != mpole.num_blocks(); ++k)
              out->element(i, j) += conj(*dat[k]++) * localJ_[k];
      }
    );
  }
  tasks.compute();

  return out;
}


shared_ptr<const ZVectorB> Box::compute_Fock_ffX(const VectorB& ocoeff_j) const {
  assert(nchild() == 0);
  auto out = make_shared<ZVectorB>(ocoeff_j.size());

  VectorB ocoeff_uj(msize_);
  int start = 0;
  for (auto& cu : coffsets_u_) {
    const int mdim = cu.second;
    const VectorB c_uj = ocoeff_j.slice(cu.first, cu.first + mdim);
    for (int idx = 0; idx != mdim; ++idx)
      ocoeff_uj(idx+start) = c_uj(idx);
    start += mdim;
  }
  assert(start == msize_);

  VectorB olm_sj_real(nsize_);
  VectorB olm_sj_imag(nsize_);
  for (int k = 0; k != nmult_; ++k) {
    olm_sj_real = *box_olm_[k]->get_real_part() * ocoeff_uj;
    olm_sj_imag = *box_olm_[k]->get_imag_part() * ocoeff_uj;
  }
  ZVectorB olm_sj(nsize_);
  for (int i = 0; i != nsize_; ++i)
    olm_sj(i) = complex<double>(olm_sj_real(i), olm_sj_imag(i));

  start = 0;
  for (auto& cs : coffsets_s_) {
    const int ndim = cs.second;
    int index = cs.first;
    for (int i = 0; i != ndim; ++i) {
      for (int k = 0; k != nmult_; ++k)
        (*out)[index] += olm_sj[i + start] * mlm_ji_[k];
      ++index;
    }
    start += ndim;
  }

  return out;
}


double Box::compute_exact_energy_ff(shared_ptr<const Matrix> density) const { //for debug

  auto jff = make_shared<Matrix>(density->ndim(), density->ndim());
  jff->zero();

  int ntask = 0;
  for (auto& inter : inter_)
    ntask += inter->sp().size();
  ntask *= sp_.size();

  TaskQueue<function<void(void)>> tasks(ntask);
  mutex jmutex;

  for (auto& v01 : sp_) {
    if (v01->schwarz() < 1e-15) continue;
    shared_ptr<const Shell> b0 = v01->shell(1);
    const int b0offset = v01->offset(1);
    const int b0size = b0->nbasis();

    shared_ptr<const Shell> b1 = v01->shell(0);
    const int b1offset = v01->offset(0);
    const int b1size = b1->nbasis();

    for (auto& inter : inter_) {
      for (auto& v23 : inter->sp()) {
        if (v23->schwarz() < 1e-15) continue;
        shared_ptr<const Shell> b2 = v23->shell(1);
        const int b2offset = v23->offset(1);
        const int b2size = b2->nbasis();

        shared_ptr<const Shell> b3 = v23->shell(0);
        const int b3offset = v23->offset(0);
        const int b3size = b3->nbasis();

        array<shared_ptr<const Shell>,4> input = {{b3, b2, b1, b0}};

        tasks.emplace_back(
          [this, &jff, &density, input, b0offset, b0size, b1offset, b1size, b2offset, b2size, b3offset, b3size, &jmutex]() {
            const double* density_data = density->data();

            ERIBatch eribatch(input, 0.0);
            eribatch.compute();
            const double* eridata = eribatch.data();

            lock_guard<mutex> lock(jmutex);
            for (int j0 = b0offset; j0 != b0offset + b0size; ++j0) {
              for (int j1 = b1offset; j1 != b1offset + b1size; ++j1) {
                const int j1n = j1 * density->ndim();
                for (int j2 = b2offset; j2 != b2offset + b2size; ++j2) {
                  const int j2n = j2 * density->ndim();
                  for (int j3 = b3offset; j3 != b3offset + b3size; ++j3, ++eridata) {
                    const double eri = *eridata;
                    jff->element(j1, j0) += density_data[j2n + j3] * eri;
                    jff->element(j3, j0) -= density_data[j1n + j2] * eri * 0.5;
                  }
                }
              }
            }

          }
        );

      }
    }
  }
  tasks.compute();
  const double out = 0.5*jff->dot_product(*density);

  return out;
}


// M2M for X
vector<shared_ptr<const ZMatrix>> Box::shift_multipolesX(const vector<shared_ptr<ZMatrix>>& oa, array<double, 3> rab) const {

  const double r = sqrt(rab[0]*rab[0] + rab[1]*rab[1] + rab[2]*rab[2]);
  const double ctheta = (r > numerical_zero__) ? rab[2]/r : 0.0;
  const double phi = atan2(rab[1], rab[0]);

  vector<shared_ptr<const ZMatrix>> ob(nmult_);

  for (int l = 0; l <= lmax_; ++l) {
    for (int m = 0; m <= 2 * l; ++m) {

      auto out = make_shared<ZMatrix>(oa[0]->ndim(), oa[0]->mdim());

      for (int j = 0; j <= lmax_; ++j) {
        for (int k = 0; k <= 2 * j; ++k) {

          const int a = l - j;
          const int b = m - l - k + j;
          if (abs(b) <= a && a >= 0) {
            double prefactor = pow(r, a) * plm.compute(a, abs(b), ctheta);
            double ft = 1.0;
            for (int i = 1; i <= a + abs(b); ++i) {
              prefactor /= ft;
              ++ft;
            }
            const double real = (b >= 0) ? (prefactor * cos(abs(b) * phi)) : (pow(-1.0, b) * prefactor * cos(abs(b) * phi));
            const double imag = (b >= 0) ? (prefactor * sin(abs(b) * phi)) : (pow(-1.0, b+1) * prefactor * sin(abs(b) * phi));
            const complex<double> Oab(real, -imag);

            out->add_block(Oab, 0, 0, oa[0]->ndim(), oa[0]->mdim(), oa[j*j+k]->data());
          }

        }
      }
      ob[l*l+m] = make_shared<const ZMatrix>(*out);
    }
  }

  return ob;
}


// L2L for X
vector<shared_ptr<const ZMatrix>> Box::shift_localLX(const vector<shared_ptr<ZMatrix>>& mr, array<double, 3> rb) const {

  const double r = sqrt(rb[0]*rb[0] + rb[1]*rb[1] + rb[2]*rb[2]);
  const double ctheta = (r > numerical_zero__) ? rb[2]/r : 0.0;
  const double phi = atan2(rb[1], rb[0]);

  vector<shared_ptr<const ZMatrix>> mrb(nmult_);

  for (int l = 0; l <= lmax_; ++l) {
    for (int m = 0; m <= 2 * l; ++m) {

      auto out = make_shared<ZMatrix>(mr[0]->ndim(), mr[0]->mdim());

      for (int j = 0; j <= lmax_; ++j) {
        for (int k = 0; k <= 2 * j; ++k) {

          const int a = j - l;
          const int b = k - j - m + l;
          if (abs(b) <= a && a >= 0) {
            double prefactor = pow(r, a) * plm.compute(a, abs(b), ctheta);
            double ft = 1.0;
            for (int i = 1; i <= a + abs(b); ++i) {
              prefactor /= ft;
              ++ft;
            }
            const double real = (b >= 0) ? (prefactor * cos(abs(b) * phi)) : (pow(-1.0, b) * prefactor * cos(abs(b) * phi));
            const double imag = (b >= 0) ? (prefactor * sin(abs(b) * phi)) : (pow(-1.0, b+1) * prefactor * sin(abs(b) * phi));
            const complex<double> Oab(real, -imag);

            out->add_block(Oab, 0, 0, mr[0]->ndim(), mr[0]->mdim(),  mr[j*j+k]->data());
          }

        }
      }

      mrb[l*l+m] = make_shared<const ZMatrix>(*out);
    }
  }

  return mrb;
}


// M2L for X
vector<shared_ptr<const ZMatrix>> Box::shift_localMX(const vector<shared_ptr<ZMatrix>>& olm, array<double, 3> r12) const {

  const double r = sqrt(r12[0]*r12[0] + r12[1]*r12[1] + r12[2]*r12[2]);
  const double ctheta = (r > numerical_zero__) ? r12[2]/r : 0.0;
  const double phi = atan2(r12[1], r12[0]);

  vector<shared_ptr<const ZMatrix>> mb(nmult_);

  for (int l = 0; l <= lmax_; ++l) {
    for (int m = 0; m <= 2 * l; ++m) {
 
      auto out = make_shared<ZMatrix>(olm[0]->ndim(), olm[0]->mdim());

      for (int j = 0; j <= lmax_; ++j) {
        for (int k = 0; k <= 2 * j; ++k) {

          const int a = l + j;
          const int b = m - l + k - j;

          double prefactor = plm.compute(a, abs(b), ctheta) / pow(r, a + 1);
          double ft = 1.0;
          for (int i = 1; i <= a - abs(b); ++i) {
            prefactor *= ft;
            ++ft;
          }
          const double real = (b >= 0) ? (prefactor * cos(abs(b) * phi)) : (pow(-1.0, b) * prefactor * cos(abs(b) * phi));
          const double imag = (b >= 0) ? (prefactor * sin(abs(b) * phi)) : (pow(-1.0, b+1) * prefactor * sin(abs(b) * phi));
          const complex<double> Mab(real, imag);

          const complex<double> coeff = pow(-1.0, l) * Mab;
          out->add_block(coeff, 0, 0, olm[0]->ndim(), olm[0]->mdim(),  olm[j*j+k]->data());
        }
      }

      mb[l*l+m] = make_shared<const ZMatrix>(*out);
    }
  }

  return mb;
}


void Box::get_offsets() {
  offsets_.resize(nsp()); //starting pos of sp in box (n x m)

  map<int, int> shells;   //key = shell index, starting pos of shell ndim or mdim
  int nishell = 0;
  int isize = 0;
  for (int isp = 0; isp != nsp(); ++isp) {
    const int key = sp_[isp]->shell_ind(0);
    map<int,int>::iterator sh = shells.find(key);
    const bool sh_found = (sh != shells.end());
    if (sh_found) {
      offsets_[isp][0] = shells.find(key)->second;
    } else {
      shells.insert(shells.end(), pair<int, int>(key, nishell));
      offsets_[isp][0] = isize;

      coffsets_s_.resize(nishell+1);
      auto ci = make_pair<int, int>(sp_[isp]->offset(0), sp_[isp]->shell(0)->nbasis());
      coffsets_s_[nishell] = ci;

      ++nishell;
      isize += sp_[isp]->shell(0)->nbasis();
    }
  }
  nsize_ = isize;

  shells.clear();
  int njshell = 0;
  int jsize = 0;
  for (int isp = 0; isp != nsp(); ++isp) {
    const int key = sp_[isp]->shell_ind(1);
    map<int,int>::iterator sh = shells.find(key);
    const bool sh_found = (sh != shells.end());
    if (sh_found) {
      offsets_[isp][1] = shells.find(key)->second;
    } else {
      shells.insert(shells.end(), pair<int, int>(key, njshell));
      offsets_[isp][0] = jsize;

      coffsets_u_.resize(njshell+1);
      auto cj = make_pair<int, int>(sp_[isp]->offset(1), sp_[isp]->shell(1)->nbasis());
      coffsets_u_[njshell] = cj;

      ++njshell;
      jsize += sp_[isp]->shell(1)->nbasis();
    }
  }
  msize_ = jsize;
}
