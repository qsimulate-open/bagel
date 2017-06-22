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
#include <src/util/math/factorial.h>

using namespace bagel;
using namespace std;

static const double pisq__ = pi__ * pi__;
const static Legendre plm;
const static Factorial f;

void Box::init() {

  nbasis0_ = 0;   nbasis1_ = 0;
  for (auto& i : sp_) {
    nbasis0_ += i->nbasis0();
    nbasis1_ += i->nbasis1();
  }

  centre_ = {{0, 0, 0}};
  extent_ = 0.0;
  nsp_ = 0;
  for (auto& i : sp_) {
    if (i->schwarz() < schwarz_thresh_) continue;
    ++nsp_;
  }

  if (nchild() == 0) {
    for (auto& i : sp_) {
      if (i->schwarz() < schwarz_thresh_) continue;
      for (int j = 0; j != 3; ++j) centre_[j] += i->centre(j);
    }
    centre_[0] /= nsp_;
    centre_[1] /= nsp_;
    centre_[2] /= nsp_;
    for (auto& i : sp_) {
      if (i->schwarz() < schwarz_thresh_) continue;
      double rad = 0;
      for (int j = 0; j != 3; ++j) rad += pow(i->centre(j)-centre_[j], 2.0);
      const double ei = sqrt(rad) + i->extent();
      assert(ei > 0);
      extent_ = max(extent_, ei);
    }
  } else {
    for (int n = 0; n != nchild(); ++n) {
      shared_ptr<const Box> i = child(n);
      for (int j = 0; j != 3; ++j)  centre_[j] += i->centre(j);
    }
    centre_[0] /= nchild();
    centre_[1] /= nchild();
    centre_[2] /= nchild();
    for (int n = 0; n != nchild(); ++n) {
      shared_ptr<const Box> i = child(n);
      double rad = 0.0;
      for (int j = 0; j != 3; ++j) rad += pow(i->centre(j)-centre_[j], 2.0);
      const double ei = sqrt(rad) + i->extent();
      assert(ei > 0);
      extent_ = max(extent_, ei);
    }
  }

  ws_ = (int) ceil(extent_/boxsize_);
  nmult_ = (lmax_ + 1) * (lmax_ + 1);
  multipole_ = make_shared<ZVectorB>(nmult_);
  localJ_ = make_shared<ZVectorB>(nmult_);
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


void Box::get_neigh(const vector<shared_ptr<Box>>& box, const double ws) {

  if ((box.empty()) || (box[0]->rank() != rank_)) return;

  neigh_.resize(box.size());
  int nn = 0;
  for (auto& b : box) {
    if (is_neigh(b, ws)) {
      neigh_[nn] = b;
      ++nn;
    }
  }
  neigh_.resize(nn);
}


bool Box::is_neigh(shared_ptr<const Box> box, const double ws) const {

#if 1
  if (nchild() == 0) {
    double rr = 0;
    for (int i = 0; i != 3; ++i) rr += pow(centre_[i] - box->centre(i), 2);
    return (sqrt(rr) <= (1.0+ws)*(extent_ + box->extent()));
  } else {
    for (int n0 = 0; n0 != nchild(); ++n0) {
      shared_ptr<const Box> i = child(n0);
      for (int n1 = 0; n1 != box->nchild(); ++n1) {
        shared_ptr<const Box> j = box->child(n1);
        if (i->is_neigh(j, ws)) return true;
      }
    }
    return false;
  }
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

#if 0
  const double ws1 = max(2.0*extent_/boxsize_, ws);
  const double ws2 = max(2.0*box->extent()/box->boxsize(), ws);
  const double sep = max(max(abs(tvec_[0]-box->tvec()[0]), abs(tvec_[1]-box->tvec()[1])), abs(tvec_[2]-box->tvec()[2]));
  const bool out = (2.0*sep <= ws1+ws2);

  //debug
  double rr = 0;
  for (int i = 0; i != 3; ++i)
    rr += pow(centre_[i] - box->centre(i), 2);

  if ((extent_+box->extent() > sqrt(rr)) && !out)
    cout << "  WARNING!!!!  " << setprecision(2) << sqrt(rr) << "  <= " << box->extent()+extent_ << endl;
  //end debug

  return out;
#endif
}

void Box::get_inter(const vector<shared_ptr<Box>>& box, const double ws) {

  if ((box.empty()) || (box[0]->rank() != rank_)) return;

  inter_.resize(box.size());
  int ni = 0;
  for (auto& b : box) {
    if ((!is_neigh(b, ws)) && (!parent_ || parent_->is_neigh(b->parent(), ws))) {
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
  assert(newsp.size() == sp_.size());
  sp_ = newsp;
}


void Box::compute_M2M(shared_ptr<const Matrix> density) {

  multipole_->fill(0.0);

  if (nchild() == 0) { // leaf
    TaskQueue<function<void(void)>> tasks(nsp_);
    mutex jmutex;
    for (auto& v : sp_) {
      if (v->schwarz() < thresh_) continue;

      tasks.emplace_back(
        [this, &v, &density, &jmutex]() {
          auto olm = make_shared<ZVectorB>(nmult_);
          MultipoleBatch mpole(v->shells(), centre_, lmax_);
          mpole.compute();
          const int dimb0 = v->shell(0)->nbasis();
          const int dimb1 = v->shell(1)->nbasis();

          for (int k = 0; k != mpole.num_blocks(); ++k) {
            size_t cnt = 0;
            for (int i = v->offset(1); i != dimb1 + v->offset(1); ++i) {
              (*olm)(k) += conj(blas::dot_product_noconj(mpole.data() + mpole.size_block()*k + cnt, dimb0, density->element_ptr(v->offset(0), i)));
              cnt += dimb0;
            }
          }

          lock_guard<mutex> lock(jmutex);
          blas::ax_plus_y_n(1.0, olm->data(), nmult_, multipole_->data());
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
      shared_ptr<const ZVectorB> smoment = shift_multipoles(c->multipole(), r12);
#if 0
      //*** debug using K
      auto cmultipole = make_shared<ZMatrix>(1, nmult_);
      copy_n(c->multipole()->data(), nmult_, cmultipole->data());
      shared_ptr<const ZMatrix> tmp = shift_multipolesX(lmax_, cmultipole, r12);
      auto smoment = make_shared<ZVectorB>(nmult_);
      copy_n(tmp->data(), nmult_, smoment->data());
      //*** end debug using K
#endif
      blas::ax_plus_y_n(1.0, smoment->data(), nmult_, multipole_->data());
    }
  }
}


//|su) C_ui C_sj
void Box::compute_M2M_X(shared_ptr<const Matrix> ocoeff_sj, shared_ptr<const Matrix> ocoeff_ui) {

  olm_ndim_ = ocoeff_sj->mdim();  olm_mdim_ = ocoeff_ui->mdim();
  const int nmult_k = (lmax_k_+1)*(lmax_k_+1);
  olm_ji_ = make_shared<ZMatrix>(olm_ndim_*olm_mdim_, nmult_k);

  if (nchild() == 0) { // leaf
    TaskQueue<function<void(void)>> tasks(nsp_);
    mutex kmutex;
    for (auto& v : sp_) {
      if (v->schwarz() < thresh_) continue;

      tasks.emplace_back(
        [this, nmult_k, &v, &ocoeff_sj, &ocoeff_ui, &kmutex]() {
          MultipoleBatch olm_su(v->shells(), centre_, lmax_k_);
          olm_su.compute();
          const int dimb0 = v->shell(0)->nbasis();
          const int dimb1 = v->shell(1)->nbasis();
          auto zsj = make_shared<const ZMatrix>(*(ocoeff_sj->cut(v->offset(0), v->offset(0)+dimb0)), 1.0);
          auto zui = make_shared<const ZMatrix>(*(ocoeff_ui->cut(v->offset(1), v->offset(1)+dimb1)), 1.0);
          auto tmp_olm_ji = make_shared<ZMatrix>(olm_ndim_*olm_mdim_, nmult_k);
          for (int k = 0; k != nmult_k; ++k) {
            ZMatrix olm_si(dimb0, ocoeff_ui->mdim());
            zgemm3m_("N", "N", dimb0, ocoeff_ui->mdim(), dimb1, 1.0, olm_su.data() + olm_su.size_block()*k, dimb0, zui->data(), dimb1, 0.0, olm_si.data(), dimb0);
            zgemm3m_("C", "N", ocoeff_sj->mdim(), ocoeff_ui->mdim(), dimb0, 1.0, zsj->data(), dimb0, olm_si.data(), dimb0, 0.0, tmp_olm_ji->data()+k*olm_ndim_*olm_mdim_, ocoeff_sj->mdim());
          }
          tmp_olm_ji = tmp_olm_ji->get_conjg()->copy();
          lock_guard<mutex> lock(kmutex);
          blas::ax_plus_y_n(1.0, tmp_olm_ji->data(), olm_ji_->size(), olm_ji_->data());
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

      shared_ptr<const ZMatrix> smoment = shift_multipolesX(lmax_k_, c->olm_ji(), r12);
      blas::ax_plus_y_n(1.0, smoment->data(), olm_ji_->size(), olm_ji_->data());
    }
  }
}


void Box::compute_L2L_X() {

  // from parent
  if (parent_) {
    array<double, 3> r12;
    r12[0] = centre_[0] - parent_->centre(0);
    r12[1] = centre_[1] - parent_->centre(1);
    r12[2] = centre_[2] - parent_->centre(2);
    shared_ptr<const ZMatrix> slocal = shift_localLX(lmax_k_, parent_->mlm_ji(), r12);
    blas::ax_plus_y_n(1.0, slocal->data(), mlm_ji_->size(), mlm_ji_->data());
  }
}


void Box::compute_L2L() {

  // from parent
  if (parent_) {
    array<double, 3> r12;
    r12[0] = centre_[0] - parent_->centre(0);
    r12[1] = centre_[1] - parent_->centre(1);
    r12[2] = centre_[2] - parent_->centre(2);
    shared_ptr<const ZVectorB> slocal = shift_localL(parent_->localJ(), r12);
#if 0
    //*** debug using K
    auto plocalJ = make_shared<ZMatrix>(1, nmult_);
    copy_n(parent_->localJ()->data(), nmult_, plocalJ->data());
    shared_ptr<const ZMatrix> tmp = shift_localLX(lmax_, plocalJ, r12);
    auto slocal = make_shared<ZVectorB>(nmult_);
    copy_n(tmp->data(), nmult_, slocal->data());
    //*** end debug using K
#endif
    blas::ax_plus_y_n(1.0, slocal->data(), nmult_, localJ_->data());
  }
}


void Box::compute_M2L_X() {

  mlm_ji_ = olm_ji_->clone();

  // from interaction list
  for (auto& it : inter_) {
    array<double, 3> r12;
    r12[0] = centre_[0] - it->centre(0);
    r12[1] = centre_[1] - it->centre(1);
    r12[2] = centre_[2] - it->centre(2);
    shared_ptr<const ZMatrix> slocal = shift_localMX(lmax_k_, it->olm_ji(), r12);
    blas::ax_plus_y_n(1.0, slocal->data(), mlm_ji_->size(), mlm_ji_->data());
  }
}


void Box::compute_M2L() {

  localJ_->fill(0.0);

  // from interaction list
  for (auto& it : inter_) {
    array<double, 3> r12;
    r12[0] = centre_[0] - it->centre(0);
    r12[1] = centre_[1] - it->centre(1);
    r12[2] = centre_[2] - it->centre(2);
    shared_ptr<const ZVectorB> slocal = shift_localM(it->multipole(), r12);
#if 0
    //*** debug using K
    auto imultipole = make_shared<ZMatrix>(1, nmult_);
    copy_n(it->multipole()->data(), nmult_, imultipole->data());
    shared_ptr<const ZMatrix> tmp = shift_localMX(lmax_, imultipole, r12);
    auto slocal = make_shared<ZVectorB>(nmult_);
    copy_n(tmp->data(), nmult_, slocal->data());
    //*** end debug using K
#endif
    blas::ax_plus_y_n(1.0, slocal->data(), nmult_, localJ_->data());
  }
}


shared_ptr<const Matrix> Box::compute_Fock_nf(shared_ptr<const Matrix> density, shared_ptr<const VectorB> max_den) const {

  assert(nchild() == 0);

  auto out = make_shared<Matrix>(density->ndim(), density->mdim());
  out->zero();

  // NF: 4c integrals
  const int shift = sizeof(int) * 4;
  const int nsh = sqrt(max_den->size());
  assert (nsh*nsh == max_den->size());

  int ntask = 0;
  for (auto& v01 : sp_) {
    if (v01->schwarz() < thresh_) continue;
    const int i0 = v01->shell_ind(1);

    const int i1 = v01->shell_ind(0);
    if (i1 < i0) continue;

    const int i01 = i0 * nsh + i1;
    const double density_01 = (*max_den)(i01) * 4.0;

    for (auto& neigh : neigh_) {
      for (auto& v23 : neigh->sp()) {
        if (v23->schwarz() < thresh_) continue;
        const int i2 = v23->shell_ind(1);
        if (i2 < i0) continue;
        const int i3 = v23->shell_ind(0);
        if (i3 < i2) continue;
        const int i23 = i2 * nsh + i3;
        if (i23 < i01) continue;

        const double density_23 = (*max_den)(i23) * 4.0;
        const double density_02 = (*max_den)(i0 * nsh + i2);
        const double density_03 = (*max_den)(i0 * nsh + i3);
        const double density_12 = (*max_den)(i1 * nsh + i2);
        const double density_13 = (*max_den)(i1 * nsh + i3);

        const double mulfactor = max(max(max(density_01, density_02),
                                         max(density_03, density_12)),
                                         max(density_13, density_23));
        const double integral_bound = mulfactor * v01->schwarz() * v23->schwarz();
        const bool skip_schwarz = integral_bound < schwarz_thresh_;
        if (skip_schwarz) continue;
        ++ntask;
      }
    }
  }

  TaskQueue<function<void(void)>> tasks(ntask);
  mutex jmutex;

  for (auto& v01 : sp_) {
    if (v01->schwarz() < thresh_) continue;
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
    const double density_01 = (*max_den)(i01) * 4.0;

    for (auto& neigh : neigh_) {
      for (auto& v23 : neigh->sp()) {
        if (v23->schwarz() < thresh_) continue;
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

        const double density_23 = (*max_den)(i23) * 4.0;
        const double density_02 = (*max_den)(i0 * nsh + i2);
        const double density_03 = (*max_den)(i0 * nsh + i3);
        const double density_12 = (*max_den)(i1 * nsh + i2);
        const double density_13 = (*max_den)(i1 * nsh + i3);

        const double mulfactor = max(max(max(density_01, density_02),
                                         max(density_03, density_12)),
                                         max(density_13, density_23));
        const double integral_bound = mulfactor * v01->schwarz() * v23->schwarz();
        const bool skip_schwarz = integral_bound < schwarz_thresh_;
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
  cout << "Box " << boxid_ << " Rank = " << rank_ << " *** nchild = " << nchild() << " *** nsp = " << sp_.size()
       << " *** nneigh = " << nneigh() << " *** ninter = " << ninter() << " *** extent = " << extent_
       << " *** centre = " << setprecision(3) << centre_[0] << "  " << centre_[1] << "  " << centre_[2] << endl;

//  cout << " *** Shell pairs at ***" << endl;
//  for (int i = 0; i != sp_.size(); ++i)
//    cout << setprecision(5) << sp(i)->centre(0) << "  " << sp(i)->centre(1) << "  " << sp(i)->centre(2) << endl;
}


// M2L: given O(a) and R(a-b) get M(b)
shared_ptr<const ZVectorB> Box::shift_localM(shared_ptr<const ZVectorB> olm, array<double, 3> r12) const {

  const double r = sqrt(r12[0]*r12[0] + r12[1]*r12[1] + r12[2]*r12[2]);
  const double ctheta = (r > numerical_zero__) ? r12[2]/r : 0.0;
  const double phi = atan2(r12[1], r12[0]);

  auto mb = make_shared<ZVectorB>(nmult_);
  unique_ptr<double[]> plm0(new double[(2*lmax_+1)*(2*lmax_+1)]);
  for (int l = 0; l != 2*lmax_+1; ++l)
    for (int m = 0; m <= 2 * l; ++m) {
      const int b = m-l;
      const double sign = (b >= 0) ? 1.0 : (1-((b&1)<<1));
      plm0[l*l+m] = sign * plm.compute(l, abs(b), ctheta);
    }

  for (int l = 0; l <= lmax_; ++l) {
    const double phase_l = (1-((l&1)<<1));
    for (int j = 0; j <= lmax_; ++j) {
      const int a = l + j;
      const double rr = 1.0 / pow(r, a + 1);
      for (int m = 0; m <= 2 * l; ++m) {
        for (int k = 0; k <= 2 * j; ++k) {
          const int b = m - l + k - j;
          double prefactor = plm0[a*a+a+b] * rr * f(a-abs(b));
          const complex<double> Mab = polar(prefactor, b*phi);
          (*mb)(l*l+m) += phase_l*Mab * (*olm)(j*j+k);
        }
      }
    }
  }

  return mb;
}


// M2M: given O(a) and R(b-a) get O(b)
shared_ptr<const ZVectorB> Box::shift_multipoles(const shared_ptr<const ZVectorB> oa, array<double, 3> rab) const {

  const double r = sqrt(rab[0]*rab[0] + rab[1]*rab[1] + rab[2]*rab[2]);
  const double ctheta = (r > numerical_zero__) ? rab[2]/r : 0.0;
  const double phi = atan2(rab[1], rab[0]);

  auto ob = make_shared<ZVectorB>(nmult_);
  unique_ptr<double[]> plm0(new double[nmult_]);
  for (int l = 0; l != lmax_+1; ++l)
    for (int m = 0; m <= 2 * l; ++m) {
      const int b = m-l;
      const double sign = (b >= 0) ? 1.0 : (1-((b&1)<<1));
      plm0[l*l+m] = sign * plm.compute(l, abs(b), ctheta);
    }

  unique_ptr<double[]> invfac(new double[2*lmax_+1]);
  fill_n(invfac.get(), 2*lmax_+1, 1.0);
  for (int i = 1; i <= 2*lmax_; ++i)
    for (int j = i; j <= 2*lmax_; ++j)
      invfac[j] /= i;

  for (int l = 0; l <= lmax_; ++l) {
    for (int j = 0; j <= lmax_; ++j) {
      const int a = l - j;
      const double rr = pow(r, a);
      if (a < 0) continue;
      for (int m = 0; m <= 2 * l; ++m) {
        const int lo = max(-a, m-l-j);
        const int hi = min( a, m-l+j);
        for (int b = lo; b <= hi; ++b) {
          const int k = m - l - b + j;
          const double prefactor = rr * plm0[a*a+a+b] * invfac[a+abs(b)];
          const complex<double> Oab = polar(prefactor, -b*phi);
          (*ob)(l*l+m) += Oab * (*oa)(j*j+k);
        }
      }
    }
  }

  return ob;
}


// L2L: given M(r) and R(b) get M(r-b)
shared_ptr<const ZVectorB> Box::shift_localL(const shared_ptr<const ZVectorB> mr, array<double, 3> rb) const {

  const double r = sqrt(rb[0]*rb[0] + rb[1]*rb[1] + rb[2]*rb[2]);
  const double ctheta = (r > numerical_zero__) ? rb[2]/r : 0.0;
  const double phi = atan2(rb[1], rb[0]);

  auto mrb = make_shared<ZVectorB>(nmult_);

  unique_ptr<double[]> plm0(new double[nmult_]);
  for (int l = 0; l != lmax_+1; ++l)
    for (int m = 0; m <= 2 * l; ++m) {
      const int b = m-l;
      const double sign = (b >= 0) ? 1.0 : (1-((b&1)<<1));
      plm0[l*l+m] = sign * plm.compute(l, abs(b), ctheta);
    }

  unique_ptr<double[]> invfac(new double[2*lmax_+1]);
  fill_n(invfac.get(), 2*lmax_+1, 1.0);
  for (int i = 1; i <= 2*lmax_; ++i)
    for (int j = i; j <= 2*lmax_; ++j)
      invfac[j] /= i;

  for (int l = 0; l <= lmax_; ++l) {
    for (int j = 0; j <= lmax_; ++j) {
      const int a = j - l;
      const double rr = pow(r, a);
      if (a < 0) continue;
      for (int m = 0; m <= 2 * l; ++m) {
        const int lo = max(-a, l-m-j);
        const int hi = min( a, l-m+j);
        for (int b = lo; b <= hi; ++b) {
          const int k = m - l + b + j;
          const double prefactor = rr * plm0[a*a+a+b] * invfac[a+abs(b)];
          const complex<double> Oab = polar(prefactor, -b*phi);
          (*mrb)(l*l+m) += Oab * (*mr)(j*j+k);
        }
      }
    }
  }

  return mrb;
}


shared_ptr<const Matrix> Box::compute_Fock_ff(shared_ptr<const Matrix> density) const {

  assert(nchild() == 0);
  auto out = make_shared<ZMatrix>(density->ndim(), density->mdim());

  // compute multipoles again - not efficient - for now
  TaskQueue<function<void(void)>> tasks(nsp_);
  mutex jmutex;
  for (auto& v : sp_) {
    if (v->schwarz() < thresh_) continue;

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
              out->element(i, j) += conj(*dat[k]++) * (*localJ_)(k);
      }
    );
  }
  tasks.compute();
  assert(out->get_imag_part()->rms() < 1e-10);

  return out->get_real_part();
}


shared_ptr<const Matrix> Box::compute_Fock_ff_K(shared_ptr<const Matrix> ocoeff_ti) const {

  assert(nchild() == 0);
  auto out = make_shared<Matrix>(ocoeff_ti->ndim(), olm_ndim_);
  assert(ocoeff_ti->mdim() == olm_mdim_);
  const int nmult_k = (lmax_k_+1)*(lmax_k_+1);
  TaskQueue<function<void(void)>> tasks(nsp_);
  mutex kmutex;
  for (auto& v : sp_) {
    if (v->schwarz() < thresh_) continue;

    tasks.emplace_back(
      [this, nmult_k, &v, &ocoeff_ti, &out, &kmutex]() {
        MultipoleBatch olm_rt(v->shells(), centre_, lmax_k_);
        olm_rt.compute();
        const int dimb0 = v->shell(0)->nbasis();
        const int dimb1 = v->shell(1)->nbasis();
        auto zti = make_shared<const ZMatrix>(*(ocoeff_ti->cut(v->offset(1), v->offset(1)+dimb1)), 1.0);

        Matrix krj(dimb0, olm_ndim_);
        for (int k = 0; k != nmult_k; ++k) {
          ZMatrix olm_ri(dimb0, olm_mdim_);
          zgemm3m_("N", "N", dimb0, olm_mdim_, dimb1, 1.0, olm_rt.data() + olm_rt.size_block()*k, dimb0, zti->data(), dimb1, 0.0, olm_ri.data(), dimb0);
          dgemm_("N", "T", dimb0, olm_ndim_, olm_mdim_, 1.0, (olm_ri.get_real_part())->data(), dimb0, 
                                              (mlm_ji_->get_real_part())->data()+k*olm_ndim_*olm_mdim_, olm_ndim_, 1.0, krj.data(), dimb0);
          dgemm_("N", "T", dimb0, olm_ndim_, olm_mdim_, 1.0, (olm_ri.get_imag_part())->data(), dimb0, 
                                              (mlm_ji_->get_imag_part())->data()+k*olm_ndim_*olm_mdim_, olm_ndim_, 1.0, krj.data(), dimb0);
        }
        lock_guard<mutex> lock(kmutex);
        out->add_block(1.0, v->offset(0), 0, dimb0, olm_ndim_, krj.data());
      }
    );
  }
  tasks.compute();

  return out;
}


double Box::compute_exact_energy_ff(shared_ptr<const Matrix> density) const { //for debug

  auto jff = make_shared<Matrix>(density->ndim(), density->mdim());
  jff->zero();

  int ntask = 0;
  for (auto& inter : inter_)
    ntask += inter->sp().size();
  ntask *= sp_.size();

  TaskQueue<function<void(void)>> tasks(ntask);
  mutex jmutex;

  for (auto& v01 : sp_) {
    if (v01->schwarz() < thresh_) continue;
    shared_ptr<const Shell> b0 = v01->shell(1);
    const int b0offset = v01->offset(1);
    const int b0size = b0->nbasis();

    shared_ptr<const Shell> b1 = v01->shell(0);
    const int b1offset = v01->offset(0);
    const int b1size = b1->nbasis();

    for (auto& inter : inter_) {
      for (auto& v23 : inter->sp()) {
        if (v23->schwarz() < thresh_) continue;
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
shared_ptr<const ZMatrix> Box::shift_multipolesX(const int lmax, shared_ptr<const ZMatrix> oa, array<double, 3> rab) const {

  const double r = sqrt(rab[0]*rab[0] + rab[1]*rab[1] + rab[2]*rab[2]);
  const double ctheta = (r > numerical_zero__) ? rab[2]/r : 0.0;
  const double phi = atan2(rab[1], rab[0]);
  const int nmult = (lmax+1)*(lmax+1);
  const int olm_size_block = oa->ndim();

  shared_ptr<ZMatrix> ob = oa->clone();
  unique_ptr<double[]> plm0(new double[nmult]);
  for (int l = 0; l != lmax+1; ++l)
    for (int m = 0; m <= 2 * l; ++m) {
      const int b = m-l;
      const double sign = (b >= 0) ? 1.0 : (1-((b&1)<<1));
      plm0[l*l+m] = sign * plm.compute(l, abs(b), ctheta);
    }

  unique_ptr<double[]> invfac(new double[2*lmax+1]);
  fill_n(invfac.get(), 2*lmax+1, 1.0);
  for (int i = 1; i <= 2*lmax; ++i)
    for (int j = i; j <= 2*lmax; ++j)
      invfac[j] /= i;

  ZMatrix lmjk(nmult, nmult);
  for (int l = 0; l <= lmax; ++l) {
    for (int j = 0; j <= lmax; ++j) {
      const int a = l - j;
      const double rr = pow(r, a);
      if (a < 0) continue;
      for (int m = 0; m <= 2 * l; ++m) {
        const int lo = max(-a, m-l-j);
        const int hi = min( a, m-l+j);
        for (int b = lo; b <= hi; ++b) {
          const int k = m - l - b + j;
          const double prefactor = rr * plm0[a*a+a+b] * invfac[a+abs(b)];
          const complex<double> Oab = polar(prefactor, -b*phi);
          lmjk.element(l*l+m, j*j+k) = Oab;
        }
      }
    }
  }
  zgemm3m_("N", "C", olm_size_block, nmult, nmult, 1.0, oa->data(), olm_size_block, lmjk.data(), nmult, 0.0, ob->data(), olm_size_block);

  return ob;
}


// L2L for X
shared_ptr<const ZMatrix> Box::shift_localLX(const int lmax, const shared_ptr<const ZMatrix> mr, array<double, 3> rb) const {

  const double r = sqrt(rb[0]*rb[0] + rb[1]*rb[1] + rb[2]*rb[2]);
  const double ctheta = (r > numerical_zero__) ? rb[2]/r : 0.0;
  const double phi = atan2(rb[1], rb[0]);
  const int nmult = (lmax+1)*(lmax+1);
  const int olm_size_block = mr->ndim();

  shared_ptr<ZMatrix> mrb  = mr->clone();
  unique_ptr<double[]> plm0(new double[nmult]);
  for (int l = 0; l != lmax+1; ++l)
    for (int m = 0; m <= 2 * l; ++m) {
      const int b = m-l;
      const double sign = (b >= 0) ? 1.0 : (1-((b&1)<<1));
      plm0[l*l+m] = sign * plm.compute(l, abs(b), ctheta);
    }

  unique_ptr<double[]> invfac(new double[2*lmax+1]);
  fill_n(invfac.get(), 2*lmax+1, 1.0);
  for (int i = 1; i <= 2*lmax; ++i)
    for (int j = i; j <= 2*lmax; ++j)
      invfac[j] /= i;

  ZMatrix lmjk(nmult, nmult);
  for (int l = 0; l <= lmax; ++l) {
    for (int j = 0; j <= lmax; ++j) {
      const int a = j - l;
      const double rr = pow(r, a);
      if (a < 0) continue;
      for (int m = 0; m <= 2 * l; ++m) {
        const int lo = max(-a, l-m-j);
        const int hi = min( a, l-m+j);
        for (int b = lo; b <= hi; ++b) {
          const int k = m - l + b + j;
          const double prefactor = rr * plm0[a*a+a+b] * invfac[a+abs(b)];
          const complex<double> Oab = polar(prefactor, -b*phi);
          lmjk.element(l*l+m, j*j+k) = Oab;        
        }
      }
    }
  }

  zgemm3m_("N", "C", olm_size_block, nmult, nmult, 1.0, mr->data(), olm_size_block, lmjk.data(), nmult, 0.0, mrb->data(), olm_size_block);
  return mrb;
}


// M2L for X
shared_ptr<const ZMatrix> Box::shift_localMX(const int lmax, shared_ptr<const ZMatrix> olm, array<double, 3> r12) const {

  const double r = sqrt(r12[0]*r12[0] + r12[1]*r12[1] + r12[2]*r12[2]);
  const double ctheta = (r > numerical_zero__) ? r12[2]/r : 0.0;
  const double phi = atan2(r12[1], r12[0]);
  const int nmult = (lmax+1)*(lmax+1);
  const int olm_size_block = olm->ndim();

  shared_ptr<ZMatrix> mb = olm->clone();
  unique_ptr<double[]> plm0(new double[(2*lmax+1)*(2*lmax+1)]);
  for (int l = 0; l != 2*lmax+1; ++l)
    for (int m = 0; m <= 2 * l; ++m) {
      const int b = m-l;
      const double sign = (b >= 0) ? 1.0 : (1-((b&1)<<1));
      plm0[l*l+m] = sign * plm.compute(l, abs(b), ctheta);
    }

  ZMatrix lmjk(nmult, nmult);
  for (int l = 0; l <= lmax; ++l) {
    const double phase_l = (1-((l&1)<<1));
    for (int j = 0; j <= lmax; ++j) {
      const int a = l + j;
      const double rr = 1.0 / pow(r, a + 1);
      for (int m = 0; m <= 2 * l; ++m) {
        for (int k = 0; k <= 2 * j; ++k) {
          const int b = m - l + k - j;
          double prefactor = plm0[a*a+a+b] * rr * f(a-abs(b));
          const complex<double> Mab = phase_l * polar(prefactor, b*phi);
          lmjk.element(l*l+m, j*j+k) = Mab;        
        }
      }
    }
  }

  zgemm3m_("N", "C", olm_size_block, nmult, nmult, 1.0, olm->data(), olm_size_block, lmjk.data(), nmult, 0.0, mb->data(), olm_size_block);
  return mb;
}


double Box::coulomb_ff() const {
 
  const complex<double> en = blas::dot_product(multipole_->data(), nmult_, localJ_->data());
  assert(abs(en.imag()) < 1e-10);
  return en.real();
}


shared_ptr<const Matrix> Box::compute_Fock_nf_J(shared_ptr<const Matrix> density, shared_ptr<const VectorB> max_den) const {

  assert(nchild() == 0);
  auto out = make_shared<Matrix>(density->ndim(), density->mdim());
  out->zero();

  // NF: 4c integrals
  const int shift = sizeof(int) * 4;
  const int nsh = sqrt(max_den->size());
  assert (nsh*nsh == max_den->size());

  int ntask = 0;
  for (auto& v01 : sp_) {
    if (v01->schwarz() < thresh_) continue;
    const int i0 = v01->shell_ind(1);

    const int i1 = v01->shell_ind(0);
    if (i1 < i0) continue;

    const int i01 = i0 * nsh + i1;
    const double density_01 = (*max_den)(i01) * 4.0;

    for (auto& neigh : neigh_) {
      for (auto& v23 : neigh->sp()) {
        if (v23->schwarz() < thresh_) continue;
        const int i2 = v23->shell_ind(1);
        if (i2 < i0) continue;
        const int i3 = v23->shell_ind(0);
        if (i3 < i2) continue;
        const int i23 = i2 * nsh + i3;
        if (i23 < i01) continue;

        const double density_23 = (*max_den)(i23) * 4.0;

        const double mulfactor = max(density_01, density_23);
        const double integral_bound = mulfactor * v01->schwarz() * v23->schwarz();
        const bool skip_schwarz = integral_bound < schwarz_thresh_;
        if (skip_schwarz) continue;
        ++ntask;
      }
    }
  }

  TaskQueue<function<void(void)>> tasks(ntask);
  mutex jmutex;

  for (auto& v01 : sp_) {
    if (v01->schwarz() < thresh_) continue;
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
    const double density_01 = (*max_den)(i01) * 4.0;

    for (auto& neigh : neigh_) {
      for (auto& v23 : neigh->sp()) {
        if (v23->schwarz() < thresh_) continue;
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

        const double density_23 = (*max_den)(i23) * 4.0;

        const double mulfactor = max(density_01, density_23);
        const double integral_bound = mulfactor * v01->schwarz() * v23->schwarz();
        const bool skip_schwarz = integral_bound < schwarz_thresh_;
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


shared_ptr<const Matrix> Box::compute_Fock_nf_K(shared_ptr<const Matrix> density, shared_ptr<const VectorB> max_den) const {

  assert(nchild() == 0);
  auto out = make_shared<Matrix>(density->ndim(), density->mdim());
  out->zero();

  // NF: 4c integrals
  const int shift = sizeof(int) * 4;
  const int nsh = sqrt(max_den->size());
  assert (nsh*nsh == max_den->size());

  int ntask = 0;
  for (auto& v01 : sp_) {
    if (v01->schwarz() < thresh_) continue;
    const int i0 = v01->shell_ind(1);

    const int i1 = v01->shell_ind(0);
    if (i1 < i0) continue;

    const int i01 = i0 * nsh + i1;
    for (auto& neigh : neigh_) {
      for (auto& v23 : neigh->sp()) {
        if (v23->schwarz() < thresh_) continue;
        const int i2 = v23->shell_ind(1);
        if (i2 < i0) continue;
        const int i3 = v23->shell_ind(0);
        if (i3 < i2) continue;
        const int i23 = i2 * nsh + i3;
        if (i23 < i01) continue;

        const double density_02 = (*max_den)(i0 * nsh + i2);
        const double density_03 = (*max_den)(i0 * nsh + i3);
        const double density_12 = (*max_den)(i1 * nsh + i2);
        const double density_13 = (*max_den)(i1 * nsh + i3);

        const double mulfactor = max(max(density_13, density_02),
                                     max(density_03, density_12));
        const double integral_bound = mulfactor * v01->schwarz() * v23->schwarz();
        const bool skip_schwarz = integral_bound < schwarz_thresh_;
        if (skip_schwarz) continue;
        ++ntask;
      }
    }
  }

  TaskQueue<function<void(void)>> tasks(ntask);
  mutex jmutex;

  for (auto& v01 : sp_) {
    if (v01->schwarz() < thresh_) continue;
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
    for (auto& neigh : neigh_) {
      for (auto& v23 : neigh->sp()) {
        if (v23->schwarz() < thresh_) continue;
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

        const double density_02 = (*max_den)(i0 * nsh + i2);
        const double density_03 = (*max_den)(i0 * nsh + i3);
        const double density_12 = (*max_den)(i1 * nsh + i2);
        const double density_13 = (*max_den)(i1 * nsh + i3);

        const double mulfactor = max(max(density_13, density_02),
                                     max(density_03, density_12));
        const double integral_bound = mulfactor * v01->schwarz() * v23->schwarz();
        const bool skip_schwarz = integral_bound < schwarz_thresh_;
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
                  for (int j3 = b3offset; j3 != b3offset + b3size; ++j3, ++eridata) {
                    if (j3 < j2) continue;
                    const unsigned int nj23 = (j2 << shift) + j3;
                    if (nj23 < nj01 && i01 == i23) continue;
                    const double scale23 = (j2 == j3) ? 0.5 : 1.0;
                    const double scale = (nj01 == nj23) ? 0.25 : 0.5;

                    const double eri = *eridata;
                    const double intval = eri * scale * scale01 * scale23;
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
