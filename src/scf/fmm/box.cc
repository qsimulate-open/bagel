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


#include <mutex>
#include <src/scf/fmm/box.h>
#include <src/integral/os/multipolebatch.h>
#include <src/integral/rys/eribatch.h>
#include <src/util/taskqueue.h>
#include <src/util/math/factorial.h>
#include <src/util/math/legendre.h>

using namespace bagel;
using namespace std;

const static Legendre plm;
const static Factorial f;

void Box::init() {

  nchild_ = child_.size();
  centre_ = {{0, 0, 0}};
  extent_ = 0.0;
  nsp_ = 0;
  for (auto& i : sp_) {
    if (i->schwarz() < schwarz_thresh_) continue;
    ++nsp_;
  }

  if (nchild_ == 0) {
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
    for (int n = 0; n != nchild_; ++n) {
      shared_ptr<const Box> i = child_[n].lock();
      for (int j = 0; j != 3; ++j)  centre_[j] += i->centre(j);
    }
    centre_[0] /= nchild_;
    centre_[1] /= nchild_;
    centre_[2] /= nchild_;
    for (int n = 0; n != nchild_; ++n) {
      shared_ptr<const Box> i = child_[n].lock();
      double rad = 0.0;
      for (int j = 0; j != 3; ++j) rad += pow(i->centre(j)-centre_[j], 2.0);
      const double ei = sqrt(rad) + i->extent();
      assert(ei > 0);
      extent_ = max(extent_, ei);
    }
  }

  ws_ = static_cast<int>(ceil(extent_/boxsize_));
  nmult_ = (lmax_ + 1) * (lmax_ + 1);
  olm_ = make_shared<ZVectorB>(nmult_);
  mlm_ = make_shared<ZVectorB>(nmult_);

  nshell0_ = 0;
  int cnt = 0;
  for (auto& i : sp_) {
    auto it = shell0_.insert(make_pair(i->shell(0), make_tuple(nshell0_, i->offset(0), cnt)));
    if (it.second) {
      nshell0_ += i->shell(0)->nbasis();
      ++cnt;
    }
  }
}


void Box::insert_sp(const vector<shared_ptr<const ShellPair>>& sp) {

  const int nsp = sp_.size();
  sp_.resize(nsp + sp.size());
  for (int i = 0; i != sp.size(); ++i)
    sp_[i + nsp] = sp[i];

}


void Box::insert_child(weak_ptr<const Box> c) {

  const int nchild = child_.size();
  shared_ptr<const Box> child = c.lock();
  if (child) {
    child_.resize(nchild + 1);
    child_[nchild] = child;
  }
}


void Box::insert_parent(weak_ptr<const Box> p) {

  parent_ = p;
}


void Box::get_neigh(const vector<weak_ptr<Box>>& box, const double ws) {

  shared_ptr<const Box> box0 = box[0].lock();
  if ((box.empty()) || (box0->rank() != rank_)) return;

  neigh_.resize(box.size());
  nonneigh_.resize(box.size());
  int nn = 0;
  int nnn = 0;
  for (int i = 0; i != box.size(); ++i) {
    shared_ptr<const Box> b = box[i].lock();
    if (is_neigh(b, ws)) {
      neigh_[nn] = b;
      ++nn;
    } else {
      nonneigh_[nnn] = b;
      ++nnn;
    }
  }
  neigh_.resize(nn);
  nonneigh_.resize(nnn);
  nneigh_ = neigh_.size();
}


bool Box::is_neigh(weak_ptr<const Box> b, const double ws) const {

  shared_ptr<const Box> box = b.lock();
  if (nchild_ == 0) {
    double rr = 0;
    for (int i = 0; i != 3; ++i) rr += pow(centre_[i] - box->centre(i), 2);
    return (sqrt(rr) <= (1.0+ws)*(extent_ + box->extent()));
  } else {
    for (int n0 = 0; n0 != nchild_; ++n0) {
      shared_ptr<const Box> i = child_[n0].lock();
      for (int n1 = 0; n1 != box->child().size(); ++n1) {
        shared_ptr<const Box> j = box->child()[n1].lock();
        if (i->is_neigh(j, ws)) return true;
      }
    }
    return false;
  }
}

void Box::get_inter(const vector<weak_ptr<Box>>& box, const double ws) {

  shared_ptr<const Box> box0 = box[0].lock();
  if ((box.empty()) || (box0->rank() != rank_)) return;

  shared_ptr<const Box> parent = parent_.lock();
  inter_.resize(box.size());
  int ni = 0;
  for (int i = 0; i != box.size(); ++i) {
    shared_ptr<const Box> b = box[i].lock();
    if ((!is_neigh(b, ws)) && (!parent || parent->is_neigh(b->parent(), ws))) {
      inter_[ni] = b;
      ++ni;
    }
  }
  inter_.resize(ni);
  ninter_ = inter_.size();
}


void Box::sort_sp() {

  if (nchild_ == 0) return;

  assert(nchild_ > 0);
  vector<shared_ptr<const ShellPair>> newsp;
  for (int n = 0; n != nchild_; ++n) {
    shared_ptr<const Box> c = child_[n].lock();
    const vector<shared_ptr<const ShellPair>>& csp = c->sp();
    newsp.insert(newsp.end(), csp.begin(), csp.end());
  }
  assert(newsp.size() == sp_.size());
  sp_ = newsp;
}


void Box::compute_M2M(shared_ptr<const Matrix> density) {

  olm_->fill(0.0);

  if (nchild_ == 0) { // leaf
    TaskQueue<function<void(void)>> tasks(nsp_);
    mutex jmutex;
    for (auto& v : sp_) {
      if (v->schwarz() < schwarz_thresh_) continue;

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
          blas::ax_plus_y_n(1.0, olm->data(), nmult_, olm_->data());
        }
      );
    }
    tasks.compute();
  } else { // shift children's multipoles
    for (int n = 0; n != nchild_; ++n) {
      shared_ptr<const Box> c = child_[n].lock();
      const array<double, 3> r12 = {{c->centre(0) - centre_[0], c->centre(1) - centre_[1], c->centre(2) - centre_[2]}};
      auto cmultipole = make_shared<ZMatrix>(1, nmult_);
      copy_n(c->olm()->data(), nmult_, cmultipole->data());
      shared_ptr<const ZMatrix> tmp = shift_multipolesX(lmax_, cmultipole, r12);
      auto smoment = make_shared<ZVectorB>(nmult_);
      copy_n(tmp->data(), nmult_, smoment->data());
      blas::ax_plus_y_n(1.0, smoment->data(), nmult_, olm_->data());
    }
  }
}


//|su) C_ui C_sj
void Box::compute_M2M_X(shared_ptr<const Matrix> ocoeff_sj, shared_ptr<const Matrix> ocoeff_ui) {

  const int nmult_k = (lmax_k_+1)*(lmax_k_+1);
  olm_ndim_ = ocoeff_sj->mdim();  olm_mdim_ = ocoeff_ui->mdim();
  olm_ji_ = make_shared<ZMatrix>(olm_ndim_*olm_mdim_, nmult_k);

  if (nchild_ == 0) { // leaf
    ZMatrix interm(nshell0_, ocoeff_ui->mdim()*nmult_k, true);
    TaskQueue<function<void(void)>> tasks(nsp_);
    vector<mutex> mmutex(shell0_.size());

    for (auto& v : sp_) {
      if (v->schwarz() < schwarz_thresh_) continue;

      tasks.emplace_back(
        [this, v, ocoeff_ui, &mmutex, nmult_k, &interm]() {
          MultipoleBatch olm_su(v->shells(), centre_, lmax_k_);
          olm_su.compute();
          const int dimb0 = v->shell(0)->nbasis();
          const int dimb1 = v->shell(1)->nbasis();
          const ZMatrix zui(*ocoeff_ui->cut(v->offset(1), v->offset(1)+dimb1), 1.0);
          unique_ptr<complex<double>[]> tmp(new complex<double>[dimb0*ocoeff_ui->mdim()*nmult_k]);
          for (int k = 0; k != nmult_k; ++k)
            zgemm_("N", "N", dimb0, ocoeff_ui->mdim(), dimb1, 1.0, olm_su.data() + olm_su.size_block()*k, dimb0, zui.data(), dimb1,
                                                              0.0, tmp.get()+dimb0*ocoeff_ui->mdim()*k, dimb0);
          lock_guard<mutex> lock(mmutex[get<2>(shell0_.at(v->shell(0)))]);
          const int off = get<0>(shell0_.at(v->shell(0)));
          interm.add_block(1.0, off, 0, dimb0, ocoeff_ui->mdim()*nmult_k, tmp.get());
        }
      );
    }
    tasks.compute();
    blas::conj_n(interm.data(), interm.size());

    ZMatrix icoeff(nshell0_, ocoeff_sj->mdim(), true);
    for (auto& i : shell0_) {
      const int local_offset = get<0>(i.second);
      const int all_offset = get<1>(i.second);
      const int size = i.first->nbasis();
      for (int m = 0; m != ocoeff_sj->mdim(); ++m)
        for (int n = 0; n != size; ++n)
          icoeff(n+local_offset, m) = ocoeff_sj->element(n+all_offset, m);
    }
    zgemm_("C", "N", icoeff.mdim(), interm.mdim(), nshell0_, 1.0, icoeff.data(), nshell0_, interm.data(), nshell0_, 0.0, olm_ji_->data(), icoeff.mdim());

  } else { // shift children's multipoles
    for (int n = 0; n != nchild_; ++n) {
      shared_ptr<const Box> c = child_[n].lock();
      const array<double, 3> r12 = {{c->centre(0) - centre_[0], c->centre(1) - centre_[1], c->centre(2) - centre_[2]}};
      shared_ptr<const ZMatrix> smoment = shift_multipolesX(lmax_k_, c->olm_ji(), r12);
      blas::ax_plus_y_n(1.0, smoment->data(), olm_ji_->size(), olm_ji_->data());
    }
  }
}


void Box::compute_L2L_X() {

  // from parent
  shared_ptr<const Box> parent = parent_.lock();
  if (parent) {
    const array<double, 3> r12 = {{centre_[0] - parent->centre(0), centre_[1] - parent->centre(1), centre_[2] - parent->centre(2)}};
    shared_ptr<const ZMatrix> slocal = shift_localLX(lmax_k_, parent->mlm_ji(), r12);
    blas::ax_plus_y_n(1.0, slocal->data(), mlm_ji_->size(), mlm_ji_->data());
  }
}


void Box::compute_L2L() {

  // from parent
  shared_ptr<const Box> parent = parent_.lock();
  if (parent) {
    const array<double, 3> r12 = {{centre_[0] - parent->centre(0), centre_[1] - parent->centre(1), centre_[2] - parent->centre(2)}};
    auto plocalJ = make_shared<ZMatrix>(1, nmult_);
    copy_n(parent->mlm()->data(), nmult_, plocalJ->data());
    shared_ptr<const ZMatrix> tmp = shift_localLX(lmax_, plocalJ, r12);
    auto slocal = make_shared<ZVectorB>(nmult_);
    copy_n(tmp->data(), nmult_, slocal->data());
    blas::ax_plus_y_n(1.0, slocal->data(), nmult_, mlm_->data());
  }
}


void Box::compute_M2L_X() {

  mlm_ji_ = olm_ji_->clone();

  // from interaction list
  for (int i = 0; i != ninter_; ++i) {
    shared_ptr<const Box> it = inter_[i].lock();
    const array<double, 3> r12 = {{centre_[0] - it->centre(0), centre_[1] - it->centre(1), centre_[2] - it->centre(2)}};
    shared_ptr<const ZMatrix> slocal = shift_localMX(lmax_k_, it->olm_ji(), r12);
    blas::ax_plus_y_n(1.0, slocal->data(), mlm_ji_->size(), mlm_ji_->data());
  }
}


void Box::compute_M2L() {

  mlm_->fill(0.0);
  // from interaction list
  for (int i = 0; i != ninter_; ++i) {
    shared_ptr<const Box> it = inter_[i].lock();
    const array<double, 3> r12 = {{centre_[0] - it->centre(0), centre_[1] - it->centre(1), centre_[2] - it->centre(2)}};
    auto imultipole = make_shared<ZMatrix>(1, nmult_);
    copy_n(it->olm()->data(), nmult_, imultipole->data());
    shared_ptr<const ZMatrix> tmp = shift_localMX(lmax_, imultipole, r12);
    auto slocal = make_shared<ZVectorB>(nmult_);
    copy_n(tmp->data(), nmult_, slocal->data());
    blas::ax_plus_y_n(1.0, slocal->data(), nmult_, mlm_->data());
  }
}


shared_ptr<const Matrix> Box::compute_Fock_nf(shared_ptr<const Matrix> density, shared_ptr<const VectorB> max_den) const {

  assert(nchild_ == 0);

  auto out = make_shared<Matrix>(density->ndim(), density->mdim());
  out->zero();

  // NF: 4c integrals
  const int shift = sizeof(int) * 4;
  const int nsh = sqrt(max_den->size());
  assert (nsh*nsh == max_den->size());

  int ntask = 0;
  for (auto& v01 : sp_) {
    if (v01->schwarz() < schwarz_thresh_) continue;
    const int i0 = v01->shell_ind(1);

    const int i1 = v01->shell_ind(0);
    if (i1 < i0) continue;

    const int i01 = i0 * nsh + i1;
    const double density_01 = (*max_den)(i01) * 4.0;

    for (int n = 0; n != nneigh_; ++n) {
      shared_ptr<const Box> neigh = neigh_[n].lock();
      for (auto& v23 : neigh->sp()) {
        if (v23->schwarz() < schwarz_thresh_) continue;
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
    if (v01->schwarz() < schwarz_thresh_) continue;
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

    for (int n = 0; n != nneigh_; ++n) {
      shared_ptr<const Box> neigh = neigh_[n].lock();
      for (auto& v23 : neigh->sp()) {
        if (v23->schwarz() < schwarz_thresh_) continue;
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
  cout << "Box " << boxid_ << " Rank = " << rank_ << " *** nchild = " << nchild_ << " *** nsp = " << sp_.size()
       << " *** nneigh = " << nneigh_ << " *** ninter = " << ninter_ << " *** extent = " << extent_
       << " *** centre = " << setprecision(3) << centre_[0] << "  " << centre_[1] << "  " << centre_[2] << endl;

//  cout << " *** Shell pairs at ***" << endl;
//  for (int i = 0; i != sp_.size(); ++i)
//    cout << setprecision(5) << sp(i)->centre(0) << "  " << sp(i)->centre(1) << "  " << sp(i)->centre(2) << endl;
}


shared_ptr<const Matrix> Box::compute_Fock_ff(shared_ptr<const Matrix> density) const {

  assert(nchild_ == 0);
  auto out = make_shared<ZMatrix>(density->ndim(), density->mdim());

  // compute multipoles again - not efficient - for now
  TaskQueue<function<void(void)>> tasks(nsp_);
  mutex jmutex;
  for (auto& v : sp_) {
    if (v->schwarz() < schwarz_thresh_) continue;

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
              out->element(i, j) += conj(*dat[k]++) * (*mlm_)(k);
      }
    );
  }
  tasks.compute();
  assert(out->get_imag_part()->rms() < 1e-10);

  return out->get_real_part();
}


shared_ptr<const Matrix> Box::compute_Fock_ff_K(shared_ptr<const Matrix> ocoeff_ti) const {

  assert(nchild_ == 0);
  auto out = make_shared<Matrix>(ocoeff_ti->ndim(), olm_ndim_);
  assert(ocoeff_ti->mdim() == olm_mdim_);
  const int nmult_k = (lmax_k_+1)*(lmax_k_+1);
  TaskQueue<function<void(void)>> tasks(nsp_);

  ZMatrix interm(nshell0_, ocoeff_ti->mdim()*nmult_k, true);
  vector<mutex> mmutex(shell0_.size());

  for (auto& v : sp_) {
    if (v->schwarz() < schwarz_thresh_) continue;

    tasks.emplace_back(
      [this, v, ocoeff_ti, &mmutex, nmult_k, &interm]() {
        MultipoleBatch olm_su(v->shells(), centre_, lmax_k_);
        olm_su.compute();
        const int dimb0 = v->shell(0)->nbasis();
        const int dimb1 = v->shell(1)->nbasis();
        const ZMatrix zti(*ocoeff_ti->cut(v->offset(1), v->offset(1)+dimb1), 1.0);
        unique_ptr<complex<double>[]> tmp(new complex<double>[dimb0*ocoeff_ti->mdim()*nmult_k]);
        for (int k = 0; k != nmult_k; ++k)
          zgemm_("N", "N", dimb0, ocoeff_ti->mdim(), dimb1, 1.0, olm_su.data() + olm_su.size_block()*k, dimb0, zti.data(), dimb1,
                                                            0.0, tmp.get()+dimb0*ocoeff_ti->mdim()*k, dimb0);
        lock_guard<mutex> lock(mmutex[get<2>(shell0_.at(v->shell(0)))]);
        const int off = get<0>(shell0_.at(v->shell(0)));
        interm.add_block(1.0, off, 0, dimb0, ocoeff_ti->mdim()*nmult_k, tmp.get());
      }
    );
  }
  tasks.compute();

  ZMatrix kjr(olm_ndim_, nshell0_);
  zgemm_("N", "C", olm_ndim_, nshell0_, olm_mdim_*nmult_k, 1.0, mlm_ji_->data(), olm_ndim_, interm.data(), nshell0_, 0.0, kjr.data(), olm_ndim_);

  auto krj = kjr.get_real_part()->transpose();

  for (auto& i : shell0_) {
    const int dimb0 = i.first->nbasis();
    const int local_offset = get<0>(i.second);
    const int all_offset = get<1>(i.second);
    for (int i = 0; i != olm_ndim_; ++i)
      blas::ax_plus_y_n(1.0, krj->element_ptr(local_offset, i), dimb0, out->element_ptr(all_offset, i));
  }

  return out;
}


shared_ptr<const Matrix> Box::compute_exact_ff(shared_ptr<const Matrix> density) const { //for debug

  auto jff = make_shared<Matrix>(density->ndim(), density->mdim());
  jff->zero();

  int ntask = 0;
  for (auto& v01 : sp_) {
    if (v01->schwarz() < schwarz_thresh_) continue;
    for (int n = 0; n != nonneigh_.size(); ++n) {
      shared_ptr<const Box> nonneigh = neigh_[n].lock();
      for (auto& v23 : nonneigh->sp()) {
        if (v23->schwarz() < schwarz_thresh_) continue;
        ++ntask;
      }
    }
  }

  TaskQueue<function<void(void)>> tasks(ntask);
  mutex jmutex;

  for (auto& v01 : sp_) {
    if (v01->schwarz() < schwarz_thresh_) continue;
    shared_ptr<const Shell> b0 = v01->shell(1);
    const int b0offset = v01->offset(1);
    const int b0size = b0->nbasis();

    shared_ptr<const Shell> b1 = v01->shell(0);
    const int b1offset = v01->offset(0);
    const int b1size = b1->nbasis();

    for (int n = 0; n != nonneigh_.size(); ++n) {
      shared_ptr<const Box> nonneigh = neigh_[n].lock();
      for (auto& v23 : nonneigh->sp()) {
        if (v23->schwarz() < schwarz_thresh_) continue;
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

  return jff;
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
  zgemm_("N", "T", olm_size_block, nmult, nmult, 1.0, oa->data(), olm_size_block, lmjk.data(), nmult, 0.0, ob->data(), olm_size_block);

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

  zgemm_("N", "T", olm_size_block, nmult, nmult, 1.0, mr->data(), olm_size_block, lmjk.data(), nmult, 0.0, mrb->data(), olm_size_block);
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

  zgemm_("N", "T", olm_size_block, nmult, nmult, 1.0, olm->data(), olm_size_block, lmjk.data(), nmult, 0.0, mb->data(), olm_size_block);
  return mb;
}


shared_ptr<const Matrix> Box::compute_Fock_nf_J(shared_ptr<const Matrix> density, shared_ptr<const VectorB> max_den) const {

  assert(nchild_ == 0);
  auto out = make_shared<Matrix>(density->ndim(), density->mdim());
  out->zero();

  // NF: 4c integrals
  const int shift = sizeof(int) * 4;
  const int nsh = sqrt(max_den->size());
  assert (nsh*nsh == max_den->size());

  int ntask = 0;
  for (auto& v01 : sp_) {
    if (v01->schwarz() < schwarz_thresh_) continue;
    const int i0 = v01->shell_ind(1);

    const int i1 = v01->shell_ind(0);
    if (i1 < i0) continue;

    const int i01 = i0 * nsh + i1;
    const double density_01 = (*max_den)(i01) * 4.0;

    for (int n = 0; n != nneigh_; ++n) {
      shared_ptr<const Box> neigh = neigh_[n].lock();
      for (auto& v23 : neigh->sp()) {
        if (v23->schwarz() < schwarz_thresh_) continue;
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
    if (v01->schwarz() < schwarz_thresh_) continue;
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

    for (int n = 0; n != nneigh_; ++n) {
      shared_ptr<const Box> neigh = neigh_[n].lock();
      for (auto& v23 : neigh->sp()) {
        if (v23->schwarz() < schwarz_thresh_) continue;
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

  assert(nchild_ == 0);
  auto out = make_shared<Matrix>(density->ndim(), density->mdim());
  out->zero();

  // NF: 4c integrals
  const int shift = sizeof(int) * 4;
  const int nsh = sqrt(max_den->size());
  assert (nsh*nsh == max_den->size());

  int ntask = 0;
  for (auto& v01 : sp_) {
    if (v01->schwarz() < schwarz_thresh_) continue;
    const int i0 = v01->shell_ind(1);

    const int i1 = v01->shell_ind(0);
    if (i1 < i0) continue;

    const int i01 = i0 * nsh + i1;
    for (int n = 0; n != nneigh_; ++n) {
      shared_ptr<const Box> neigh = neigh_[n].lock();
      for (auto& v23 : neigh->sp()) {
        if (v23->schwarz() < schwarz_thresh_) continue;
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
    if (v01->schwarz() < schwarz_thresh_) continue;
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
    for (int n = 0; n != nneigh_; ++n) {
      shared_ptr<const Box> neigh = neigh_[n].lock();
      for (auto& v23 : neigh->sp()) {
        if (v23->schwarz() < schwarz_thresh_) continue;
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
