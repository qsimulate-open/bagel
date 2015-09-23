//
// BAGEL - Parallel electron correlation program.
// Filename: pfmm.cc
// Copyright (C) 2015 Toru Shiozaki
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


#include <src/periodic/pfmm.h>

using namespace std;
using namespace bagel;

const static Legendre plm;
const static Factorial f;
const static Gamma_scaled sgamma;

const static double beta__ = sqrt(pi__); // convergence parameter

PFMM::PFMM(shared_ptr<const SimulationCell> scell, const int lmax, const int ws, const int extent, const double thresh, shared_ptr<StackMem> stack)
  : scell_(scell), auxfile_(scell->geom()->auxfile()), lmax_(lmax), ws_(ws), extent_sum_(extent), thresh_(thresh) {

  if (stack == nullptr) {
    stack_ = resources__->get();
    allocated_here_ = true;
  } else {
    stack_ = stack;
    allocated_here_ = false;
  }

  assert(lmax_ <= 10);
  ndim_ = scell->ndim();
  if (ndim_ > 3 || ndim_ < 1)
    throw runtime_error("System must be periodic in 1-, 2-, or 3-D");

  num_multipoles_ = (2*lmax_ + 1) * (2*lmax_ + 1);
  mlm_.resize(num_multipoles_);
  max_rank_ = (lmax_ * 2) + 1;

  // should be provided from input
  max_height_ = 21; // tree construction 21 is absolute max
  do_contract_ = true;

  compute_mlm();
  stack_->release(size_allocated_, buff_);
  resources__->release(stack_);
}


bool PFMM::is_in_cff(array<double, 3> L) {

  const double extent = scell_->extent();

  const double rsq = L[0]*L[0] + L[1]*L[1] + L[2]*L[2];
  const bool out = (rsq > 2.0 * (1 + ws_) *  extent) ? true : false;

  return out;
}


void PFMM::compute_mlm() { // rectangular scell for now

  const double pibeta = pi__ * pi__ / (beta__ * beta__);
  const int nvec = pow(2*extent_sum_+1, ndim_);
  allocate_arrays(nvec);
  vector<array<int, 3>> vidx(nvec);

  int cnt = 0;
  if (ndim_ == 3) {
    for (int n3 = -extent_sum_; n3 <= extent_sum_; ++n3)
      for (int n2 = -extent_sum_; n2 <= extent_sum_; ++n2)
        for (int n1 = -extent_sum_; n1 <= extent_sum_; ++n1, ++cnt)
          vidx[cnt] = {{n1, n2, n3}};
  } else if (ndim_ == 2) {
    for (int n2 = -extent_sum_; n2 <= extent_sum_; ++n2)
      for (int n1 = -extent_sum_; n1 <= extent_sum_; ++n1, ++cnt)
      vidx[cnt] = {{n1, n2, 0}};
  } else if (ndim_ == 1) {
    for (int n1 = -extent_sum_; n1 <= extent_sum_; ++n1, ++cnt)
      vidx[cnt] = {{n1, 0, 0}};
  }
  assert(cnt == nvec);

  std::sort(vidx.begin(), vidx.end(), sort_vector); // sort to sum spherically

  vector<array<double, 3>> primvecs(ndim_);
  for (int i = 0; i != ndim_; ++i)
    primvecs[i] = scell_->primitive_vectors(i);

  for (int ivec = 0; ivec != nvec; ++ivec) {
    const int pos = ivec * 3;
    array<int, 3> idx = vidx[ivec];
    rvec_[pos    ] = idx[0] * primvecs[0][0] + idx[1] * primvecs[1][0] + idx[2] * primvecs[2][0];
    rvec_[pos + 1] = idx[0] * primvecs[0][1] + idx[1] * primvecs[1][1] + idx[2] * primvecs[2][1];
    rvec_[pos + 2] = idx[0] * primvecs[0][2] + idx[1] * primvecs[1][2] + idx[2] * primvecs[2][2];
    Rsq_[ivec] = rvec_[pos]*rvec_[pos] + rvec_[pos+1]*rvec_[pos+1] + rvec_[pos+2]*rvec_[pos+2];
    T_[ivec] = Rsq_[ivec] * beta__ * beta__;
  }

  for (int l = 0; l < max_rank_; ++l) {
    root_weight(l, nvec);
    const int rank = l + 1;

    for (int ivec = 1; ivec != nvec; ++ivec) {
      array<int, 3> idx = vidx[ivec];

      const int pos = ivec * 3;
      const double rsq = Rsq_[ivec];
      const double r = sqrt(rsq);
      const double ctheta = (rsq > numerical_zero__) ? rvec_[pos+2]/r : 0.0;
      const double phi = atan2(rvec_[pos+1], rvec_[pos]);
      const double* croots = roots_ + ivec * rank;
      const double* cweights = weights_ + ivec * rank;

      double glower = 0.0;
      if (l == 0) {
        for (int i = 0; i != rank; ++i)
          glower += cweights[i];
      } else {
        for (int i = 0; i != rank; ++i)
          glower += cweights[i] * pow(croots[i], l);
      }

      const double coeff = 2.0 * pow(beta__, 2*l+1) * sgamma(l, r);
      const double gupper = 1.0 / pow(r, l+1.0) - glower * coeff;

      for (int mm = 0; mm <= 2 * l; ++mm) {
        const int m = mm - l;
        const int am = abs(m);
        const int imul = l * l + mm;

        double plm_tilde = plm.compute(l, am, ctheta);
        double ft = 1.0;
        for (int i = 1; i <= l - am; ++i) {
          plm_tilde *= ft;
          ft += 1.0;
        }
        const double sign = (m >=0) ? (cos(am * phi)) : (-1.0 * cos(am * phi));

        if (abs(idx[0]) > ws_ && abs(idx[1]) > ws_ && abs(idx[2]) > ws_) {
          // real term
          const double real = gupper * sign * plm_tilde;
          const double imag = gupper * sin(am * phi) * plm_tilde;
          mlm_[imul] += complex<double>(real, imag);
        }
        if (abs(idx[0]) <= ws_ && abs(idx[1]) <= ws_ && abs(idx[2]) <= ws_) {
          // substract smooth part within ws_
          const double real = coeff * glower * sign * plm_tilde;
          const double imag = coeff * glower * sin(am * phi) * plm_tilde;
          mlm_[imul] -= complex<double>(real, imag);
        }
      }
    }
  }

  vector<array<double, 3>> primkvecs(ndim_);
  if (ndim_ == 1) {
    const double a1sq = sqrt(primvecs[0][0]*primvecs[0][1] + primvecs[0][1]*primvecs[0][1] + primvecs[0][2]*primvecs[0][2]);
    for (int i = 0; i != 3; ++i)
      primkvecs[0][i] = primvecs[0][i] / a1sq;
  } else if (ndim_ == 2) {
    array<double, 3> a12 = cross(primvecs[0], primvecs[1]);
    const double scale = 1.0 / sqrt(dot(a12, a12));
    primkvecs[0] = cross(primvecs[1], a12, scale);
    primkvecs[1] = cross(a12, primvecs[0], scale);
  } else {
    array<double, 3> a23 = cross(primvecs[1], primvecs[2]);
    const double scale = 1.0 / dot(primvecs[0], a23);
    for (int i = 0; i != 3; ++i) {
      const int j = (i+1 < ndim_) ? i+1 : i+1-ndim_;
      const int k = (i+2 < ndim_) ? i+2 : i+2-ndim_;
      primkvecs[i] = cross(primvecs[j], primvecs[k], scale);
    }
  }

  for (int ivec = 0; ivec != nvec; ++ivec) {
    const int pos = ivec * 3;
    array<int, 3> idx = vidx[ivec];
    kvec_[pos    ] = idx[0] * primkvecs[0][0] + idx[1] * primkvecs[1][0] + idx[2] * primkvecs[2][0];
    kvec_[pos + 1] = idx[0] * primkvecs[0][1] + idx[1] * primkvecs[1][1] + idx[2] * primkvecs[2][1];
    kvec_[pos + 2] = idx[0] * primkvecs[0][2] + idx[1] * primkvecs[1][2] + idx[2] * primkvecs[2][2];
    Rsq_[ivec] = kvec_[pos]*kvec_[pos] + kvec_[pos+1]*kvec_[pos+1] + kvec_[pos+2]*kvec_[pos+2];
    T_[ivec] = Rsq_[ivec] * beta__ * beta__;
  }

  fill_n(roots_, max_rank_ * nvec, 0.0);
  fill_n(weights_, max_rank_ * nvec, 0.0);
  for (int l = 0; l < max_rank_; ++l) {
    const complex<double> coeffl = std::pow(complex<double>(0.0, 1.0), l) * pow(pi__, l-0.5);
    root_weight(l, nvec);
    const int rank = l + 1;


    for (int ivec = 1; ivec != nvec; ++ivec) {
      const int pos = ivec * 3;
      const double rsq = Rsq_[ivec];
      const double r = sqrt(rsq);
      const double ctheta = (rsq > numerical_zero__) ? kvec_[pos+2]/r : 0.0;
      const double phi = atan2(kvec_[pos+1], kvec_[pos]);
      const double* croots = roots_ + ivec * rank;
      const double* cweights = weights_ + ivec * rank;

      double glower = 0.0;
      if (l == 0) {
        for (int i = 0; i != rank; ++i)
          glower += cweights[i];
      } else {
        for (int i = 0; i != rank; ++i)
          glower += cweights[i] * pow(croots[i], l);
      }

      for (int mm = 0; mm <= 2 * l; ++mm) {
        const int m = mm - l;
        const int am = abs(m);
        const int imul = l * l + mm;

        double plm_tilde = plm.compute(l, am, ctheta);
        double ft = 1.0;
        for (int i = 1; i <= l - am; ++i) {
          plm_tilde *= ft;
          ft += 1.0;
        }

        const double sign = (m >=0) ? (cos(am * phi)) : (-1.0 * cos(am * phi));

        // smooth term
        const double coeffm = plm_tilde * sgamma(l, r) * exp(-rsq * pibeta) / rsq;
        double real = coeffm * sign;
        double imag = coeffm * sin(am * phi);
        mlm_[imul] += coeffl * complex<double>(real, imag);
      }
    }
  }
}


void PFMM::root_weight(const int l, const int size) {
  if (l == 0) {
    for (int i = 0; i != size; ++i) {
      if (abs(T_[i]) < thresh_) {
        weights_[i] = 1.0;
      } else {
        const double sqrtt = sqrt(T_[i]);
        const double erfsqt = erf(sqrtt);
        weights_[i] = erfsqt * sqrt(pi__) * 0.5 / sqrtt;
      }
    }
  } else {
    const int rank = l + 1;
    eriroot__.root(rank, T_, roots_, weights_, size);
  }
}


void PFMM::allocate_arrays(const size_t ps) {

  size_allocated_ = (max_rank_ * 2 + 8) * ps;
  buff_ = stack_->get(size_allocated_);
  double* pointer = buff_;

  rvec_  = pointer;   pointer += ps * 3;
  kvec_  = pointer;   pointer += ps * 3;
  Rsq_ = pointer;       pointer += ps;
  T_ = pointer;       pointer += ps;
  roots_ = pointer;   pointer += max_rank_ * ps;
  weights_ = pointer; pointer += max_rank_ * ps;
}


void PFMM::compute_Slm() {

  const size_t nbasis1 = scell_->multipoles().front()->ndim();
  const size_t nbasis0 = scell_->multipoles().front()->mdim();
  for (int l = 0; l <= lmax_; ++l) {
    for (int m = 0; m <= 2*l; ++m) {
      const int im1 = l * l + m;

      ZMatrix local(nbasis1, nbasis0);
      for (int j = 0; j <= lmax_; ++j) {
        for (int k = 0; k <= 2*j; ++k) {
          const int im2 = j * j + k;
          const int a = l + j;
          const int b = m + k - l - j;
          const int im = a * a + m + k;
          zaxpy_(nbasis0 * nbasis1, mlm_[im], scell_->multipoles().at(im2)->data(), 1, local.data(), 1);
        }
      }
      slm_[im1] = make_shared<const ZMatrix>(local);

    }
  }
}


shared_ptr<const PData> PFMM::compute_cfmm(shared_ptr<const PData> density) const { // within ws

  Timer time;
  const int nvec = pow(2*ws_+1, ndim_);
  vector<array<int, 3>> vidx(nvec);

  int cnt = 0;
  if (ndim_ == 3) {
    for (int n3 = -ws_; n3 <= ws_; ++n3)
      for (int n2 = -ws_; n2 <= ws_; ++n2)
        for (int n1 = -ws_; n1 <= ws_; ++n1, ++cnt)
          vidx[cnt] = {{n1, n2, n3}};
  } else if (ndim_ == 2) {
    for (int n2 = -ws_; n2 <= ws_; ++n2)
      for (int n1 = -ws_; n1 <= ws_; ++n1, ++cnt)
      vidx[cnt] = {{n1, n2, 0}};
  } else if (ndim_ == 1) {
    for (int n1 = -ws_; n1 <= ws_; ++n1, ++cnt)
      vidx[cnt] = {{n1, 0, 0}};
  }
  assert(cnt == nvec);

  vector<array<double, 3>> primvecs(3);
  for (int i = 0; i != ndim_; ++i)
    primvecs[i] = scell_->primitive_vectors(i);
  for (int i = ndim_; i != 3 - ndim_; ++i)
    primvecs[i] = {{0.0, 0.0, 0.0}};

  // concatenate all cells within ws into one supercell
  vector<shared_ptr<const Geometry>> geoms;
  for (int ivec = 0; ivec != nvec; ++ivec) {
    array<int, 3> idx = vidx[ivec];
    array<double, 3> disp;
    disp[0] = idx[0] * primvecs[0][0] + idx[1] * primvecs[1][0] + idx[2] * primvecs[2][0];
    disp[1] = idx[0] * primvecs[0][1] + idx[1] * primvecs[1][1] + idx[2] * primvecs[2][1];
    disp[2] = idx[0] * primvecs[0][2] + idx[1] * primvecs[1][2] + idx[2] * primvecs[2][2];
    geoms.push_back(make_shared<const Geometry>(*scell_->geom(), disp));
  }
  auto supergeom = make_shared<const Geometry>(geoms, true/*nodf*/);
  time.tick_print("  Construct a supercell for crystal near-field");

  // get density for supergeom: ncell_ in lattice should be set to ws
  const size_t nbas = nvec * scell_->nbasis();
  auto superden = make_shared<Matrix>(nbas, nbas);
  int offset = 0;
  for (int i = 0; i != nvec; ++i, offset += scell_->nbasis()) {
    array<int, 3> idx = vidx[i];
    int block = -1;
    switch(ndim_) {
      case 1 :
        {
          block = idx[0] + ws_;
          break;
        }
      case 2 :
        {
          block = idx[0] + ws_ + (2*ws_+1) * (idx[1] + ws_);
          break;
        }
      case 3 :
        {
          block = idx[0] + ws_ + (2*ws_+1) * ((idx[1] + ws_) + (2*ws_+1) * (idx[2] + ws_));
          break;
        }
    }
    assert (block >= 0);
    superden->copy_block(offset, offset, scell_->nbasis(), scell_->nbasis(), *density->pdata(block)->get_real_part());
  }
  time.tick_print("  Construct superdensity");

  // construct a tree from the super-geometry
  Tree fmm_tree(supergeom, max_height_, do_contract_);
  fmm_tree.fmm(lmax_, superden, true /*dodf*/, auxfile_);
  shared_ptr<const ZMatrix> coulomb = fmm_tree.coulomb();

  vector<shared_ptr<const ZMatrix>> out(nvec);
  offset = 0;
  for (int i = 0; i != nvec; ++i, offset += scell_->nbasis()) {
    auto tmp = coulomb->get_submatrix(offset, offset, scell_->nbasis(), scell_->nbasis());
    out[i] = make_shared<const ZMatrix>(*tmp);
  }

  return make_shared<const PData>(out);
}


shared_ptr<const PData> PFMM::pcompute_Jop(shared_ptr<const PData> density) const {

  shared_ptr<PData> out;

  // call functions to compute jOp here

  return out;
}
