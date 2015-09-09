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
  : scell_(scell), lmax_(lmax), ws_(ws), extent_sum_(extent), thresh_(thresh) {

  if (stack == nullptr) {
    stack_ = resources__->get();
    allocated_here_ = true;
  } else {
    stack_ = stack;
    allocated_here_ = false;
  }

  assert(lmax_ < RYS_MAX);
  ndim_ = scell->ndim();
  if (ndim_ > 3 || ndim_ < 1)
    throw runtime_error("System must be periodic in 1-, 2-, or 3-D");

  num_multipoles_ = (lmax_ + 1) * (lmax_ + 1);
  mlm_.resize(num_multipoles_);
  max_rank_ = lmax_ + 1;
  compute_mlm();
}


bool PFMM::is_in_cff(array<double, 3> L) {

  const double extent = scell_->extent();

  const double rsq = L[0]*L[0] + L[1]*L[1] + L[2]*L[2];
  const bool out = (rsq > 2.0 * (1 + ws_) *  extent) ? true : false;

  return out;
}


void PFMM::compute_mlm() {

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
    T_[ivec] = Rsq_[ivec] * beta_ * beta_;
  }

  for (int l = 0; l <= lmax_; ++l) {
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

      const double coeff = 2.0 * pow(beta_, 2*l+1) * sgamma(l, r);
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

#if 0
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
    T_[ivec] = Rsq_[ivec] * beta_ * beta_;
  }

  fill_n(roots_, max_rank_ * nvec, 0.0);
  fill_n(weights_, max_rank_ * nvec, 0.0);
  for (int l = 0; l <= lmax_; ++l) {
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
#endif
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


#if 0
void PFMM::compute_Sn(const int max_iter) { // S(n+1) = U_M[S(n)] O* + M*

  for (int iter = 0; iter != max_iter; ++i) {

    const double error = 0.0;
    if (error < thresh) {
      cout << "  * Sn converged." << endl << endl;
      break;
    } else if (iter == max_iter-1) {
      cout << "  * Max iteration reached when in compute_Sn." << endl << endl;
      break;
    }
  }
}


shared_ptr<PData> PFMM::compute_Jop(shared_ptr<const PData> density) {

}
#endif
