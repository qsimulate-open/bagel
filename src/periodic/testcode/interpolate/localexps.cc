//
// Author  : Hai-Anh Le <anh@u.northwestern.edu>
// Date    : July 2015
//


#include <src/util/math/legendre.h>
#include <src/util/math/gamma.h>
#include <src/integral/rys/erirootlist.h>
#include "localexps.h"

using namespace std;
using namespace bagel;

const static Legendre plm;

LocalExps::LocalExps(const int ws, const int l, const int n, const double thr)
  : ws_(ws), lmax_(l), limit_(n), thresh_(thr) {

  stack_ = resources__->get();
  allocated_here_ = true;

  beta_ = sqrt(pi__);
  rank_ = lmax_ + 1;
  compute_mlm();
}

bool LocalExps::is_in_cff(const int ws_, const int n0, const int n1, const int n2) {

  const bool out = (abs(n0) > ws_ && abs(n1) > ws_ && abs(n2) > ws_) ? true : false;
  return out;
};


void LocalExps::compute_mlm() {

  const size_t ndim = 3;
  const size_t num_multipoles = (lmax_ + 1) * (lmax_ + 1);
  mlm_.resize(num_multipoles);

  double pibeta = (pi__ / beta_) * (pi__ / beta_);
  const int nvec = std::pow(2*limit_+1, ndim);
  allocate_arrays(nvec);
  vector<array<int, 3>> vindex(nvec);

  // generate lattice vectors - 3D for now
  vector<array<double, 3>> primitive_vectors(3);
  primitive_vectors[0] = {{1.0, 0.0, 0.0}};
  primitive_vectors[1] = {{0.0, 1.0, 0.0}};
  primitive_vectors[2] = {{0.0, 0.0, 1.0}};


  int cnt = 0;
  for (int n3 = -limit_; n3 <= limit_; ++n3) {
    for (int n2 = -limit_; n2 <= limit_; ++n2) {
      for (int n1 = -limit_; n1 <= limit_; ++n1, ++cnt) {
        const int pos = cnt * 3;
        rvec_[pos    ] = n1 * primitive_vectors[0][0] + n2 * primitive_vectors[1][0] + n3 * primitive_vectors[2][0];
        rvec_[pos + 1] = n1 * primitive_vectors[0][1] + n2 * primitive_vectors[1][1] + n3 * primitive_vectors[2][1];
        rvec_[pos + 2] = n1 * primitive_vectors[0][2] + n2 * primitive_vectors[1][2] + n3 * primitive_vectors[2][2];
        Rsq_[cnt] = rvec_[pos]*rvec_[pos] + rvec_[pos+1]*rvec_[pos+1] * rvec_[pos+2]*rvec_[pos+2];
        T_[cnt] = Rsq_[cnt] * beta_ * beta_;
        vindex[cnt] = {{n1, n2, n3}};
      }
    }
  }

  array<double, 3> a23 = cross(primitive_vectors[1], primitive_vectors[2], 1.0);
  const double scale = 2.0 * pi__ / dot(primitive_vectors[0], a23);
  vector<array<double, 3>> primitive_kvectors(3);
  primitive_kvectors[0] = cross(primitive_vectors[1], primitive_vectors[2], scale);
  primitive_kvectors[1] = cross(primitive_vectors[2], primitive_vectors[0], scale);
  primitive_kvectors[2] = cross(primitive_vectors[0], primitive_vectors[1], scale);

  cnt = 0;
  for (int n3 = -limit_; n3 <= limit_; ++n3) {
    for (int n2 = -limit_; n2 <= limit_; ++n2) {
      for (int n1 = -limit_; n1 <= limit_; ++n1, ++cnt) {
        const int pos = cnt * 3;
        kvec_[pos    ] = n1 * primitive_kvectors[0][0] + n2 * primitive_kvectors[1][0] + n3 * primitive_kvectors[2][0];
        kvec_[pos + 1] = n1 * primitive_kvectors[0][1] + n2 * primitive_kvectors[1][1] + n3 * primitive_kvectors[2][1];
        kvec_[pos + 2] = n1 * primitive_kvectors[0][2] + n2 * primitive_kvectors[1][2] + n3 * primitive_kvectors[2][2];
        Rsq_[nvec + cnt] = kvec_[pos]*kvec_[pos] + kvec_[pos+1]*kvec_[pos+1] * kvec_[pos+2]*kvec_[pos+2];
        T_[nvec + cnt] = Rsq_[nvec + cnt] * beta_ * beta_;
      }
    }
  }

  root_weight(nvec*2);

  for (int ivec = 0; ivec != nvec; ++ivec) {
    array<int, 3> id = vindex[ivec];

    const int pos = ivec * 3;
    const double rsq = Rsq_[ivec];
    const double r = sqrt(rsq);
    const double ctheta = (rsq > numerical_zero__) ? rvec_[pos+2]/r : 0.0;
    const double phi = atan2(rvec_[pos+1], rvec_[pos]);
    const double b2r2 = T_[ivec];
    const double* croots = roots_ + ivec * rank_;
    const double* cweights = weights_ + ivec * rank_;

    int cnt = 0;
    for (int l = 0; l <= lmax_; ++l) {
      Gamma gamma_func;
      const double gamma = gamma_func(2*l+1);
      const double glower = croots[l];
      Gamma_lower_scaled bagel_gs;
      const double tmp = bagel_gs(l, b2r2, beta_);
      cout << setprecision(9) << glower << " wt = " << cweights[l] << " bagel = " << tmp << endl;
//      const double gupper = gamma / pow(r, l+1.0) - gupper;

#if 0
      for (int mm = 0; mm <= 2 * l; ++mm, ++cnt) {
        const int m = mm - l;
        const int am = abs(m);

        double plm_tilde = plm.compute(l, abs(m), ctheta);
        double ft = 1.0;
        for (int i = 1; i <= l - abs(m); ++i) {
          plm_tilde *= ft;
          ft += 1.0;
        }
        const double sign = (m >=0) ? (cos(am * phi)) : (-1.0 * cos(am * phi));

        if (is_in_cff(ws_, id[0], id[1], id[2])) {
          // real term
          const double real = gupper * sign * plm_tilde / gamma;
          const double imag = gupper * sin(am * phi) * plm_tilde;
          mlm[cnt] += complex<double>(real, imag);
        }
      }
#endif
    }

  }

#if 0
  ivec = 0;
  for (auto& v : kvec_) {
    array<int, 3> id = vindex[ivec];

    const double rsq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    const double r = sqrt(rsq);
    const double ctheta = (rsq > numerical_zero__) ? v[2]/r : 0.0;
    const double phi = atan2(v[1], v[0]);
    const double b2r2 = beta_ * beta_ * rsq;

    int cnt = 0;
    for (int l = 0; l <= lmax_; ++l) {
      const complex<double> coeffl = std::pow(complex<double>(0.0, 1.0), l) * pow(pi__, l-0.5);
      Gamma gamma_func;
      const double gamma  = gamma_func(2*l+1);
      Gamma_lower_scaled lower_gamma;
      const double glower = lower_gamma(l, b2r2, beta_);
      for (int mm = 0; mm <= 2 * l; ++mm, ++cnt) {
        const int m = mm - l;
        const int am = abs(m);

        double plm_tilde = plm.compute(l, abs(m), ctheta);
        double ft = 1.0;
        for (int i = 1; i <= l - abs(m); ++i) {
          plm_tilde *= ft;
          ft += 1.0;
        }

        const double sign = (m >=0) ? (cos(am * phi)) : (-1.0 * cos(am * phi));

        if (!is_in_cff(ws_, id[0], id[1], id[2])) {
          // substract smooth part within ws_
          const double real = glower * sign * plm_tilde / gamma;
          const double imag = glower * sin(am * phi) * plm_tilde;
          mlm_[cnt] -= complex<double>(real, imag);
        }
        // smooth term
        const double coeffm = plm_tilde * pow(r, l-2) *  exp(-rsq * pibeta);
        double real = coeffm * sign / gamma;
        double imag = coeffm * sin(am * phi);
        mlm_[cnt] += coeffl * complex<double>(real, imag);
      }
    }

    ++ivec;
  }
#endif
};


void LocalExps::root_weight(const int size) {
  if (lmax_ == 0) {
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
    eriroot__.root(lmax_, T_, roots_, weights_, size);
  }
}


void LocalExps::allocate_arrays(const size_t ps) {

  buff_ = stack_->get(size_allocated_);
  double* pointer = buff_;

  rvec_  = pointer;   pointer += ps * 3;
  kvec_  = pointer;   pointer += ps * 3;
  Rsq_ = pointer;       pointer += ps * 2;
  T_ = pointer;       pointer += ps * 2;
  coeff_ = pointer;   pointer += ps;
  roots_ = pointer;   pointer += rank_ * ps;
  weights_ = pointer; pointer += rank_ * ps;
}
