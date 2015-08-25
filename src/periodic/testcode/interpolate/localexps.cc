//
// Author  : Hai-Anh Le <anh@u.northwestern.edu>
// Date    : July 2015
//


#include <src/util/math/legendre.h>
#include <src/util/math/sphharmonics.h>
#include <src/util/math/gamma.h>
#include <src/integral/rys/erirootlist.h>
#include "localexps.h"
#include <boost/math/special_functions/erf.hpp>

using namespace std;
using namespace bagel;
using namespace test;

const static Legendre plm;
const static Gamma_scaled gamma_sc;

LocalExps::LocalExps(const int ws, const int l, const int n, const double thr)
  : ws_(ws), lmax_(l), limit_(n), thresh_(thr) {

  stack_ = resources__->get();
  allocated_here_ = true;

  beta_ = sqrt(pi__);
  max_rank_ = lmax_ + 1;
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
  int cnt = 0;
  for (int n3 = -limit_; n3 <= limit_; ++n3)
    for (int n2 = -limit_; n2 <= limit_; ++n2)
      for (int n1 = -limit_; n1 <= limit_; ++n1, ++cnt)
        vindex[cnt] = {{n1, n2, n3}};
  assert(cnt == nvec);

  std::sort(vindex.begin(), vindex.end(), sort_vector); // sort to sum spherically

  vector<array<double, 3>> primitive_vectors(3);
  primitive_vectors[0] = {{1.0, 0.0, 0.0}};
  primitive_vectors[1] = {{0.0, 1.0, 0.0}};
  primitive_vectors[2] = {{0.0, 0.0, 1.0}};


  for (int ivec = 0; ivec != nvec; ++ivec) {
    const int pos = ivec * 3;
    array<int, 3> idx = vindex[ivec];
    rvec_[pos    ] = idx[0] * primitive_vectors[0][0] + idx[1] * primitive_vectors[1][0] + idx[2] * primitive_vectors[2][0];
    rvec_[pos + 1] = idx[0] * primitive_vectors[0][1] + idx[1] * primitive_vectors[1][1] + idx[2] * primitive_vectors[2][1];
    rvec_[pos + 2] = idx[0] * primitive_vectors[0][2] + idx[1] * primitive_vectors[1][2] + idx[2] * primitive_vectors[2][2];
    Rsq_[ivec] = rvec_[pos]*rvec_[pos] + rvec_[pos+1]*rvec_[pos+1] + rvec_[pos+2]*rvec_[pos+2];
    T_[ivec] = Rsq_[ivec] * beta_ * beta_;
  }

  for (int l = 0; l <= lmax_; ++l) {
    root_weight(l, nvec);
    const int rank = l + 1;

    for (int ivec = 0; ivec != nvec; ++ivec) {
      array<int, 3> id = vindex[ivec];

      const int pos = ivec * 3;
      const double rsq = Rsq_[ivec];
      const double r = sqrt(rsq);
      const double ctheta = (rsq > numerical_zero__) ? rvec_[pos+2]/r : 0.0;
      const double phi = atan2(rvec_[pos+1], rvec_[pos]);
      const double b2r2 = T_[ivec];
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

      const double coeff = 2.0 * pow(beta_, 2*l+1) * gamma_sc(l, r);
#if 0 ////DEBUG
      const double boost_gamma_sc = pow(r, l) / boost::math::tgamma(l+0.5);
      const double bagel_gamma_sc = gamma_sc(l, r);
      const double error1 = boost_gamma_sc - bagel_gamma_sc;
      if (abs(error1) > 1e-14) {
        cout << setprecision(9) << rvec_[pos] << "  " << rvec_[pos+1] << "  " << rvec_[pos+2] << "   ";
        cout << "(l,x) = (" << l << ", " << r << ")  " << setprecision(9) << bagel_gamma_sc << " ***  boost = " << boost_gamma_sc << endl;
      }

      const double boost_gamma = 0.5 * boost::math::tgamma_lower(l+0.5, b2r2)/std::pow(b2r2, l+0.5);
      const double error = boost_gamma - glower;
      if (abs(error) > 1e-14) {
        cout << setprecision(9) << rvec_[pos] << "  " << rvec_[pos+1] << "  " << rvec_[pos+2] << "   ";
        cout << "(l,x) = (" << l << ", " << b2r2 << ")  " << setprecision(9) << glower << " ***  boost = " << boost_gamma << endl;
      }
#endif ///END OF DEBUG
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

        if (is_in_cff(ws_, id[0], id[1], id[2])) {
          // real term
          const double real = gupper * sign * plm_tilde;
          const double imag = gupper * sin(am * phi) * plm_tilde;
          mlm_[imul] += complex<double>(real, imag);
        }
      }
    }
  }

  array<double, 3> a23 = cross(primitive_vectors[1], primitive_vectors[2], 1.0);
  const double scale = 1.0 / dot(primitive_vectors[0], a23);
  /////const double scale = 2.0 * pi__ / dot(primitive_vectors[0], a23);
  vector<array<double, 3>> primitive_kvectors(3);
  primitive_kvectors[0] = cross(primitive_vectors[1], primitive_vectors[2], scale);
  primitive_kvectors[1] = cross(primitive_vectors[2], primitive_vectors[0], scale);
  primitive_kvectors[2] = cross(primitive_vectors[0], primitive_vectors[1], scale);

  for (int ivec = 0; ivec != nvec; ++ivec) {
    const int pos = ivec * 3;
    array<int, 3> idx = vindex[ivec];
    kvec_[pos    ] = idx[0] * primitive_kvectors[0][0] + idx[1] * primitive_kvectors[1][0] + idx[2] * primitive_kvectors[2][0];
    kvec_[pos + 1] = idx[0] * primitive_kvectors[0][1] + idx[1] * primitive_kvectors[1][1] + idx[2] * primitive_kvectors[2][1];
    kvec_[pos + 2] = idx[0] * primitive_kvectors[0][2] + idx[1] * primitive_kvectors[1][2] + idx[2] * primitive_kvectors[2][2];
    Rsq_[ivec] = kvec_[pos]*kvec_[pos] + kvec_[pos+1]*kvec_[pos+1] + kvec_[pos+2]*kvec_[pos+2];
    T_[ivec] = Rsq_[ivec] * beta_ * beta_;
  }

#if 1
  fill_n(roots_, max_rank_ * nvec, 0.0);
  fill_n(weights_, max_rank_ * nvec, 0.0);
  for (int l = 0; l <= lmax_; ++l) {
    const complex<double> coeffl = std::pow(complex<double>(0.0, 1.0), l) * pow(pi__, l-0.5);
    root_weight(l, nvec);
    const int rank = l + 1;


    for (int ivec = 0; ivec != nvec; ++ivec) {
      array<int, 3> id = vindex[ivec];

      const int pos = ivec * 3;
      const double rsq = Rsq_[ivec];
      const double r = sqrt(rsq);
      const double ctheta = (rsq > numerical_zero__) ? kvec_[pos+2]/r : 0.0;
      const double phi = atan2(kvec_[pos+1], kvec_[pos]);
      const double b2r2 = T_[ivec];
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

      const double gamma_coeff = gamma_sc(l, r);
      const double coeff = 2.0 * pow(beta_, 2*l+1) * gamma_coeff;

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

        if (abs(id[0]) <= ws_ && abs(id[1]) <= ws_ && abs(id[2]) <= ws_) {
          // substract smooth part within ws_
          const double real = coeff * glower * sign * plm_tilde;
          const double imag = coeff * glower * sin(am * phi) * plm_tilde;
          mlm_[imul] -= complex<double>(real, imag);
        }
        // smooth term
        const double coeffm = (r > numerical_zero__) ? plm_tilde * gamma_coeff *  exp(-rsq * pibeta) / rsq : 0.0;
        double real = coeffm * sign;
        double imag = coeffm * sin(am * phi);
        mlm_[imul] += coeffl * complex<double>(real, imag);
      }
    }
  }
#endif
};


void LocalExps::root_weight(const int l, const int size) {
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


void LocalExps::allocate_arrays(const size_t ps) {

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
