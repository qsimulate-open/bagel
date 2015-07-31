//
// Author  : Hai-Anh Le <anh@u.northwestern.edu>
// Date    : July 2015
//


#include <boost/math/special_functions/gamma.hpp>
#include <src/util/math/legendre.h>

using namespace std;
using namespace bagel;

const static Legendre plm;
const static Factorial f;

const static double beta__ = sqrt(pi__); // convergence parameter

bool is_in_cff(const int ws, const int n0, const int n1, const int n2) {

  const bool out = (abs(n0) > ws && abs(n1) > ws and abs(n2) > ws) ? true : false;
  return out;
};


vector<complex<double>> compute_mlm(const int ws, const int lmax, const int limit, const double thresh) {

  const size_t ndim = 3;
  const size_t num_multipoles = (lmax + 1) * (lmax + 1);
  vector<complex<double>> out(num_multipoles);

  const double pibeta = pi__ * pi__ / (beta__ * beta__);

  vector<array<double, 3>> primitive_vectors(3);
  primitive_vectors[0] = {{1.0, 0.0, 0.0}};
  primitive_vectors[1] = {{0.0, 1.0, 0.0}};
  primitive_vectors[2] = {{0.0, 0.0, 1.0}};

  array<double, 3> v = {{0.0, 0.0, 0.0}};
  for (int n0 = -limit; n0 != limit; ++n0) {
    array<double, 3> prim0 = primitive_vectors[0];
    v[0] += n0 * prim0[0];
    v[1] += n0 * prim0[1];
    v[2] += n0 * prim0[2];
    for (int n1 = -limit; n1 != limit; ++n1) {
      if (ndim > 1) {
        array<double, 3> prim1 = primitive_vectors[1];
        v[0] += n1 * prim1[0];
        v[1] += n1 * prim1[1];
        v[2] += n1 * prim1[2];
      }
      for (int n2 = -limit; n2 != limit; ++n2) {
        if (ndim > 2) {
          array<double, 3> prim2 = primitive_vectors[2];
          v[0] += n2 * prim2[0];
          v[1] += n2 * prim2[1];
          v[2] += n2 * prim2[2];
        }

        const double rsq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
        const double ctheta = (rsq > numerical_zero__) ? v[2]/sqrt(rsq) : 0.0;
        const double phi = atan2(v[1], v[0]);
        const double b2r2 = beta__ * beta__ * rsq;

        if (is_in_cff(ws, n0, n1, n2)) {
          // real term
          int count = 1;
          out[0] += 1.0;
          for (int l = 1; l <= lmax; ++l) {
            for (int mm = 0; mm <= 2 * l; ++mm, ++count) {
              const int m = mm - l;
              const int am = abs(m);

              double coeff = plm.compute(l, abs(m), ctheta) / pow(rsq, (l + 1)/2) * boost::math::gamma_q(l+0.5, b2r2);
              double ft = 1.0;
              for (int i = 1; i <= l - abs(m); ++i) {
                coeff *= ft;
                ++ft;
              }

              const double real = (m >=0) ? (coeff * cos(am * phi)) : (-1.0 * coeff * cos(am * phi));
              const double imag = coeff * sin(am * phi);
              out[count] += complex<double>(real, imag);

            }
          }
        } else {
          // substract smooth part within ws
          int count = 1;
          for (int l = 1; l <= lmax; ++l) {
            for (int mm = 0; mm <= 2 * l; ++mm, ++count) {
              const int m = mm - l;
              const int am = abs(m);

              double coeff = plm.compute(l, abs(m), ctheta) / pow(rsq, (l + 1)/2) * boost::math::gamma_p(l+0.5, b2r2);
              double ft = 1.0;
              for (int i = 1; i <= l - abs(m); ++i) {
                coeff *= ft;
                ++ft;
              }

              const double real = (m >=0) ? (coeff * cos(am * phi)) : (-1.0 * coeff * cos(am * phi));
              const double imag = coeff * sin(am * phi);
              out[count] -= complex<double>(real, imag);
            }
          }
        }

        // smooth term
        int count = 1;
        for (int l = 1; l <= lmax; ++l) {
          const complex<double> coeffl = pow(complex<double>(0.0, 1.0), l) * pow(pi__, l-0.5) / boost::math::tgamma(l+0.5);
          for (int mm = 0; mm <= 2 * l; ++mm, ++count) {
            const int m = mm - l;
            const int am = abs(m);

            double coeffm = plm.compute(l, abs(m), ctheta) * pow(rsq, (l - 2)/2) * exp(-rsq * pibeta);
            double ft = 1.0;
            for (int i = 1; i <= l - abs(m); ++i) {
              coeffm *= ft;
              ++ft;
            }

            const double real = (m >=0) ? (coeffm * cos(am * phi)) : (-1.0 * coeffm * cos(am * phi));
            const double imag = coeffm * sin(am * phi);
            out[count] += coeffl * complex<double>(real, imag);
          }
        }

      } //n2
    } //n1
  } //n0

  return out;
};
