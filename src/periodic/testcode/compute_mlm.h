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

  const bool out = (abs(n0) > ws && abs(n1) > ws && abs(n2) > ws) ? true : false;
  return out;
};


vector<complex<double>> compute_mlm(const int ws, const int lmax, const int limit, const double thresh) {

  const size_t ndim = 3;
  const size_t num_multipoles = (lmax + 1) * (lmax + 1);
  vector<complex<double>> out(num_multipoles);

  const double pibeta = pi__ * pi__ / (beta__ * beta__);

  // generate lattice vectors - 3D for now
  vector<array<double, 3>> primitive_vectors(3);
  primitive_vectors[0] = {{1.0, 0.0, 0.0}};
  primitive_vectors[1] = {{0.0, 1.0, 0.0}};
  primitive_vectors[2] = {{0.0, 0.0, 1.0}};

  const int nvec = pow(2*limit+1, ndim);
  vector<array<double, 3>> rvec(nvec);
  vector<array<int, 3>> vindex(nvec);

  array<double, 3> v = {{0.0, 0.0, 0.0}};
  int cnt = 0;
  for (int n3 = -limit; n3 <= limit; ++n3) {
    for (int n2 = -limit; n2 <= limit; ++n2) {
      for (int n1 = -limit; n1 <= limit; ++n1, ++cnt) {
        v[0] = n1 * primitive_vectors[0][0] + n2 * primitive_vectors[1][0] + n3 * primitive_vectors[2][0];
        v[1] = n1 * primitive_vectors[0][1] + n2 * primitive_vectors[1][1] + n3 * primitive_vectors[2][1];
        v[2] = n1 * primitive_vectors[0][2] + n2 * primitive_vectors[1][2] + n3 * primitive_vectors[2][2];
        vindex[cnt] = {{n1, n2, n3}};
        rvec[cnt] = v;
      }
    }
  }

  cnt = 0;
  for (auto& v : rvec) {
    array<int, 3> id = vindex[cnt];

    const double rsq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    const double r = sqrt(rsq);
    const double ctheta = (rsq > numerical_zero__) ? v[2]/r : 0.0;
    const double phi = atan2(v[1], v[0]);
    const double b2r2 = beta__ * beta__ * rsq;

    if (is_in_cff(ws, id[0], id[1], id[2])) {
      // real term
      int count = 0;
      for (int l = 0; l <= lmax; ++l) {
        for (int mm = 0; mm <= 2 * l; ++mm, ++count) {
          const int m = mm - l;
          const int am = abs(m);

          double coeff = plm.compute(l, abs(m), ctheta) / pow(r, l + 1) * boost::math::gamma_q(l+0.5, b2r2);
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
    } else if(abs(id[0]) <= ws && abs(id[1]) <= ws && abs(id[2]) <= ws) {
      // substract smooth part within ws
      int count = 0;
      for (int l = 0; l <= lmax; ++l) {
        for (int mm = 0; mm <= 2 * l; ++mm, ++count) {
          const int m = mm - l;
          const int am = abs(m);

          double coeff = plm.compute(l, abs(m), ctheta) / pow(r, l + 1) * boost::math::gamma_p(l+0.5, b2r2);
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
    int count = 0;
    for (int l = 0; l <= lmax; ++l) {
      const complex<double> coeffl = pow(complex<double>(0.0, 1.0), l) * pow(pi__, l-0.5) / boost::math::tgamma(l+0.5);
      for (int mm = 0; mm <= 2 * l; ++mm, ++count) {
        const int m = mm - l;
        const int am = abs(m);

        double coeffm = plm.compute(l, abs(m), ctheta) * pow(r, l - 2) * exp(-rsq * pibeta);
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

    ++cnt;
  }

  return out;
};
