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

double dot(array<double, 3> b, array<double, 3> c) { return b[0] * c[0] + b[1] * c[1] + b[2] * c[2]; }


array<double, 3> cross(array<double, 3> b, array<double, 3> c, double s) {

  array<double, 3> out;
  out[0] = (b[1] * c[2] - b[2] * c[1]) * s;
  out[1] = (b[2] * c[0] - b[0] * c[2]) * s;
  out[2] = (b[0] * c[1] - b[1] * c[0]) * s;

  return out;
}


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
  vector<array<double, 3>> kvec(nvec);
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

  int ivec = 0;
  for (auto& v : rvec) {
    array<int, 3> id = vindex[ivec];

    const double rsq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    const double r = sqrt(rsq);
    const double ctheta = (rsq > numerical_zero__) ? v[2]/r : 0.0;
    const double phi = atan2(v[1], v[0]);
    const double b2r2 = beta__ * beta__ * rsq;

    int cnt = 0;
    for (int l = 0; l <= lmax; ++l) {
      assert(boost::math::gamma_q(l+0.5, b2r2) - gsl_sf_gamma_inc_Q(l+0.5, b2r2) < 1e-10);

      const double gupper = boost::math::tgamma(l+0.5, b2r2) / pow(r, l+1.0);
      for (int mm = 0; mm <= 2 * l; ++mm, ++cnt) {
        const int m = mm - l;
        const int am = abs(m);

        double plm_tilde = plm.compute(l, abs(m), ctheta);
        double ft = 1.0;
        for (int i = 1; i <= l - abs(m); ++i) {
          plm_tilde *= ft;
          ++ft;
        }

        const double sign = (m >=0) ? (cos(am * phi)) : (-1.0 * cos(am * phi));

        if (is_in_cff(ws, id[0], id[1], id[2])) {
          // real term
          const double real = gupper * sign * plm_tilde;
          const double imag = gupper * sin(am * phi) * plm_tilde;
          out[cnt] += complex<double>(real, imag) / boost::math::tgamma(l+0.5);
        }
      }
    }

    ++ivec;
  }

  array<double, 3> a23 = cross(primitive_vectors[1], primitive_vectors[2], 1.0);
  const double scale = 1.0 / dot(primitive_vectors[0], a23);
  //const double scale = 2.0 * pi__ / dot(primitive_vectors[0], a23);
  vector<array<double, 3>> primitive_kvectors(3);
  primitive_kvectors[0] = cross(primitive_vectors[1], primitive_vectors[2], scale);
  primitive_kvectors[1] = cross(primitive_vectors[2], primitive_vectors[0], scale);
  primitive_kvectors[2] = cross(primitive_vectors[0], primitive_vectors[1], scale);
//  cout << "[" << setprecision(9) << primitive_kvectors[0][0] << ", " << primitive_kvectors[0][1] << ", " << primitive_kvectors[0][2] << "]" << endl;
//  cout << "[" << setprecision(9) << primitive_kvectors[1][0] << ", " << primitive_kvectors[1][1] << ", " << primitive_kvectors[1][2] << "]" << endl;
//  cout << "[" << setprecision(9) << primitive_kvectors[2][0] << ", " << primitive_kvectors[2][1] << ", " << primitive_kvectors[2][2] << "]" << endl;

  cnt = 0;
  for (int n3 = -limit; n3 <= limit; ++n3) {
    for (int n2 = -limit; n2 <= limit; ++n2) {
      for (int n1 = -limit; n1 <= limit; ++n1, ++cnt) {
        v[0] = n1 * primitive_kvectors[0][0] + n2 * primitive_kvectors[1][0] + n3 * primitive_kvectors[2][0];
        v[1] = n1 * primitive_kvectors[0][1] + n2 * primitive_kvectors[1][1] + n3 * primitive_kvectors[2][1];
        v[2] = n1 * primitive_kvectors[0][2] + n2 * primitive_kvectors[1][2] + n3 * primitive_kvectors[2][2];
        kvec[cnt] = v;
      }
    }
  }

  ivec = 0;
  for (auto& v : kvec) {
    array<int, 3> id = vindex[ivec];

    const double rsq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    const double r = sqrt(rsq);
    const double ctheta = (rsq > numerical_zero__) ? v[2]/r : 0.0;
    const double phi = atan2(v[1], v[0]);
    const double b2r2 = beta__ * beta__ * rsq;

    int cnt = 0;
    for (int l = 0; l <= lmax; ++l) {
      assert(boost::math::gamma_p(l+0.5, b2r2) - gsl_sf_gamma_inc_P(l+0.5, b2r2) < 1e-10);
      const complex<double> coeffl = pow(complex<double>(0.0, 1.0), l) * pow(pi__, l-0.5);
      const double glower = boost::math::tgamma_lower(l+0.5, b2r2) / pow(r, l+1.0);
      cout << "l = " << l << " b2r2 = " << b2r2 << endl;
      bagel::Gamma_lower_scaled gamma_s;
      const double bagel_glower = gamma_s(l, b2r2, beta__);
      if (abs(bagel_glower - glower) > 1e-14) {
        cout << "(l, z) = (" << l << ", " << setprecision(5) << b2r2 << ")" << endl;
        cout << "boost      = " << setw(20) << setprecision(12) << glower << endl;
        cout << "bagel      = " << setw(20) << setprecision(12) << bagel_glower << endl;
      }
      for (int mm = 0; mm <= 2 * l; ++mm, ++cnt) {
        const int m = mm - l;
        const int am = abs(m);

        double plm_tilde = plm.compute(l, abs(m), ctheta);
        double ft = 1.0;
        for (int i = 1; i <= l - abs(m); ++i) {
          plm_tilde *= ft;
          ++ft;
        }

        const double sign = (m >=0) ? (cos(am * phi)) : (-1.0 * cos(am * phi));

        if (!is_in_cff(ws, id[0], id[1], id[2])) {
          // substract smooth part within ws
          const double real = glower * sign * plm_tilde;
          const double imag = glower * sin(am * phi) * plm_tilde;
          out[cnt] -= complex<double>(real, imag) / boost::math::tgamma(l+0.5);
        }
        // smooth term
        double coeffm = plm_tilde * pow(r, l - 2) * exp(-rsq * pibeta);
        const double real = coeffm * sign;
        const double imag = coeffm * sin(am * phi);
        out[cnt] += coeffl * complex<double>(real, imag) / boost::math::tgamma(l+0.5);
//        cout << "(lm) = (" << l << m << ")  =   " << setprecision(9) << out[cnt].real() << endl;
      }
    }

    ++ivec;
  }

  return out;
};
