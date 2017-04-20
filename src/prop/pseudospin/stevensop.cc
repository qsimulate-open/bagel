//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: stevensop.cc
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Ryan D. Reynolds <rreynolds2018@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

// A function to compute the coefficients of Extended Stevens Operators, for the pseudospin Hamiltonian
// Notation follows I. D. Ryabov, Appl. Magn. Reson. (2009) 35, 481-494.
// Some equations also come from I. D. Ryabov, J. Magn. Reson. (1999) 140, 141-145.

#include <src/prop/pseudospin/pseudospin.h>
#include <src/util/math/factorial.h>
#include <src/util/math/comb.h>
#include <src/util/math/matop.h>


using namespace std;
using namespace bagel;

// Some helper functions only used here
namespace {

const Factorial fact;
const Comb comb;

// Odd k:  alpha = 1.0
// Even k: alpha = 1.0 or 0.5 for even or odd q, respectively
double compute_alpha(const int k, const int q) {
  double out;
  if (k % 2 == 0)
    out = (q % 2 == 0) ? 1.0 : 0.5;
  else
    out = 1.0;
  return out;
}

// Compute Nkq for q = 0
double compute_Nk0(const int k) {
  assert (k >= 0);
  const double twok = pow(2.0,k);
  const double nkk_denomenator = twok * fact(k);
  const double out = 1.0 / nkk_denomenator;
  return out;
}

// a(k, q; m, i) in Ryabov's notation
// Uses a downward recurrence relation:  the input vector akq1m contains the output for k, q+1, and all values of m, i
// Output corresponds to a specific k, q; m but all values of i
vector<long long> compute_akqmi(const int k, const int q, const int m, const int nspin, const vector<vector<long long>>& akq1mi) {
  assert (k >= 1 && q >= 0 && m >= 0 && q < k && m <= k - q);

  const int nval = (k - q - m) / 2 + 1;
  vector<long long> out(nval);

  for (int i=0; i!=nval; ++i) {

    out[i] = m > 0 ? (2 * q + m + 1) * akq1mi[m - 1][i] : 0;

    if (m <= k - q - 1 && i <= (k - q - m - 1) / 2)
      out[i] += (q * (q + 1) - m * (m + 1) / 2) * akq1mi[m][i];

    for (int n = 1; n <= k - q - m - 1; ++n) {
      const long long sign = (n % 2 == 0) ? 1 : -1;
      const long long coeff1 = comb(m + n, m);
      const long long coeff2 = m > 0 ? comb(m + n, m - 1) : 0;
      const long long coeff3 = m > 1 ? comb(m + n, m - 2) : 0;

      if (i <= (k - q - m - n + 1) / 2 && i > 0)
        out[i] += sign * coeff1 * akq1mi[m + n][i-1];

      if (i <= (k - q - m - n - 1) / 2)
        out[i] -= sign * (coeff2 + coeff3) * akq1mi[m + n][i];
    }
  }
  return out;
}

long long greatest_common_factor(const long long a, const long long b) {
  return b == 0 ? a : greatest_common_factor(b, a % b);
}

// fkq is the greatest common factor of all a(k, q; m, i) with the same k, q
long long compute_fkq(const vector<vector<long long>>& input) {
  long long out = abs(input[0][0]);
  for(int m = 0; m != input.size(); ++m)
    for(int i = 0; i != input[m].size(); ++i)
      out = greatest_common_factor(out, abs(input[m][i]));
  return out;
}

} // end of anonymous namespace


// And now the driver
vector<Stevens_Operator> Pseudospin::build_extended_stevens_operators(const vector<int> ranks) const {
  if (ranks.size()) {
    cout << "    Using extended Stevens operators for rank" << ((ranks.size() > 1) ? "s " : " ");
    for (int i = 0; i != ranks.size(); ++i)
      cout << ranks[i] << " ";
    cout << endl;
  }
  const double ss1 = nspin_ * (nspin_ + 2.0) / 4.0; // S(S+1)

  vector<Stevens_Operator> stevensop = {};

  for (auto& k : ranks) {

    // Requires factorial of k
    // Above k = 12, coefficients become so large that long long fails to capture them.
    // Potentially we could tabulate a(k, q; m) / fkq or something, but for now this is fine.
    const int kmax = min(fact.max(), 13);
    if (k >= kmax)
      throw runtime_error("Sorry, numerical issues currently limit us to Stevens operators of order " + to_string(kmax - 1) + " and lower");
    if (k < 0)
      throw runtime_error("Ranks of Extended Stevens Operators must be whole numbers.");

    vector<double> alpha(k + 1);
    vector<double> nkq(k + 1);  // positive q
    vector<double> nk_q(k + 1); // negative q

    // a(k, q; m) in Ryabov's notation
    // access is akqm[q][m]
    vector<vector<double>> akqm(k + 1);

    // same but for negative q
    vector<vector<double>> ak_qm(k + 1);

    // The coefficients that go into computation of akqm
    vector<vector<vector<long long>>> akqmi(k + 1);

    alpha[0] = 1.0;
    nkq[0] = compute_Nk0(k);

    for (int q = 1; q <= k; ++q) {
      alpha[q] = compute_alpha(k, q);
      nkq[q] = -1.0 * nkq[q-1] * sqrt((k + q) * (k - q + 1.0));
    }

    for (int q = k; q >= 0; --q) {

      vector<double> akqm_current(k - q + 1);  // positive q
      vector<double> ak_qm_current(k - q + 1); // negative q

      vector<vector<long long>> akqmi_current(k - q + 1);

      nk_q[q] = nkq[q] * ((k % 2 == 0) ? 1.0 : -1.0);

      for (int m = 0; m <= k - q; ++m) {
        if (q == k)
          akqmi_current[m] = vector<long long>(1, 1);
        else
          akqmi_current[m] = compute_akqmi(k, q, m, nspin_, akqmi[q + 1]);

        akqm_current[m] = 0.0;
        for (int i=0; i <= (k-q-m)/2; ++i)
          akqm_current[m] += pow(ss1, i) * akqmi_current[m][i];

        ak_qm_current[m] = akqm_current[m] * (m % 2 == 0 ? 1.0 : -1.0);
      }

      const long long fkq = compute_fkq(akqmi_current);

      akqm[q] = akqm_current;
      akqmi[q] = akqmi_current;
      ak_qm[q] = ak_qm_current;

      shared_ptr<ZMatrix> tkq = spin_plus()->clone();
      shared_ptr<ZMatrix> tk_q = spin_plus()->clone();
      const double sign1 = (q % 2 == 0) ? 1.0 : -1.0;

      for (int m = 0; m <= k - q; ++m) {
        // Need spin-z matrix to power of m and spin-plus to power of q
        shared_ptr<ZMatrix> szm = spin_xyz(2)->clone();
        shared_ptr<ZMatrix> spq = spin_xyz(2)->clone();
        szm->unit();
        spq->unit();
        for (int i = 0; i != m; ++i)
          *szm *= *spin_xyz(2);
        for (int i = 0; i != q; ++i)
          *spq *= *spin_plus();

        const double sign2 = ((k - m) % 2 == 0) ? 1.0 : -1.0;
        const double coeff = sign1 * sign2 * nkq[q] * akqm[q][m];
        *tkq += (coeff * *szm * *spq);
      }

      *tk_q = sign1 * *tkq->transpose_conjg();

      const double ckq = alpha[q] / (nkq[q] * fkq);

      auto ocos_kq = make_shared<const ZMatrix>(complex<double>(0.5,  0.0) * ckq * (*tkq + *tkq->transpose_conjg()));
      auto osin_kq = make_shared<const ZMatrix>(complex<double>(0.0, -0.5) * ckq * (*tkq - *tkq->transpose_conjg()));

      stevensop.emplace_back(ocos_kq, k, q);
      if (q != 0)
        stevensop.emplace_back(osin_kq, k, -q);
    }
  }
  return stevensop;
}


