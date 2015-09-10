//
// BAGEL - Parallel electron correlation program.
// Filename: stevensop.cc
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Ryan D. Reynolds <rreynolds2018@u.northwestern.edu>
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

// A function to compute the coefficients of Extended Stevens Operators, for the pseudospin Hamiltonian
// Notation follows I. D. Ryabov, Appl. Magn. Reson. (2009) 35, 481-494.
// Some equations also come from I. D. Ryabov, J. Magn. Reson. (1999) 140, 141-145.

#include <src/ci/zfci/zharrison.h>
#include <src/util/math/factorial.h>


using namespace std;
using namespace bagel;

// Some helper functions only used here
namespace {

const Factorial fact;

// Even k:  alpha = 1.0
// Odd k:   alpha = 1.0 or 0.5 for even or odd q, respectively
double compute_alpha(const int k, const int q) {
  double out;
  if (k % 2 == 0)
    out = 1.0;
  else
    out = (q % 2 == 0) ? 1.0 : 0.5;
  return out;
}

double compute_Nkk(const int k) {
  const double sign = (k % 2 == 0) ? 1.0 : -1.0;
  const double twok = std::pow(2.0,k);
  const double denomenator = twok * fact(k);
  const double f2k = fact(2 * k);
  const double numerator = std::sqrt(f2k);
  const double out = sign * numerator / denomenator;
  return out;
}

double compute_Nkq(const int k, const int q, const double Nkk) {
  double out;
  assert (k >= 0 && q >= 0 && q <= k);
  const double sign = (k - q) % 2 == 0 ? 1.0 : -1.0;
  const double numerator = fact(k + q);
  const double denomenator = (fact(k - q) * fact(2 * k));
  out = sign * Nkk * std::sqrt(numerator / denomenator);
  return out;
}

// a(k, q; m) in Ryabov's notation
double compute_akqm(const int k, const int q, const int m) {
  assert (k >= 0 && q >= 0 && m >= 0 && q <= k && m <= k - q);
  double out = 0.0;
  // TODO
  return out;
}

} // end of anonymous namespace


// And now the driver
void ZHarrison::compute_extended_stevens_operators() const {

  /***************************************************/
  /* TODO This section is repeated from zharrison.cc */
  /***************************************************/
  int nspin = idata_->get<int>("aniso_spin", states_.size()-1);
  int nspin1 = nspin + 1;

  // S_x, S_y, and S_z operators in pseudospin basis
  array<shared_ptr<ZMatrix>,3> pspinmat;
  for (int i = 0; i != 3; ++i)
    pspinmat[i] = make_shared<ZMatrix>(nspin1, nspin1);

  auto spin_plus = make_shared<ZMatrix>(nspin1, nspin1);
  auto spin_minus = make_shared<ZMatrix>(nspin1, nspin1);
  const double sval = nspin/2.0;
  const double ssp1 = sval*(sval+1.0);

  for (int i = 0; i != nspin1; ++i) {
    const double ml1 = sval - i;
    pspinmat[2]->element(i,i) = ml1;
    if (i < nspin) {
      const double ml1m = ml1 - 1.0;
      spin_plus->element(i,i+1) = std::sqrt(ssp1 - ml1*(ml1m));
    }
    if (i > 0) {
      const double ml1p = ml1 + 1.0;
      spin_minus->element(i,i-1) = std::sqrt(ssp1 - ml1*(ml1p));
    }
  }

  pspinmat[0]->add_block( 0.5, 0, 0, nspin1, nspin1, spin_plus);
  pspinmat[0]->add_block( 0.5, 0, 0, nspin1, nspin1, spin_minus);
  pspinmat[1]->add_block( complex<double>( 0.0, -0.5), 0, 0, nspin1, nspin1, spin_plus);
  pspinmat[1]->add_block( complex<double>( 0.0,  0.5), 0, 0, nspin1, nspin1, spin_minus);

  /***************************************************/
  /*                    fin                          */
  /***************************************************/

  cout << "    Computing extended stevens operators for S = " << nspin/2 << (nspin % 2 == 0 ? "" : " 1/2") << endl;

  const int kmax = idata_->get<bool>("aniso_extrastevens", false) ? 8 : nspin;

  // Requires (2k)!, and factorials are tabulated up to 20
  if (kmax > 10)
    throw runtime_error("Sorry, numerical issues currently limit us to Stevens operators of 10th order and lower");

  cout << fixed << setprecision(6);
  for (int k = 0; k <= kmax; ++k) {

    const double Nkk = compute_Nkk(k);
    vector<double> alpha(k + 1);
    vector<double> Nkq(k + 1);  // positive q
    vector<double> Nk_q(k + 1); // negative q

    // a(k, q; m) in Ryabov's notation
    // access is akqm[q][m]
    vector<vector<double>> akqm(k + 1);

    // same but for negative q
    vector<vector<double>> ak_qm(k + 1);

    for (int q = k; q >= 0; --q) {

      vector<double> akqm_current(k - q + 1);  // positive q
      vector<double> ak_qm_current(k - q + 1); // negative q

      alpha[q] = compute_alpha(k, q);
      Nkq[q] = compute_Nkq(k, q, Nkk);
      Nk_q[q] = Nkq[q] * (k % 2 == 0) ? 1.0 : -1.0;


      for (int m = 0; m <= k - q; ++m) {
        // TODO
        akqm_current[m] = compute_akqm(k, q, m);
        ak_qm_current[m] = 0.0;

        cout << "       k = " << setw(4) << k << ", q = " << setw(4) << q << ", m = " << setw(4) << m << ", alpha = " << setw(12) << alpha[q] << ", Nkk = " << setw(12) << Nkk << ", Nkq = " << setw(12) << Nkq[q] << ", Nk_q = " << setw(12) << Nk_q[q] << ", Nkq/Nkk = " << setw(12) << Nkq[q] / Nkk << "  acoeff = " << setw(12) << akqm_current[m] << endl;
      }

      assert(q == k ? true : (std::abs(Nkq[q] + (Nkq[q+1] / std::sqrt((k + q + 1)*(k - q))))) < 1.0e-12); // check against a recurrence relation

      akqm[q] = akqm_current;
      ak_qm[q] = ak_qm_current;

    }
    cout << endl;
  }
}

