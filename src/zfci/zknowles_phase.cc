//
// BAGEL - Parallel electron correlation program.
// Filename: zknowles_phase.cc
// Copyright (C) 2013 Michael Caldwell
//
// Author: Michael Caldwell <caldwell@u.northwestern.edu>>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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

#include <cmath>
#include <random>
#include <src/zfci/zknowles.h>

using namespace std;
using namespace bagel;

void ZKnowlesHandy::mult_phase_factor() {

  // (1) implement ZMOFile and replace MOFile in zfci.h with ZMOFile (and check)
  // (2) these two lines should be in complex
  // (3) think about compression
  // (4) make a matrix like below


  ZMatrix phase(norb_, norb_);
  for (int i = 0; i != norb_; ++i) {
    const double ran = rand();
    const complex<double> fac(cos(ran), sin(ran));
    phase(i,i) = fac;
  }

  auto mo1e = make_shared<Matrix>(norb_, norb_);
  unique_ptr<double[]> mo2e(new double[nij()]);
#if 0
  ZMatrix complex_mo2e(ij, ij);
//generating complex matrix from mofile real matrix
  if (nelec == 2) {
  //temporary fix to allow for non-complex mo2e.
    cout << endl << endl << "printing complex_m02e" << endl << endl;
    for (int j = 0; j < ij; j++) {
      cout << j << setw(1);
      for (int i = 0; i < ij; i++) {
        complex_moe(i, j) = complex<double>(jop->mo2e(i,j), 0.0);
      }
    }
    complex_moe.print();
  } else if (nelec==1) {
  //temporary fix to allow for non-complex mo1e.
    cout << endl << endl << "printing complex_m01e" << endl << endl;
    for (int j = 0; j < ij;j++) {
      cout << j << setw(1);
      for (int i =0; i < ij; i++) {
        complex_moe[i] = complex<double>(jop->mo1e(i), 0.0);
        cout << setprecision(1) << complex_moe[i+j*ij] << setw(1);
      }
      cout << endl;
    }
  }
//introducing phase factor to 2e integrals to check for bugs in complex portion of zfci
//random numbers for use in phase factor
  complex<double> random[ij];
  for (int j=0;j<ij;j++) {
    double a = rand() % 10 + 1;
    random[j] = complex<double>(0,a);
  }
//declaring and filling intermediate w/ 0.0
  complex<double> intermediate[ij*ij];
  for (int j=0; j<(ij*ij); j++) intermediate[j] =0.0;

//generating diagonal phase factor matrix
  complex<double> phase_factor[ij*ij];
  cout << endl << endl << "printing phase factor" << endl << endl;
  for (int j=0;j<(ij);j++) {
    cout << j << setw(1);
    for (int i=0; i<(ij); i++) {
      if (i==j) phase_factor[i+j*ij] = random[i];
      else phase_factor[i+j*ij] = 0.0;
      assert(phase_factor[i+j*ij].real()<1e-8);
      cout << setprecision(1) << phase_factor[i+j*ij] << setw(1);
    }
    cout << endl;
  }
// (A*U)
  zgemm3m_("n", "n", ij, ij, ij, 1.0, complex_moe, ij, phase_factor, ij, 0.0, intermediate, ij);

  cout << endl << endl << "printing intermediate" << endl << endl;
  for (int j=0; j<(ij); j++) {
    cout << j << setw(1);
    for (int i=0; i<(ij); i++) {
      cout << setprecision(1) << intermediate[i+j*ij] << setw(1);
    }
    cout << endl;
  }
  cout << endl;

// Udagger*(A*U)
  zgemm3m_("c", "n", ij, ij, ij, 1.0, phase_factor, ij, intermediate, ij, 0.0, complex_moe, ij);

//end of phase factor
//there has to be a better way to do this, but this is temporary anyway
  unique_ptr<complex<double>[]> out(new complex<double>[ij*ij]);
  for (int i=0;i<(ij*ij);i++) out[i] = complex_moe[i];
#endif

  jop_ =  make_shared<Htilde>(ref_, 0, 0, mo1e, move(mo2e));

}
