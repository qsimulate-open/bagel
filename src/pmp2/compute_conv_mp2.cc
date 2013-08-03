//
// BAGEL - Parallel electron correlation program.
// Filename: compute_conv_mp2.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
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


#include <iostream>
#include <src/pmp2/pmp2.h>

using namespace std;
using namespace bagel;

#ifdef HAVE_LIBSLATER

void PMP2::compute_conv_mp2() {
  cout << "  Now computing mp2 energy" << endl << endl;
  const int K = eri_ii_pp_->K();
  const int maxK1 = max(K, 1);
  const int KK = K + K;

  complex<double>* workoovv1 = new complex<double>[noovv_];
  complex<double>* workoovv2 = new complex<double>[noovv_];

  double energy = 0.0;
  int nbja = 0;
  for (int kb = -K; kb < maxK1; ++kb) {
    vector<double>::const_iterator eb = eig_.begin() + (kb + K) * nbasis_;

    for (int kj = -K; kj != maxK1; ++kj) {
      vector<double>::const_iterator ej = eig_.begin() + (kj + K) * nbasis_;

      for (int ka = -K; ka != maxK1; ++ka, ++nbja) {
        vector<double>::const_iterator ea = eig_.begin() + (ka + K) * nbasis_;

        int ki = ka + kb - kj;
        if (ki < -K) ki += K * 2;
        else if (ki >= K) ki -= K * 2;
        vector<double>::const_iterator ei = eig_.begin() + (ki + K) * nbasis_;

        assert(K == 0 || (ki < K && ki >= -K));

        {
          const long tmp2 = kb + K + KK * (kj + K + KK * (ka + K));
          eri_ii_pp_->get_block(noovv_ * nbja, noovv_, workoovv1);
          eri_ii_pp_->get_block(noovv_ * tmp2, noovv_, workoovv2);
        }

//      #pragma omp parallel for reduction(+:energy)
        for (int i = nfrc_; i < nocc_; ++i) {
          const int xi = i - nfrc_;
          for (int j = nfrc_, xj = 0; j != nocc_; ++j, ++xj) {
            const int ijoffset = nbasis_ * (xj + nocc_act_ * xi);
            for (int a = nocc_; a != nbasis_; ++a) {
              for (int b = nocc_; b != nbasis_; ++b) {
                int cnt1 = b + nbasis_ * (a + ijoffset);
                int cnt2 = a + nbasis_ * (b + ijoffset);
                const complex<double> v1 = workoovv1[cnt1];
                const complex<double> v2 = workoovv2[cnt2];
                energy += (real(v1 * conj(v1 + v1 - v2))) / (ei[i] + ej[j] - ea[a] - eb[b]);
              }
            }
          }
        }
      }
    }
  }
  cout << "  MP2 energy: " << setprecision(10) << energy / (max(pow(2.0 * K, 3.0), 1.0))  << endl << endl;
  delete[] workoovv1;
  delete[] workoovv2;

}

#endif
