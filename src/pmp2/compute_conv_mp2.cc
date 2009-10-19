/*
 * compute_conv_m2.cc
 *
 *  Created on: Oct 19, 2009
 *      Author: shiozaki
 */

#include <cstring>
#include <iostream>
#include <algorithm>
#include <src/pmp2/pmp2.h>
#include <src/pscf/poverlap.h>
#include <src/pscf/ptildex.h>
#include <src/macros.h>
#include <src/scf/scf_macros.h>
#include <src/slater/slaterbatch.h>
#include <src/util/paircompfile.h>
#include <src/util/pcompcabsfile.h>
#include <src/util/pmofile.h>

using namespace std;

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
    const double* eb = eig_ + (kb + K) * nbasis_;

    for (int kj = -K; kj != maxK1; ++kj) {
      const double* ej = eig_ + (kj + K) * nbasis_;

      for (int ka = -K; ka != maxK1; ++ka, ++nbja) {
        const double* ea = eig_ + (ka + K) * nbasis_;

        int ki = ka + kb - kj;
        if (ki < -K) ki += K * 2;
        else if (ki >= K) ki -= K * 2;
        const double* ei = eig_ + (ki + K) * nbasis_;

        assert(K == 0 || (ki < K && ki >= -K));

        {
          const long tmp2 = kb + K + KK * (kj + K + KK * (ka + K));
          eri_ii_pp_->get_block(noovv_ * nbja, noovv_, workoovv1);
          eri_ii_pp_->get_block(noovv_ * tmp2, noovv_, workoovv2);
        }

        #pragma omp parallel for reduction(+:energy)
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

