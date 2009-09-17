//
// Author : Toru Shiozaki
// Date   : August 2009
//

#include <src/pmp2/pmp2.h>
#include <src/macros.h>
#include <src/scf/scf_macros.h>
#include <cstring>
#include <algorithm>
#include <src/slater/slaterbatch.h>
#include <src/util/paircompfile.h>

typedef boost::shared_ptr<Atom> RefAtom;
typedef boost::shared_ptr<PGeometry> RefGeom;
typedef boost::shared_ptr<PCoeff> RefPCoeff;
typedef boost::shared_ptr<Shell> RefShell;

using namespace std;
using namespace boost;

PMP2::PMP2(const RefGeom g, const RefPCoeff co, const double* eg, const shared_ptr<PCompFile<ERIBatch> > fl)
 : geom_(g), coeff_(co), eig_(eg), ao_eri_(fl) {

  nfrc_ = geom_->nfrc() / 2;
  nocc_ = geom_->nocc() / 2;
  nocc_act_ = nocc_ - nfrc_;
  nbasis_ = geom_->nbasis();
  nvir_ = nbasis_ - nocc_;
  noovv_ = nocc_act_ * nocc_act_ * nvir_ * nvir_;

}


PMP2::~PMP2() {

}


void PMP2::compute() {

  // AO ERI has been computed in the SCF class.

  // Fully transform aa/ii integrals and dump them to disk (... forcus is on MP2-R12).
  eri_aa_ii_ = ao_eri_->mo_transform(coeff_,
                                     nfrc_, nocc_,
                                     nfrc_, nocc_,
                                     nocc_, nbasis_,
                                     nocc_, nbasis_);

  // Compute the conventional MP2 contribution
  compute_conv_mp2();

  // Calculate Yukawa potential integrals
//shared_ptr<PCompFile<SlaterBatch> > slater(new PCompFile<SlaterBatch>(geom_));
  shared_ptr<PairCompFile<SlaterBatch> > slater_and_yukawa(new PairCompFile<SlaterBatch>(geom_));

}


void PMP2::compute_conv_mp2() {
  cout << "  Now computing mp2 energy" << endl << endl;
  const int K = eri_aa_ii_->K();
  const int maxK1 = max(K, 1);
  const int KK = K + K;

  complex<double>* workoovv1 = new complex<double>[noovv_];
  complex<double>* workoovv2 = new complex<double>[noovv_];

  double energy = 0.0;
  for (int kb = -K, nbja = 0; kb != maxK1; ++kb) {
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
          eri_aa_ii_->get_block(noovv_ * nbja, noovv_, workoovv1);
          eri_aa_ii_->get_block(noovv_ * tmp2, noovv_, workoovv2);
        }

        int cnt1 = 0;
        const int ov = nocc_act_ * nvir_;
        for (int i = nfrc_, xi = 0; i != nocc_; ++i, ++xi) {
          for (int a = nocc_, xa = 0; a != nbasis_; ++a, ++xa) {
            for (int j = nfrc_, xj = 0; j != nocc_; ++j, ++xj) {
              int cnt2 = xa + nvir_ * (xj + ov * xi);

              for (int b = nocc_; b != nbasis_; ++b, ++cnt1, cnt2 += ov) {
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

