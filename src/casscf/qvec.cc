//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#include <src/casscf/qvec.h>

using namespace std;

Qvec::Qvec(const int n, const int m, shared_ptr<DensityFit> df, shared_ptr<Coeff> coeff, const size_t nclosed, shared_ptr<FCI> fci)
 : QFile(n,m) {

  const int nbasis_ = df->nbasis();
  // target dimensions
  const int nbas = n;
  const int nact_ = m;
  const int nx_ = df->naux();

  shared_ptr<RDM<2> > rdm = fci->rdm2_av();

  // three index integral without compression.
  const double* const buf = df->data_3index();
  const double* const jm = df->data_2index();
  // area for the full transformed quantity.
  unique_ptr<double[]> full1(new double[nact_*nact_*nx_]);
  unique_ptr<double[]> full2(new double[nact_*nact_*nx_]);

  // half-transformed integrals
  const double* const half = fci->jop()->mo2e_1ext_ptr(); 

  // second transform third index: (D|ii) = (D|ix) Cxi
  dgemm_("N", "N", nx_*nact_, nact_, nbasis_, 1.0, half, nx_*nact_, coeff->data()+nclosed*nbasis_, nbasis_, 0.0, full1.get(), nx_*nact_); 

  // then multiply S^-1/2 twice from the front. S^-1/2 (D|ii)
  dgemm_("N", "N", nx_, nact_*nact_, nx_, 1.0, jm, nx_, full1.get(), nx_, 0.0, full2.get(), nx_); 
  dgemm_("T", "N", nx_, nact_*nact_, nx_, 1.0, jm, nx_, full2.get(), nx_, 0.0, full1.get(), nx_); 
  // then multiply 2rdm from the back
  dgemm_("N", "N", nx_, nact_*nact_, nact_*nact_, 1.0, full1.get(), nx_, rdm->data(), nact_*nact_, 0.0, full2.get(), nx_);

  // assemble to full1
  dgemm_("T", "N", nbasis_, nact_, nx_*nact_, 1.0, half, nx_*nact_, full2.get(), nx_*nact_, 0.0, full1.get(), nbasis_); 
  // index transform
  dgemm_("T", "N", nbas, nact_, nbasis_, 1.0, coeff->data(), nbasis_, full1.get(), nbasis_, 0.0, data(), nbas);

}
