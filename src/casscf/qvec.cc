//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#include <src/casscf/qvec.h>

using namespace std;

Qvec::Qvec(const int n, const int m, shared_ptr<DensityFit> df, shared_ptr<Coeff> coeff, const size_t nclosed, shared_ptr<RDM<2> > rdm)
 : QFile(n,m) {

  const int nbasis_ = n;
  assert(nbasis_ == df->nbasis());
  const int nact_ = m;
  const int nx_ = df->naux();

  // three index integral without compression.
  const double* const buf = df->data_3index();
  const double* const jm = df->data_2index();
  // area for the first half transformation.
  unique_ptr<double[]> half(new double[nact_*nbasis_*nx_]);
  // area for the full transformed quantity.
  unique_ptr<double[]> full1(new double[nact_*nact_*nx_]);
  unique_ptr<double[]> full2(new double[nact_*nact_*nx_]);

  // first transform second index: (D|ix) = (D|yx) Cyi
  for (size_t i = 0; i != nbasis_; ++i)
    dgemm_("N", "N", nx_, nact_, nbasis_, 1.0, buf+i*nbasis_*nx_, nx_, coeff->data()+nclosed*nbasis_, nbasis_, 0.0, half.get()+i*nx_*nact_, nx_);

  // second transform third index: (D|ii) = (D|ix) Cxi
  dgemm_("N", "N", nx_*nact_, nact_, nbasis_, 1.0, half.get(), nx_*nact_, coeff->data()+nclosed*nbasis_, nbasis_, 0.0, full1.get(), nx_*nact_); 

  // then multiply S^-1/2 twice from the front. S^-1/2 (D|ii)
  dgemm_("N", "N", nx_, nact_*nact_, nx_, 1.0, jm, nx_, full1.get(), nx_, 0.0, full2.get(), nx_); 
  dgemm_("T", "N", nx_, nact_*nact_, nx_, 1.0, jm, nx_, full2.get(), nx_, 0.0, full1.get(), nx_); 
  // then multiply 2rdm from the back
  dgemm_("N", "N", nx_, nact_*nact_, nact_*nact_, 1.0, full1.get(), nx_, rdm->data(), nact_*nact_, 0.0, full2.get(), nx_);

  // assemble to data 
  dgemm_("T", "N", nbasis_, nact_, nx_*nact_, 1.0, half.get(), nx_*nact_, full2.get(), nx_*nact_, 0.0, data(), nbasis_); 

}
