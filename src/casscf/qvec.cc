//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#include <src/casscf/qvec.h>

using namespace std;

Qvec::Qvec(const int n, const int m, shared_ptr<DensityFit> df, shared_ptr<Coeff> coeff, const size_t nclosed, shared_ptr<FCI> fci)
 : QFile(n,m) {

  const int nbasis = df->nbasis();

  shared_ptr<DF_Half> half = fci->jop()->mo2e_1ext();

  shared_ptr<DF_Full> full = half->compute_second_transform(coeff->data()+nclosed*nbasis, m)->apply_J()->apply_J(); 

  shared_ptr<DF_Full> prdm = full->apply_2rdm(fci->rdm2_av()->data());

  unique_ptr<double[]> tmp(new double[nbasis*m]);
  half->form_2index(tmp, prdm);
  dgemm_("T", "N", n, m, nbasis, 1.0, coeff->data(), nbasis, tmp.get(), nbasis, 0.0, data(), n);

}
