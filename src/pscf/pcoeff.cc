//
// Author : Toru Shiozaki
// Date   : July 2009
//

#include <src/pscf/pcoeff.h>
#include <src/pscf/f77.h>

typedef std::complex<double> Complex;

using namespace std;

PCoeff::PCoeff(const PMatrix1e& inp) : PMatrix1e(inp.geom()) {

  const int unit = 1;
  const int ndim_ = inp.ndim();
  const int mdim_ = inp.mdim();
  zcopy_(&totalsize_, inp.data()->front(), &unit, data()->front(), &unit); 

}

PCoeff::~PCoeff() {

}

PMatrix1e PCoeff::form_density_rhf(const bool return_ao) {
  const int nocc = geom_->nocc() / 2;
  assert(geom_->nocc() % 2 == 0);
  const Complex one(1.0, 0.0);
  const Complex zero(0.0, 0.0);

  // first, form the density matrix in k space. 
  PMatrix1e k_density(geom_);
  int kcount = 0;
  for (int k = -K(); k <= K(); ++k, ++kcount) { 
    const int koffset = kcount * blocksize_;
    zgemm_("N", "C", &nbasis_, &nbasis_, &nocc, &one, data()->pointer(koffset), &nbasis_, 
                                                      data()->pointer(koffset), &nbasis_,
                                     &zero, k_density.data()->pointer(koffset), &nbasis_); 
  }

  // back Fourier transform
  if (return_ao)
      return k_density.bft();
  else
      return k_density;
}
