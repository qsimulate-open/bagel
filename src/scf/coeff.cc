//
// Author: Toru Shiozaki
// Date  : May 2009
//

#include <src/scf/coeff.h>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <src/scf/f77.h>


Coeff::Coeff(const Matrix1e& inp) : Matrix1e(inp.geom()) {

  const int unit = 1;
  const int size = nbasis_ * nbasis_;
  const int ndim_ = inp.ndim();
  const int mdim_ = inp.mdim();
  dcopy_(&size, inp.data(), &unit, data_, &unit); 

}


Coeff::~Coeff() {

}


Matrix1e Coeff::form_density_rhf() const {
  const int nocc = geom_->nocc() / 2;
  assert(geom_->nocc() % 2 == 0);
  const double one = 1.0;
  const double zero = 0.0;

  Matrix1e out(geom_);
  double* out_data = out.data();

  dgemm_("N", "T", &nbasis_, &nbasis_, &nocc, &one, data_, &nbasis_, data_, &nbasis_, &zero, out_data, &nbasis_); 

  return out;
}
 

