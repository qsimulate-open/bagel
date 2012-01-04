//
// Author: Toru Shiozaki
// Date  : May 2009
//

#include <src/scf/coeff.h>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <src/util/f77.h>


Coeff::Coeff(const Matrix1e& inp) : Matrix1e(inp.geom()) {

  const int ndim_ = inp.ndim();
  const int mdim_ = inp.mdim();
  dcopy_(nbasis_*nbasis_, inp.data(), 1, data(), 1); 

}


Coeff::~Coeff() {

}


Matrix1e Coeff::form_density_rhf() const {
  const int nocc = geom_->nocc() / 2;
  assert(geom_->nocc() % 2 == 0);

  Matrix1e out(geom_);
  double* out_data = out.data();

  dgemm_("N", "T", nbasis_, nbasis_, nocc, 1.0, data(), nbasis_, data(), nbasis_, 0.0, out_data, nbasis_); 

  return out;
}
 

Matrix1e Coeff::form_core_density_rhf() const {
  const int nocc = geom_->nfrc() / 2;
  assert(geom_->nfrc() % 2 == 0);

  Matrix1e out(geom_);
  double* out_data = out.data();

  dgemm_("N", "T", nbasis_, nbasis_, nocc, 1.0, data(), nbasis_, data(), nbasis_, 0.0, out_data, nbasis_); 

  return out;
}
 

