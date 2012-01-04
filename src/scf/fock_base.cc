//
// Author: Toru Shiozaki
// Date  : May 2009
//

#include <src/scf/fock.h>
#include <src/util/f77.h>
#include <src/scf/symmat.h>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <algorithm>

typedef std::shared_ptr<Geometry> RefGeometry;
typedef std::shared_ptr<Hcore> RefHcore;
typedef std::shared_ptr<Matrix1e> RefAODensity;
typedef std::shared_ptr<Shell> RefShell;
typedef std::shared_ptr<Fock_base> RefFock_base;

using namespace std;

Fock_base::Fock_base(const RefGeometry geom, const RefFock_base previous, const RefAODensity den, const vector<double>& schwarz)
 : Matrix1e(geom), previous_(previous), density_(den), schwarz_(schwarz) {

  schwarz_thresh_ = geom->schwarz_thresh();

  init(); // zero here
}


void Fock_base::fock_one_electron_part() {

  const int nirrep = geom_->nirrep();
  const int unit = 1;
  const double one = 1.0;
  const int size = nbasis_ * nbasis_;
  if (nirrep != 1) {
    Matrix1e intermediate(geom_);
    for (int i = 1; i != nirrep; ++i) {
      SymMat symm(geom_, i);
      Matrix1e tmp = symm % (*this) * symm;
      intermediate += tmp;
    }
    double* idata = intermediate.data();
    daxpy_(size, 1.0, idata, 1, data(), 1);
    dscal_(size, 1.0/nirrep, data(), 1);
  }
  const double* previous_data = previous_->data();
  daxpy_(size, 1.0, previous_data, 1, data(), 1); 

  symmetrize();
}


Fock_base::Fock_base(const RefGeometry geom, const RefHcore hcore)
 : Matrix1e(geom) {

  init(); // zero here
  daxpy_(nbasis_*nbasis_, 1.0, hcore->data(), 1, data(), 1); 

  symmetrize();
}


Fock_base::~Fock_base() {

}


void Fock_base::computebatch(const vector<RefShell>& input, const int offsetb0, const int offsetb1, const int ndim) {

  // input = [b1, b0]
  assert(input.size() == 2);
  const int dimb0 = input[1]->nbasis(); 
  const int dimb1 = input[0]->nbasis(); 

  for (int i = offsetb0; i != dimb0 + offsetb0; ++i) {
    for (int j = offsetb1; j != dimb1 + offsetb1; ++j) {
      data_[i * ndim + j] = 0.0; 
    }
  }
}


