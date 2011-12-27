//
// Author: Toru Shiozaki
// Date  : May 2009
//

#include <src/scf/fock.h>
#include <src/scf/f77.h>
#include <src/scf/symmat.h>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <algorithm>

typedef std::shared_ptr<Geometry> RefGeometry;
typedef std::shared_ptr<Hcore> RefHcore;
typedef std::shared_ptr<Matrix1e> RefAODensity;
typedef std::shared_ptr<Shell> RefShell;
typedef std::shared_ptr<Fock> RefFock;

using namespace std;

Fock::Fock(const RefGeometry geom, const RefFock previous, const RefAODensity den, const vector<double>& schwarz)
 : Matrix1e(geom), previous_(previous), density_(den), schwarz_(schwarz) {

  schwarz_thresh_ = geom->schwarz_thresh();

  init(); // zero here

  fock_two_electron_part();

  for (int i = 0; i != ndim_; ++i) {
    const int diagpos = i * nbasis_ + i;
    data_[diagpos] *= 2.0;
  }
  symmetrize();

  const int nirrep = geom_->nirrep();
  const int unit = 1;
  const double one = 1.0;
  const int size = nbasis_ * nbasis_;
  if (nirrep != 1) {
    Matrix1e intermediate(geom);
    for (int i = 1; i != nirrep; ++i) {
      SymMat symm(geom_, i);
      Matrix1e tmp = symm % (*this) * symm;
      intermediate += tmp;
    }
    double* idata = intermediate.data();
    daxpy_(&size, &one, idata, &unit, data_, &unit);
    const double scale = 1.0 / nirrep;
    dscal_(&size, &scale, data_, &unit);
  }
  const double* previous_data = previous_->data();
  daxpy_(&size, &one, previous_data, &unit, data_, &unit); 

  symmetrize();
}


Fock::Fock(const RefGeometry geom, const RefHcore hcore)
 : Matrix1e(geom) {

  init(); // zero here
  const int unit = 1;
  const int size = nbasis_ * nbasis_; 
  const double one = 1.0;
  const double* hcore_data = hcore->data();
  daxpy_(&size, &one, hcore_data, &unit, data_, &unit); 

  symmetrize();
}


Fock::~Fock() {

}


void Fock::computebatch(const vector<RefShell>& input, const int offsetb0, const int offsetb1, const int ndim) {

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


