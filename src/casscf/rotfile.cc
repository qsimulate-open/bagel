//
// Author : Toru Shiozaki
// Date   : Dec 2011
//

#include <src/casscf/rotfile.h>

using namespace std;


shared_ptr<Matrix1e> RotFile::unpack(shared_ptr<Geometry> geom) const {

  const int nocc_ = nclosed_ + nact_;
  const int nbasis_ = nclosed_ + nact_ + nvirt_; 
  shared_ptr<Matrix1e> out(new Matrix1e(geom, nbasis_, nbasis_));

  for (int i = 0; i != nact_; ++i) {
    for (int j = 0; j != nvirt_;   ++j) {
      out->element(j+nocc_, i+nclosed_) = ele_va(j, i);
    }
    for (int j = 0; j != nclosed_; ++j) {
      out->element(i+nclosed_, j) = ele_ca(j, i);
    }
  }
  for (int i = 0; i != nclosed_; ++i) {
    for (int j = 0; j != nvirt_;   ++j) {
      out->element(j+nocc_, i) = ele_vc(j, i);
    }
  }
  for (int i = 0; i != nbasis_; ++i) {
    for (int j = 0; j <= i; ++j) {
      out->element(j, i) = -out->element(i, j);
    }
  }
  return out;
}

void RotFile::print() const {
  if (nact_ && nclosed_) {
    cout << " printing closed-active block" << endl;
    for (int i = 0; i != nact_; ++i) {
      for (int j = 0; j != nclosed_; ++j) {
        cout << setw(10) << setprecision(6) << ele_ca(j,i); 
      }
      cout << endl;
    }
  }
  if (nact_ && nvirt_) {
    cout << " printing virtual-active block" << endl;
    for (int i = 0; i != nact_; ++i) {
      for (int j = 0; j != nvirt_; ++j) {
        cout << setw(10) << setprecision(6) << ele_va(j,i); 
      }
      cout << endl;
    }
  }
  if (nclosed_ && nvirt_) {
    cout << " printing virtual-closed block" << endl;
    for (int i = 0; i != nclosed_; ++i) {
      for (int j = 0; j != nvirt_; ++j) {
        cout << setw(10) << setprecision(6) << ele_vc(j,i); 
      }
      cout << endl;
    }
  }
}
