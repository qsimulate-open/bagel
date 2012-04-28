//
// Author : Toru Shiozaki
// Date   : April 2012
//

#include <src/mp2/f12mat.h>

using namespace std;

void F12Mat::symmetrize(const bool braket) {
  // first electron 1-2 symmetry
  for (int i = 0; i != nocc_; ++i) {
    for (int j = 0; j != nocc_; ++j) {
      for (int k = 0; k != nocc_; ++k) {
        for (int l = 0; l != nocc_; ++l) {
          if ((i<<8) + j > (k<<8) + l) continue; 
          const double ijkl = (data(l,k,j,i) + data(j,i,l,k)) * 0.5;
          data(l,k,j,i) = ijkl;
          data(j,i,l,k) = ijkl;
        }
      }
    }
  }
  // bra-ket symmetry 
  if (braket) {
    for (int i = 0; i != nocc_; ++i) {
      for (int j = 0; j != nocc_; ++j) {
        for (int k = 0; k != nocc_; ++k) {
          for (int l = 0; l != nocc_; ++l) {
            if ((i<<8) + j > (k<<8) + l) continue; 
            const double ijkl = (data(l,k,j,i) + data(j,i,l,k)) * 0.5;
            data(l,k,j,i) = ijkl;
            data(k,l,i,j) = ijkl;
          }
        }
      }
    }
  }

}
