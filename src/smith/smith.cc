//
// Author : Toru Shiozaki
// Date   : Feb 2012
//

#include <src/smith/smith.h>

using namespace std;
using namespace SMITH;

SMITH_info::SMITH_info() {
  maxiter_ = 5;
  thresh_residual_ = 1.0e-8;
}


SMITH_info::~SMITH_info() {

}
