//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#include <src/wfn/reference.h>

using namespace std;

Reference::Reference(SCF& a) : geom_(a.geom()), coeff_(a.coeff()), hcore_(a.hcore()), schwarz_(a.schwarz()) {

}


