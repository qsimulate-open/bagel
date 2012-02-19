//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#include <src/wfn/reference.h>

using namespace std;

Reference::Reference(shared_ptr<Geometry> g, shared_ptr<Coeff> c,  shared_ptr<Hcore> h, const vector<double>& s,
                     const int& ncl, const int& nac, const int& nvi)
 : geom_(g), coeff_(c), hcore_(h), schwarz_(s), nclosed_(ncl), nact_(nac), nvirt_(nvi) {

}

