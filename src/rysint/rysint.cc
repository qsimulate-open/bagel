//
// Author: Toru Shiozaki
// Date  : May 2009
//

#include <src/rysint/rysint.h>

using namespace std;

RysInt::RysInt(const vector<boost::shared_ptr<Shell> > info)
 : basisinfo_(info), spherical_(info.front()->spherical()), sort_(info.front()->spherical()) {

}


RysInt::~RysInt() {

}


