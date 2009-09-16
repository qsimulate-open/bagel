//
// Author : Toru Shiozaki
// Date   : June 2009
//

#include <src/slater/f77.h>
#include <src/slater/srootlist.h>

SRootList::SRootList() {

  srfunc[1] = &root1_;
  srfunc[2] = &root2_;
  srfunc[3] = &root3_;
  srfunc[4] = &root4_;
  srfunc[5] = &root5_;
  srfunc[6] = &root6_;
  srfunc[7] = &root7_;
  srfunc[8] = &root8_;
  srfunc[9] = &root9_;
}


SRootList::~SRootList() {

}
