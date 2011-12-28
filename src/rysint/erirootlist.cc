//
// Author : Toru Shiozaki
// Date   : June 2009
//

#include <src/rysint/erirootlist.h>
#include <src/rysint/f77.h>

ERIRootList::ERIRootList() {

  rfunc[1] = &eriroot1_;
  rfunc[2] = &eriroot2_;
  rfunc[3] = &eriroot3_;
  rfunc[4] = &eriroot4_;
  rfunc[5] = &eriroot5_;
  rfunc[6] = &eriroot6_;
  rfunc[7] = &eriroot7_;
  rfunc[8] = &eriroot8_;
  rfunc[9] = &eriroot9_;
}


ERIRootList::~ERIRootList() {

}
