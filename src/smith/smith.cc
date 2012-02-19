//
// Author : Toru Shiozaki
// Date   : Feb 2012
//

#include <src/smith/storage.h>
#include <src/smith/tensor.h>
#include <iostream>
#include <src/util/f77.h>
#include <src/smith/moint.h>
#include <src/wfn/reference.h>

using namespace SMITH;
using namespace std;

void a(shared_ptr<Reference> r){
  const int max = 7;
  IndexRange closed(10, max);
  IndexRange acc(7, max);
  IndexRange virt(20, max);

#if 0
  closed.print();
  acc.print();
  virt.print();
#endif

  vector<IndexRange> o;
  o.push_back(closed);
  o.push_back(virt);
  o.push_back(closed);
  o.push_back(virt);

  MOInt<Storage_Incore> a(r, o);
}
