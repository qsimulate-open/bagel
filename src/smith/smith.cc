//
// Author : Toru Shiozaki
// Date   : Feb 2012
//

#include <src/smith/storage.h>
#include <src/smith/tensor.h>
#include <iostream>

using namespace SMITH;
using namespace std;

void a() {
  const int max = 7;
  IndexRange closed(10, max);
  IndexRange acc(7, max);
  IndexRange virt(20, max);

#if 0
  closed.print();
  acc.print();
  virt.print();
#endif

  list<IndexRange> o;
  o.push_back(closed);
  o.push_back(closed);
  o.push_back(virt);
  o.push_back(virt);

  Tensor<Storage_Incore> t(o);
}
