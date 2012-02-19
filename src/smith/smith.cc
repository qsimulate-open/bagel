//
// Author : Toru Shiozaki
// Date   : Feb 2012
//

#include <src/smith/storage.h>
#include <src/smith/tensor.h>
#include <iostream>
#include <src/util/f77.h>

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

  vector<size_t> p(4,0);
  unique_ptr<double[]> tmp = t.get_block(p);

  const size_t size = t.get_size(p);

  fill(&tmp[0], &tmp[size], 1.0);
  t.put_block(p,tmp);

  cout << ddot_(size, tmp, 1, tmp, 1) << endl;
}
