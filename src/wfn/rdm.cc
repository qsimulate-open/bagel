//
// Author : Toru Shiozaki
// Date   : October 2012
//

#include <src/wfn/rdm.h>

using namespace bagel;
using namespace std;

RDM_base::RDM_base(const int n, const int rank) : norb_(n), rank_(rank) {
  assert(rank > 0);
  dim_ = 1;
  for (int i = 0; i != rank; ++i) dim_ *= n;
  data_ = unique_ptr<double[]>(new double[dim_*dim_]);
}


RDM_base::RDM_base(const RDM_base& o) : norb_(o.norb_), dim_(o.dim_), rank_(o.rank_) {
  data_ = unique_ptr<double[]>(new double[dim_*dim_]);
  copy(o.data(), o.data()+dim_*dim_, data());
}
