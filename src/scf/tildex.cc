//
// Author: Toru Shiozaki
// Date  : May 2009
//

#include <src/scf/tildex.h>
#include <src/scf/f77.h>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <vector>

using namespace std;

#define USE_CANONICAL

TildeX::TildeX(const std::shared_ptr<Overlap> olp, const double thresh) : Matrix1e(olp->geom())  {

  // Use canonical orthogonalization (Szabo pp.144)
  nbasis_ = geom_->nbasis();
  ndim_ = nbasis_;
  const int size = nbasis_ * nbasis_;

  double* eig = new double[nbasis_];
  const int lwork = 5 * nbasis_;
  double* work = new double[lwork];
  const double* S = olp->data();

  const int unit = 1;
  dcopy_(&size, S, &unit, data_, &unit);

  int info;
  dsyev_("V", "L", &ndim_, data_, &ndim_, eig, work, &lwork, &info); 
  assert(info == 0);
  delete[] work;
  const double largest = fabs(eig[ndim_ - 1]);

  // counting how many orbital must be deleted owing to the linear dependency
  int cnt = 0;
  for (int i = 0; i != ndim_; ++i) {
    if (fabs(eig[i]) < largest * thresh) ++cnt;
    else break;
  }
  if (cnt != 0) 
    cout << "  Caution: ignored " << cnt << " orbital" << (cnt == 1 ? "" : "s") << " in the orthogonalization." << endl << endl;

  for (int i = cnt; i != ndim_; ++i) {
#ifdef USE_CANONICAL
    const double scale = 1.0 / ::sqrt(eig[i]);
#else
    const double scale = 1.0 / ::sqrt(::sqrt(eig[i]));
#endif
    const int offset = i * ndim_;
    for (int j = 0; j != ndim_; ++j) {
      data_[j + offset] *= scale; 
    }
  }
  delete[] eig;
#ifdef USE_CANONICAL
  mdim_ = ndim_ - cnt;
  if (cnt != 0) { 
    for (int i = 0; i != mdim_; ++i) {
      dcopy_(&ndim_, &data_[(i + cnt) * ndim_], &unit, &data_[i * ndim_], &unit); 
    }
  }
#else
  mdim_ = ndim_;
  double* tmp = new double[size];
  const int m = ndim_ - cnt;
  const double one = 1.0;
  const double zero = 0.0;
  dgemm_("N", "T", &ndim_, &ndim_, &m, &one, &data_[cnt * ndim_], &ndim_, &data_[cnt * ndim_], &ndim_, &zero, tmp, &ndim_); 
  dcopy_(&size, tmp, &unit, data_, &unit);
  delete[] tmp;
#endif
  

}


TildeX::~TildeX() {

}


