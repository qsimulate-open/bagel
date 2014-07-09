#ifndef __BTAS_TARRAY_H
#define __BTAS_TARRAY_H 1

#include <btas/tensor.h>

namespace btas {

  /// Fixed-rank version of TArray
  template<typename _T,
           unsigned long _N,
           CBLAS_ORDER _Order = CblasRowMajor,
           class _Container = DEFAULT::storage<_T>>
  using TArray = Tensor<_T, RangeNd<_Order, std::array<long, _N> >, _Container >;

} // namespace btas

#endif // __BTAS_TARRAY_H
