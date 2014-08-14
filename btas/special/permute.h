#ifndef __BTAS_TARRAY_PERMUTE_H
#define __BTAS_TARRAY_PERMUTE_H 1

#include <array>
#include <type_traits>

#include <btas/types.h>
#include <btas/tarray.h>

#include <btas/generic/permute.h>
#include <btas/special/reindex.h>

namespace btas {

template<typename _T, size_type _N>
std::array<_T, _N> __permute_index (const std::array<_T, _N>& x, const std::array<size_type, _N>& index)
{
   std::array<_T, _N> y;
   for (size_type i = 0; i < _N; ++i)
   {
      y[i] = x[index[i]];
   }
   return y;
}

template<typename _T, size_type _N>
void permute (const TArray<_T, _N>& x, const std::array<size_type, _N>& index, TArray<_T, _N>& y)
{
   y.resize(__permute_index(x.shape(), index));
   Reindex(x.data(), y.data(), __permute_index(x.stride(), index), y.shape());
}

} // namespace btas

#endif // __BTAS_TARRAY_PERMUTE_H
