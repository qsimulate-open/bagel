#ifndef __BTAS_UTIL_DOT_H
#define __BTAS_UTIL_DOT_H 1

#include <cassert>
#include <array>
#include <type_traits>

namespace btas {

template<class _Vector>
typename _Vector::value_type dot (const _Vector& x, const _Vector& y)
{
   assert(x.size() == y.size());
   auto ix = x.begin();
   auto iy = y.begin();
   typename _Vector::value_type value = (*ix) * (*iy);
   ++ix;
   ++iy;
   for(; ix != x.end(); ++ix, ++iy) {
      value += (*ix) * (*iy);
   }
   return value;
}

//
// small overhead version for std::array
//

template<typename _T, unsigned long _I, unsigned long _N>
struct __dot_helper
{
   static _T multiply (const std::array<_T, _N>& x, const std::array<_T, _N>& y)
   {
      return x[_I-1]*y[_I-1]+__dot_helper<_T, _I+1, _N>::multiply(x, y);
   }
};

template<typename _T, unsigned long _N>
struct __dot_helper<_T, _N, _N>
{
   static _T multiply (const std::array<_T, _N>& x, const std::array<_T, _N>& y)
   {
      return x[_N-1]*y[_N-1];
   }
};

template<typename _T, unsigned long _N, class = typename std::enable_if<(_N > 0)>::type>
_T dot (const std::array<_T, _N>& x, const std::array<_T, _N>& y)
{
   return __dot_helper<_T, 1, _N>::multiply(x, y);
}

} // namespace btas

#endif // __BTAS_UTIL_DOT_H
