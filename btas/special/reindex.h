#ifndef __BTAS_TARRAY_REINDEX_H
#define __BTAS_TARRAY_REINDEX_H 1

#include <array>
#include <type_traits>

namespace btas {

/// NDloop class for Reindex
template<size_type _I, size_type _N>
struct __NDloop_reindex
{
   /// loop upon construction
   /// NOTE: pX and pY are passed as a reference of pointer to the next loop
   /// NOTE: on the other hand, addrX is passed as a value so that offset position (by addrX) is kept in this scope
   template<typename _T, class = typename std::enable_if<(_I < _N)>::type>
   __NDloop_reindex (const _T*& pX, _T*& pY, size_type addrX, const std::array<size_type, _N>& strX, const std::array<size_type, _N>& shapeY)
   {
      for (size_type i = 0; i < shapeY[_I-1]; ++i)
      {
         __NDloop_reindex<_I+1, _N> loop(pX, pY, addrX+i*strX[_I-1], strX, shapeY);
      }
   }
};

/// NDloop class for Reindex, specialized for the last index
template<size_type _N>
struct __NDloop_reindex<_N, _N>
{
   /// loop upon construction
   template<typename _T>
   __NDloop_reindex (const _T*& pX, _T*& pY, size_type addrX, const std::array<size_type, _N>& strX, const std::array<size_type, _N>& shapeY)
   {
      for (size_type i = 0; i < shapeY[_N-1]; ++i, ++pY)
      {
         *pY = pX[addrX+i*strX[_N-1]];
      }
   }
};

/// reindex (i.e. permute) for "any-rank" tensor
/// multiple loop is expanded at compile time
/// FIXME: how slower than explicit looping?
/// if considerably slower, should be specialized for small ranks (_N = 1 ~ 8?)
template<typename _T, size_type _N>
void Reindex (const _T* pX, _T* pY, const std::array<size_type, _N>& strX, const std::array<size_type, _N>& shapeY)
{
   __NDloop_reindex<1, _N> loop(pX, pY, 0, strX, shapeY);
}

} // namespace btas

#endif // __BTAS_TARRAY_REINDEX_H
