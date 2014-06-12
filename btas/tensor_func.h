/*
 * tensor_func.h
 *
 *  Created on: Dec 30, 2013
 *      Author: evaleev
 */

#ifndef BTAS_TENSOR_FUNC_H_
#define BTAS_TENSOR_FUNC_H_


namespace btas {

  // Helper template for TensorViewOf
  template<typename _T>
  using Nref = typename std::remove_reference<_T>::type;

  // Maps Tensor     -> TensorView,
  //      TensorView -> TensorView
  // appropriately transferring constness of the storage, that is,
  // if _T is const, uses const _T::storage_type, otherwise just _T::storage_type
  template<typename _T>
  using TensorViewOf = TensorView<typename Nref<_T>::value_type,
                                  typename Nref<_T>::range_type,
                                  typename std::conditional<std::is_const<Nref<_T>>::value,
                                                            const typename Nref<_T>::storage_type,
                                                            typename Nref<_T>::storage_type
                                                           >::type>;

  template<typename _T,
           typename _Permutation>
  TensorViewOf<_T>
  permute( _T&& t,
           _Permutation p) {
      return TensorViewOf<_T>( permute(t.range(), p), t.storage() );
  }

  template<typename _T,
           typename _U>
  TensorViewOf<_T>
  permute( _T&& t,
           std::initializer_list<_U> p) {
      return TensorViewOf<_T>( permute(t.range(), p), t.storage() );
  }

  template <typename _T>
  TensorViewOf<_T>
  diag(_T&& T)
    {
    return TensorViewOf<_T>(diag(T.range()),T.storage());
    }

  template <typename _T, typename ArrayType>
  TensorViewOf<_T>
  tieIndex(_T&& T,
           const ArrayType& inds)
    {
    return TensorViewOf<_T>(tieIndex(T.range(),inds),T.storage());
    }

  template <typename _T, typename... _args>
  TensorViewOf<_T>
  tieIndex(_T&& T,
           size_t i0,
           const _args&... rest)
    {
    const auto size = 1 + sizeof...(rest);
    std::array<size_t,size> inds = { i0, static_cast<size_t>(rest)...};
    return TensorViewOf<_T>(tieIndex(T.range(),inds),T.storage());
    }

  template <typename _T>
  TensorViewOf<_T>
  group(_T&& T,
        size_t istart,
        size_t iend)
    {
    return TensorViewOf<_T>(group(T.range(),istart,iend),T.storage());
    }

  template <typename _T>
  TensorViewOf<_T>
  flatten(_T&& T)
    {
    return TensorViewOf<_T>(flatten(T.range()),T.storage());
    }

} // namespace btas


#endif /* BTAS_TENSOR_FUNC_H_ */
