#ifndef __BTAS_TENSOR_ITERATOR_WRAPPER_H
#define __BTAS_TENSOR_ITERATOR_WRAPPER_H 1

#include <type_traits>

#include <btas/tensor_traits.h>

namespace btas {

template<bool _HasData>
struct tensor_iterator_wrapper
{
   template<class _Tensor>
   static auto begin(_Tensor& x) -> decltype(x.begin()) { return x.begin(); }

   template<class _Tensor>
   static auto end  (_Tensor& x) -> decltype(x.end  ()) { return x.end  (); }
};

template<>
struct tensor_iterator_wrapper<true>
{
   template<class _Tensor>
   static auto begin(_Tensor& x) -> decltype(x.data()) { return x.data(); }

   template<class _Tensor>
   static auto end  (_Tensor& x) -> decltype(x.data()) { return x.data()+x.size(); }
};

/// \return wrapped iterator of Tensor to the first
template<class _Tensor, class = typename std::enable_if<is_tensor<_Tensor>::value>::type>
auto tbegin (_Tensor& x) -> decltype(tensor_iterator_wrapper<has_data<_Tensor>::value>::begin(x))
{
   return tensor_iterator_wrapper<has_data<_Tensor>::value>::begin(x);
}

/// \return wrapped iterator of Tensor to the last
template<class _Tensor, class = typename std::enable_if<is_tensor<_Tensor>::value>::type>
auto tend (_Tensor& x) -> decltype(tensor_iterator_wrapper<has_data<_Tensor>::value>::end(x))
{
   return tensor_iterator_wrapper<has_data<_Tensor>::value>::end(x);
}

} // namespace btas

#endif // __BTAS_TENSOR_ITERATOR_WRAPPER_H
