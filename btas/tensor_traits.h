#ifndef __BTAS_TENSOR_TRAITS_H
#define __BTAS_TENSOR_TRAITS_H 1

#include <iterator>
#include <type_traits>

#include <btas/index_traits.h>
#include <btas/range_traits.h>

namespace btas {

/// test T has data() member
/// this will be used to detect whether or not the storage is consecutive
template<class T>
class has_data {
   /// true case
   template<class U>
   static auto __test(U* p) -> decltype(p->data(), std::true_type());
   /// false case
   template<class>
   static std::false_type __test(...);
public:
   static constexpr const bool value = std::is_same<std::true_type, decltype(__test<T>(0))>::value;
};

/// test T has range_type
template<class T>
class has_range_type {
   /// true case
   template<class U>
   static std::true_type __test(typename U::range_type*);
   /// false case
   template<class>
   static std::false_type __test(...);
public:
   static constexpr const bool value = std::is_same<std::true_type, decltype(__test<T>(0))>::value;
};

/// test T has range_type && is_boxrange<range_type>::value is true
template<class T>
class has_boxrange_range_type {
   /// true case
   template<class U>
   static std::true_type __test(typename U::range_type*,
                                typename std::enable_if<is_boxrange<typename U::range_type>::value, void*>::type = 0);
   /// false case
   template<class>
   static std::false_type __test(...);
public:
   static constexpr const bool value = std::is_same<std::true_type, decltype(__test<T>(0))>::value;
};

/// test T has storage_type
template<class T>
class has_storage_type {
   /// true case
   template<class U>
   static std::true_type __test(typename U::storage_type*);
   /// false case
   template<class>
   static std::false_type __test(...);
public:
   static constexpr const bool value = std::is_same<std::true_type, decltype(__test<T>(0))>::value;
};

/// checks _Tensor meets the TWG.Tensor concept requirements
/// checks only value_type, range_type, storage_type, and rank() member TODO check the rest
template<class _Tensor>
class is_tensor {
public:
   static constexpr const bool
   value = has_value_type<_Tensor>::value & has_range_type<_Tensor>::value &
           has_storage_type<_Tensor>::value & has_rank<_Tensor>::value;
};

/// checks _Tensor meets the TWG.BoxTensor concept requirements
template<class _Tensor>
class is_boxtensor {
public:
   static constexpr const bool
   value = is_tensor<_Tensor>::value && has_boxrange_range_type<_Tensor>::value;
};

/// checks _Tensor meets the TWG.BoxTensor concept requirements
template<class _Tensor>
class boxtensor_storage_order {
public:
   enum {row_major = boxrange_iteration_order<typename _Tensor::range_type>::row_major,
         other = boxrange_iteration_order<typename _Tensor::range_type>::other,
         column_major = boxrange_iteration_order<typename _Tensor::range_type>::column_major
   };
   static constexpr const int
   value = boxrange_iteration_order<typename _Tensor::range_type>::value;
};

} // namespace btas

#endif // __BTAS_TENSOR_TRAITS_H
