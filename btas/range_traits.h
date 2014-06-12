#ifndef __BTAS_RANGE_TRAITS_H
#define __BTAS_RANGE_TRAITS_H 1

#include <iterator>
#include <type_traits>
#include <btas/index_traits.h>

namespace btas {

/// test T has rank() member
template<class T>
class has_rank {
   /// true case
   template<class U>
   static auto __test(U* p) -> decltype(p->rank(), std::true_type());
   /// false case
   template<class>
   static std::false_type __test(...);
public:
   static constexpr const bool value = std::is_same<std::true_type, decltype(__test<T>(0))>::value;
};

/// test T has index_type
template<class T>
class has_index_type {
   /// true case
   template<class U>
   static std::true_type __test(typename U::index_type*);
   /// false case
   template<class>
   static std::false_type __test(...);
public:
   static constexpr const bool value = std::is_same<std::true_type, decltype(__test<T>(0))>::value;
};

/// test T has ordinal_type
template<class T>
class has_ordinal_type {
   /// true case
   template<class U>
   static std::true_type __test(typename U::ordinal_type*);
   /// false case
   template<class>
   static std::false_type __test(...);
public:
   static constexpr const bool value = std::is_same<std::true_type, decltype(__test<T>(0))>::value;
};

/// test _Range conforms the TWG.Range concept
/// check only index_type, ordinal_type, and rank() member
template<class _Range>
class is_range {
public:
   static constexpr const bool
   value = has_index_type<_Range>::value & has_ordinal_type<_Range>::value &
           has_rank<_Range>::value & has_begin<_Range>::value & has_end<_Range>::value;
};

/// test T has extents() member
template<class T>
class has_extent {
   /// true case
   template<class U>
   static auto __test(U* p) -> decltype(p->extent(), std::true_type());
   /// false case
   template<class>
   static std::false_type __test(...);
public:
   static constexpr const bool value = std::is_same<std::true_type, decltype(__test<T>(0))>::value;
};

/// test T has range_size_type
template<class T>
class has_extent_type {
   /// true case
   template<class U>
   static std::true_type __test(typename U::extent_type*);
   /// false case
   template<class>
   static std::false_type __test(...);
public:
   static constexpr const bool value = std::is_same<std::true_type, decltype(__test<T>(0))>::value;
};

/// test _Range conforms the TWG.BoxRange concept
/// in addition to Range, check extent() member and extent_type
template<class _Range>
class is_boxrange {
public:
   static constexpr const bool
   value = is_range<_Range>::value &
           has_extent<_Range>::value &
           has_extent_type<_Range>::value;
};

template <class _Range>
class boxrange_iteration_order {
public:
  enum {row_major = -1, other = 0, column_major = 1};
  // must specialize this trait
  static constexpr const int
  value = other;
};

/// Range Traits
template <typename Range>
struct range_traits;

} // namespace btas

#endif // __BTAS_RANGE_TRAITS_H
