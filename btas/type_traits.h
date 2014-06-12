#ifndef __BTAS_TYPE_TRAITS_H
#define __BTAS_TYPE_TRAITS_H 1

#include <type_traits>

namespace btas {

  /// extends std::common_type to yield a signed integer type if one of the arguments is a signed type
  template <typename I0, typename I1>
  struct common_signed_type {
      typedef typename std::common_type<I0,I1>::type common_type;
      typedef typename std::conditional<
          std::is_signed<I0>::value || std::is_signed<I1>::value,
          typename std::make_signed<common_type>::type,
          common_type
        >::type type;
  }; // common_signed_type

  /// test T has begin() member
  template<class T>
  class has_begin {
      /// true case
      template<class U>
      static auto __test(U* p) -> decltype(p->begin(), std::true_type());
      /// false case
      template<class >
      static std::false_type __test(...);
    public:
      static constexpr const bool value = std::is_same<std::true_type,
          decltype(__test<T>(0))>::value;
  };

  /// test T has end() member
  template<class T>
  class has_end {
      /// true case
      template<class U>
      static auto __test(U* p) -> decltype(p->end(), std::true_type());
      /// false case
      template<class >
      static std::false_type __test(...);
    public:
      static constexpr const bool value = std::is_same<std::true_type,
          decltype(__test<T>(0))>::value;
  };

  /// test T has value_type
  template<class T>
  class has_value_type {
      /// true case
      template<class U>
      static std::true_type __test(typename U::value_type*);
      /// false case
      template<class >
      static std::false_type __test(...);
    public:
      static constexpr const bool value = std::is_same<std::true_type,
          decltype(__test<T>(0))>::value;
  };

  /// test _C conforms to the standard Container concept; basic tests only
  template<class _C>
  class is_container {
    public:
      static constexpr const bool value = has_value_type<_C>::value
          & has_begin<_C>::value & has_end<_C>::value;
  };

  /// test T has operator[] member
  template<class T>
  class has_squarebraket {
      /// true case
      template<class U>
      static auto __test(
          U* p, std::size_t i) -> decltype(p->operator[](i), std::true_type());
      /// false case
      template<class >
      static std::false_type __test(...);
    public:
      static constexpr const bool value = std::is_same<std::true_type,
          decltype(__test<T>(0,std::size_t(0)))>::value;
  };

} // namespace btas

#endif // __BTAS_TYPE_TRAITS_H
