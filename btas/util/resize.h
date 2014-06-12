#ifndef __BTAS_RESIZE_H
#define __BTAS_RESIZE_H 1

#include <type_traits>

//
// This provides generic wrapper to treat variable- and fixed-size vectors within the same interface
// e.g. std::vector vs. std::array
//

namespace btas {

/// test T has resize(size_type) member
/// T is assumed to be a vector
template<class T>
class is_resizable {
   /// true case
   template<class U>
   static auto __test(U* p) -> decltype(p->resize(0), std::true_type());
   /// false case
   template<class>
   static std::false_type __test(...);
public:
   static constexpr const bool value = std::is_same<std::true_type, decltype(__test<T>(0))>::value;
};

/// decl.
template<bool _IsResizable> struct __resize_wrapper { };

/// wrapper for variable size vector
template<>
struct __resize_wrapper<true>
{
   /// resize x by n
   template<class _Vector>
   static void resize (_Vector& x, const typename _Vector::size_type& n) { x.resize(n); }
};

/// wrapper for fixed-size vector
template<>
struct __resize_wrapper<false>
{
   /// nothing to do
   template<class _Vector>
   static void resize (_Vector& x, const typename _Vector::size_type& n) { }
};

/// resize vector x
/// if x is resizable, resize x by n, otherwise, do nothing
template<class _Vector>
void resize (_Vector& x, const typename _Vector::size_type& n)
{
   __resize_wrapper<is_resizable<_Vector>::value>::resize(x, n);
}

} // namespace btas

#endif // __BTAS_RESIZE_H
