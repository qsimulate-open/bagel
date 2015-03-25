#ifndef __BTAS_SERIALIZATION_H
#define __BTAS_SERIALIZATION_H 1

#include <array>
#include <boost/serialization/is_bitwise_serializable.hpp>
#include <boost/serialization/array.hpp>

namespace boost { namespace serialization {
  // this is needed to serialize  efficiently corner cases, like std::vector<std::array<std::complex<T>>>.
  // since bitwise serialization is not portable anyway, this is OK in the context of btas
  template <typename T, size_t N>
  struct is_bitwise_serializable<std::array<T,N> > : is_bitwise_serializable<T> { };
}}

#endif
