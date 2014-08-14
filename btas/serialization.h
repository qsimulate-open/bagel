#ifndef __BTAS_SERIALIZATION_H
#define __BTAS_SERIALIZATION_H 1

#include <complex>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/array.hpp>

namespace btas {

  template <typename T>
  auto make_array(T* data, const size_t n) -> decltype(boost::serialization::make_array(data, n)) {
    return boost::serialization::make_array(data, n);
  }
  // specialization for complex
  template <typename T>
  auto make_array(std::complex<T>* data, const size_t n) -> decltype(boost::serialization::make_array(reinterpret_cast<T*>(data), 2*n)) {
    return boost::serialization::make_array(reinterpret_cast<T*>(data), 2*n);
  }
  template <typename T>
  auto make_array(const std::complex<T>* data, const size_t n) -> decltype(boost::serialization::make_array(reinterpret_cast<const T*>(data), 2*n)) {
    return boost::serialization::make_array(reinterpret_cast<const T*>(data), 2*n);
  }

}

#endif
