#ifndef __BTAS_ARRAYADAPTOR_H_
#define __BTAS_ARRAYADAPTOR_H_

// adaptors for std "array" containers

#include <vector>
#include <array>
#include <cassert>
#include <btas/tensor_traits.h>
#include <btas/generic/numeric_type.h>
#include <btas/varray/varray.h>

namespace btas {

  template <typename Array> struct array_adaptor;

  /// Adaptor from std::array
  template <typename T, size_t N>
  struct array_adaptor< std::array<T, N> > {
      typedef std::array<T, N> array;
      typedef typename array::value_type value_type;

      static array construct(std::size_t n) {
        assert(n <= N);
        array result;
        return result;
      }
      static array construct(std::size_t n, T value) {
        assert(n <= N);
        array result;
        std::fill_n(result.begin(), n, value);
        return result;
      }

      static void resize(array& x, std::size_t n) {
        assert(x.size() == N);
        assert(x.size() >= n);
      }

      static void print(const array& a, std::ostream& os) {
        os << "{";
        for(std::size_t i = 0; i != N; ++i) {
          os << a[i];
          if (i != (N - 1))
            os << ",";
        }
        os << "}";
      }
  };

  template <typename T, size_t N>
  std::size_t rank(const std::array<T, N>& x) {
    return N;
  }

  /// Adaptor from const-size array
  template <typename T, size_t N>
  struct array_adaptor< T[N] > {
      typedef T (array)[N];
      typedef T value_type;

      static void print(const array& a, std::ostream& os) {
        os << "{";
        for(std::size_t i = 0; i != N; ++i) {
          os << a[i];
          if (i != (N - 1))
            os << ",";
        }
        os << "}";
      }
  };

  template <typename T, size_t N>
  std::size_t rank(const T (&x)[N]) {
    return N;
  }

  template <typename T, size_t N>
  std::ostream& operator<<(std::ostream& os, const std::array<T, N>& x) {
    array_adaptor<std::array<T, N> >::print(x,os);
    return os;
  }

  /// Adaptors for sequence container, e.g. std::vector, btas::varray, and std::initializer_list

  template <typename Array, class = typename std::enable_if<not btas::is_tensor<Array>::value>::type>
  std::size_t rank(const Array& x) {
    return x.size();
  }

  template <typename Array>
  struct array_adaptor {
      typedef Array array;
      typedef typename Array::value_type value_type;

      static array construct(std::size_t N) {
        return array(N);
      }
      static array construct(std::size_t N,
                             value_type value) {
        return array(N, value);
      }
      static void resize(array& x, std::size_t N) {
        x.resize(N);
      }
      static void print(const array& a, std::ostream& os) {
        std::size_t n = rank(a);
        os << "{";
        for(std::size_t i = 0; i != n; ++i) {
          os << a[i];
          if (i != (n - 1))
            os << ",";
        }
        os << "}";
      }
  };

  template <typename T>
  std::ostream& operator<<(std::ostream& os, const btas::varray<T>& x) {
    array_adaptor<btas::varray<T> >::print(x,os);
    return os;
  }

  template <typename T>
  std::ostream& operator<<(std::ostream& os, const std::vector<T>& x) {
    array_adaptor<std::vector<T> >::print(x,os);
    return os;
  }

  template <typename T>
  std::ostream& operator<<(std::ostream& os, const std::initializer_list<T>& x) {
    array_adaptor<std::vector<T> >::print(x,os);
    return os;
  }
}

namespace std {

  template <typename T, size_t N>
  auto cbegin(const std::array<T,N>& x) -> decltype(x.cbegin()) {
    return x.cbegin();
  }
  template <typename T, size_t N>
  auto cend(const std::array<T,N>& x) -> decltype(x.cend()) {
    return x.cend();
  }

  template <typename T>
  auto cbegin(const btas::varray<T>& x) -> decltype(x.cbegin()) {
    return x.cbegin();
  }
  template <typename T>
  auto cend(const btas::varray<T>& x) -> decltype(x.cend()) {
    return x.cend();
  }

  template <typename T>
  auto cbegin(const std::vector<T>& x) -> decltype(x.cbegin()) {
    return x.cbegin();
  }
  template <typename T>
  auto cend(const std::vector<T>& x) -> decltype(x.cend()) {
    return x.cend();
  }

  template <typename T>
  struct make_unsigned<std::vector<T> > {
      typedef std::vector<typename make_unsigned<T>::type > type;
  };
  template <typename T>
  struct make_unsigned<std::initializer_list<T> > {
      typedef std::initializer_list<typename make_unsigned<T>::type > type;
  };
  template <typename T, size_t N>
  struct make_unsigned<std::array<T, N> > {
      typedef std::array<typename make_unsigned<T>::type, N> type;
  };
  template <typename T>
  struct make_unsigned<btas::varray<T> > {
      typedef btas::varray<typename make_unsigned<T>::type > type;
  };
  template <typename T, size_t N>
  struct make_unsigned<T[N]> {
      typedef typename make_unsigned<T>::type uT;
      typedef uT (type)[N];
  };

}

namespace btas {
  template <typename Array, typename T>
  struct replace_value_type;

  template <typename T, typename U>
  struct replace_value_type<std::vector<T>, U> {
      typedef std::vector<U> type;
  };
  template <typename T, typename U>
  struct replace_value_type<std::initializer_list<T>,U> {
      typedef std::initializer_list<U> type;
  };
  template <typename T, size_t N, typename U>
  struct replace_value_type<std::array<T, N>,U> {
      typedef std::array<U, N> type;
  };
  template <typename T, typename U>
  struct replace_value_type<btas::varray<T>,U> {
      typedef btas::varray<U> type;
  };
  template <typename T, size_t N, typename U>
  struct replace_value_type<T[N],U> {
      typedef U (type)[N];
  };
}

#endif /* __BTAS_ARRAYADAPTOR_H_ */
