//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: serialization.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

#ifndef __SRC_UTIL_SERIALIZATION_H
#define __SRC_UTIL_SERIALIZATION_H

#include <array>
#include <memory>
#include <tuple>
#include <unordered_map>
#include <type_traits>
#include <boost/version.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/bitset.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/array.hpp>
#include <boost/property_tree/ptree_serialization.hpp>

#if BOOST_VERSION <= 105600
#include <boost/archive/shared_ptr_helper.hpp>
#endif

#if BOOST_VERSION >= 105600
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/weak_ptr.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/unordered_set.hpp>
#endif

// to avoid duplicate serialize functions, we need to include this here
#include <src/util/math/btas_interface.h>

namespace bagel {
  // default implementation of bagel::base_of
  template <class T, typename = void>
  struct base_of {
    // always return non-const type
    typedef typename std::remove_cv<T>::type type;
  };

  // make_array
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

namespace boost {
  namespace serialization {

// as of 1.57 there is a problem with serialization of std::shared_ptr<const T>
#if BOOST_VERSION >= 105600
    template<class Archive, class T>
    inline void serialize(Archive& ar, std::shared_ptr<const T>& t, const unsigned int version) {
      BOOST_STATIC_ASSERT(tracking_level<T>::value != track_never);
      split_free(ar, t, version);
    }

    template<class Archive, class T>
    inline void save(Archive& ar, const std::shared_ptr<const T>& t, const unsigned int version) {
      auto tt =  std::const_pointer_cast<T>(t);
      save(ar, tt, version);
    }

    template<class Archive, class T>
    inline void load(Archive& ar, std::shared_ptr<const T>& t, const unsigned int version) {
      std::shared_ptr<T> tt;
      load(ar, tt, version);
      t = tt;
    }
#else

    template<class Archive, class T>
    inline void serialize(Archive& ar, std::shared_ptr<T>& t, const unsigned int version) {
      BOOST_STATIC_ASSERT(tracking_level<T>::value != track_never);
      split_free(ar, t, version);
    }

    // base, non-const classes
    template<class Archive, class T, typename = void>
    struct save_impl {
      static void save(Archive& ar, T*& t) {
        ar << t;
      }
    };

    // non-base non-const classes
    template<class Archive, class T>
    struct save_impl<Archive, T,
                     typename std::enable_if<not std::is_same<T, typename bagel::base_of<T>::type>::value
                                         and not std::is_const<T>::value
                                            >::type
                    > {
      static void save(Archive& ar, T*& t) {
        typedef typename bagel::base_of<T>::type base;
        base* u = static_cast<base*>(t);
        save_impl<Archive, base>::save(ar, u);
      }
    };

    // const classes
    template<class Archive, class T>
    struct save_impl<Archive, T,
                     typename std::enable_if<std::is_const<T>::value>::type
                    > {
      static void save(Archive& ar, T*& t) {
        typedef typename std::remove_cv<T>::type non_const;
        non_const* u = const_cast<non_const*>(t);
        save_impl<Archive, non_const>::save(ar, u);
      }
    };

    template<class Archive, class T>
    inline void save(Archive& ar, const std::shared_ptr<T>& t, const unsigned int) {
      T* u = t.get();
      save_impl<Archive, T>::save(ar, u);
    }

    // base, non-const classes
    template<class Archive, class T, typename = void>
    struct load_impl {
      static void load(Archive& ar, std::shared_ptr<T>& t) {
        BOOST_STATIC_ASSERT((tracking_level<T>::value != track_never));
        T* ptr;
        ar >> ptr;

        static std::unordered_map<T*, std::weak_ptr<T>> hash;
        if (hash[ptr].expired()) {
          t = std::shared_ptr<T>(ptr);
          hash[ptr] = t;
        } else {
          t = hash[ptr].lock();
        }
      }
    };

    // non-base non-const classes
    template<class Archive, class T>
    struct load_impl<Archive, T,
                     typename std::enable_if<not std::is_same<T, typename bagel::base_of<T>::type>::value
                                         and not std::is_const<T>::value
                                            >::type
                    > {
      static void load(Archive& ar, std::shared_ptr<T>& t) {
        typedef typename bagel::base_of<T>::type base;
        std::shared_ptr<base> u;
        load_impl<Archive, base>::load(ar, u);
        t = dynamic_pointer_cast<T>(u);
      }
    };

    // const classes
    template<class Archive, class T>
    struct load_impl<Archive, T,
                     typename std::enable_if<std::is_const<T>::value>::type
                    > {
      static void load(Archive& ar, std::shared_ptr<T>& t) {
        typedef typename std::remove_cv<T>::type non_const;
        std::shared_ptr<non_const> u;
        load_impl<Archive, non_const>::load(ar, u);
        t = u;
      }
    };

    template<class Archive, class T>
    inline void load(Archive& ar, std::shared_ptr<T>& t, const unsigned int) {
      load_impl<Archive, T>::load(ar, t);
    }

    // serialization of weak_ptr
    template<class Archive, class T>
    void serialize(Archive& ar, std::weak_ptr<T>& t, const unsigned int version) {
      BOOST_STATIC_ASSERT(tracking_level<T>::value != track_never);
      split_free(ar, t, version);
    }

    template<class Archive, class T>
    void save(Archive& ar, const std::weak_ptr<T>& t, const unsigned int version) {
      std::shared_ptr<T> u = t.lock();
      ar << u;
    }

    template<class Archive, class T>
    void load(Archive& ar, std::weak_ptr<T>& t, const unsigned int version) {
      std::shared_ptr<T> u;
      ar >> u;
      t = u;
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////

    // serialization of unordered_map
    template<class Archive, class T, class U>
    void serialize(Archive& ar, std::unordered_map<T,U>& t, const unsigned int version) {
      split_free(ar, t, version);
    }

    template<class Archive, class T, class U>
    void save(Archive& ar, const std::unordered_map<T,U>& t, const unsigned int version) {
      int size = t.size();
      ar << size;
      for (auto& i : t) ar << i.first << i.second;
    }

    template<class Archive, class T, class U>
    void load(Archive& ar, std::unordered_map<T,U>& t, const unsigned int version) {
      int size;
      ar >> size;
      for (int i = 0; i != size; ++i) {
        typename std::remove_cv<T>::type a;
        typename std::remove_cv<U>::type b;
        ar >> a >> b;
        t.emplace(a, b);
      }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    // serialization of std::array
#ifndef BOOST_SERIALIZATION_STD_ARRAY
    template<class Archive, typename T, size_t N>
    void serialize(Archive& ar, std::array<T,N>& t, const unsigned int) {
      ar & bagel::make_array(t.data(), N);
    }
#endif
#endif

    // serialization of tuple
    template <size_t N>
    struct Serialize {
      template<class Archive, typename... Args>
      static void serialize(Archive& ar, std::tuple<Args...>& t, const unsigned int version) {
        ar & std::get<N-1>(t);
        Serialize<N-1>::serialize(ar, t, version);
      }
    };

    template<>
    struct Serialize<0LLU> {
      template<class Archive, typename... Args>
      static void serialize(Archive&, std::tuple<Args...>&, const unsigned int) { }
    };

    template<class Archive, typename... Args>
    void serialize(Archive& ar, std::tuple<Args...>& t, const unsigned int version) {
      Serialize<sizeof...(Args)>::serialize(ar, t, version);
    }

  }
}

// If serialization is diabled, we reset the macros
#ifdef DISABLE_SERIALIZATION
  #undef BOOST_CLASS_EXPORT_KEY
  #define BOOST_CLASS_EXPORT_KEY(x)
  #undef BOOST_CLASS_EXPORT_IMPLEMENT
  #define BOOST_CLASS_EXPORT_IMPLEMENT(x)
#endif

#endif
