//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: citraits.h
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

#ifndef __SRC_CIUTIL_CITRAITS_H
#define __SRC_CIUTIL_CITRAITS_H

#include <type_traits>
#include <src/ci/ciutil/cistring.h>

namespace bagel {
namespace {

// Traits for strings (not really robust)
template <class T>
struct is_cistring {
  private:
    const static std::bitset<nbit__> dum;
    template<class U> static auto __lexical_zero(U* p) -> decltype(p->lexical_zero(dum), std::true_type());
    template<class  > static std::false_type __lexical_zero(...);
    template<class U> static auto __offset(U* p) -> decltype(p->offset(), std::true_type());
    template<class  > static std::false_type __offset(...);
    template<class U> static auto __phi(U* p) -> decltype(p->phi(), std::true_type());
    template<class  > static std::false_type __phi(...);
  public:
    static constexpr const bool value = std::is_base_of<CIString_base, T>::value
                                    and std::is_same<std::true_type, decltype(__offset<T>(0))>::value
                                    and std::is_same<std::true_type, decltype(__lexical_zero<T>(0))>::value
                                    and std::is_same<std::false_type, decltype(__phi<T>(0))>::value; // to rule out StringSet
};


// Traits for stringsets
template <class T>
struct is_cistringset {
  private:
    template<class U> static auto __phi(U* p) -> decltype(p->phi(), std::true_type());
    template<class  > static std::false_type __phi(...);
    template<class U> static auto __find_string(U* p) -> decltype(p->find_string(0,0), std::true_type());
    template<class  > static std::false_type __find_string(...);
  public:
    static constexpr const bool value = std::is_same<std::true_type, decltype(__phi<T>(0))>::value
                                    and std::is_same<std::true_type, decltype(__find_string<T>(0))>::value;
};


// Traits for blocks
template <class T>
struct has_string_type {
  private:
    template<class U> static std::true_type __test(typename U::string_type*);
    template<class  > static std::false_type __test(...);
  public:
    static constexpr const bool value = std::is_same<std::true_type, decltype(__test<T>(0))>::value;
};

template <class T>
struct has_data_type {
  private:
    template<class U> static std::true_type __test(typename U::data_type*);
    template<class  > static std::false_type __test(...);
  public:
    static constexpr const bool value = std::is_same<std::true_type, decltype(__test<T>(0))>::value;
};


template <class T>
struct is_ciblockinfo {
  private:
    const static std::bitset<nbit__> dum;
    template<class U> static auto __test(U* p) -> decltype(p->stringsa(), std::true_type());
    template<class  > static std::false_type __test(...);
    template<class U> static auto __bit(U* p) -> decltype(p->string_bits_a(0), std::true_type());
    template<class  > static std::false_type __bit(...);
  public:
    static constexpr const bool value = std::is_same<std::true_type, decltype(__test<T>(0))>::value
                                    and std::is_same<std::true_type, decltype(__bit<T>(0))>::value
                                    and is_cistring<typename std::conditional<has_string_type<T>::value, typename T::string_type, void>::type>::value;
};

template <class T>
struct is_ciblock {
  public:
    static constexpr const bool value = is_ciblockinfo<T>::value
                                    and has_data_type<T>::value;
};

}
}

#endif
