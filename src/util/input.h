//
// BAGEL - Parallel electron correlation program.
// Filename: input.h
// Copyright (C) 2011 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


// new interface to the input

#ifndef __NEWINT_UTIL_INPUT_H
#define __NEWINT_UTIL_INPUT_H

#include <map>
#include <list>
#include <string>
#include <vector>
#include <stdexcept>
#include <src/util/lexical_cast.h>

namespace bagel {

class InputData {
  protected:
    std::list<std::pair<std::string, std::multimap<std::string, std::string>>> data_;
    const std::string inputfile_;

  public:
    InputData(const std::string filename);

    std::multimap<std::string, std::string> get_input(const std::string t) const {
      auto iter = data_.begin();
      for (; iter != data_.end(); ++iter) if (iter->first == t) break;
      if (iter == data_.end())
        throw std::runtime_error(t + " does not appear to be present in your input");
      return iter->second;
    }

    bool exist (const std::string t) const {
      auto iter = data_.begin();
      for (; iter != data_.end(); ++iter) if (iter->first == t) break;
      return data_.end() != iter;
    }

    std::list<std::pair<std::string, std::multimap<std::string, std::string>>> data() { return data_; }

};

static std::vector<std::string> split(const std::string& s, const std::string& delimiters) {
  std::vector<std::string> out;
  size_t current;
  size_t next = -1;
  do {
    current = next + 1;
    next = s.find_first_of(delimiters, current);
    const std::string now = s.substr(current, next - current);
    out.push_back(now);
  } while (next != std::string::npos);
  return out;
}
static std::string &ltrim(std::string &s) { s.erase(s.begin(), find_if(s.begin(), s.end(), not1(std::ptr_fun<int, int>(isspace)))); return s; }
static std::string &rtrim(std::string &s) { s.erase(find_if(s.rbegin(), s.rend(), not1(std::ptr_fun<int, int>(isspace))).base(), s.end()); return s; }
static std::string &trim(std::string &s) { return ltrim(rtrim(s)); }

template <typename T> T read_input(const std::multimap<std::string, std::string> idat, const std::string key, const T defvalue) {
  T out = defvalue;
  auto iter = idat.find(key);
  if (iter != idat.end()) out = lexical_cast<T>(iter->second);
  return out;
}

template<typename T> std::vector<std::vector<T>> read_input_vector_range(const std::multimap<std::string, std::string> idat, const std::string key, const std::string defvalue) {
  std::vector<std::vector<T>> out;

  auto range = idat.equal_range(key);
  for (auto iter = range.first; iter != range.second; ++iter) {
    std::vector<std::string> vecstring = split(iter->second, ",");
    std::vector<T> tmp;
    for ( auto& i : vecstring ) {
      tmp.push_back(lexical_cast<T>(i));
    }
    out.push_back(tmp);
  }

  return out;
}

}

#endif
