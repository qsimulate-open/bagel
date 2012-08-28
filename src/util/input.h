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
#include <stdexcept>
#include <boost/lexical_cast.hpp>

namespace bagel {

class InputData {
  protected:
    std::list<std::pair<std::string, std::multimap<std::string, std::string> > >data_;
    const std::string inputfile_;

  public:
    InputData(const std::string filename);
    ~InputData() {};

    std::multimap<std::string, std::string> get_input(const std::string t) const {
      auto iter = data_.begin();
      for (; iter != data_.end(); ++iter) if (iter->first == t) break;
      if (iter == data_.end())
        throw std::runtime_error(t + " does not appear to be present in your input");
      return iter->second;
    };

    bool exist (const std::string t) const {
      auto iter = data_.begin();
      for (; iter != data_.end(); ++iter) if (iter->first == t) break;
      return data_.end() != iter;
    };

    std::list<std::pair<std::string, std::multimap<std::string, std::string> > > data() { return data_; };

};

template <typename T> T read_input(const std::multimap<std::string, std::string> idat, const std::string key, const T defvalue) {
  T out = defvalue;
  auto iter = idat.find(key);
  if (iter != idat.end()) out = boost::lexical_cast<T>(iter->second);
  return out;
};

}

#endif
