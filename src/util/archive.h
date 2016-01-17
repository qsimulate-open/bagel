//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: archive.h
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

#ifndef __SRC_UTIL_ARCHIVE_H
#define __SRC_UTIL_ARCHIVE_H

#include <string>
#include <fstream>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/export.hpp>

namespace bagel {

class OArchive {
  protected:
    std::string filename_;
    std::ofstream os_;

    using Ostream = boost::archive::binary_oarchive;
    Ostream archive_;

  public:
    OArchive(std::string name) : filename_(name+".archive"), os_(filename_), archive_(os_) {
    }

    template<typename T>
    OArchive& operator<<(const T& val) {
      archive_ << val;
      return *this;
    }

    template<typename T>
    OArchive& operator>>(T& val) = delete;
};


class IArchive {
  protected:
    std::string filename_;
    std::ifstream is_;

    using Istream = boost::archive::binary_iarchive;
    Istream archive_;

  public:
    IArchive(std::string name) : filename_(name+".archive"), is_(filename_), archive_(is_) {
      if (!is_.is_open())
        throw std::runtime_error(name+".archive not found");
    }

    template<typename T>
    IArchive& operator>>(T& val) {
      archive_ >> val;
      return *this;
    }

    template<typename T>
    IArchive& operator<<(T& val) = delete;
};

}

#endif
