//
// BAGEL - Parallel electron correlation program.
// Filename: archive.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
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
