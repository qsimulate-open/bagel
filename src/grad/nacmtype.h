//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: nacmtype.h
// Copyright (C) 2017 Toru Shiozaki
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


#ifndef __SRC_GRAD_NACMTYPE_H
#define __SRC_GRAD_NACMTYPE_H

namespace bagel {

enum Nacms { full, interstate, etf, noweight };

class NacmType {
  protected:
    Nacms type_;

  public:
    NacmType() : type_(Nacms::full) { }
    NacmType(std::string input) {
      if (input == "full") {
        type_ = Nacms::full;
      } else if (input == "interstate") {
        type_ = Nacms::interstate;
      } else if (input == "etf") {
        type_ = Nacms::etf;
      } else if (input == "noweight") {
        type_ = Nacms::noweight;
      } else {
        throw std::logic_error ("Available nacmtype: \"full\", \"interstate\", \"etf\", or \"noweight\".");
      }
    }

    bool is_full() const { return type_ == Nacms::full; }
    bool is_interstate() const { return type_ == Nacms::interstate; }
    bool is_etf() const { return type_ == Nacms::etf; }
    bool is_noweight() const { return type_ == Nacms::noweight; }
};

}

#endif
