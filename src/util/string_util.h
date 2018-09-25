//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: string_util.h
// Copyright (C) 2013 Toru Shiozaki
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

// Degraded version of lexical_cast.

#ifndef __SRC_UTIL_STRING_UTIL_H
#define __SRC_UTIL_STRING_UTIL_H

#include <sstream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

namespace bagel {
namespace {

template<typename T, typename U> T lexical_cast(U in) { return boost::lexical_cast<T>(in); }

std::string to_lower(const std::string& in) { std::string tmp(in); boost::algorithm::to_lower(tmp); return tmp; }
std::string to_upper(const std::string& in) { std::string tmp(in); boost::algorithm::to_upper(tmp); return tmp; }

}}

#endif

