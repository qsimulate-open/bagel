//
// BAGEL - Parallel electron correlation program.
// Filename: unpack.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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

// Small helper class which is needed to store variadic template arguments as a tuple and reuse later
// Currently only used in <src/integral/smallints1e_london.h>
// See http://stackoverflow.com/questions/7858817/

#ifndef __SRC_UTIL_UNPACK_H
#define __SRC_UTIL_UNPACK_H

namespace {

  template<int ...>
  struct seq { };

  template<int N, int ...S>
  struct gens : gens<N-1, N-1, S...> { };

  template<int ...S>
  struct gens<0, S...> {
    typedef seq<S...> type;
  };

}

#endif
