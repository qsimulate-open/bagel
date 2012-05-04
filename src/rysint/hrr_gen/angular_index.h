//
// Newint - Parallel electron correlation program.
// Filename: angular_index.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include <boost/tuple/tuple.hpp>
#include <string>

class Angular_Index {
  protected:
    boost::tuple<int, int, int> index_;
    int rank_;

  public:
    Angular_Index() {}; 
    Angular_Index(int i, int j, int k) : index_(boost::make_tuple(i, j, k)), rank_(i + j + k) {};
    ~Angular_Index() {};

    const std::string show() const;

    const bool operator==(const Angular_Index& o) const { return x() == o.x() && y() == o.y() && z() == o.z(); };

    const int x() const { return boost::get<0>(index_); };
    const int y() const { return boost::get<1>(index_); };
    const int z() const { return boost::get<2>(index_); };

};


class Angular_Pair {
  protected:
    std::pair<Angular_Index, Angular_Index> indices_;

  public:
    Angular_Pair() {};
    Angular_Pair(std::pair<Angular_Index, Angular_Index> a) : indices_(a) {};

    ~Angular_Pair() {};

    const bool operator==(const Angular_Pair& o) const 
    { return indices().first == o.indices().first && indices().second == o.indices().second; }; 

    boost::tuple<Angular_Pair, Angular_Pair, int> hrr_formula() const; 

    const std::pair<Angular_Index, Angular_Index> indices() const { return indices_; };

    const std::string show() const;
};


