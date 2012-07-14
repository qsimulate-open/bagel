//
// Newint - Parallel electron correlation program.
// Filename: atommap.cc
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


#include <src/scf/atommap.h>
#include <stdexcept>

using namespace std;

AtomMap::AtomMap () {
  atommap.insert(make_pair("h" ,  1));
  atommap.insert(make_pair("he",  2));
  atommap.insert(make_pair("li",  3));
  atommap.insert(make_pair("be",  4));
  atommap.insert(make_pair("b" ,  5));
  atommap.insert(make_pair("c" ,  6));
  atommap.insert(make_pair("n" ,  7));
  atommap.insert(make_pair("o" ,  8));
  atommap.insert(make_pair("f" ,  9));
  atommap.insert(make_pair("ne", 10));

  angmap.insert(make_pair("s", 0));
  angmap.insert(make_pair("p", 1));
  angmap.insert(make_pair("d", 2));
  angmap.insert(make_pair("f", 3));
  angmap.insert(make_pair("g", 4));
  angmap.insert(make_pair("h", 5));
  angmap.insert(make_pair("i", 6));
// Since they are not implemented yet
//angmap.insert(make_pair("j", 7));
//angmap.insert(make_pair("k", 8));
//angmap.insert(make_pair("l", 9));

}


AtomMap::~AtomMap () {

}


const int AtomMap::angular_number(const string input) const {
  auto miter = angmap.find(input);  
  if (miter == angmap.end()) throw runtime_error("Unknown angular number in a basis set file.");
  return miter->second;
}


const int AtomMap::atom_number(const string input) const {
  auto miter = atommap.find(input);  
  if (miter == atommap.end()) throw runtime_error("Unknown Atom number in a basis set file.");
  return miter->second;
}
