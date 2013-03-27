//
// BAGEL - Parallel electron correlation program.
// Filename: atommap.cc
// Copyright (C) 2009 Toru Shiozaki
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


#include <src/util/atommap.h>
#include <src/util/constants.h>
#include <sstream>
#include <stdexcept>

using namespace std;
using namespace bagel;

AtomMap::AtomMap () {
  atommap.insert(make_pair("q", 0));
  atommap.insert(make_pair("h", 1));
  atommap.insert(make_pair("he", 2));
  atommap.insert(make_pair("li", 3));
  atommap.insert(make_pair("be", 4));
  atommap.insert(make_pair("b", 5));
  atommap.insert(make_pair("c", 6));
  atommap.insert(make_pair("n", 7));
  atommap.insert(make_pair("o", 8));
  atommap.insert(make_pair("f", 9));
  atommap.insert(make_pair("ne", 10));
  atommap.insert(make_pair("na", 11));
  atommap.insert(make_pair("mg", 12));
  atommap.insert(make_pair("al", 13));
  atommap.insert(make_pair("si", 14));
  atommap.insert(make_pair("p", 15));
  atommap.insert(make_pair("s", 16));
  atommap.insert(make_pair("cl", 17));
  atommap.insert(make_pair("ar", 18));
  atommap.insert(make_pair("k", 19));
  atommap.insert(make_pair("ca", 20));
  atommap.insert(make_pair("sc", 21));
  atommap.insert(make_pair("ti", 22));
  atommap.insert(make_pair("v", 23));
  atommap.insert(make_pair("cr", 24));
  atommap.insert(make_pair("mn", 25));
  atommap.insert(make_pair("fe", 26));
  atommap.insert(make_pair("co", 27));
  atommap.insert(make_pair("ni", 28));
  atommap.insert(make_pair("cu", 29));
  atommap.insert(make_pair("zn", 30));
  atommap.insert(make_pair("ga", 31));
  atommap.insert(make_pair("ge", 32));
  atommap.insert(make_pair("as", 33));
  atommap.insert(make_pair("se", 34));
  atommap.insert(make_pair("br", 35));
  atommap.insert(make_pair("kr", 36));
  atommap.insert(make_pair("rb", 37));
  atommap.insert(make_pair("sr", 38));
  atommap.insert(make_pair("y", 39));
  atommap.insert(make_pair("zr", 40));
  atommap.insert(make_pair("nb", 41));
  atommap.insert(make_pair("mo", 42));
  atommap.insert(make_pair("tc", 43));
  atommap.insert(make_pair("ru", 44));
  atommap.insert(make_pair("rh", 45));
  atommap.insert(make_pair("pd", 46));
  atommap.insert(make_pair("ag", 47));
  atommap.insert(make_pair("cd", 48));
  atommap.insert(make_pair("in", 49));
  atommap.insert(make_pair("sn", 50));
  atommap.insert(make_pair("sb", 51));
  atommap.insert(make_pair("te", 52));
  atommap.insert(make_pair("i", 53));
  atommap.insert(make_pair("xe", 54));
  atommap.insert(make_pair("cs", 55));
  atommap.insert(make_pair("ba", 56));
  atommap.insert(make_pair("la", 57));
  atommap.insert(make_pair("ce", 58));
  atommap.insert(make_pair("pr", 59));
  atommap.insert(make_pair("nd", 60));
  atommap.insert(make_pair("pm", 61));
  atommap.insert(make_pair("sm", 62));
  atommap.insert(make_pair("eu", 63));
  atommap.insert(make_pair("gd", 64));
  atommap.insert(make_pair("tb", 65));
  atommap.insert(make_pair("dy", 66));
  atommap.insert(make_pair("ho", 67));
  atommap.insert(make_pair("er", 68));
  atommap.insert(make_pair("tm", 69));
  atommap.insert(make_pair("yb", 70));
  atommap.insert(make_pair("lu", 71));
  atommap.insert(make_pair("hf", 72));
  atommap.insert(make_pair("ta", 73));
  atommap.insert(make_pair("w", 74));
  atommap.insert(make_pair("re", 75));
  atommap.insert(make_pair("os", 76));
  atommap.insert(make_pair("ir", 77));
  atommap.insert(make_pair("pt", 78));
  atommap.insert(make_pair("au", 79));
  atommap.insert(make_pair("hg", 80));
  atommap.insert(make_pair("tl", 81));
  atommap.insert(make_pair("pb", 82));
  atommap.insert(make_pair("bi", 83));
  atommap.insert(make_pair("po", 84));
  atommap.insert(make_pair("at", 85));
  atommap.insert(make_pair("rn", 86));
  atommap.insert(make_pair("fr", 87));
  atommap.insert(make_pair("ra", 88));
  atommap.insert(make_pair("ac", 89));
  atommap.insert(make_pair("th", 90));
  atommap.insert(make_pair("pa", 91));
  atommap.insert(make_pair("u", 92));
  atommap.insert(make_pair("np", 93));
  atommap.insert(make_pair("pu", 94));
  atommap.insert(make_pair("am", 95));
  atommap.insert(make_pair("cm", 96));
  atommap.insert(make_pair("bk", 97));
  atommap.insert(make_pair("cf", 98));
  atommap.insert(make_pair("es", 99));
  atommap.insert(make_pair("fm", 100));
  atommap.insert(make_pair("md", 101));
  atommap.insert(make_pair("no", 102));
  atommap.insert(make_pair("lr", 103));
  atommap.insert(make_pair("rf", 104));
  atommap.insert(make_pair("db", 105));
  atommap.insert(make_pair("sg", 106));
  atommap.insert(make_pair("bh", 107));
  atommap.insert(make_pair("hs", 108));
  atommap.insert(make_pair("mt", 109));
  atommap.insert(make_pair("ds", 110));
  atommap.insert(make_pair("rg", 111));
  atommap.insert(make_pair("cn", 112));
  atommap.insert(make_pair("uut", 113));
  atommap.insert(make_pair("fl", 114));
  atommap.insert(make_pair("uup", 115));
  atommap.insert(make_pair("lv", 116));
  atommap.insert(make_pair("uus", 117));
  atommap.insert(make_pair("uuo", 118));

  // atom sizes (Bragg-Slater radii)
  bsradii.insert(make_pair("h", 0.25));
  bsradii.insert(make_pair("he", 0.25));
  bsradii.insert(make_pair("li", 1.45));
  bsradii.insert(make_pair("be", 1.05));
  bsradii.insert(make_pair("b", 0.85));
  bsradii.insert(make_pair("c", 0.7));
  bsradii.insert(make_pair("n", 0.65));
  bsradii.insert(make_pair("o", 0.6));
  bsradii.insert(make_pair("f", 0.5));
  bsradii.insert(make_pair("ne", 0.45));
  bsradii.insert(make_pair("na", 1.8));
  bsradii.insert(make_pair("mg", 1.5));
  bsradii.insert(make_pair("al", 1.25));
  bsradii.insert(make_pair("si", 1.1));
  bsradii.insert(make_pair("p", 1.0));
  bsradii.insert(make_pair("s", 1.0));
  bsradii.insert(make_pair("cl", 1.0));
  bsradii.insert(make_pair("ar", 1.0));
  bsradii.insert(make_pair("k", 2.2));
  bsradii.insert(make_pair("ca", 1.8));
  bsradii.insert(make_pair("sc", 1.6));
  bsradii.insert(make_pair("ti", 1.4));
  bsradii.insert(make_pair("v", 1.35));
  bsradii.insert(make_pair("cr", 1.4));
  bsradii.insert(make_pair("mn", 1.4));
  bsradii.insert(make_pair("fe", 1.4));
  bsradii.insert(make_pair("co", 1.35));
  bsradii.insert(make_pair("ni", 1.35));
  bsradii.insert(make_pair("cu", 1.35));
  bsradii.insert(make_pair("zn", 1.35));
  bsradii.insert(make_pair("ga", 1.3));
  bsradii.insert(make_pair("ge", 1.25));
  bsradii.insert(make_pair("as", 1.15));
  bsradii.insert(make_pair("se", 1.15));
  bsradii.insert(make_pair("br", 1.15));
  bsradii.insert(make_pair("kr", 1.15));
  bsradii.insert(make_pair("rb", 2.35));
  bsradii.insert(make_pair("sr", 2.0));
  bsradii.insert(make_pair("y", 1.8));
  bsradii.insert(make_pair("zr", 1.55));
  bsradii.insert(make_pair("nb", 1.45));
  bsradii.insert(make_pair("mo", 1.45));
  bsradii.insert(make_pair("tc", 1.35));
  bsradii.insert(make_pair("ru", 1.3));
  bsradii.insert(make_pair("rh", 1.35));
  bsradii.insert(make_pair("pd", 1.4));
  bsradii.insert(make_pair("ag", 1.6));
  bsradii.insert(make_pair("cd", 1.55));
  bsradii.insert(make_pair("in", 1.55));
  bsradii.insert(make_pair("sn", 1.45));
  bsradii.insert(make_pair("sb", 1.45));
  bsradii.insert(make_pair("te", 1.4));
  bsradii.insert(make_pair("i", 1.4));
  bsradii.insert(make_pair("xe", 1.4));
  bsradii.insert(make_pair("cs", 2.6));
  bsradii.insert(make_pair("ba", 2.15));
  bsradii.insert(make_pair("la", 1.95));
  bsradii.insert(make_pair("ce", 1.85));
  bsradii.insert(make_pair("pr", 1.85));
  bsradii.insert(make_pair("nd", 1.85));
  bsradii.insert(make_pair("pm", 1.85));
  bsradii.insert(make_pair("sm", 1.85));
  bsradii.insert(make_pair("eu", 1.85));
  bsradii.insert(make_pair("gd", 1.8));
  bsradii.insert(make_pair("tb", 1.75));
  bsradii.insert(make_pair("dy", 1.75));
  bsradii.insert(make_pair("ho", 1.75));
  bsradii.insert(make_pair("er", 1.75));
  bsradii.insert(make_pair("tm", 1.75));
  bsradii.insert(make_pair("yb", 1.75));
  bsradii.insert(make_pair("lu", 1.75));
  bsradii.insert(make_pair("hf", 1.55));
  bsradii.insert(make_pair("ta", 1.45));
  bsradii.insert(make_pair("w", 1.35));
  bsradii.insert(make_pair("re", 1.35));
  bsradii.insert(make_pair("os", 1.3));
  bsradii.insert(make_pair("ir", 1.35));
  bsradii.insert(make_pair("pt", 1.35));
  bsradii.insert(make_pair("au", 1.35));
  bsradii.insert(make_pair("hg", 1.5));
  bsradii.insert(make_pair("tl", 1.9));
  bsradii.insert(make_pair("pb", 1.8));
  bsradii.insert(make_pair("bi", 1.6));
  bsradii.insert(make_pair("po", 1.9));
  bsradii.insert(make_pair("at", 1.9));
  bsradii.insert(make_pair("rn", 1.9));
  bsradii.insert(make_pair("fr", 2.85));
  bsradii.insert(make_pair("ra", 2.15));
  bsradii.insert(make_pair("ac", 1.95));
  bsradii.insert(make_pair("th", 1.8));
  bsradii.insert(make_pair("pa", 1.8));
  bsradii.insert(make_pair("u", 1.75));
  bsradii.insert(make_pair("np", 1.75));
  bsradii.insert(make_pair("pu", 1.75));
  bsradii.insert(make_pair("am", 1.75));
  bsradii.insert(make_pair("cm", 1.75));
  bsradii.insert(make_pair("bk", 1.75));
  bsradii.insert(make_pair("cf", 1.75));
  bsradii.insert(make_pair("es", 1.75));
  bsradii.insert(make_pair("fm", 1.75));
  bsradii.insert(make_pair("md", 1.75));
  bsradii.insert(make_pair("no", 1.75));
  bsradii.insert(make_pair("lr", 1.75));
#if 0
  bsradii.insert(make_pair("rf", 104));
  bsradii.insert(make_pair("db", 105));
  bsradii.insert(make_pair("sg", 106));
  bsradii.insert(make_pair("bh", 107));
  bsradii.insert(make_pair("hs", 108));
  bsradii.insert(make_pair("mt", 109));
  bsradii.insert(make_pair("ds", 110));
  bsradii.insert(make_pair("rg", 111));
  bsradii.insert(make_pair("cn", 112));
  bsradii.insert(make_pair("uut", 113));
  bsradii.insert(make_pair("fl", 114));
  bsradii.insert(make_pair("uup", 115));
  bsradii.insert(make_pair("lv", 116));
  bsradii.insert(make_pair("uus", 117));
  bsradii.insert(make_pair("uuo", 118));
#endif

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

int AtomMap::angular_number(const string input) const {
  auto miter = angmap.find(input);
  stringstream ss; ss << "Unknown angular number in a basis set file. Requested: " << input << endl;
  if (miter == angmap.end()) throw runtime_error(ss.str());
  return miter->second;
}


double AtomMap::radius(const string input) const {
  auto miter = bsradii.find(input);
  if (miter == bsradii.end()) throw runtime_error("Unknown atom (Bragg-Slater radii).");
  return miter->second*ang2bohr__;
}


int AtomMap::atom_number(const string input) const {
  auto miter = atommap.find(input);
  stringstream ss; ss << "Unknown Atom number in a basis set file. Requested: " << input << endl;
  if (miter == atommap.end()) throw runtime_error(ss.str());
  return miter->second;
}

const string AtomMap::angular_string(const int input) {
  for(auto& m : angmap) {
    if(m.second == input) { return m.first; }
  }
  return "X";
}
