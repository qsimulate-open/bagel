#!/usr/bin/python

l = [ "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Uut", "Fl", "Uup", "Lv", "Uus", "Uuo" ]

ss = "\
//\n\
// BAGEL - Parallel electron correlation program.\n\
// Filename: atommap.cc\n\
// Copyright (C) 2009 Toru Shiozaki\n\
//\n\
// Author: Toru Shiozaki <shiozaki@northwestern.edu>\n\
// Maintainer: Shiozaki group\n\
//\n\
// This file is part of the BAGEL package.\n\
//\n\
// The BAGEL package is free software; you can redistribute it and\/or modify\n\
// it under the terms of the GNU Library General Public License as published by\n\
// the Free Software Foundation; either version 2, or (at your option)\n\
// any later version.\n\
//\n\
// The BAGEL package is distributed in the hope that it will be useful,\n\
// but WITHOUT ANY WARRANTY; without even the implied warranty of\n\
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n\
// GNU Library General Public License for more details.\n\
//\n\
// You should have received a copy of the GNU Library General Public License\n\
// along with the BAGEL package; see COPYING.  If not, write to\n\
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.\n\
//\n\
\n\
\n\
#include <src/util/atommap.h>\n\
#include <stdexcept>\n\
\n\
using namespace std;\n\
using namespace bagel;\n\
\n\
AtomMap::AtomMap () {\n"
ss += "  atommap.insert(make_pair(\"q\", 0));\n"
cnt = 1
for i in l:
  ss += "  atommap.insert(make_pair(\"" + i.lower() + "\", " + str(cnt) + "));\n"
  cnt += 1

ss += "\
  angmap.insert(make_pair(\"s\", 0));\n\
  angmap.insert(make_pair(\"p\", 1));\n\
  angmap.insert(make_pair(\"d\", 2));\n\
  angmap.insert(make_pair(\"f\", 3));\n\
  angmap.insert(make_pair(\"g\", 4));\n\
  angmap.insert(make_pair(\"h\", 5));\n\
  angmap.insert(make_pair(\"i\", 6));\n\
// Since they are not implemented yet\n\
//angmap.insert(make_pair(\"j\", 7));\n\
//angmap.insert(make_pair(\"k\", 8));\n\
//angmap.insert(make_pair(\"l\", 9));\n\
\n\
}\n\
\n\
AtomMap::~AtomMap () {\n\
\n\
}\n\
\n\
int AtomMap::angular_number(const string input) const {\n\
  auto miter = angmap.find(input);\n\
  if (miter == angmap.end()) throw runtime_error(\"Unknown angular number in a basis set file.\");\n\
  return miter->second;\n\
}\n\
\n\
\n\
int AtomMap::atom_number(const string input) const {\n\
  auto miter = atommap.find(input);\n\
  if (miter == atommap.end()) throw runtime_error(\"Unknown Atom number in a basis set file.\");\n\
  return miter->second;\n\
}\n\
\n\
const string AtomMap::angular_string(const int input) {\n\
  for(auto& m : angmap) {\n\
    if(m.second == input) { return m.first; }\n\
  }\n\
  return \"X\";\n\
}\n"

fl = open("atommap.cc", "w")
fl.write(ss)
fl.close()
