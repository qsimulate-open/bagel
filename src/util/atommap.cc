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


#include <src/util/atommap.h>
#include <src/util/constants.h>
#include <sstream>
#include <stdexcept>

using namespace std;
using namespace bagel;

AtomMap::AtomMap () {
  atommap.emplace("q", 0);
  atommap.emplace("h", 1);
  atommap.emplace("he", 2);
  atommap.emplace("li", 3);
  atommap.emplace("be", 4);
  atommap.emplace("b", 5);
  atommap.emplace("c", 6);
  atommap.emplace("n", 7);
  atommap.emplace("o", 8);
  atommap.emplace("f", 9);
  atommap.emplace("ne", 10);
  atommap.emplace("na", 11);
  atommap.emplace("mg", 12);
  atommap.emplace("al", 13);
  atommap.emplace("si", 14);
  atommap.emplace("p", 15);
  atommap.emplace("s", 16);
  atommap.emplace("cl", 17);
  atommap.emplace("ar", 18);
  atommap.emplace("k", 19);
  atommap.emplace("ca", 20);
  atommap.emplace("sc", 21);
  atommap.emplace("ti", 22);
  atommap.emplace("v", 23);
  atommap.emplace("cr", 24);
  atommap.emplace("mn", 25);
  atommap.emplace("fe", 26);
  atommap.emplace("co", 27);
  atommap.emplace("ni", 28);
  atommap.emplace("cu", 29);
  atommap.emplace("zn", 30);
  atommap.emplace("ga", 31);
  atommap.emplace("ge", 32);
  atommap.emplace("as", 33);
  atommap.emplace("se", 34);
  atommap.emplace("br", 35);
  atommap.emplace("kr", 36);
  atommap.emplace("rb", 37);
  atommap.emplace("sr", 38);
  atommap.emplace("y", 39);
  atommap.emplace("zr", 40);
  atommap.emplace("nb", 41);
  atommap.emplace("mo", 42);
  atommap.emplace("tc", 43);
  atommap.emplace("ru", 44);
  atommap.emplace("rh", 45);
  atommap.emplace("pd", 46);
  atommap.emplace("ag", 47);
  atommap.emplace("cd", 48);
  atommap.emplace("in", 49);
  atommap.emplace("sn", 50);
  atommap.emplace("sb", 51);
  atommap.emplace("te", 52);
  atommap.emplace("i", 53);
  atommap.emplace("xe", 54);
  atommap.emplace("cs", 55);
  atommap.emplace("ba", 56);
  atommap.emplace("la", 57);
  atommap.emplace("ce", 58);
  atommap.emplace("pr", 59);
  atommap.emplace("nd", 60);
  atommap.emplace("pm", 61);
  atommap.emplace("sm", 62);
  atommap.emplace("eu", 63);
  atommap.emplace("gd", 64);
  atommap.emplace("tb", 65);
  atommap.emplace("dy", 66);
  atommap.emplace("ho", 67);
  atommap.emplace("er", 68);
  atommap.emplace("tm", 69);
  atommap.emplace("yb", 70);
  atommap.emplace("lu", 71);
  atommap.emplace("hf", 72);
  atommap.emplace("ta", 73);
  atommap.emplace("w", 74);
  atommap.emplace("re", 75);
  atommap.emplace("os", 76);
  atommap.emplace("ir", 77);
  atommap.emplace("pt", 78);
  atommap.emplace("au", 79);
  atommap.emplace("hg", 80);
  atommap.emplace("tl", 81);
  atommap.emplace("pb", 82);
  atommap.emplace("bi", 83);
  atommap.emplace("po", 84);
  atommap.emplace("at", 85);
  atommap.emplace("rn", 86);
  atommap.emplace("fr", 87);
  atommap.emplace("ra", 88);
  atommap.emplace("ac", 89);
  atommap.emplace("th", 90);
  atommap.emplace("pa", 91);
  atommap.emplace("u", 92);
  atommap.emplace("np", 93);
  atommap.emplace("pu", 94);
  atommap.emplace("am", 95);
  atommap.emplace("cm", 96);
  atommap.emplace("bk", 97);
  atommap.emplace("cf", 98);
  atommap.emplace("es", 99);
  atommap.emplace("fm", 100);
  atommap.emplace("md", 101);
  atommap.emplace("no", 102);
  atommap.emplace("lr", 103);
  atommap.emplace("rf", 104);
  atommap.emplace("db", 105);
  atommap.emplace("sg", 106);
  atommap.emplace("bh", 107);
  atommap.emplace("hs", 108);
  atommap.emplace("mt", 109);
  atommap.emplace("ds", 110);
  atommap.emplace("rg", 111);
  atommap.emplace("cn", 112);
  atommap.emplace("uut", 113);
  atommap.emplace("fl", 114);
  atommap.emplace("uup", 115);
  atommap.emplace("lv", 116);
  atommap.emplace("uus", 117);
  atommap.emplace("uuo", 118);

  cov_radii.emplace("h",  0.31);
  cov_radii.emplace("he", 0.28);
  cov_radii.emplace("li", 1.28);
  cov_radii.emplace("be", 0.96);
  cov_radii.emplace("b",  0.84);
  cov_radii.emplace("c",  0.76);
  cov_radii.emplace("n",  0.71);
  cov_radii.emplace("o",  0.66);
  cov_radii.emplace("f",  0.57);
  cov_radii.emplace("ne", 0.58);
  cov_radii.emplace("na", 1.66);
  cov_radii.emplace("mg", 1.41);
  cov_radii.emplace("al", 1.21);
  cov_radii.emplace("si", 1.11);
  cov_radii.emplace("p",  1.07);
  cov_radii.emplace("s",  1.05);
  cov_radii.emplace("cl", 1.02);
  cov_radii.emplace("ar", 1.06);
  cov_radii.emplace("k",  2.03);
  cov_radii.emplace("ca", 1.76);
  cov_radii.emplace("sc", 1.70);
  cov_radii.emplace("ti", 1.60);
  cov_radii.emplace("v",  1.53);
  cov_radii.emplace("cr", 1.39);
  cov_radii.emplace("mn", 1.39);
  cov_radii.emplace("fe", 1.32);
  cov_radii.emplace("co", 1.26);
  cov_radii.emplace("ni", 1.24);
  cov_radii.emplace("cu", 1.32);
  cov_radii.emplace("zn", 1.22);
  cov_radii.emplace("ga", 1.22);
  cov_radii.emplace("ge", 1.20);
  cov_radii.emplace("as", 1.19);
  cov_radii.emplace("se", 1.20);
  cov_radii.emplace("br", 1.20);
  cov_radii.emplace("kr", 1.16);
  cov_radii.emplace("rb", 2.20);
  cov_radii.emplace("sr", 1.95);
  cov_radii.emplace("y",  1.90);
  cov_radii.emplace("zr", 1.75);
  cov_radii.emplace("nb", 1.64);
  cov_radii.emplace("mo", 1.54);
  cov_radii.emplace("tc", 1.47);
  cov_radii.emplace("ru", 1.46);
  cov_radii.emplace("rh", 1.42);
  cov_radii.emplace("pd", 1.39);
  cov_radii.emplace("ag", 1.45);
  cov_radii.emplace("cd", 1.44);
  cov_radii.emplace("in", 1.42);
  cov_radii.emplace("sn", 1.39);
  cov_radii.emplace("sb", 1.39);
  cov_radii.emplace("te", 1.38);
  cov_radii.emplace("i",  1.39);
  cov_radii.emplace("xe", 1.40);
  cov_radii.emplace("cs", 2.44);
  cov_radii.emplace("ba", 2.15);
  cov_radii.emplace("la", 2.07);
  cov_radii.emplace("ce", 2.04);
  cov_radii.emplace("pr", 2.03);
  cov_radii.emplace("nd", 2.01);
  cov_radii.emplace("pm", 1.99);
  cov_radii.emplace("sm", 1.98);
  cov_radii.emplace("eu", 1.98);
  cov_radii.emplace("gd", 1.96);
  cov_radii.emplace("tb", 1.94);
  cov_radii.emplace("dy", 1.92);
  cov_radii.emplace("ho", 1.92);
  cov_radii.emplace("er", 1.89);
  cov_radii.emplace("tm", 1.90);
  cov_radii.emplace("yb", 1.87);
  cov_radii.emplace("lu", 1.87);
  cov_radii.emplace("hf", 1.75);
  cov_radii.emplace("ta", 1.70);
  cov_radii.emplace("w",  1.62);
  cov_radii.emplace("re", 1.51);
  cov_radii.emplace("os", 1.44);
  cov_radii.emplace("ir", 1.41);
  cov_radii.emplace("pt", 1.36);
  cov_radii.emplace("au", 1.36);
  cov_radii.emplace("hg", 1.32);
  cov_radii.emplace("tl", 1.45);
  cov_radii.emplace("pb", 1.46);
  cov_radii.emplace("bi", 1.48);
  cov_radii.emplace("po", 1.40);
  cov_radii.emplace("at", 1.50);
  cov_radii.emplace("rn", 1.50);
  cov_radii.emplace("fr", 2.60);
  cov_radii.emplace("ra", 2.21);
  cov_radii.emplace("ac", 2.15);
  cov_radii.emplace("th", 2.06);
  cov_radii.emplace("pa", 2.00);
  cov_radii.emplace("u",  1.96);
  cov_radii.emplace("np", 1.90);
  cov_radii.emplace("pu", 1.87);
  cov_radii.emplace("am", 1.80);
  cov_radii.emplace("cm", 1.69);

  // atom sizes (Bragg-Slater radii)
  bsradii.emplace("h", 0.25);
  bsradii.emplace("he", 0.25);
  bsradii.emplace("li", 1.45);
  bsradii.emplace("be", 1.05);
  bsradii.emplace("b", 0.85);
  bsradii.emplace("c", 0.7);
  bsradii.emplace("n", 0.65);
  bsradii.emplace("o", 0.6);
  bsradii.emplace("f", 0.5);
  bsradii.emplace("ne", 0.45);
  bsradii.emplace("na", 1.8);
  bsradii.emplace("mg", 1.5);
  bsradii.emplace("al", 1.25);
  bsradii.emplace("si", 1.1);
  bsradii.emplace("p", 1.0);
  bsradii.emplace("s", 1.0);
  bsradii.emplace("cl", 1.0);
  bsradii.emplace("ar", 1.0);
  bsradii.emplace("k", 2.2);
  bsradii.emplace("ca", 1.8);
  bsradii.emplace("sc", 1.6);
  bsradii.emplace("ti", 1.4);
  bsradii.emplace("v", 1.35);
  bsradii.emplace("cr", 1.4);
  bsradii.emplace("mn", 1.4);
  bsradii.emplace("fe", 1.4);
  bsradii.emplace("co", 1.35);
  bsradii.emplace("ni", 1.35);
  bsradii.emplace("cu", 1.35);
  bsradii.emplace("zn", 1.35);
  bsradii.emplace("ga", 1.3);
  bsradii.emplace("ge", 1.25);
  bsradii.emplace("as", 1.15);
  bsradii.emplace("se", 1.15);
  bsradii.emplace("br", 1.15);
  bsradii.emplace("kr", 1.15);
  bsradii.emplace("rb", 2.35);
  bsradii.emplace("sr", 2.0);
  bsradii.emplace("y", 1.8);
  bsradii.emplace("zr", 1.55);
  bsradii.emplace("nb", 1.45);
  bsradii.emplace("mo", 1.45);
  bsradii.emplace("tc", 1.35);
  bsradii.emplace("ru", 1.3);
  bsradii.emplace("rh", 1.35);
  bsradii.emplace("pd", 1.4);
  bsradii.emplace("ag", 1.6);
  bsradii.emplace("cd", 1.55);
  bsradii.emplace("in", 1.55);
  bsradii.emplace("sn", 1.45);
  bsradii.emplace("sb", 1.45);
  bsradii.emplace("te", 1.4);
  bsradii.emplace("i", 1.4);
  bsradii.emplace("xe", 1.4);
  bsradii.emplace("cs", 2.6);
  bsradii.emplace("ba", 2.15);
  bsradii.emplace("la", 1.95);
  bsradii.emplace("ce", 1.85);
  bsradii.emplace("pr", 1.85);
  bsradii.emplace("nd", 1.85);
  bsradii.emplace("pm", 1.85);
  bsradii.emplace("sm", 1.85);
  bsradii.emplace("eu", 1.85);
  bsradii.emplace("gd", 1.8);
  bsradii.emplace("tb", 1.75);
  bsradii.emplace("dy", 1.75);
  bsradii.emplace("ho", 1.75);
  bsradii.emplace("er", 1.75);
  bsradii.emplace("tm", 1.75);
  bsradii.emplace("yb", 1.75);
  bsradii.emplace("lu", 1.75);
  bsradii.emplace("hf", 1.55);
  bsradii.emplace("ta", 1.45);
  bsradii.emplace("w", 1.35);
  bsradii.emplace("re", 1.35);
  bsradii.emplace("os", 1.3);
  bsradii.emplace("ir", 1.35);
  bsradii.emplace("pt", 1.35);
  bsradii.emplace("au", 1.35);
  bsradii.emplace("hg", 1.5);
  bsradii.emplace("tl", 1.9);
  bsradii.emplace("pb", 1.8);
  bsradii.emplace("bi", 1.6);
  bsradii.emplace("po", 1.9);
  bsradii.emplace("at", 1.9);
  bsradii.emplace("rn", 1.9);
  bsradii.emplace("fr", 2.85);
  bsradii.emplace("ra", 2.15);
  bsradii.emplace("ac", 1.95);
  bsradii.emplace("th", 1.8);
  bsradii.emplace("pa", 1.8);
  bsradii.emplace("u", 1.75);
  bsradii.emplace("np", 1.75);
  bsradii.emplace("pu", 1.75);
  bsradii.emplace("am", 1.75);
  bsradii.emplace("cm", 1.75);
  bsradii.emplace("bk", 1.75);
  bsradii.emplace("cf", 1.75);
  bsradii.emplace("es", 1.75);
  bsradii.emplace("fm", 1.75);
  bsradii.emplace("md", 1.75);
  bsradii.emplace("no", 1.75);
  bsradii.emplace("lr", 1.75);
#if 0
  bsradii.emplace("rf", 104);
  bsradii.emplace("db", 105);
  bsradii.emplace("sg", 106);
  bsradii.emplace("bh", 107);
  bsradii.emplace("hs", 108);
  bsradii.emplace("mt", 109);
  bsradii.emplace("ds", 110);
  bsradii.emplace("rg", 111);
  bsradii.emplace("cn", 112);
  bsradii.emplace("uut", 113);
  bsradii.emplace("fl", 114);
  bsradii.emplace("uup", 115);
  bsradii.emplace("lv", 116);
  bsradii.emplace("uus", 117);
  bsradii.emplace("uuo", 118);
#endif

  angmap.emplace("s", 0);
  angmap.emplace("p", 1);
  angmap.emplace("d", 2);
  angmap.emplace("f", 3);
  angmap.emplace("g", 4);
  angmap.emplace("h", 5);
  angmap.emplace("i", 6);
// Since they are not implemented yet
//angmap.emplace("j", 7);
//angmap.emplace("k", 8);
//angmap.emplace("l", 9);

  nclosed.emplace("h",  make_tuple(0,0,0,0));
  nclosed.emplace("he", make_tuple(2,0,0,0));
  nclosed.emplace("li", make_tuple(2,0,0,0));
  nclosed.emplace("be", make_tuple(4,0,0,0));
  nclosed.emplace("b",  make_tuple(4,0,0,0));
  nclosed.emplace("c",  make_tuple(4,0,0,0));
  nclosed.emplace("n",  make_tuple(4,0,0,0));
  nclosed.emplace("o",  make_tuple(4,0,0,0));
  nclosed.emplace("f",  make_tuple(4,0,0,0));
  nclosed.emplace("ne", make_tuple(4,6,0,0));
  nclosed.emplace("na", make_tuple(4,6,0,0));
  nclosed.emplace("mg", make_tuple(6,6,0,0));
  nclosed.emplace("al", make_tuple(6,6,0,0));
  nclosed.emplace("si", make_tuple(6,6,0,0));
  nclosed.emplace("p",  make_tuple(6,6,0,0));
  nclosed.emplace("s",  make_tuple(6,6,0,0));
  nclosed.emplace("cl", make_tuple(6,6,0,0));
  nclosed.emplace("ar", make_tuple(6,12,0,0));
  nclosed.emplace("k",  make_tuple(6,12,0,0));
  nclosed.emplace("ca", make_tuple(8,12,0,0));
  nclosed.emplace("sc", make_tuple(8,12,0,0));
  nclosed.emplace("ti", make_tuple(8,12,0,0));
  nclosed.emplace("v",  make_tuple(8,12,0,0));
  nclosed.emplace("cr", make_tuple(8,12,0,0));
  nclosed.emplace("mn", make_tuple(8,12,0,0));
  nclosed.emplace("fe", make_tuple(8,12,0,0));
  nclosed.emplace("co", make_tuple(8,12,0,0));
  nclosed.emplace("ni", make_tuple(8,12,0,0));
  nclosed.emplace("cu", make_tuple(8,12,0,0));
  nclosed.emplace("zn", make_tuple(8,12,10,0));
  nclosed.emplace("ga", make_tuple(8,12,10,0));
  nclosed.emplace("ge", make_tuple(8,12,10,0));
  nclosed.emplace("as", make_tuple(8,12,10,0));
  nclosed.emplace("se", make_tuple(8,12,10,0));
  nclosed.emplace("br", make_tuple(8,12,10,0));
  nclosed.emplace("kr", make_tuple(8,18,10,0));
  nclosed.emplace("rb", make_tuple(8,18,10,0));
  nclosed.emplace("sr", make_tuple(10,18,10,0));
  nclosed.emplace("y",  make_tuple(10,18,10,0));
  nclosed.emplace("zr", make_tuple(10,18,10,0));
  nclosed.emplace("nb", make_tuple(10,18,10,0));
  nclosed.emplace("mo", make_tuple(10,18,10,0));
  nclosed.emplace("tc", make_tuple(10,18,10,0));
  nclosed.emplace("ru", make_tuple(10,18,10,0));
  nclosed.emplace("rh", make_tuple(10,18,10,0));
  nclosed.emplace("pd", make_tuple(10,18,10,0));
  nclosed.emplace("ag", make_tuple(10,18,10,0));
  nclosed.emplace("cd", make_tuple(10,18,20,0));
  nclosed.emplace("in", make_tuple(10,18,20,0));
  nclosed.emplace("sn", make_tuple(10,18,20,0));
  nclosed.emplace("sb", make_tuple(10,18,20,0));
  nclosed.emplace("te", make_tuple(10,18,20,0));
  nclosed.emplace("i",  make_tuple(10,18,20,0));
  nclosed.emplace("xe", make_tuple(10,24,20,0));
  nclosed.emplace("cs", make_tuple(10,24,20,0));
  nclosed.emplace("ba", make_tuple(10,24,20,0));
  nclosed.emplace("la", make_tuple(10,24,20,0));
  nclosed.emplace("ce", make_tuple(10,24,20,0));
  nclosed.emplace("pr", make_tuple(10,24,20,0));
  nclosed.emplace("nd", make_tuple(10,24,20,0));
  nclosed.emplace("pm", make_tuple(10,24,20,0));
  nclosed.emplace("sm", make_tuple(12,24,20,0));
  nclosed.emplace("eu", make_tuple(12,24,20,0));
  nclosed.emplace("gd", make_tuple(12,24,20,0));
  nclosed.emplace("tb", make_tuple(12,24,20,0));
  nclosed.emplace("dy", make_tuple(12,24,20,0));
  nclosed.emplace("ho", make_tuple(12,24,20,0));
  nclosed.emplace("er", make_tuple(12,24,20,0));
  nclosed.emplace("tm", make_tuple(12,24,20,0));
  nclosed.emplace("yb", make_tuple(12,24,20,14));
  nclosed.emplace("lu", make_tuple(12,24,20,14));
  nclosed.emplace("hf", make_tuple(12,24,20,14));
  nclosed.emplace("ta", make_tuple(12,24,20,14));
  nclosed.emplace("w",  make_tuple(12,24,20,14));
  nclosed.emplace("re", make_tuple(12,24,20,14));
  nclosed.emplace("os", make_tuple(12,24,20,14));
  nclosed.emplace("ir", make_tuple(12,24,20,14));
  nclosed.emplace("pt", make_tuple(10,24,20,14));
  nclosed.emplace("au", make_tuple(10,24,30,14));
  nclosed.emplace("hg", make_tuple(12,24,30,14));
  nclosed.emplace("tl", make_tuple(12,24,30,14));
  nclosed.emplace("pb", make_tuple(12,24,30,14));
  nclosed.emplace("bi", make_tuple(12,24,30,14));
  nclosed.emplace("po", make_tuple(12,24,30,14));
  nclosed.emplace("at", make_tuple(12,24,30,14));
  nclosed.emplace("rn", make_tuple(12,30,30,14));
  nclosed.emplace("fr", make_tuple(12,30,30,14));
  nclosed.emplace("ra", make_tuple(14,30,30,14));
  nclosed.emplace("ac", make_tuple(14,30,30,14));
  nclosed.emplace("th", make_tuple(14,30,30,14));
  nclosed.emplace("pa", make_tuple(14,30,30,14));
  nclosed.emplace("u",  make_tuple(14,30,30,14));
  nclosed.emplace("np", make_tuple(14,30,30,14));
  nclosed.emplace("pu", make_tuple(14,30,30,14));
  nclosed.emplace("am", make_tuple(14,30,30,14));
  nclosed.emplace("cm", make_tuple(14,30,30,14));
  nclosed.emplace("bk", make_tuple(14,30,30,14));
  nclosed.emplace("cf", make_tuple(14,30,30,14));
  nclosed.emplace("es", make_tuple(14,30,30,14));
  nclosed.emplace("fm", make_tuple(14,30,30,14));
  nclosed.emplace("md", make_tuple(14,30,30,14));
  nclosed.emplace("no", make_tuple(14,30,30,14));
  nclosed.emplace("lr", make_tuple(14,30,30,28));
  nclosed.emplace("rf", make_tuple(14,30,30,28));
  nclosed.emplace("db", make_tuple(14,30,30,28));
  nclosed.emplace("sg", make_tuple(14,30,30,28));
  nclosed.emplace("bh", make_tuple(14,30,30,28));
  nclosed.emplace("hs", make_tuple(14,30,30,28));
  nclosed.emplace("mt", make_tuple(14,30,30,28));
  nclosed.emplace("ds", make_tuple(14,30,30,28));
  nclosed.emplace("rg", make_tuple(14,30,30,28));
  nclosed.emplace("cn", make_tuple(14,30,40,28));
  nclosed.emplace("uut", make_tuple(14,30,40,28));
  nclosed.emplace("fl",  make_tuple(14,30,40,28));
  nclosed.emplace("uup", make_tuple(14,30,40,28));
  nclosed.emplace("lv",  make_tuple(14,30,40,28));
  nclosed.emplace("uus", make_tuple(14,30,40,28));
  nclosed.emplace("uuo", make_tuple(14,36,40,28));

  nopen.emplace("h",  make_tuple(1,0,0,0));
  nopen.emplace("he", make_tuple(0,0,0,0));
  nopen.emplace("li", make_tuple(1,0,0,0));
  nopen.emplace("be", make_tuple(0,0,0,0));
  nopen.emplace("b",  make_tuple(0,1,0,0));
  nopen.emplace("c",  make_tuple(0,2,0,0));
  nopen.emplace("n",  make_tuple(0,3,0,0));
  nopen.emplace("o",  make_tuple(0,4,0,0));
  nopen.emplace("f",  make_tuple(0,5,0,0));
  nopen.emplace("ne", make_tuple(0,0,0,0));
  nopen.emplace("na", make_tuple(1,0,0,0));
  nopen.emplace("mg", make_tuple(0,0,0,0));
  nopen.emplace("al", make_tuple(0,1,0,0));
  nopen.emplace("si", make_tuple(0,2,0,0));
  nopen.emplace("p",  make_tuple(0,3,0,0));
  nopen.emplace("s",  make_tuple(0,4,0,0));
  nopen.emplace("cl", make_tuple(0,5,0,0));
  nopen.emplace("ar", make_tuple(0,0,0,0));
  nopen.emplace("k",  make_tuple(1,0,0,0));
  nopen.emplace("ca", make_tuple(0,0,0,0));
  nopen.emplace("sc", make_tuple(0,0,1,0));
  nopen.emplace("ti", make_tuple(0,0,2,0));
  nopen.emplace("v",  make_tuple(0,0,3,0));
  nopen.emplace("cr", make_tuple(0,0,4,0));
  nopen.emplace("mn", make_tuple(0,0,5,0));
  nopen.emplace("fe", make_tuple(0,0,6,0));
  nopen.emplace("co", make_tuple(0,0,7,0));
  nopen.emplace("ni", make_tuple(0,0,8,0));
  nopen.emplace("cu", make_tuple(0,0,9,0));
  nopen.emplace("zn", make_tuple(0,0,0,0));
  nopen.emplace("ga", make_tuple(0,1,0,0));
  nopen.emplace("ge", make_tuple(0,2,0,0));
  nopen.emplace("as", make_tuple(0,3,0,0));
  nopen.emplace("se", make_tuple(0,4,0,0));
  nopen.emplace("br", make_tuple(0,5,0,0));
  nopen.emplace("kr", make_tuple(0,0,0,0));
  nopen.emplace("rb", make_tuple(1,0,0,0));
  nopen.emplace("sr", make_tuple(0,0,0,0));
  nopen.emplace("y",  make_tuple(0,0,1,0));
  nopen.emplace("zr", make_tuple(0,0,2,0));
  nopen.emplace("nb", make_tuple(0,0,3,0));
  nopen.emplace("mo", make_tuple(0,0,4,0));
  nopen.emplace("tc", make_tuple(0,0,5,0));
  nopen.emplace("ru", make_tuple(0,0,6,0));
  nopen.emplace("rh", make_tuple(0,0,7,0));
  nopen.emplace("pd", make_tuple(0,0,8,0));
  nopen.emplace("ag", make_tuple(0,0,9,0));
  nopen.emplace("cd", make_tuple(0,0,0,0));
  nopen.emplace("in", make_tuple(0,1,0,0));
  nopen.emplace("sn", make_tuple(0,2,0,0));
  nopen.emplace("sb", make_tuple(0,3,0,0));
  nopen.emplace("te", make_tuple(0,4,0,0));
  nopen.emplace("i",  make_tuple(0,5,0,0));
  nopen.emplace("xe", make_tuple(0,0,0,0));
  nopen.emplace("cs", make_tuple(1,0,0,0));
  nopen.emplace("ba", make_tuple(0,0,0,0));
  nopen.emplace("la", make_tuple(0,0,1,0));
  nopen.emplace("ce", make_tuple(0,0,1,1));
  nopen.emplace("pr", make_tuple(0,0,0,3));
  nopen.emplace("nd", make_tuple(0,0,0,4));
  nopen.emplace("pm", make_tuple(0,0,0,5));
  nopen.emplace("sm", make_tuple(0,0,0,6));
  nopen.emplace("eu", make_tuple(0,0,0,7));
  nopen.emplace("gd", make_tuple(0,0,1,7));
  nopen.emplace("tb", make_tuple(0,0,0,9));
  nopen.emplace("dy", make_tuple(0,0,0,10));
  nopen.emplace("ho", make_tuple(0,0,0,11));
  nopen.emplace("er", make_tuple(0,0,0,12));
  nopen.emplace("tm", make_tuple(0,0,0,13));
  nopen.emplace("yb", make_tuple(0,0,0,0));
  nopen.emplace("lu", make_tuple(0,0,1,0));
  nopen.emplace("hf", make_tuple(0,0,2,0));
  nopen.emplace("ta", make_tuple(0,0,3,0));
  nopen.emplace("w",  make_tuple(0,0,4,0));
  nopen.emplace("re", make_tuple(0,0,5,0));
  nopen.emplace("os", make_tuple(0,0,6,0));
  nopen.emplace("ir", make_tuple(0,0,7,0));
  nopen.emplace("pt", make_tuple(1,0,9,0));
  nopen.emplace("au", make_tuple(1,0,0,0));
  nopen.emplace("hg", make_tuple(0,0,0,0));
  nopen.emplace("tl", make_tuple(0,1,0,0));
  nopen.emplace("pb", make_tuple(0,2,0,0));
  nopen.emplace("bi", make_tuple(0,3,0,0));
  nopen.emplace("po", make_tuple(0,4,0,0));
  nopen.emplace("at", make_tuple(0,5,0,0));
  nopen.emplace("rn", make_tuple(0,0,0,0));
  nopen.emplace("fr", make_tuple(1,0,0,0));
  nopen.emplace("ra", make_tuple(0,0,0,0));
  nopen.emplace("ac", make_tuple(0,0,1,0));
  nopen.emplace("th", make_tuple(0,0,2,0));
  nopen.emplace("pa", make_tuple(0,0,1,2));
  nopen.emplace("u",  make_tuple(0,0,1,3));
  nopen.emplace("np", make_tuple(0,0,1,4));
  nopen.emplace("pu", make_tuple(0,0,0,6));
  nopen.emplace("am", make_tuple(0,0,0,7));
  nopen.emplace("cm", make_tuple(0,0,1,7));
  nopen.emplace("bk", make_tuple(0,0,0,9));
  nopen.emplace("cf", make_tuple(0,0,0,10));
  nopen.emplace("es", make_tuple(0,0,0,11));
  nopen.emplace("fm", make_tuple(0,0,0,12));
  nopen.emplace("md", make_tuple(0,0,0,13));
  nopen.emplace("no", make_tuple(0,0,0,0));
  nopen.emplace("lr", make_tuple(0,1,0,0));
  nopen.emplace("rf", make_tuple(0,0,2,0));
  nopen.emplace("db", make_tuple(0,0,3,0));
  nopen.emplace("sg", make_tuple(0,0,4,0));
  nopen.emplace("bh", make_tuple(0,0,5,0));
  nopen.emplace("hs", make_tuple(0,0,6,0));
  nopen.emplace("mt", make_tuple(0,0,7,0));
  nopen.emplace("ds", make_tuple(0,0,8,0));
  nopen.emplace("rg", make_tuple(0,0,9,0));
  nopen.emplace("cn", make_tuple(0,0,0,0));
  nopen.emplace("uut", make_tuple(0,1,0,0));
  nopen.emplace("fl",  make_tuple(0,2,0,0));
  nopen.emplace("uup", make_tuple(0,3,0,0));
  nopen.emplace("lv",  make_tuple(0,4,0,0));
  nopen.emplace("uus", make_tuple(0,5,0,0));
  nopen.emplace("uuo", make_tuple(0,0,0,0));

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
  return miter->second/au2angstrom__;
}


double AtomMap::cov_radius(const string input) const {
  auto miter = cov_radii.find(input);
  if (miter == cov_radii.end()) throw runtime_error("Unknown atom (Covalent radii).");
  return miter->second/au2angstrom__;
}


tuple<int,int,int,int> AtomMap::num_closed(const string input) const {
  auto miter = nclosed.find(input);
  if (miter == nclosed.end()) throw runtime_error("Unknown atom (#closed shell).");
  return miter->second;
}


tuple<int,int,int,int> AtomMap::num_open(const string input) const {
  auto miter = nopen.find(input);
  if (miter == nopen.end()) throw runtime_error("Unknown atom (#open shell).");
  return miter->second;
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
