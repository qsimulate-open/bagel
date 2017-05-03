//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: atommap.cc
// Copyright (C) 2009 Toru Shiozaki
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

  // finite nuclear exponents, from L. Visscher and K. G. Dyall, At. Data Nucl. Data Tables 67, 207 (1997)
  nuclear_exponents.emplace("q",  0.0000000000E+00);
  nuclear_exponents.emplace("h",  2.1248239171E+09);
  nuclear_exponents.emplace("he", 1.1671538870E+09);
  nuclear_exponents.emplace("li", 8.9266848806E+08);
  nuclear_exponents.emplace("be", 7.8788802914E+08);
  nuclear_exponents.emplace("b",  7.1178709563E+08);
  nuclear_exponents.emplace("c",  6.8077502929E+08);
  nuclear_exponents.emplace("n",  6.2865615725E+08);
  nuclear_exponents.emplace("o",  5.8631436655E+08);
  nuclear_exponents.emplace("f",  5.3546911034E+08);
  nuclear_exponents.emplace("ne", 5.2105715255E+08);
  nuclear_exponents.emplace("na", 4.8349721509E+08);
  nuclear_exponents.emplace("mg", 4.7254270882E+08);
  nuclear_exponents.emplace("al", 4.4335984491E+08);
  nuclear_exponents.emplace("si", 4.3467748823E+08);
  nuclear_exponents.emplace("p",  4.1117553148E+08);
  nuclear_exponents.emplace("s",  4.0407992047E+08);
  nuclear_exponents.emplace("cl", 3.8463852873E+08);
  nuclear_exponents.emplace("ar", 3.5722217300E+08);
  nuclear_exponents.emplace("k",  3.6228128110E+08);
  nuclear_exponents.emplace("ca", 3.5722217300E+08);
  nuclear_exponents.emplace("sc", 3.3451324570E+08);
  nuclear_exponents.emplace("ti", 3.2263108827E+08);
  nuclear_exponents.emplace("v",  3.1181925878E+08);
  nuclear_exponents.emplace("cr", 3.0842641793E+08);
  nuclear_exponents.emplace("mn", 2.9881373610E+08);
  nuclear_exponents.emplace("fe", 2.9578406371E+08);
  nuclear_exponents.emplace("co", 2.8716667270E+08);
  nuclear_exponents.emplace("ni", 2.8996391416E+08);
  nuclear_exponents.emplace("cu", 2.7665979354E+08);
  nuclear_exponents.emplace("zn", 2.7419021043E+08);
  nuclear_exponents.emplace("ga", 2.6267002737E+08);
  nuclear_exponents.emplace("ge", 2.5235613399E+08);
  nuclear_exponents.emplace("as", 2.5042024280E+08);
  nuclear_exponents.emplace("se", 2.4130163719E+08);
  nuclear_exponents.emplace("br", 2.4305454351E+08);
  nuclear_exponents.emplace("kr", 2.3461213272E+08);
  nuclear_exponents.emplace("rb", 2.3301551109E+08);
  nuclear_exponents.emplace("sr", 2.2839354730E+08);
  nuclear_exponents.emplace("y",  2.2690621893E+08);
  nuclear_exponents.emplace("zr", 2.2544431039E+08);
  nuclear_exponents.emplace("nb", 2.2120420724E+08);
  nuclear_exponents.emplace("mo", 2.1458511597E+08);
  nuclear_exponents.emplace("tc", 2.1458511597E+08);
  nuclear_exponents.emplace("ru", 2.0965270287E+08);
  nuclear_exponents.emplace("rh", 2.0846586999E+08);
  nuclear_exponents.emplace("pd", 2.0500935221E+08);
  nuclear_exponents.emplace("ag", 2.0389047621E+08);
  nuclear_exponents.emplace("cd", 1.9648639618E+08);
  nuclear_exponents.emplace("in", 1.9548577691E+08);
  nuclear_exponents.emplace("sn", 1.9067718154E+08);
  nuclear_exponents.emplace("sb", 1.8975246242E+08);
  nuclear_exponents.emplace("te", 1.8193056289E+08);
  nuclear_exponents.emplace("i",  1.8444240538E+08);
  nuclear_exponents.emplace("xe", 1.8030529331E+08);
  nuclear_exponents.emplace("cs", 1.7950688281E+08);
  nuclear_exponents.emplace("ba", 1.7565009043E+08);
  nuclear_exponents.emplace("la", 1.7490463170E+08);
  nuclear_exponents.emplace("ce", 1.7416744147E+08);
  nuclear_exponents.emplace("pr", 1.7343837120E+08);
  nuclear_exponents.emplace("nd", 1.7129844956E+08);
  nuclear_exponents.emplace("pm", 1.7060044589E+08);
  nuclear_exponents.emplace("sm", 1.6591550422E+08);
  nuclear_exponents.emplace("eu", 1.6527352089E+08);
  nuclear_exponents.emplace("gd", 1.6215880671E+08);
  nuclear_exponents.emplace("tb", 1.6155419421E+08);
  nuclear_exponents.emplace("dy", 1.5977529080E+08);
  nuclear_exponents.emplace("ho", 1.5977529080E+08);
  nuclear_exponents.emplace("er", 1.5636673634E+08);
  nuclear_exponents.emplace("tm", 1.5581702004E+08);
  nuclear_exponents.emplace("yb", 1.5314257850E+08);
  nuclear_exponents.emplace("lu", 1.5262201512E+08);
  nuclear_exponents.emplace("hf", 1.5008710340E+08);
  nuclear_exponents.emplace("ta", 1.4959325643E+08);
  nuclear_exponents.emplace("w",  1.4813689532E+08);
  nuclear_exponents.emplace("re", 1.4671710337E+08);
  nuclear_exponents.emplace("os", 1.4442808782E+08);
  nuclear_exponents.emplace("ir", 1.4398142103E+08);
  nuclear_exponents.emplace("pt", 1.4309883584E+08);
  nuclear_exponents.emplace("au", 1.4223027307E+08);
  nuclear_exponents.emplace("hg", 1.4011788914E+08);
  nuclear_exponents.emplace("tl", 1.3888925203E+08);
  nuclear_exponents.emplace("pb", 1.3768840081E+08);
  nuclear_exponents.emplace("bi", 1.3729411599E+08);
  nuclear_exponents.emplace("po", 1.3729411599E+08);
  nuclear_exponents.emplace("at", 1.3690277000E+08);
  nuclear_exponents.emplace("rn", 1.3242350205E+08);
  nuclear_exponents.emplace("fr", 1.3206733609E+08);
  nuclear_exponents.emplace("ra", 1.3101367628E+08);
  nuclear_exponents.emplace("ac", 1.3066730974E+08);
  nuclear_exponents.emplace("th", 1.2897067480E+08);
  nuclear_exponents.emplace("pa", 1.2930539512E+08);
  nuclear_exponents.emplace("u",  1.2700881714E+08);
  nuclear_exponents.emplace("np", 1.2733038109E+08);
  nuclear_exponents.emplace("pu", 1.2512299012E+08);
  nuclear_exponents.emplace("am", 1.2543221826E+08);
  nuclear_exponents.emplace("cm", 1.2420711085E+08);
  nuclear_exponents.emplace("bk", 1.2420711085E+08);
  nuclear_exponents.emplace("cf", 1.2301273547E+08);
  nuclear_exponents.emplace("es", 1.2271879740E+08);
  nuclear_exponents.emplace("fm", 1.2127611477E+08);
  nuclear_exponents.emplace("md", 1.2099285491E+08);
  nuclear_exponents.emplace("no", 1.2071131346E+08);
  nuclear_exponents.emplace("lr", 1.1987683191E+08);
  nuclear_exponents.emplace("db", 1.2015331850E+08);
  nuclear_exponents.emplace("jl", 1.1987683191E+08);
  nuclear_exponents.emplace("rf", 1.1960199758E+08);
  nuclear_exponents.emplace("bh", 1.1987683191E+08);
  nuclear_exponents.emplace("hn", 1.1905722195E+08);
  nuclear_exponents.emplace("mt", 1.1878724932E+08);

  //Atomic masses taken from CRC Handbook of Chemistry and Physics
  //average isotope mass
  averaged_masses.emplace("q", 0.0);
  averaged_masses.emplace("h", 1.00794);
  averaged_masses.emplace("he", 4.00260);
  averaged_masses.emplace("li", 6.941);
  averaged_masses.emplace("be", 9.01218);
  averaged_masses.emplace("b", 10.81);
  averaged_masses.emplace("c", 12.011);
  averaged_masses.emplace("n", 14.0067);
  averaged_masses.emplace("o", 15.9994);
  averaged_masses.emplace("f", 18.998403);
  averaged_masses.emplace("ne", 20.179);
  averaged_masses.emplace("na", 22.989877);
  averaged_masses.emplace("mg", 24.305);
  averaged_masses.emplace("al", 26.98154);
  averaged_masses.emplace("si", 28.0955);
  averaged_masses.emplace("p", 30.97376);
  averaged_masses.emplace("s", 32.06);
  averaged_masses.emplace("cl", 35.453);
  averaged_masses.emplace("ar", 39.948);
  averaged_masses.emplace("k", 39.0983);
  averaged_masses.emplace("ca", 40.08);
  averaged_masses.emplace("sc", 44.9559);
  averaged_masses.emplace("ti", 47.88);
  averaged_masses.emplace("v", 50.9415);
  averaged_masses.emplace("cr", 51.996);
  averaged_masses.emplace("mn", 54.9380);
  averaged_masses.emplace("fe", 55.847);
  averaged_masses.emplace("co", 58.9332);
  averaged_masses.emplace("ni", 58.69);
  averaged_masses.emplace("cu", 63.546);
  averaged_masses.emplace("zn", 65.39);
  averaged_masses.emplace("ga", 69.72);
  averaged_masses.emplace("ge", 72.59);
  averaged_masses.emplace("as", 74.9216);
  averaged_masses.emplace("se", 78.96);
  averaged_masses.emplace("br", 79.904);
  averaged_masses.emplace("kr", 83.80);
  averaged_masses.emplace("rb", 85.4678);
  averaged_masses.emplace("sr", 87.62);
  averaged_masses.emplace("y", 88.9059);
  averaged_masses.emplace("zr", 91.22);
  averaged_masses.emplace("nb", 92.9064);
  averaged_masses.emplace("mo", 95.94);
  averaged_masses.emplace("tc", 98.0);
  averaged_masses.emplace("ru", 101.07);
  averaged_masses.emplace("rh", 102.9055);
  averaged_masses.emplace("pd", 106.42);
  averaged_masses.emplace("ag", 107.8682);
  averaged_masses.emplace("cd", 112.41);
  averaged_masses.emplace("in", 114.82);
  averaged_masses.emplace("sn", 118.69);
  averaged_masses.emplace("sb", 121.75);
  averaged_masses.emplace("te", 127.60);
  averaged_masses.emplace("i", 126.9054);
  averaged_masses.emplace("xe", 131.29);
  averaged_masses.emplace("cs", 132.9054);
  averaged_masses.emplace("ba", 137.33);
  averaged_masses.emplace("la", 138.9055);
  averaged_masses.emplace("ce", 140.12);
  averaged_masses.emplace("pr", 140.9077);
  averaged_masses.emplace("nd", 144.24);
  averaged_masses.emplace("pm", 145.0);
  averaged_masses.emplace("sm", 150.36);
  averaged_masses.emplace("eu", 151.96);
  averaged_masses.emplace("gd", 157.25);
  averaged_masses.emplace("tb", 158.9254);
  averaged_masses.emplace("dy", 162.5);
  averaged_masses.emplace("ho", 164.9304);
  averaged_masses.emplace("er", 167.26);
  averaged_masses.emplace("tm", 168.9342);
  averaged_masses.emplace("yb", 173.04);
  averaged_masses.emplace("lu", 174.976);
  averaged_masses.emplace("hf", 178.49);
  averaged_masses.emplace("ta", 180.9479);
  averaged_masses.emplace("w", 183.85);
  averaged_masses.emplace("re", 186.207);
  averaged_masses.emplace("os", 190.2);
  averaged_masses.emplace("ir", 192.22);
  averaged_masses.emplace("pt", 195.08);
  averaged_masses.emplace("au", 196.9665);
  averaged_masses.emplace("hg", 200.59);
  averaged_masses.emplace("tl", 204.383);
  averaged_masses.emplace("pb", 207.2);
  averaged_masses.emplace("bi", 208.9804);
  averaged_masses.emplace("po", 208.9825);
  averaged_masses.emplace("at", 210.9875);
  averaged_masses.emplace("rn", 222.0175);
  averaged_masses.emplace("fr", 223.0198);
  averaged_masses.emplace("ra", 226.0254);
  averaged_masses.emplace("ac", 227.0278);
  averaged_masses.emplace("th", 232.0381);
  averaged_masses.emplace("pa", 231.0359);
  averaged_masses.emplace("u", 238.0289);
  averaged_masses.emplace("np", 237.0482);
  averaged_masses.emplace("pu", 244.0);
  averaged_masses.emplace("am", 243.0);
  averaged_masses.emplace("cm", 247.0);

  angmap.emplace("s", 0);
  angmap.emplace("p", 1);
  angmap.emplace("d", 2);
  angmap.emplace("f", 3);
  angmap.emplace("g", 4);
  angmap.emplace("h", 5);
  angmap.emplace("i", 6);
  angmap.emplace("j", 7);
// Since they are not implemented yet
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
  nclosed.emplace("ba", make_tuple(12,24,20,0));
  nclosed.emplace("la", make_tuple(12,24,20,0));
  nclosed.emplace("ce", make_tuple(12,24,20,0));
  nclosed.emplace("pr", make_tuple(12,24,20,0));
  nclosed.emplace("nd", make_tuple(12,24,20,0));
  nclosed.emplace("pm", make_tuple(12,24,20,0));
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

  // P in the hyperfine coupling formula, copied from Takeshi's ORZ package
  // which in turn originates from Orca's output according to Takeshi
  // P = g_e * g_N * beta_e * beta_N :
  //   g_e = 2.000... (g factor of free electron)
  //   g_N = g factor of nucleus
  //   beta_e = electron magneton = alpha(fine structure constant) / 2.0
  //   beta_N = nuclear magneton
  hfccp.emplace("h",   533.5514);
  hfccp.emplace("he",    0.0000);
  hfccp.emplace("li",  207.3726);
  hfccp.emplace("be",    0.0000);
  hfccp.emplace("b",   171.2143);
  hfccp.emplace("c",   134.1900);
  hfccp.emplace("n",    38.5677);
  hfccp.emplace("o",   -72.3588);
  hfccp.emplace("f",   502.2244);
  hfccp.emplace("ne",    0.0000);
  hfccp.emplace("na",  141.2175);
  hfccp.emplace("mg",    0.0000);
  hfccp.emplace("al",  139.1361);
  hfccp.emplace("si",    0.0000);
  hfccp.emplace("p",   216.1834);
  hfccp.emplace("s",     0.0000);
  hfccp.emplace("cl",    0.0000);
  hfccp.emplace("ar",    0.0000);
  hfccp.emplace("k",    24.9301);
  hfccp.emplace("ca",    0.0000);
  hfccp.emplace("sc",    0.0000);
  hfccp.emplace("ti",  -30.1264);
  hfccp.emplace("v",   140.2594);
  hfccp.emplace("cr",  -30.0605);
  hfccp.emplace("mn",  132.0006);
  hfccp.emplace("fe",    0.0000);
  hfccp.emplace("co",    0.0000);
  hfccp.emplace("ni",    0.0000);
  hfccp.emplace("cu",  141.7533);
  hfccp.emplace("zn",   33.4622);
  hfccp.emplace("ga",    0.0000);
  hfccp.emplace("ge",    0.0000);
  hfccp.emplace("as",    0.0000);
  hfccp.emplace("se",    0.0000);
  hfccp.emplace("br",    0.0000);
  hfccp.emplace("kr",    0.0000);
//hfccp.emplace("rb",    0.0000);
//hfccp.emplace("sr",    0.0000);
//hfccp.emplace("y",     0.0000);
//hfccp.emplace("zr",    0.0000);
//hfccp.emplace("nb",    0.0000);
//hfccp.emplace("mo",    0.0000);
//hfccp.emplace("tc",    0.0000);
//hfccp.emplace("ru",    0.0000);
  hfccp.emplace("rh",   -16.8881);
  hfccp.emplace("pd",   -24.4534);
  hfccp.emplace("ag",   -21.7071);
  hfccp.emplace("cd",  -113.7112);
//hfccp.emplace("in",    0.0000);
//hfccp.emplace("sn",    0.0000);
//hfccp.emplace("sb",    0.0000);
//hfccp.emplace("te",    0.0000);
//hfccp.emplace("i",     0.0000);
//hfccp.emplace("xe",    0.0000);
//hfccp.emplace("cs",    0.0000);
//hfccp.emplace("ba",    0.0000);
//hfccp.emplace("la",    0.0000);
//hfccp.emplace("ce",    0.0000);
//hfccp.emplace("pr",    0.0000);
//hfccp.emplace("nd",    0.0000);
//hfccp.emplace("pm",    0.0000);
//hfccp.emplace("sm",    0.0000);
//hfccp.emplace("eu",    0.0000);
//hfccp.emplace("gd",    0.0000);
//hfccp.emplace("tb",    0.0000);
//hfccp.emplace("dy",    0.0000);
//hfccp.emplace("ho",    0.0000);
//hfccp.emplace("er",    0.0000);
//hfccp.emplace("tm",    0.0000);
//hfccp.emplace("yb",    0.0000);
//hfccp.emplace("lu",    0.0000);
//hfccp.emplace("hf",    0.0000);
//hfccp.emplace("ta",    0.0000);
//hfccp.emplace("w",     0.0000);
//hfccp.emplace("re",    0.0000);
//hfccp.emplace("os",    0.0000);
//hfccp.emplace("ir",    0.0000);
//hfccp.emplace("pt",    0.0000);
  hfccp.emplace("au",    9.3580);
//hfccp.emplace("hg",    0.0000);
//hfccp.emplace("tl",    0.0000);
//hfccp.emplace("pb",    0.0000);
//hfccp.emplace("bi",    0.0000);
//hfccp.emplace("po",    0.0000);
//hfccp.emplace("at",    0.0000);
//hfccp.emplace("rn",    0.0000);
//hfccp.emplace("fr",    0.0000);
//hfccp.emplace("ra",    0.0000);
//hfccp.emplace("ac",    0.0000);
//hfccp.emplace("th",    0.0000);
//hfccp.emplace("pa",    0.0000);
//hfccp.emplace("u",     0.0000);
//hfccp.emplace("np",    0.0000);
//hfccp.emplace("pu",    0.0000);
//hfccp.emplace("am",    0.0000);
//hfccp.emplace("cm",    0.0000);
//hfccp.emplace("bk",    0.0000);
//hfccp.emplace("cf",    0.0000);
//hfccp.emplace("es",    0.0000);
//hfccp.emplace("fm",    0.0000);
//hfccp.emplace("md",    0.0000);
//hfccp.emplace("no",    0.0000);
//hfccp.emplace("lr",    0.0000);
//hfccp.emplace("rf",    0.0000);

}

int AtomMap::angular_number(const string input) const {
  auto miter = angmap.find(input);
  stringstream ss; ss << "Unknown angular number in a basis set file. Requested: " << input << endl;
  if (miter == angmap.end()) throw runtime_error(ss.str());
#ifndef COMPILE_J_ORB
  if (input == "j")
    throw runtime_error("j-orbitals requested in a basis set file.  BAGEL must be compiled with the -DCOMPILE_J_ORB flag to use this.");
#endif
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


double AtomMap::nuclear_exponent(const string input) const {
  auto miter = nuclear_exponents.find(input);
  if (miter == nuclear_exponents.end()) throw runtime_error("Unknown atom (Finite nucleus exponent).");
  return miter->second;
}

double AtomMap::averaged_mass(const string input) const {
  auto miter = averaged_masses.find(input);
  if (miter == averaged_masses.end()) throw runtime_error("Unknown atom (Averaged atomic mass).");
  return miter->second;
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


bool AtomMap::hfcc_exists(const string input) const {
  return hfccp.find(input) != hfccp.end();
}


double AtomMap::hfcc_pfac(const string input) const {
  auto miter = hfccp.find(input);
  if (miter == hfccp.end()) throw runtime_error("Unknown atom (HPCC P factor).");
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
