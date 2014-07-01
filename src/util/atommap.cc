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
  atommap.insert({"q", 0});
  atommap.insert({"h", 1});
  atommap.insert({"he", 2});
  atommap.insert({"li", 3});
  atommap.insert({"be", 4});
  atommap.insert({"b", 5});
  atommap.insert({"c", 6});
  atommap.insert({"n", 7});
  atommap.insert({"o", 8});
  atommap.insert({"f", 9});
  atommap.insert({"ne", 10});
  atommap.insert({"na", 11});
  atommap.insert({"mg", 12});
  atommap.insert({"al", 13});
  atommap.insert({"si", 14});
  atommap.insert({"p", 15});
  atommap.insert({"s", 16});
  atommap.insert({"cl", 17});
  atommap.insert({"ar", 18});
  atommap.insert({"k", 19});
  atommap.insert({"ca", 20});
  atommap.insert({"sc", 21});
  atommap.insert({"ti", 22});
  atommap.insert({"v", 23});
  atommap.insert({"cr", 24});
  atommap.insert({"mn", 25});
  atommap.insert({"fe", 26});
  atommap.insert({"co", 27});
  atommap.insert({"ni", 28});
  atommap.insert({"cu", 29});
  atommap.insert({"zn", 30});
  atommap.insert({"ga", 31});
  atommap.insert({"ge", 32});
  atommap.insert({"as", 33});
  atommap.insert({"se", 34});
  atommap.insert({"br", 35});
  atommap.insert({"kr", 36});
  atommap.insert({"rb", 37});
  atommap.insert({"sr", 38});
  atommap.insert({"y", 39});
  atommap.insert({"zr", 40});
  atommap.insert({"nb", 41});
  atommap.insert({"mo", 42});
  atommap.insert({"tc", 43});
  atommap.insert({"ru", 44});
  atommap.insert({"rh", 45});
  atommap.insert({"pd", 46});
  atommap.insert({"ag", 47});
  atommap.insert({"cd", 48});
  atommap.insert({"in", 49});
  atommap.insert({"sn", 50});
  atommap.insert({"sb", 51});
  atommap.insert({"te", 52});
  atommap.insert({"i", 53});
  atommap.insert({"xe", 54});
  atommap.insert({"cs", 55});
  atommap.insert({"ba", 56});
  atommap.insert({"la", 57});
  atommap.insert({"ce", 58});
  atommap.insert({"pr", 59});
  atommap.insert({"nd", 60});
  atommap.insert({"pm", 61});
  atommap.insert({"sm", 62});
  atommap.insert({"eu", 63});
  atommap.insert({"gd", 64});
  atommap.insert({"tb", 65});
  atommap.insert({"dy", 66});
  atommap.insert({"ho", 67});
  atommap.insert({"er", 68});
  atommap.insert({"tm", 69});
  atommap.insert({"yb", 70});
  atommap.insert({"lu", 71});
  atommap.insert({"hf", 72});
  atommap.insert({"ta", 73});
  atommap.insert({"w", 74});
  atommap.insert({"re", 75});
  atommap.insert({"os", 76});
  atommap.insert({"ir", 77});
  atommap.insert({"pt", 78});
  atommap.insert({"au", 79});
  atommap.insert({"hg", 80});
  atommap.insert({"tl", 81});
  atommap.insert({"pb", 82});
  atommap.insert({"bi", 83});
  atommap.insert({"po", 84});
  atommap.insert({"at", 85});
  atommap.insert({"rn", 86});
  atommap.insert({"fr", 87});
  atommap.insert({"ra", 88});
  atommap.insert({"ac", 89});
  atommap.insert({"th", 90});
  atommap.insert({"pa", 91});
  atommap.insert({"u", 92});
  atommap.insert({"np", 93});
  atommap.insert({"pu", 94});
  atommap.insert({"am", 95});
  atommap.insert({"cm", 96});
  atommap.insert({"bk", 97});
  atommap.insert({"cf", 98});
  atommap.insert({"es", 99});
  atommap.insert({"fm", 100});
  atommap.insert({"md", 101});
  atommap.insert({"no", 102});
  atommap.insert({"lr", 103});
  atommap.insert({"rf", 104});
  atommap.insert({"db", 105});
  atommap.insert({"sg", 106});
  atommap.insert({"bh", 107});
  atommap.insert({"hs", 108});
  atommap.insert({"mt", 109});
  atommap.insert({"ds", 110});
  atommap.insert({"rg", 111});
  atommap.insert({"cn", 112});
  atommap.insert({"uut", 113});
  atommap.insert({"fl", 114});
  atommap.insert({"uup", 115});
  atommap.insert({"lv", 116});
  atommap.insert({"uus", 117});
  atommap.insert({"uuo", 118});

  cov_radii.insert({"h",  0.31});
  cov_radii.insert({"he", 0.28});
  cov_radii.insert({"li", 1.28});
  cov_radii.insert({"be", 0.96});
  cov_radii.insert({"b",  0.84});
  cov_radii.insert({"c",  0.76});
  cov_radii.insert({"n",  0.71});
  cov_radii.insert({"o",  0.66});
  cov_radii.insert({"f",  0.57});
  cov_radii.insert({"ne", 0.58});
  cov_radii.insert({"na", 1.66});
  cov_radii.insert({"mg", 1.41});
  cov_radii.insert({"al", 1.21});
  cov_radii.insert({"si", 1.11});
  cov_radii.insert({"p",  1.07});
  cov_radii.insert({"s",  1.05});
  cov_radii.insert({"cl", 1.02});
  cov_radii.insert({"ar", 1.06});
  cov_radii.insert({"k",  2.03});
  cov_radii.insert({"ca", 1.76});
  cov_radii.insert({"sc", 1.70});
  cov_radii.insert({"ti", 1.60});
  cov_radii.insert({"v",  1.53});
  cov_radii.insert({"cr", 1.39});
  cov_radii.insert({"mn", 1.39});
  cov_radii.insert({"fe", 1.32});
  cov_radii.insert({"co", 1.26});
  cov_radii.insert({"ni", 1.24});
  cov_radii.insert({"cu", 1.32});
  cov_radii.insert({"zn", 1.22});
  cov_radii.insert({"ga", 1.22});
  cov_radii.insert({"ge", 1.20});
  cov_radii.insert({"as", 1.19});
  cov_radii.insert({"se", 1.20});
  cov_radii.insert({"br", 1.20});
  cov_radii.insert({"kr", 1.16});
  cov_radii.insert({"rb", 2.20});
  cov_radii.insert({"sr", 1.95});
  cov_radii.insert({"y",  1.90});
  cov_radii.insert({"zr", 1.75});
  cov_radii.insert({"nb", 1.64});
  cov_radii.insert({"mo", 1.54});
  cov_radii.insert({"tc", 1.47});
  cov_radii.insert({"ru", 1.46});
  cov_radii.insert({"rh", 1.42});
  cov_radii.insert({"pd", 1.39});
  cov_radii.insert({"ag", 1.45});
  cov_radii.insert({"cd", 1.44});
  cov_radii.insert({"in", 1.42});
  cov_radii.insert({"sn", 1.39});
  cov_radii.insert({"sb", 1.39});
  cov_radii.insert({"te", 1.38});
  cov_radii.insert({"i",  1.39});
  cov_radii.insert({"xe", 1.40});
  cov_radii.insert({"cs", 2.44});
  cov_radii.insert({"ba", 2.15});
  cov_radii.insert({"la", 2.07});
  cov_radii.insert({"ce", 2.04});
  cov_radii.insert({"pr", 2.03});
  cov_radii.insert({"nd", 2.01});
  cov_radii.insert({"pm", 1.99});
  cov_radii.insert({"sm", 1.98});
  cov_radii.insert({"eu", 1.98});
  cov_radii.insert({"gd", 1.96});
  cov_radii.insert({"tb", 1.94});
  cov_radii.insert({"dy", 1.92});
  cov_radii.insert({"ho", 1.92});
  cov_radii.insert({"er", 1.89});
  cov_radii.insert({"tm", 1.90});
  cov_radii.insert({"yb", 1.87});
  cov_radii.insert({"lu", 1.87});
  cov_radii.insert({"hf", 1.75});
  cov_radii.insert({"ta", 1.70});
  cov_radii.insert({"w",  1.62});
  cov_radii.insert({"re", 1.51});
  cov_radii.insert({"os", 1.44});
  cov_radii.insert({"ir", 1.41});
  cov_radii.insert({"pt", 1.36});
  cov_radii.insert({"au", 1.36});
  cov_radii.insert({"hg", 1.32});
  cov_radii.insert({"tl", 1.45});
  cov_radii.insert({"pb", 1.46});
  cov_radii.insert({"bi", 1.48});
  cov_radii.insert({"po", 1.40});
  cov_radii.insert({"at", 1.50});
  cov_radii.insert({"rn", 1.50});
  cov_radii.insert({"fr", 2.60});
  cov_radii.insert({"ra", 2.21});
  cov_radii.insert({"ac", 2.15});
  cov_radii.insert({"th", 2.06});
  cov_radii.insert({"pa", 2.00});
  cov_radii.insert({"u",  1.96});
  cov_radii.insert({"np", 1.90});
  cov_radii.insert({"pu", 1.87});
  cov_radii.insert({"am", 1.80});
  cov_radii.insert({"cm", 1.69});

  // atom sizes (Bragg-Slater radii)
  bsradii.insert({"h", 0.25});
  bsradii.insert({"he", 0.25});
  bsradii.insert({"li", 1.45});
  bsradii.insert({"be", 1.05});
  bsradii.insert({"b", 0.85});
  bsradii.insert({"c", 0.7});
  bsradii.insert({"n", 0.65});
  bsradii.insert({"o", 0.6});
  bsradii.insert({"f", 0.5});
  bsradii.insert({"ne", 0.45});
  bsradii.insert({"na", 1.8});
  bsradii.insert({"mg", 1.5});
  bsradii.insert({"al", 1.25});
  bsradii.insert({"si", 1.1});
  bsradii.insert({"p", 1.0});
  bsradii.insert({"s", 1.0});
  bsradii.insert({"cl", 1.0});
  bsradii.insert({"ar", 1.0});
  bsradii.insert({"k", 2.2});
  bsradii.insert({"ca", 1.8});
  bsradii.insert({"sc", 1.6});
  bsradii.insert({"ti", 1.4});
  bsradii.insert({"v", 1.35});
  bsradii.insert({"cr", 1.4});
  bsradii.insert({"mn", 1.4});
  bsradii.insert({"fe", 1.4});
  bsradii.insert({"co", 1.35});
  bsradii.insert({"ni", 1.35});
  bsradii.insert({"cu", 1.35});
  bsradii.insert({"zn", 1.35});
  bsradii.insert({"ga", 1.3});
  bsradii.insert({"ge", 1.25});
  bsradii.insert({"as", 1.15});
  bsradii.insert({"se", 1.15});
  bsradii.insert({"br", 1.15});
  bsradii.insert({"kr", 1.15});
  bsradii.insert({"rb", 2.35});
  bsradii.insert({"sr", 2.0});
  bsradii.insert({"y", 1.8});
  bsradii.insert({"zr", 1.55});
  bsradii.insert({"nb", 1.45});
  bsradii.insert({"mo", 1.45});
  bsradii.insert({"tc", 1.35});
  bsradii.insert({"ru", 1.3});
  bsradii.insert({"rh", 1.35});
  bsradii.insert({"pd", 1.4});
  bsradii.insert({"ag", 1.6});
  bsradii.insert({"cd", 1.55});
  bsradii.insert({"in", 1.55});
  bsradii.insert({"sn", 1.45});
  bsradii.insert({"sb", 1.45});
  bsradii.insert({"te", 1.4});
  bsradii.insert({"i", 1.4});
  bsradii.insert({"xe", 1.4});
  bsradii.insert({"cs", 2.6});
  bsradii.insert({"ba", 2.15});
  bsradii.insert({"la", 1.95});
  bsradii.insert({"ce", 1.85});
  bsradii.insert({"pr", 1.85});
  bsradii.insert({"nd", 1.85});
  bsradii.insert({"pm", 1.85});
  bsradii.insert({"sm", 1.85});
  bsradii.insert({"eu", 1.85});
  bsradii.insert({"gd", 1.8});
  bsradii.insert({"tb", 1.75});
  bsradii.insert({"dy", 1.75});
  bsradii.insert({"ho", 1.75});
  bsradii.insert({"er", 1.75});
  bsradii.insert({"tm", 1.75});
  bsradii.insert({"yb", 1.75});
  bsradii.insert({"lu", 1.75});
  bsradii.insert({"hf", 1.55});
  bsradii.insert({"ta", 1.45});
  bsradii.insert({"w", 1.35});
  bsradii.insert({"re", 1.35});
  bsradii.insert({"os", 1.3});
  bsradii.insert({"ir", 1.35});
  bsradii.insert({"pt", 1.35});
  bsradii.insert({"au", 1.35});
  bsradii.insert({"hg", 1.5});
  bsradii.insert({"tl", 1.9});
  bsradii.insert({"pb", 1.8});
  bsradii.insert({"bi", 1.6});
  bsradii.insert({"po", 1.9});
  bsradii.insert({"at", 1.9});
  bsradii.insert({"rn", 1.9});
  bsradii.insert({"fr", 2.85});
  bsradii.insert({"ra", 2.15});
  bsradii.insert({"ac", 1.95});
  bsradii.insert({"th", 1.8});
  bsradii.insert({"pa", 1.8});
  bsradii.insert({"u", 1.75});
  bsradii.insert({"np", 1.75});
  bsradii.insert({"pu", 1.75});
  bsradii.insert({"am", 1.75});
  bsradii.insert({"cm", 1.75});
  bsradii.insert({"bk", 1.75});
  bsradii.insert({"cf", 1.75});
  bsradii.insert({"es", 1.75});
  bsradii.insert({"fm", 1.75});
  bsradii.insert({"md", 1.75});
  bsradii.insert({"no", 1.75});
  bsradii.insert({"lr", 1.75});
#if 0
  bsradii.insert({"rf", 104});
  bsradii.insert({"db", 105});
  bsradii.insert({"sg", 106});
  bsradii.insert({"bh", 107});
  bsradii.insert({"hs", 108});
  bsradii.insert({"mt", 109});
  bsradii.insert({"ds", 110});
  bsradii.insert({"rg", 111});
  bsradii.insert({"cn", 112});
  bsradii.insert({"uut", 113});
  bsradii.insert({"fl", 114});
  bsradii.insert({"uup", 115});
  bsradii.insert({"lv", 116});
  bsradii.insert({"uus", 117});
  bsradii.insert({"uuo", 118});
#endif

  angmap.insert({"s", 0});
  angmap.insert({"p", 1});
  angmap.insert({"d", 2});
  angmap.insert({"f", 3});
  angmap.insert({"g", 4});
  angmap.insert({"h", 5});
  angmap.insert({"i", 6});
// Since they are not implemented yet
//angmap.insert({"j", 7});
//angmap.insert({"k", 8});
//angmap.insert({"l", 9});

  nclosed.insert({"h",  make_tuple(0,0,0,0)});
  nclosed.insert({"he", make_tuple(2,0,0,0)});
  nclosed.insert({"li", make_tuple(2,0,0,0)});
  nclosed.insert({"be", make_tuple(4,0,0,0)});
  nclosed.insert({"b",  make_tuple(4,0,0,0)});
  nclosed.insert({"c",  make_tuple(4,0,0,0)});
  nclosed.insert({"n",  make_tuple(4,0,0,0)});
  nclosed.insert({"o",  make_tuple(4,0,0,0)});
  nclosed.insert({"f",  make_tuple(4,0,0,0)});
  nclosed.insert({"ne", make_tuple(4,6,0,0)});
  nclosed.insert({"na", make_tuple(4,6,0,0)});
  nclosed.insert({"mg", make_tuple(6,6,0,0)});
  nclosed.insert({"al", make_tuple(6,6,0,0)});
  nclosed.insert({"si", make_tuple(6,6,0,0)});
  nclosed.insert({"p",  make_tuple(6,6,0,0)});
  nclosed.insert({"s",  make_tuple(6,6,0,0)});
  nclosed.insert({"cl", make_tuple(6,6,0,0)});
  nclosed.insert({"ar", make_tuple(6,12,0,0)});
  nclosed.insert({"k",  make_tuple(6,12,0,0)});
  nclosed.insert({"ca", make_tuple(8,12,0,0)});
  nclosed.insert({"sc", make_tuple(8,12,0,0)});
  nclosed.insert({"ti", make_tuple(8,12,0,0)});
  nclosed.insert({"v",  make_tuple(8,12,0,0)});
  nclosed.insert({"cr", make_tuple(8,12,0,0)});
  nclosed.insert({"mn", make_tuple(8,12,0,0)});
  nclosed.insert({"fe", make_tuple(8,12,0,0)});
  nclosed.insert({"co", make_tuple(8,12,0,0)});
  nclosed.insert({"ni", make_tuple(8,12,0,0)});
  nclosed.insert({"cu", make_tuple(8,12,0,0)});
  nclosed.insert({"zn", make_tuple(8,12,10,0)});
  nclosed.insert({"ga", make_tuple(8,12,10,0)});
  nclosed.insert({"ge", make_tuple(8,12,10,0)});
  nclosed.insert({"as", make_tuple(8,12,10,0)});
  nclosed.insert({"se", make_tuple(8,12,10,0)});
  nclosed.insert({"br", make_tuple(8,12,10,0)});
  nclosed.insert({"kr", make_tuple(8,18,10,0)});
  nclosed.insert({"rb", make_tuple(8,18,10,0)});
  nclosed.insert({"sr", make_tuple(10,18,10,0)});
  nclosed.insert({"y",  make_tuple(10,18,10,0)});
  nclosed.insert({"zr", make_tuple(10,18,10,0)});
  nclosed.insert({"nb", make_tuple(10,18,10,0)});
  nclosed.insert({"mo", make_tuple(10,18,10,0)});
  nclosed.insert({"tc", make_tuple(10,18,10,0)});
  nclosed.insert({"ru", make_tuple(10,18,10,0)});
  nclosed.insert({"rh", make_tuple(10,18,10,0)});
  nclosed.insert({"pd", make_tuple(10,18,10,0)});
  nclosed.insert({"ag", make_tuple(10,18,10,0)});
  nclosed.insert({"cd", make_tuple(10,18,20,0)});
  nclosed.insert({"in", make_tuple(10,18,20,0)});
  nclosed.insert({"sn", make_tuple(10,18,20,0)});
  nclosed.insert({"sb", make_tuple(10,18,20,0)});
  nclosed.insert({"te", make_tuple(10,18,20,0)});
  nclosed.insert({"i",  make_tuple(10,18,20,0)});
  nclosed.insert({"xe", make_tuple(10,24,20,0)});
  nclosed.insert({"cs", make_tuple(10,24,20,0)});
  nclosed.insert({"ba", make_tuple(10,24,20,0)});
  nclosed.insert({"la", make_tuple(10,24,20,0)});
  nclosed.insert({"ce", make_tuple(10,24,20,0)});
  nclosed.insert({"pr", make_tuple(10,24,20,0)});
  nclosed.insert({"nd", make_tuple(10,24,20,0)});
  nclosed.insert({"pm", make_tuple(10,24,20,0)});
  nclosed.insert({"sm", make_tuple(12,24,20,0)});
  nclosed.insert({"eu", make_tuple(12,24,20,0)});
  nclosed.insert({"gd", make_tuple(12,24,20,0)});
  nclosed.insert({"tb", make_tuple(12,24,20,0)});
  nclosed.insert({"dy", make_tuple(12,24,20,0)});
  nclosed.insert({"ho", make_tuple(12,24,20,0)});
  nclosed.insert({"er", make_tuple(12,24,20,0)});
  nclosed.insert({"tm", make_tuple(12,24,20,0)});
  nclosed.insert({"yb", make_tuple(12,24,20,14)});
  nclosed.insert({"lu", make_tuple(12,24,20,14)});
  nclosed.insert({"hf", make_tuple(12,24,20,14)});
  nclosed.insert({"ta", make_tuple(12,24,20,14)});
  nclosed.insert({"w",  make_tuple(12,24,20,14)});
  nclosed.insert({"re", make_tuple(12,24,20,14)});
  nclosed.insert({"os", make_tuple(12,24,20,14)});
  nclosed.insert({"ir", make_tuple(12,24,20,14)});
  nclosed.insert({"pt", make_tuple(10,24,20,14)});
  nclosed.insert({"au", make_tuple(10,24,30,14)});
  nclosed.insert({"hg", make_tuple(12,24,30,14)});
  nclosed.insert({"tl", make_tuple(12,24,30,14)});
  nclosed.insert({"pb", make_tuple(12,24,30,14)});
  nclosed.insert({"bi", make_tuple(12,24,30,14)});
  nclosed.insert({"po", make_tuple(12,24,30,14)});
  nclosed.insert({"at", make_tuple(12,24,30,14)});
  nclosed.insert({"rn", make_tuple(12,30,30,14)});
  nclosed.insert({"fr", make_tuple(12,30,30,14)});
  nclosed.insert({"ra", make_tuple(14,30,30,14)});
  nclosed.insert({"ac", make_tuple(14,30,30,14)});
  nclosed.insert({"th", make_tuple(14,30,30,14)});
  nclosed.insert({"pa", make_tuple(14,30,30,14)});
  nclosed.insert({"u",  make_tuple(14,30,30,14)});
  nclosed.insert({"np", make_tuple(14,30,30,14)});
  nclosed.insert({"pu", make_tuple(14,30,30,14)});
  nclosed.insert({"am", make_tuple(14,30,30,14)});
  nclosed.insert({"cm", make_tuple(14,30,30,14)});
  nclosed.insert({"bk", make_tuple(14,30,30,14)});
  nclosed.insert({"cf", make_tuple(14,30,30,14)});
  nclosed.insert({"es", make_tuple(14,30,30,14)});
  nclosed.insert({"fm", make_tuple(14,30,30,14)});
  nclosed.insert({"md", make_tuple(14,30,30,14)});
  nclosed.insert({"no", make_tuple(14,30,30,14)});
  nclosed.insert({"lr", make_tuple(14,30,30,28)});
  nclosed.insert({"rf", make_tuple(14,30,30,28)});
  nclosed.insert({"db", make_tuple(14,30,30,28)});
  nclosed.insert({"sg", make_tuple(14,30,30,28)});
  nclosed.insert({"bh", make_tuple(14,30,30,28)});
  nclosed.insert({"hs", make_tuple(14,30,30,28)});
  nclosed.insert({"mt", make_tuple(14,30,30,28)});
  nclosed.insert({"ds", make_tuple(14,30,30,28)});
  nclosed.insert({"rg", make_tuple(14,30,30,28)});
  nclosed.insert({"cn", make_tuple(14,30,40,28)});
  nclosed.insert({"uut", make_tuple(14,30,40,28)});
  nclosed.insert({"fl",  make_tuple(14,30,40,28)});
  nclosed.insert({"uup", make_tuple(14,30,40,28)});
  nclosed.insert({"lv",  make_tuple(14,30,40,28)});
  nclosed.insert({"uus", make_tuple(14,30,40,28)});
  nclosed.insert({"uuo", make_tuple(14,36,40,28)});

  nopen.insert({"h",  make_tuple(1,0,0,0)});
  nopen.insert({"he", make_tuple(0,0,0,0)});
  nopen.insert({"li", make_tuple(1,0,0,0)});
  nopen.insert({"be", make_tuple(0,0,0,0)});
  nopen.insert({"b",  make_tuple(0,1,0,0)});
  nopen.insert({"c",  make_tuple(0,2,0,0)});
  nopen.insert({"n",  make_tuple(0,3,0,0)});
  nopen.insert({"o",  make_tuple(0,4,0,0)});
  nopen.insert({"f",  make_tuple(0,5,0,0)});
  nopen.insert({"ne", make_tuple(0,0,0,0)});
  nopen.insert({"na", make_tuple(1,0,0,0)});
  nopen.insert({"mg", make_tuple(0,0,0,0)});
  nopen.insert({"al", make_tuple(0,1,0,0)});
  nopen.insert({"si", make_tuple(0,2,0,0)});
  nopen.insert({"p",  make_tuple(0,3,0,0)});
  nopen.insert({"s",  make_tuple(0,4,0,0)});
  nopen.insert({"cl", make_tuple(0,5,0,0)});
  nopen.insert({"ar", make_tuple(0,0,0,0)});
  nopen.insert({"k",  make_tuple(1,0,0,0)});
  nopen.insert({"ca", make_tuple(0,0,0,0)});
  nopen.insert({"sc", make_tuple(0,0,1,0)});
  nopen.insert({"ti", make_tuple(0,0,2,0)});
  nopen.insert({"v",  make_tuple(0,0,3,0)});
  nopen.insert({"cr", make_tuple(0,0,4,0)});
  nopen.insert({"mn", make_tuple(0,0,5,0)});
  nopen.insert({"fe", make_tuple(0,0,6,0)});
  nopen.insert({"co", make_tuple(0,0,7,0)});
  nopen.insert({"ni", make_tuple(0,0,8,0)});
  nopen.insert({"cu", make_tuple(0,0,9,0)});
  nopen.insert({"zn", make_tuple(0,0,0,0)});
  nopen.insert({"ga", make_tuple(0,1,0,0)});
  nopen.insert({"ge", make_tuple(0,2,0,0)});
  nopen.insert({"as", make_tuple(0,3,0,0)});
  nopen.insert({"se", make_tuple(0,4,0,0)});
  nopen.insert({"br", make_tuple(0,5,0,0)});
  nopen.insert({"kr", make_tuple(0,0,0,0)});
  nopen.insert({"rb", make_tuple(1,0,0,0)});
  nopen.insert({"sr", make_tuple(0,0,0,0)});
  nopen.insert({"y",  make_tuple(0,0,1,0)});
  nopen.insert({"zr", make_tuple(0,0,2,0)});
  nopen.insert({"nb", make_tuple(0,0,3,0)});
  nopen.insert({"mo", make_tuple(0,0,4,0)});
  nopen.insert({"tc", make_tuple(0,0,5,0)});
  nopen.insert({"ru", make_tuple(0,0,6,0)});
  nopen.insert({"rh", make_tuple(0,0,7,0)});
  nopen.insert({"pd", make_tuple(0,0,8,0)});
  nopen.insert({"ag", make_tuple(0,0,9,0)});
  nopen.insert({"cd", make_tuple(0,0,0,0)});
  nopen.insert({"in", make_tuple(0,1,0,0)});
  nopen.insert({"sn", make_tuple(0,2,0,0)});
  nopen.insert({"sb", make_tuple(0,3,0,0)});
  nopen.insert({"te", make_tuple(0,4,0,0)});
  nopen.insert({"i",  make_tuple(0,5,0,0)});
  nopen.insert({"xe", make_tuple(0,0,0,0)});
  nopen.insert({"cs", make_tuple(1,0,0,0)});
  nopen.insert({"ba", make_tuple(0,0,0,0)});
  nopen.insert({"la", make_tuple(0,0,1,0)});
  nopen.insert({"ce", make_tuple(0,0,1,1)});
  nopen.insert({"pr", make_tuple(0,0,0,3)});
  nopen.insert({"nd", make_tuple(0,0,0,4)});
  nopen.insert({"pm", make_tuple(0,0,0,5)});
  nopen.insert({"sm", make_tuple(0,0,0,6)});
  nopen.insert({"eu", make_tuple(0,0,0,7)});
  nopen.insert({"gd", make_tuple(0,0,1,7)});
  nopen.insert({"tb", make_tuple(0,0,0,9)});
  nopen.insert({"dy", make_tuple(0,0,0,10)});
  nopen.insert({"ho", make_tuple(0,0,0,11)});
  nopen.insert({"er", make_tuple(0,0,0,12)});
  nopen.insert({"tm", make_tuple(0,0,0,13)});
  nopen.insert({"yb", make_tuple(0,0,0,0)});
  nopen.insert({"lu", make_tuple(0,0,1,0)});
  nopen.insert({"hf", make_tuple(0,0,2,0)});
  nopen.insert({"ta", make_tuple(0,0,3,0)});
  nopen.insert({"w",  make_tuple(0,0,4,0)});
  nopen.insert({"re", make_tuple(0,0,5,0)});
  nopen.insert({"os", make_tuple(0,0,6,0)});
  nopen.insert({"ir", make_tuple(0,0,7,0)});
  nopen.insert({"pt", make_tuple(1,0,9,0)});
  nopen.insert({"au", make_tuple(1,0,0,0)});
  nopen.insert({"hg", make_tuple(0,0,0,0)});
  nopen.insert({"tl", make_tuple(0,1,0,0)});
  nopen.insert({"pb", make_tuple(0,2,0,0)});
  nopen.insert({"bi", make_tuple(0,3,0,0)});
  nopen.insert({"po", make_tuple(0,4,0,0)});
  nopen.insert({"at", make_tuple(0,5,0,0)});
  nopen.insert({"rn", make_tuple(0,0,0,0)});
  nopen.insert({"fr", make_tuple(1,0,0,0)});
  nopen.insert({"ra", make_tuple(0,0,0,0)});
  nopen.insert({"ac", make_tuple(0,0,1,0)});
  nopen.insert({"th", make_tuple(0,0,2,0)});
  nopen.insert({"pa", make_tuple(0,0,1,2)});
  nopen.insert({"u",  make_tuple(0,0,1,3)});
  nopen.insert({"np", make_tuple(0,0,1,4)});
  nopen.insert({"pu", make_tuple(0,0,0,6)});
  nopen.insert({"am", make_tuple(0,0,0,7)});
  nopen.insert({"cm", make_tuple(0,0,1,7)});
  nopen.insert({"bk", make_tuple(0,0,0,9)});
  nopen.insert({"cf", make_tuple(0,0,0,10)});
  nopen.insert({"es", make_tuple(0,0,0,11)});
  nopen.insert({"fm", make_tuple(0,0,0,12)});
  nopen.insert({"md", make_tuple(0,0,0,13)});
  nopen.insert({"no", make_tuple(0,0,0,0)});
  nopen.insert({"lr", make_tuple(0,1,0,0)});
  nopen.insert({"rf", make_tuple(0,0,2,0)});
  nopen.insert({"db", make_tuple(0,0,3,0)});
  nopen.insert({"sg", make_tuple(0,0,4,0)});
  nopen.insert({"bh", make_tuple(0,0,5,0)});
  nopen.insert({"hs", make_tuple(0,0,6,0)});
  nopen.insert({"mt", make_tuple(0,0,7,0)});
  nopen.insert({"ds", make_tuple(0,0,8,0)});
  nopen.insert({"rg", make_tuple(0,0,9,0)});
  nopen.insert({"cn", make_tuple(0,0,0,0)});
  nopen.insert({"uut", make_tuple(0,1,0,0)});
  nopen.insert({"fl",  make_tuple(0,2,0,0)});
  nopen.insert({"uup", make_tuple(0,3,0,0)});
  nopen.insert({"lv",  make_tuple(0,4,0,0)});
  nopen.insert({"uus", make_tuple(0,5,0,0)});
  nopen.insert({"uuo", make_tuple(0,0,0,0)});

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
