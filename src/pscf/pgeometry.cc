//
// BAGEL - Parallel electron correlation program.
// Filename: pgeometry.cc
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


#include <src/pscf/pgeometry.h>
#include <fstream>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <regex>
#include <cmath>
#include <src/util/string.h>

using namespace std;
using namespace bagel;

PGeometry::PGeometry(const string fil, const int levl) : Geometry(fil) {
  ifstream ifs(fil);
  if(!ifs.is_open()) throw runtime_error("input file could not be opened.");

  smatch what;
  regex reg("Periodic");
  int out = 0;
  bool found = false;
  while(!ifs.eof()) {
    string sline;
    getline(ifs, sline);
    string::const_iterator start = sline.begin();
    string::const_iterator end   = sline.end();
    if (regex_search(start, end, what, reg)) {
      found = true;
      break;
    }
  }

  assert(found);
  regex regL("L");
  regex regS("S");
  regex regK("K");
  regex rega("Period");
  regex regnumber("[0-9]+");
  regex regdouble("[0-9e\\.\\+-]+");
  bool foundL = false;
  bool foundS = false;
  bool foundK = false;
  bool founda = false;
  while(!ifs.eof()) {
    string sline;
    getline(ifs, sline);
    if (sline.empty()) continue;

    string::const_iterator start = sline.begin();
    string::const_iterator end   = sline.end();

    if (regex_search(start, end, what, regL)) {
      start = what[0].second;
      foundL = regex_search(start, end, what, regnumber);
      if (foundL) {
        string tmp_str(what[0].first, what[0].second);
        L_ = lexical_cast<int>(tmp_str);
      }
    } else if (regex_search(start, end, what, regS)) {
      start = what[0].second;
      foundS = regex_search(start, end, what, regnumber);
      if (foundS) {
        string tmp_str(what[0].first, what[0].second);
        S_ = lexical_cast<int>(tmp_str);
      }
    } else if (regex_search(start, end, what, regK)) {
      start = what[0].second;
      foundK = regex_search(start, end, what, regnumber);
      if (foundK) {
        string tmp_str(what[0].first, what[0].second);
        K_ = lexical_cast<int>(tmp_str);
      }
    } else if (regex_search(start, end, what, rega)) {
      start = what[0].second;
      founda = regex_search(start, end, what, regdouble);
      if (founda) {
        string tmp_str(what[0].first, what[0].second);
        A_ = lexical_cast<double>(tmp_str);
      }
    } else {
      break;
    }
  }
  ifs.close();
  assert(founda || (A_ > 0.0) );

  if (!foundK) K_ = 10;
  if (!foundL) L_ = K_;
  if (!foundS) S_ = K_;

  cout << "  " << "* Performs crystalline orbital calculation (1D)" << endl;
  cout << endl;
  cout << "  " << "    Parameters : K = " << K() << "  (L = " << L() << " and S = " << S() << ")" << endl;
  cout << "  " << "    Cell size  : " << setprecision(10) << setw(15) << A() << endl;
  cout << endl;

  pnuclear_repulsion_ = compute_pnuclear_repulsion();
}


double PGeometry::compute_pnuclear_repulsion() const {
  double out = 0.0;
  typedef std::shared_ptr<const Atom> RefAtom;
  for (auto iter = atoms_.begin(); iter != atoms_.end(); ++iter) {
    const array<double,3> tmp = (*iter)->position();
    const double c = (*iter)->atom_charge();
    for (auto titer = atoms_.begin(); titer != atoms_.end(); ++titer) {
      const RefAtom target = *titer;
      const double charge = c * target->atom_charge();
      for (int el = -L(); el <= L(); ++el) {
        if (el == 0 && titer == iter) continue;
        const double disp = el * A();
        const double x = target->position(0) - tmp[0];
        const double y = target->position(1) - tmp[1];
        const double z = target->position(2) - tmp[2] + disp;
        const double dist = sqrt(x * x + y * y + z * z);
        if (!(*iter)->dummy() || !(*titer)->dummy())
          out += charge / dist * 0.5;
      }
    }
  }
  return out;
}


