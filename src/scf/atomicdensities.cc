//
// BAGEL - Parallel electron correlation program.
// Filename: atomicdensities.cc
// Copyright (C) 2013 Toru Shiozaki
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

#include <src/scf/atomicdensities.h>
#include <src/util/atommap.h>
#include <src/scf/rohf.h>

using namespace std;
using namespace bagel;

const static AtomMap atommap_;

AtomicDensities::AtomicDensities(std::shared_ptr<const Geometry> g) : Matrix(g->nbasis(), g->nbasis()), geom_(g) {
  // first make a list of unique atoms
  string basis = geom_->basisfile();
  map<string, shared_ptr<const Matrix>> atoms;

  for (auto& i : geom_->atoms())
    if (atoms.find(i->name()) == atoms.end()) {
      // dummy buffer to suppress the output
      stringstream ss;
      std::streambuf* cout_orig = cout.rdbuf();
      cout.rdbuf(ss.rdbuf());

      shared_ptr<const Atom> atom(new Atom(i->spherical(), i->name(), {{0.0,0.0,0.0}}, basis));

      multimap<string,string> geomop;
      geomop.insert(make_pair("basis", basis));
      geomop.insert(make_pair("df_basis", basis));
      shared_ptr<const Geometry> ga(new Geometry({atom}, geomop));

      multimap<string, string> options;
      options.insert(make_pair("thresh", "1.0e-6"));
      shared_ptr<ROHF> hfa(new ROHF(options, ga)); 
      hfa->compute();
      atoms.insert(make_pair(i->name(), hfa->aodensity()));

      // restore cout
      cout.rdbuf(cout_orig);
    }

}
