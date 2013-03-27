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
#include <src/util/diis.h>
#include <src/scf/rohf.h>

using namespace std;
using namespace bagel;

const static AtomMap atommap_;

AtomicDensities::AtomicDensities(std::shared_ptr<const Geometry> g) : Matrix(g->nbasis(), g->nbasis()), geom_(g) {
  // first make a list of unique atoms
  const string basis = geom_->basisfile();
  const string dfbasis = geom_->auxfile();
  map<string, shared_ptr<const Matrix>> atoms;
  unique_ptr<double[]> eig(new double[geom_->nbasis()]);

  int offset = 0;
  for (auto& i : geom_->atoms()) {
    if (atoms.find(i->name()) == atoms.end()) {
      // dummy buffer to suppress the output
//    stringstream ss;
//    std::streambuf* cout_orig = cout.rdbuf();
//    cout.rdbuf(ss.rdbuf());

      shared_ptr<const Atom> atom(new Atom(i->spherical(), i->name(), {{0.0,0.0,0.0}}, basis));

      multimap<string,string> geomop;
      geomop.insert(make_pair("basis", basis));
      geomop.insert(make_pair("df_basis", dfbasis.empty() ? basis : dfbasis));
      shared_ptr<const Geometry> ga(new Geometry({atom}, geomop));
      atoms.insert(make_pair(i->name(), compute_atomic(ga)));

      // restore cout
//    cout.rdbuf(cout_orig);
    }
    copy_block(offset, offset, i->nbasis(), i->nbasis(), atoms[i->name()]->data());
    offset += i->nbasis();
  }

}


shared_ptr<const Matrix> AtomicDensities::compute_atomic(shared_ptr<const Geometry> ga) const {
  shared_ptr<const Overlap> overlap(new Overlap(ga));
  TildeX tildex(overlap, 1.0e-8);
  Hcore hcore(ga);
  DIIS<Matrix> diis(5);
  shared_ptr<const Matrix> coeff;
  unique_ptr<double[]> eig(new double[ga->nbasis()]);
  {
    Matrix ints = tildex % hcore * tildex;
    ints.diagonalize(eig.get());
    coeff = shared_ptr<const Matrix>(new Matrix(tildex * ints));
  }

  tuple<int,int,int,int> nclosed = atommap_.num_closed(ga->atoms().front()->name());
  tuple<int,int,int,int> nopen   = atommap_.num_open(ga->atoms().front()->name());
  if (get<0>(nclosed)+get<1>(nclosed)+get<2>(nclosed)+get<3>(nclosed)+
      get<0>(nopen)+get<1>(nopen)+get<2>(nopen)+get<3>(nopen) != ga->nele()) throw logic_error("Inconsistent nclosed and nopen. See AtomMap"); 

  int iter = 0;
  for (; iter != 20; ++iter) {
//  Fock<1> fock(ga, hcore, std::shared_ptr<const Matrix>(), coeff_->slice(0, nocc), true);
  }
//if (iter == 20) throw runtime_error("spin-averaged atomic HF did not converge");

  return shared_ptr<const Matrix>(new Matrix(ga->nbasis(), ga->nbasis()));
}
