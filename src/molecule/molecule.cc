//
// BAGEL - Parallel electron correlation program.
// Filename: molecule.cc
// Copyright (C) 2013 Toru Shiozaki
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

#include <src/molecule/molecule.h>

using namespace std;
using namespace bagel;

double Molecule::compute_nuclear_repulsion() {
  double out = 0.0;
  for (auto iter = atoms_.begin(); iter != atoms_.end(); ++iter) {
    const double c = (*iter)->atom_charge();
    for (auto titer = iter + 1; titer != atoms_.end(); ++titer) {
      // nuclear repulsion between dummy atoms are not computed here (as in Molpro)
      if (!(*iter)->dummy() || !(*titer)->dummy()) {
        const double dist = (*iter)->distance(*titer);
        const double charge = c * (*titer)->atom_charge();
// it turned out that the deviation is numerically zero in double precision if the exponent is more than 1.0e-6. I will leave them for a while
#ifdef BEYOND_DOUBLE
        if (!(*iter)->finite_nucleus() && !(*titer)->finite_nucleus()) { // both point charges
#endif
          out += charge / dist;
#ifdef BEYOND_DOUBLE
        } else if ((*iter)->finite_nucleus() && (*titer)->finite_nucleus()) { // both gaussian charges
          const double gamma0 = (*iter)->finite_nucleus();
          const double gamma1 = (*titer)->finite_nucleus();
          out += charge * erf(sqrt(gamma0 * gamma1 / (gamma0 + gamma1)) * dist) / dist;
        } else { // one point charge, the other gaussian
          const double gamma = (*iter)->finite_nucleus() ? (*iter)->atom_exponent() : (*titer)->atom_exponent();
          out += charge * erf(sqrt(gamma) * dist) / dist;
        }
#endif
      }
    }
  }
  return out;
}


void Molecule::print_atoms() const {
  cout << "  *** Geometry ***" << endl << endl;
  cout << "  Symmetry: " << symmetry() << endl;
  cout << endl;
  for (auto i : atoms_) i->print();
  cout << endl;
}


array<double,3> Molecule::charge_center() const {
  array<double,3> out{{0.0, 0.0, 0.0}};
  double sum = 0.0;
  for (auto& i : atoms_) {
    out[0] += i->atom_charge() * i->position(0);
    out[1] += i->atom_charge() * i->position(1);
    out[2] += i->atom_charge() * i->position(2);
    sum += i->atom_charge();
  }
  out[0] /= sum;
  out[1] /= sum;
  out[2] /= sum;
  return out;
}


array<double,6> Molecule::quadrupole() const {
  array<double,6> out;
  array<double,3> c = charge_center();
  for (auto& i : atoms_) {
    out[0] += i->atom_charge() * pow(i->position(0) - c[0], 2);
    out[1] += i->atom_charge() * (i->position(0) - c[0]) * (i->position(1) - c[1]);
    out[2] += i->atom_charge() * (i->position(0) - c[0]) * (i->position(2) - c[2]);
    out[3] += i->atom_charge() * pow(i->position(1) - c[1], 2);
    out[4] += i->atom_charge() * (i->position(1) - c[1]) * (i->position(2) - c[2]);
    out[5] += i->atom_charge() * pow(i->position(2) - c[2], 2);
  }
  return out;
}


bool Molecule::has_finite_nucleus() const {
  return any_of(atoms_.begin(), atoms_.end(), [](shared_ptr<const Atom> a) { return a->finite_nucleus(); });
}


shared_ptr<const XYZFile> Molecule::xyz() const {
  auto out = make_shared<XYZFile>(natom());
  int iat = 0;
  for (auto& i : atoms_) {
    out->element(0, iat) = i->position(0);
    out->element(1, iat) = i->position(1);
    out->element(2, iat) = i->position(2);
    ++iat;
  }
  return out;
}
