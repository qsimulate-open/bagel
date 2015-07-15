//
// BAGEL - Parallel electron correlation program.
// Filename: simulationcell.cc
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#include <src/periodic/simulationcell.h>

using namespace std;
using namespace bagel;

SimulationCell::SimulationCell(const shared_ptr<const Geometry> g, const int n) : atoms_(g->atoms()), ndim_(n) { }


SimulationCell::SimulationCell(vector<shared_ptr<const Atom>> a, const int n) : atoms_(a), ndim_(n) { }


SimulationCell::init() {

  radius_ = 0.0;
  nbasis_ = 0;
  const array<double, 3> charge_centre = charge_centre();
  for (auto& atom : atoms_) {
    nbasis_ += atom->nbasis();
    array<double, 3> r;
    r[0] = atom->position(0) - charge_centre[0];
    r[1] = atom->position(1) - charge_centre[1];
    r[2] = atom->position(2) - charge_centre[2];
    const double rsq = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
    if (rsq > out)
      radius_ = rsq;
  }
  radius_ = sqrt(rsq);
}


array<double,3> SimulationCell::charge_center() const {

  array<double,3> out{{0.0, 0.0, 0.0}};
  double sum = 0.0;
  for (auto& atom : atoms_) {
    out[0] += atom->atom_charge() * atom->position(0);
    out[1] += atom->atom_charge() * atom->position(1);
    out[2] += atom->atom_charge() * atom->position(2);
    sum += atom->atom_charge();
  }
  out[0] /= sum;
  out[1] /= sum;
  out[2] /= sum;

  return out;
}


void SimulationCell::compute_extent(const double thresh) { // extent at charge centre

  vector<shared_ptr<const Shell>> shells;
  for (auto& atom : atoms_)
    shells.insert(shells.end(), atom->shells().begin(), atom->shells().end());

  extent_ = 0.0;
  for (auto& ish : shells)
    for (auto& jsh : shells) {
      const vector<double> exp0 = ish->exponents();
      const vector<double> exp1 = jsh->exponents();
      array<double, 3> AB;
      AB[0] = ish->position(0) - jsh->position(0);
      AB[1] = ish->position(1) - jsh->position(1);
      AB[2] = ish->position(2) - jsh->position(2);
      const double rsq = AB[0] * AB[0] + AB[1] * AB[1] + AB[2] * AB[2];
      const double lnthresh = log(thresh);

      for (auto& expi0 : exp0) {
        for (auto& expi1 : exp1) {
          const double cxp_inv = 1.0 / (expi0 + expi1);
          const double expi01 = expi0 * expi1;
          const double lda_kl = sqrt((- lnthresh - expi01 * rsq * cxp_inv + 0.75 * log(4.0 * expi01 / (pi__ * pi__))) * cxp_inv);

          array<double, 3> tmp;
          tmp[0] = (ish->position(0) * expi0 + jsh->position(0) * expi1) * cxp_inv - centre(0);
          tmp[1] = (ish->position(1) * expi0 + jsh->position(1) * expi1) * cxp_inv - centre(1);
          tmp[2] = (ish->position(2) * expi0 + jsh->position(2) * expi1) * cxp_inv - centre(2);

          const double extent0 = sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2]) + lda_kl;
          if (extent0 > extent_) extent_ = extent0;
        }
      }
    }
}


void SimulationCell::compute_multipoles(const int lmax) {

  const int nmultipole = (lmax + 1) * (lmax + 1);
  multipoles_.resize(nmultipole);
  vector<shared_ptr<ZMatrix>> multipoles(nmultipole);
  for (int i = 0; i != nmultipole; ++i)
    multipoles[i] = make_shared<ZMatrix>(nbasis_, nbasis_);

  bool skip = false;
  size_t ob0 = 0;
  for (auto& atom0 : atoms_) {
    for (auto& b0 : atom0->shells()) {
      size_t ob1 = 0;
      for (auto& atom1 : atoms_) {
        for (auto& b1 : atom1->shells()) {
          MultipoleBatch mpole(array<shared_ptr<const Shell>, 2>{{b1, b0}}, position_, lmax);
          mpole.compute();
          for (int i = 0; i != nmultipole; ++i)
            ltipoles[i]->copy_block(ob1, ob0, b1->nbasis(), b0->nbasis(), mpole.data(i));

          ob1 += b1->nbasis();
        }
      }
    }
    ob0 += b0->nbasis();
  }

  for (int i = 0; i != nmultipole; ++i)
    multipoles_[i] = multipoles[i];
}


#if 0
void SimulationCell::compute_Sn(const double thresh, const int max_iter) { // S(n+1) = U_M[S(n)] O* + M*

  for (int iter = 0; i != max_iter; ++i) {






    const double error = 0.0;
    if (error < thresh) {
      cout << "  * Sn converged." << endl << endl;
      break;
    } else if (iter == max_iter-1) {
      cout << "  * Max iteration reached when in compute_Sn." << endl << endl;
      break;
    }
  }
}


vector<const array<double, 3>> SimulationCell::get_intermediate_layer(const int imax) {

  ws = pow(2*ncell_+1, ndim_);

  vector<const array<double, 3>> out;

  int cnt = 0;
  for (auto& L : lattice_vectors_) {

  }
}
#endif
