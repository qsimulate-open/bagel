//
// BAGEL - Parallel electron correlation program.
// Filename: current.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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


#include <src/london/current.h>
#include <src/rel/relreference.h>
#include <src/prop/momentum_london.h>
#include <src/prop/momentum_point.h>

using namespace std;
using namespace bagel;


Current::Current(const std::shared_ptr<const PTree> idata, const std::shared_ptr<const Geometry> geom,
                 const std::shared_ptr<const Reference> re) : Method(idata, geom, re) {

  // Need a GIAO-based Reference object
  assert(geom_->magnetism());
  auto newref = dynamic_pointer_cast<const RelReference>(ref_);
  if (!newref)
    throw runtime_error("Charge currents are only available when using the result of a GIAO calculation.");
  relativistic_ = newref->rel();

  // Determine type of current desired
  const string type = to_lower(idata->get<string>("current_type", "all"));
  if (type == "all" || type == "total") {
    paramagnetic_ = true;
    diamagnetic_ = true;
  } else if (type == "diamagnetic" || type == "dia") {
    paramagnetic_ = false;
    diamagnetic_ = true;
    // TODO modify integral to return just para or dia contribution
    throw runtime_error("Not yet implemented");
  } else if (type == "paramagnetic" || type == "para") {
    paramagnetic_ = true;
    diamagnetic_ = false;
    throw runtime_error("Not yet implemented");
  } else {
    throw runtime_error("Current type not understood - should be para, dia, or all (default)");
  }

  // Determine coordinates where current will be computed
  const bool angstrom = idata->get<bool>("angstrom", false);
  array<double,3> start_pos = idata->get_array<double,3>("start_pos", {{0.0, 0.0, 0.0}});
  array<double,3> inc_size = idata->get_array<double,3>("inc_size", {{0.5, 0.5, 0.5}});
  const array<int,3> ngrid = idata->get_array<int,3>("ngrid", {{1, 1, 1}});

  if (angstrom) {
    for (int i=0; i!=3; ++i) {
      start_pos[i] /= au2angstrom__;
      inc_size[i] /= au2angstrom__;
    }
  }

  for (int i=0; i!=ngrid[0]; ++i) {
    for (int j=0; j!=ngrid[1]; ++j) {
      for (int k=0; k!=ngrid[2]; ++k) {
        const array<double,3> tmp = {{ start_pos[0]+i*inc_size[0], start_pos[1]+j*inc_size[1], start_pos[2]+k*inc_size[2] }};
        coords_.push_back(tmp);
      }
    }
  }

  // Form density matrix
  const double scale = relativistic_ ? 1.0 : 2.0;
  density_ = newref->relcoeff()->form_density_rhf(newref->nclosed(), 0, scale);

  const string mtype = relativistic_ ? "Dirac-Fock" : "RHF";
  const string ctype = paramagnetic_ ? (diamagnetic_ ? "total" : "paramagnetic") : "diamagnetic";
  cout << "Computing "  << mtype << " " << ctype << " current at " << coords_.size() << " gridpoint" << ((coords_.size() > 1) ? "s" : "") << ". " << endl;

  if (relativistic_ && (!paramagnetic_ || !diamagnetic_))
    cout << "CAUTION: The diamagnetic/paramagnetic separation is not well-founded for relativistic methods.  Recommend using total current." << endl;

  cout << endl;

}


void Current::compute() {

  currents_.resize(coords_.size());

  // TODO thread, parallel
  for (int i=0; i!=coords_.size(); ++i) {
    auto mom = make_shared<Momentum_Point>(geom_, coords_[i]);
    array<shared_ptr<ZMatrix>,3> pi = mom->compute();
    currents_[i] = computepoint(pi);
  }

  {
    auto mom = make_shared<Momentum_London>(geom_);
    array<shared_ptr<ZMatrix>,3> pi = mom->compute();
    total_current_ = computepoint(pi);
  }

  print();

}


array<double,3> Current::computepoint(array<shared_ptr<ZMatrix>,3> pi) {
  array<double,3> out;
  array<shared_ptr<ZMatrix>,3> ao_current;

  // First build the current matrix in AO basis
  if (relativistic_) {
    const int n = geom_->nbasis();
    ao_current[0] = make_shared<ZMatrix>(4*n, 4*n);
    ao_current[1] = make_shared<ZMatrix>(4*n, 4*n);
    ao_current[2] = make_shared<ZMatrix>(4*n, 4*n);

    const complex<double> re(-0.5,  0.0);
    const complex<double> im( 0.0, -0.5);

    // Assumes RMB basis
    // TODO Probably there's some elegant way to build this using 2n by 2n block matrices
    ao_current[0]->add_block( re, 0*n, 2*n, n, n, *pi[0]);
    ao_current[0]->add_block( im, 0*n, 2*n, n, n, *pi[1]);
    ao_current[0]->add_block(-re, 0*n, 3*n, n, n, *pi[2]);
    ao_current[0]->add_block( re, 1*n, 2*n, n, n, *pi[2]);
    ao_current[0]->add_block( re, 1*n, 3*n, n, n, *pi[0]);
    ao_current[0]->add_block(-im, 1*n, 3*n, n, n, *pi[1]);

    ao_current[0]->add_block( re, 2*n, 0*n, n, n, *pi[0]);
    ao_current[0]->add_block(-im, 2*n, 0*n, n, n, *pi[1]);
    ao_current[0]->add_block( re, 2*n, 1*n, n, n, *pi[2]);
    ao_current[0]->add_block(-re, 3*n, 0*n, n, n, *pi[2]);
    ao_current[0]->add_block( re, 3*n, 1*n, n, n, *pi[0]);
    ao_current[0]->add_block( im, 3*n, 1*n, n, n, *pi[1]);

    ao_current[1]->add_block(-im, 0*n, 2*n, n, n, *pi[0]);
    ao_current[1]->add_block( re, 0*n, 2*n, n, n, *pi[1]);
    ao_current[1]->add_block( im, 0*n, 3*n, n, n, *pi[2]);
    ao_current[1]->add_block( im, 1*n, 2*n, n, n, *pi[2]);
    ao_current[1]->add_block( im, 1*n, 3*n, n, n, *pi[0]);
    ao_current[1]->add_block( re, 1*n, 3*n, n, n, *pi[1]);

    ao_current[1]->add_block( im, 2*n, 0*n, n, n, *pi[0]);
    ao_current[1]->add_block( re, 2*n, 0*n, n, n, *pi[1]);
    ao_current[1]->add_block(-im, 2*n, 1*n, n, n, *pi[2]);
    ao_current[1]->add_block(-im, 3*n, 0*n, n, n, *pi[2]);
    ao_current[1]->add_block(-im, 3*n, 1*n, n, n, *pi[0]);
    ao_current[1]->add_block( re, 3*n, 1*n, n, n, *pi[1]);

    ao_current[2]->add_block( re, 0*n, 2*n, n, n, *pi[2]);
    ao_current[2]->add_block( re, 0*n, 3*n, n, n, *pi[0]);
    ao_current[2]->add_block(-im, 0*n, 3*n, n, n, *pi[1]);
    ao_current[2]->add_block(-re, 1*n, 2*n, n, n, *pi[0]);
    ao_current[2]->add_block(-im, 1*n, 2*n, n, n, *pi[1]);
    ao_current[2]->add_block( re, 1*n, 3*n, n, n, *pi[2]);

    ao_current[2]->add_block( re, 2*n, 0*n, n, n, *pi[2]);
    ao_current[2]->add_block(-re, 2*n, 1*n, n, n, *pi[0]);
    ao_current[2]->add_block( im, 2*n, 1*n, n, n, *pi[1]);
    ao_current[2]->add_block( re, 3*n, 0*n, n, n, *pi[0]);
    ao_current[2]->add_block( im, 3*n, 0*n, n, n, *pi[1]);
    ao_current[2]->add_block( re, 3*n, 1*n, n, n, *pi[2]);
  } else {
    for (int i=0; i!=3; ++i) {
      ao_current[i] = make_shared<ZMatrix>(-1.0**(pi[i]));
    }
  }

  // Now compute total current using AO contributions (discarding imag. part)
  for (int i=0; i!=3; ++i) {
    out[i] = std::real(density_->dot_product(*ao_current[i]));
  }

  return out;
}


void Current::print() {
  cout << fixed << setprecision(10);
  cout << "   x-coord        y-coord        z-coord        x-current      y-current      z-current" << endl;
  for (int i=0; i!=coords_.size(); ++i) {
    cout << ((coords_[i][0] < 0) ? "" : " ") << coords_[i][0] << "  "
         << ((coords_[i][1] < 0) ? "" : " ") << coords_[i][1] << "  "
         << ((coords_[i][2] < 0) ? "" : " ") << coords_[i][2] << "  "
         << ((currents_[i][0] < 0) ? "" : " ") << currents_[i][0] << "  "
         << ((currents_[i][1] < 0) ? "" : " ") << currents_[i][1] << "  "
         << ((currents_[i][2] < 0) ? "" : " ") << currents_[i][2] << "  "
         << endl;
  }

  cout << endl << "Total integrated current = ( " << total_current_[0] << ", " << total_current_[1] << ", " << total_current_[2] << " ). " << endl << endl;;

}
