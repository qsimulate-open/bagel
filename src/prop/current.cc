//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: current.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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


#include <src/prop/current.h>
#include <src/wfn/relreference.h>
#include <src/wfn/zreference.h>
#include <src/prop/momentum_london.h>
#include <src/prop/momentum_point.h>

using namespace std;
using namespace bagel;


Current::Current(const shared_ptr<const PTree> idata, const shared_ptr<const Geometry> geom,
                 const shared_ptr<const Reference> re) : Method(idata, geom, re) {

  // Need a GIAO-based Reference object
  auto ref_rel = dynamic_pointer_cast<const RelReference>(ref_);
  auto ref_nr  = dynamic_pointer_cast<const ZReference>(ref_);
  if (!ref_rel && !ref_nr)
    throw runtime_error("Charge currents are only available when using the result of a GIAO calculation.");
  if (ref_->nact() != 0)
    throw runtime_error("Charge currents have only been implemented for closed-shell Hartree--Fock methods.");
  assert(geom_->magnetism());
  assert(!ref_rel || !ref_nr);
  relativistic_ = ref_rel ? true : false;

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
  inc_size_ = idata->get_array<double,3>("inc_size", {{0.5, 0.5, 0.5}});
  ngrid_dim_ = idata->get_array<size_t,3>("ngrid", {{1, 1, 1}});
  ngrid_ = ngrid_dim_[0]*ngrid_dim_[1]*ngrid_dim_[2];

  if (angstrom) {
    for (int i=0; i!=3; ++i) {
      start_pos[i] /= au2angstrom__;
      inc_size_[i] /= au2angstrom__;
    }
  }

  for (int i=0; i!=ngrid_dim_[0]; ++i) {
    for (int j=0; j!=ngrid_dim_[1]; ++j) {
      for (int k=0; k!=ngrid_dim_[2]; ++k) {
        coords_.push_back(start_pos[0]+i*inc_size_[0]);
        coords_.push_back(start_pos[1]+j*inc_size_[1]);
        coords_.push_back(start_pos[2]+k*inc_size_[2]);
      }
    }
  }

  assert(ngrid_*3 == coords_.size());

  // Form density matrix
  if (relativistic_)
    density_ = ref_rel->relcoeff()->form_density_rhf(2*ref_rel->nclosed() + ref_rel->nact(), 0, 1.0);
  else
    density_ = ref_nr->zcoeff()->form_density_rhf(ref_nr->nclosed(), 0, 2.0);

  const string mtype = relativistic_ ? "Dirac-Fock" : "RHF";
  const string ctype = paramagnetic_ ? (diamagnetic_ ? "total" : "paramagnetic") : "diamagnetic";
  cout << "Computing "  << mtype << " " << ctype << " current at " << ngrid_ << " gridpoint" << ((ngrid_ > 1) ? "s" : "") << ". " << endl;

  if (relativistic_ && (!paramagnetic_ || !diamagnetic_))
    cout << "CAUTION: The diamagnetic/paramagnetic separation is not well-founded for relativistic methods.  Recommend using total current." << endl;

  cout << endl;

}


namespace bagel {
  class CurrentTask {
    protected:
      Current* parent_;
      size_t pos_;

    public:
      CurrentTask(size_t pos, Current* par)
        : parent_(par), pos_(pos) { }
      void compute() const { parent_->computepoint(pos_); }
  };
}


void Current::compute() {

  // The last Task will compute integrated currents
  TaskQueue<CurrentTask> task(ngrid_+1);
  currents_.resize(3*(ngrid_+1), 0.0);

  for (int i=0; i<=ngrid_; ++i)
    if (i % mpi__->size() == mpi__->rank())
      task.emplace_back(i, this);

  task.compute();

//  for (int i=-1; i!=ngrid_; ++i)
//    mpi__->broadcast(i->data(), 3, mpi__->rank());

  mpi__->allreduce(currents_.data(), 3*(ngrid_+1));
  print();

}


void Current::computepoint(const size_t pos) {
  array<complex<double>,3> out;
  array<shared_ptr<ZMatrix>,3> ao_current;
  array<shared_ptr<ZMatrix>,3> pi;

  if (pos == ngrid_) {
    auto mom = make_shared<Momentum_London>(geom_);
    pi = mom->compute();
  } else {
    array<double,3> tmp = {{ coords_[3*pos], coords_[3*pos+1], coords_[3*pos+2] }};
    auto mom = make_shared<Momentum_Point>(geom_, tmp);
    pi = mom->compute();
  }

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

  // Now compute total current using AO contributions
  for (int i=0; i!=3; ++i) {
    out[i] = density_->dot_product(*ao_current[i]);
  }

  currents_[3*pos+0] = out[0];
  currents_[3*pos+1] = out[1];
  currents_[3*pos+2] = out[2];

}


void Current::print() const {
  array<complex<double>,3> current_sum = {{ 0.0, 0.0, 0.0 }};
  cout << fixed << setprecision(10);
  cout << "   x-coord        y-coord        z-coord           Re(x-current)  Re(y-current)  Re(z-current)       Im(x-current)  Im(y-current)  Im(z-current)" << endl;
  for (int i=0; i!=ngrid_; ++i) {
    cout << ((coords_[3*i+0] < 0) ? "" : " ") << coords_[3*i+0] << "  "
         << ((coords_[3*i+1] < 0) ? "" : " ") << coords_[3*i+1] << "  "
         << ((coords_[3*i+2] < 0) ? "" : " ") << coords_[3*i+2] << "       "
         << ((real(currents_[3*i+0]) < 0) ? "" : " ") << real(currents_[3*i+0]) << "  "
         << ((real(currents_[3*i+1]) < 0) ? "" : " ") << real(currents_[3*i+1]) << "  "
         << ((real(currents_[3*i+2]) < 0) ? "" : " ") << real(currents_[3*i+2]) << "       "
         << ((imag(currents_[3*i+0]) < 0) ? "" : " ") << imag(currents_[3*i+0]) << "  "
         << ((imag(currents_[3*i+1]) < 0) ? "" : " ") << imag(currents_[3*i+1]) << "  "
         << ((imag(currents_[3*i+2]) < 0) ? "" : " ") << imag(currents_[3*i+2])
         << endl;
     current_sum[0] += currents_[3*i+0];
     current_sum[1] += currents_[3*i+1];
     current_sum[2] += currents_[3*i+2];
  }

  cout << endl << "Sum of all gridpoints = ( " << current_sum[0] << ", " << current_sum[1] << ", " << current_sum[2] << " ). " << endl << endl;;

  // TODO some of this should be moved out of the ``print'' function
  // determine if the gridpoints form a flat plane
  array<bool,3> single;
  for (int i=0; i!=3; ++i) single[i] = (ngrid_dim_[i] == 1);
  int direction = -1;
  if (single[0] && !single[1] && !single[2]) direction = 0;
  if (!single[0] && single[1] && !single[2]) direction = 1;
  if (!single[0] && !single[1] && single[2]) direction = 2;

  // report integrated current through that plane
  if (direction != -1) {
    double area = 1.0;
    for (int i=0; i!=3; ++i)
      if (i != direction)
        area *= inc_size_[i];
    const array<string,3> plane = {{ "yz", "xz", "xy" }};
    const complex<double> int_current = current_sum[direction]*area*au2coulomb__/au2second__*1.0e9;
    cout << endl << "Integrated current through the selected slice of the " << plane[direction] << " plane = " << int_current << " nA." << endl << endl;
  }

  cout << endl << "Current integrated over all space = ( " << currents_[3*ngrid_+0] << ", " << currents_[3*ngrid_+1] << ", " << currents_[3*ngrid_+2] << " ). " << endl << endl;;

}

