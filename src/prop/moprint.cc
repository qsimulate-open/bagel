//
// BAGEL - Parallel electron correlation program.
// Filename: moprint.cc
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


#include <src/prop/moprint.h>
#include <src/wfn/relreference.h>
#include <src/mat1e/overlap.h>
#include <src/mat1e/giao/zoverlap.h>
#include <src/prop/overlap_point.h>

using namespace std;
using namespace bagel;


MOPrint::MOPrint(const std::shared_ptr<const PTree> idata, const std::shared_ptr<const Geometry> geom,
                 const std::shared_ptr<const Reference> re) : Method(idata, geom, re) {

  // For now, disable for GIAO wavefunctions
  assert(!geom_->magnetism());
  auto newref = dynamic_pointer_cast<const RelReference>(ref_);
  if (!newref) {
    relativistic_ = newref->rel();
    assert(relativistic_);
  } else {
    relativistic_ = false;
  }

  // Determine which MOs to get
  vector<int> orbitals_;
  const shared_ptr<const PTree> iorb = idata_->get_child_optional("orbitals");

  if (iorb)
    // Subtracting one so that orbitals are input in 1-based format but are stored in C format (0-based)
    for (auto& i : *iorb)
      active_indices.push_back(lexical_cast<int>(i->data()) - 1);
  else
    for (int i=0; i!=geom_->nbasis(); ++i)
      active_indices.push_back(i);

  // Determine coordinates where current will be computed
  const bool angstrom = idata->get<bool>("angstrom", false);
  array<double,3> start_pos = idata->get_array<double,3>("start_pos", {{-10.0, -10.0, -10.0}});
  inc_size_ = idata->get_array<double,3>("inc_size", {{0.25, 0.25, 0.25}});
  ngrid_dim_ = idata->get_array<size_t,3>("ngrid", {{80, 80, 80}});
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

  // Form density matrices
  const double scale = relativistic_ ? 1.0 : 2.0;
  for (int i=0; i!=orbitals_.size(); ++i)
    density_.back() = newref->relcoeff()->form_density_rhf(1, orbitals[i], scale);
  density_.push_back(newref->relcoeff()->form_density_rhf(newref->nclosed(), 0, scale));

  // TODO NONREL version, also make sure it is okay with striped vs. block coefficients...
  // TODO Should combine the spin-up and spin-down components of relativistic MOs...

  const string mtype = relativistic_ ? "relativistic" : "non-relativistic";
  cout << "Printing " << mtype << " MO densities at " << ngrid_ << " gridpoint" << ((ngrid_ > 1) ? "s" : "") << ". " << endl;

  cout << endl;

}


namespace bagel {
  class MOPrintTask {
    protected:
      MOPrint* parent_;
      size_t pos_;

    public:
      MOPrintTask(size_t pos, MOPrint* par)
        : parent_(par), pos_(pos) { }
      void compute() const { parent_->computepoint(pos_); }
  };
}


void MOPrint::compute() {

  // The last Task will compute integrated total charge
  TaskQueue<MOPrintTask> task(ngrid_+1);
  points_.resize(ngrid_*(orbitals_.size()+1), 0.0);

  for (int i=0; i<=ngrid_; ++i)
    if (i % mpi__->size() == mpi__->rank())
      task.emplace_back(i, this);

  task.compute();

  mpi__->allreduce(points_.data(), points_.size());
  print();

}


void MOPrint::computepoint(const size_t pos) {
  // TODO avoid overhead?  This is repeated for each point
  vector<complex<double>> out = {};
  out.resize(orbitals_.size() + i);

  shared_ptr<ZMatrix> ao_density;
  shared_ptr<ZMatrix> input_ovlp;

  if (pos == ngrid_) {
    auto ovlp = make_shared<Overlap>(geom_);
    input_ovlp = make_shared<ZMatrix>(ovlp->compute(), 1.0);
  } else {
    array<double,3> tmp = {{ coords_[3*pos], coords_[3*pos+1], coords_[3*pos+2] }};
    auto ovlp = make_shared<Overlap_Point>(geom_, tmp);
    input_ovlp = make_shared<ZMatrix>(ovlp->compute(), 1.0);
  }

  // First build the current matrix in AO basis
  if (relativistic_) {
    const int n = geom_->nbasis();
    ao_density[0] = make_shared<ZMatrix>(4*n, 4*n);

    const complex<double> re( 0.1,  0.0);
    const complex<double> im( 0.0,  1.0);

    // Assumes RMB basis
    ao_density->add_block( re, 0*n, 0*n, n, n, *input_ovlp);
    ao_density->add_block( re, 1*n, 1*n, n, n, *input_ovlp);
    assert(false);
    // TODO Need to fill in the small component blocks too - requires point_kinetic, I guess?  Could also use Small1e<Overlap_point>...
  } else {
    ao_density = make_shared<ZMatrix>(input_ovlp);
  }

  // Now compute total MO density using AO contributions
  for (int i=0; i<=density_.size(); ++i)
    out[i] = density_[i]->dot_product(*ao_density);

  points_[(orbitals_.size()+1)*pos+i] = out[i];

}


void MOPrint::print() const {
  array<complex<double>,3> current_sum = {{ 0.0, 0.0, 0.0 }};
  cout << fixed << setprecision(10);
  cout << "   x-coord        y-coord        z-coord           Re(x-current)  Re(y-current)  Re(z-current)       Im(x-current)  Im(y-current)  Im(z-current)" << endl;
  for (int i=0; i!=ngrid_; ++i) {
    cout << ((coords_[3*i+0] < 0) ? "" : " ") << coords_[3*i+0] << "  "
         << ((coords_[3*i+1] < 0) ? "" : " ") << coords_[3*i+1] << "  "
         << ((coords_[3*i+2] < 0) ? "" : " ") << coords_[3*i+2] << "       "
         << ((real(points_[3*i+0]) < 0) ? "" : " ") << real(points_[3*i+0]) << "  "
         << ((imag(points_[3*i+0]) < 0) ? "" : " ") << imag(points_[3*i+0])
         //<< ((real(points_[3*i+1]) < 0) ? "" : " ") << real(points_[3*i+1]) << "  "
         //<< ((real(points_[3*i+2]) < 0) ? "" : " ") << real(pionts_[3*i+2]) << "       "
         //<< ((imag(points_[3*i+0]) < 0) ? "" : " ") << imag(points_[3*i+0]) << "  "
         //<< ((imag(points_[3*i+1]) < 0) ? "" : " ") << imag(points_[3*i+1]) << "  "
         << endl;
     current_sum[0] += points_[3*i+0];
     //current_sum[1] += points_[3*i+1];
     //current_sum[2] += points_[3*i+2];
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

  cout << endl << "Orbital density integrated over all space = ( " << points_[3*ngrid_+0] << " ). " << endl << endl;;

}

