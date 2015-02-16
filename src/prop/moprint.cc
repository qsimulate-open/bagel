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
  //assert(!geom_->magnetism());
  auto newref = dynamic_pointer_cast<const RelReference>(ref_);
  if (newref) {
    relativistic_ = newref->rel();
    //assert(relativistic_);
  } else {
    relativistic_ = false;
  }

  // Determine which MOs to get
  const shared_ptr<const PTree> iorb = idata_->get_child_optional("orbitals");

  if (iorb)
    // Subtracting one so that orbitals are input in 1-based format but are stored in C format (0-based)
    for (auto& i : *iorb)
      orbitals_.push_back(lexical_cast<int>(i->data()) - 1);
  else
    for (int i=0; i!=geom_->nbasis(); ++i)
      orbitals_.push_back(i);

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
  if (newref) {
    const double scale = relativistic_ ? 1.0 : 2.0;
    for (int i=0; i!=orbitals_.size(); ++i) {
      density_.push_back(newref->relcoeff()->form_density_rhf(1, orbitals_[i], scale));
    }
    density_.push_back(newref->relcoeff()->form_density_rhf(newref->nclosed(), 0, scale));
  } else {
    // TODO Optimize - We shouldn't be storing ZMatrices with the imaginary parts all zero (let alone re-allocating them...)
    for (int i=0; i!=orbitals_.size(); ++i) {
      density_.push_back(make_shared<ZMatrix>(*ref_->coeff()->form_density_rhf(1, orbitals_[i]), 1.0));
    }
    density_.push_back(make_shared<ZMatrix>(*ref_->coeff()->form_density_rhf(ref_->nclosed(), 0), 1.0));
  }


  // TODO NONREL version, also make sure it is okay with striped vs. block coefficients...
  // TODO Should combine the spin-up and spin-down components of relativistic MOs...

  const string mtype = relativistic_ ? "relativistic" : "non-relativistic";
  cout << "Printing " << mtype << " MO densities at " << ngrid_ << " gridpoint" << ((ngrid_ > 1) ? "s" : "") << ". " << endl;

  if (relativistic_) cout << "Caution:  Currently orbital printing ignores the small components of relativistic MOs."  << endl;

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

  assert(density_.size() == orbitals_.size()+1);
  // The last Task will compute integrated total charge
  TaskQueue<MOPrintTask> task(ngrid_+1);
  points_.resize(ngrid_*(orbitals_.size()+1), 0.0);

  for (int i=0; i<=ngrid_; ++i)
    if (i % mpi__->size() == mpi__->rank())
      task.emplace_back(i, this);

  task.compute();

  mpi__->allreduce(points_.data(), points_.size());
  print();
  assert(density_.size() == orbitals_.size()+1);

}


void MOPrint::computepoint(const size_t pos) {
  assert(density_.size() == orbitals_.size()+1);
  // TODO avoid overhead?  This is repeated for each point

  shared_ptr<ZMatrix> ao_density;
  shared_ptr<ZMatrix> input_ovlp;

  if (pos == ngrid_) {
    auto ovlp = make_shared<Overlap>(geom_);
    input_ovlp = make_shared<ZMatrix>(*ovlp, 1.0);
  } else {
    array<double,3> tmp = {{ coords_[3*pos], coords_[3*pos+1], coords_[3*pos+2] }};
    auto ovlp = make_shared<Overlap_Point>(geom_, tmp);
    input_ovlp = make_shared<ZMatrix>(*ovlp->compute(), 1.0);
  }

  // First build the current matrix in AO basis
  if (relativistic_) {
    const int n = geom_->nbasis();
    ao_density = make_shared<ZMatrix>(4*n, 4*n);

    const complex<double> re( 0.1,  0.0);
    const complex<double> im( 0.0,  1.0);

    // Assumes RMB basis
    ao_density->add_block( re, 0*n, 0*n, n, n, *input_ovlp);
    ao_density->add_block( re, 1*n, 1*n, n, n, *input_ovlp);
    // TODO Need to fill in the small component blocks too - requires point_kinetic, I guess?  Could also use Small1e<Overlap_point>...
  } else {
    ao_density = make_shared<ZMatrix>(*input_ovlp);
  }

  // Now compute total MO density using AO contributions
  for (int i=0; i!=density_.size(); ++i) {
    const complex<double> out = density_[i]->dot_product(*ao_density);
    assert(std::abs(std::imag(out)) < 1.0e-8);
    points_[(orbitals_.size()+1)*pos+i] = std::real(out);
  }

  assert(density_.size() == orbitals_.size()+1);
}


void MOPrint::print() const {
  assert(density_.size() == orbitals_.size()+1);
  vector<double> density_sum = {};
  density_sum.resize(density_.size());
  cout << fixed << setprecision(10);
  cout << "   x-coord        y-coord        z-coord           Orbital 1      Orbital 2      Orbital 3           Orbital 4      Orbital 5      Orbital 6    " << endl;
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
     density_sum[0] += points_[3*i+0];
     //density_sum[1] += points_[3*i+1];
     //density_sum[2] += points_[3*i+2];
  }

  cout << endl << "Sum of all gridpoints = ( " << density_sum[0] << ", " << density_sum[1] << ", " << density_sum[2] << " ). " << endl << endl;;

  cout << endl << "Orbital density integrated over all space = ( " << points_[3*ngrid_+0] << " ). " << endl << endl;;

}

