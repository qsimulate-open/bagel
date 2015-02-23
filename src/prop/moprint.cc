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

/*
// simple local function to help things align properly
namespace {
  int digits(int x) {
    int out = (x < 10 ? 1 :
              (x < 100 ? 2 :
              (x < 1000 ? 3 :
              (x < 10000 ? 4 :
              (x < 100000 ? 5 :
              (x < 1000000 ? 6 :
              (x < 10000000 ? 7 :
              (x < 100000000 ? 8 :
              (x < 1000000000 ? 9 : 10)))))))));
    if (x < 0) out++;
    return out;
  }
}
*/

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
  norb_ = orbitals_.size();


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
    for (int i=0; i!=norb_; ++i) {
      density_.push_back(newref->relcoeff()->form_density_rhf(1, orbitals_[i], scale));
    }
    density_.push_back(newref->relcoeff()->form_density_rhf(newref->nclosed(), 0, scale));
  } else {
    // TODO Optimize - We shouldn't be storing ZMatrices with the imaginary parts all zero (let alone re-allocating them...)
    for (int i=0; i!=norb_; ++i) {
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

  assert(density_.size() == orbitals_.size()+1 && density_.size() == norb_+1);
  // The last Task will compute integrated total charge
  TaskQueue<MOPrintTask> task(ngrid_+1);
  points_.resize((ngrid_+1)*(norb_+1), 0.0);

  for (int i=0; i<=ngrid_; ++i)
    if (i % mpi__->size() == mpi__->rank())
      task.emplace_back(i, this);

  task.compute();

  mpi__->allreduce(points_.data(), points_.size());
  print();

}


void MOPrint::computepoint(const size_t pos) {
  // TODO avoid overhead?  This is repeated for each point

  shared_ptr<ZMatrix> ao_density;
  shared_ptr<ZMatrix> input_ovlp;

  if (pos != ngrid_) {
    // for each point in space
    array<double,3> tmp = {{ coords_[3*pos], coords_[3*pos+1], coords_[3*pos+2] }};
    auto ovlp = make_shared<Overlap_Point>(geom_, tmp);
    input_ovlp = make_shared<ZMatrix>(*ovlp->compute(), 1.0);
  } else {
    // total integrated overlap is also stored
    auto ovlp = make_shared<Overlap>(geom_);
    input_ovlp = make_shared<ZMatrix>(*ovlp, 1.0);
  }

  /*****/
  /*
  {
    assert(input_ovlp->ndim() == input_ovlp->mdim());
    for (int i=0; i!=input_ovlp->ndim(); ++i) {
      if (std::real(input_ovlp->element(i, i)) < -1.0e-5) {
        input_ovlp->print("AO matrix of overlaps", 40);
        cout << "pos = " << pos << ", i = " << i << ", value = " << scientific << input_ovlp->element(i, i) << endl;
      }
      assert(std::real(input_ovlp->element(i, i)) > -1.0e-5);
      assert(std::abs(std::imag(input_ovlp->element(i, i))) < 1.0e-9);
    }
  }
  */
  /*****/

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
  for (int i=0; i!=norb_+1; ++i) {
    const complex<double> out = density_[i]->dot_product(*ao_density);

    /*****/
    /*
    // density matrix has only real values along the diagonal - debug check passed
    assert(density_[i]->ndim() == density_[i]->mdim());
    for (int j=0; j!=density_[i]->ndim(); ++j) {
      cout << "j = " << j << ", value = " << scientific << density_[i]->element(j, j) << endl;
      assert(std::real(density_[i]->element(j, j)) > -1.0e-9);
      assert(std::abs(std::imag(density_[i]->element(j, j))) < 1.0e-9);
    }
    */
    /*****/

    /*
    if (std::real(out) < -1.0e-5) {
      density_[i]->print("Density matrix", 40);
      ao_density->print("AO density contributions", 40);
      cout << "orbital " << i << ", out = " << scientific << out << endl;
    }
    assert(std::real(out) > -1.0e-5);
    */
    assert(std::abs(std::imag(out)) < 1.0e-8);
    points_[(norb_+1)*pos+i] = std::real(out);
  }

}


void MOPrint::print() const {
  vector<double> density_sum = {};
  density_sum.resize(norb_+1);
  cout << fixed << setprecision(10);

  // TODO set with input
  const bool cube_format = false;

  if (cube_format) {
    throw runtime_error("Not yet implemented");
  } else {

    std::string heading = "   x-coord        y-coord        z-coord     ";
    for (int i=0; i!=norb_; ++i) {
      heading += "      Orbital " + to_string(orbitals_[i]+1);
    }
    heading += "        Total density";

    cout << heading << endl;

    for (int i=0; i!=ngrid_; ++i) {
      string line = "";
      for (int j=0; j!=3; ++j)
        line += ((coords_[3*i+j] < 0) ? "" : " ") + to_string(coords_[3*i+j]) + "  ";

      for (int j=0; j<=norb_; ++j) {
        line += ((points_[(norb_+1)*i+j] < 0) ? "" : " ") + to_string(points_[(norb_+1)*i+j]) + "  ";
        density_sum[j] += points_[(norb_+1)*i+j];
      }

      cout << line << endl;
    }


    const double scale = inc_size_[0] * inc_size_[1] * inc_size_[2];
    for (int j=0; j!=norb_; ++j)
      cout << "Sum of all gridpoints for orbital " << orbitals_[j]+1 << " = " << density_sum[j]*scale << ".  Integrated orbital density = " << points_[(norb_+1)*ngrid_+j] << "." << endl;
    cout << "Sum of all gridpoints for total density = " << density_sum.back()*scale << ".  Total integrated density = " << points_.back() << "." << endl;

    cout << "density_sum.back() = " << density_sum.back() << endl;
    cout << "density_sum[norb_] = " << density_sum[norb_] << endl;
    cout << "points_.back() = " << points_.back() << endl;
    cout << "points_[(norb_+1)*ngrid_+norb_] = " << points_[(norb_+1)*ngrid_+norb_] << endl;

    assert(density_sum.back() == density_sum[norb_]);
    assert(points_.back() == points_[(norb_+1)*ngrid_+norb_]);

  }

}

