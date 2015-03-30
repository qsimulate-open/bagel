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
#include <src/util/muffle.h>
#include <src/wfn/relreference.h>
#include <src/mat1e/overlap.h>
#include <src/mat1e/giao/zoverlap.h>
#include <src/prop/overlap_point.h>

using namespace std;
using namespace bagel;


MOPrint::MOPrint(const std::shared_ptr<const PTree> idata, const std::shared_ptr<const Geometry> geom,
                 const std::shared_ptr<const Reference> re) : Method(idata, geom, re) {

  // Assumes striped (or non-relativistic) coefficient matrix
  auto newref = dynamic_pointer_cast<const RelReference>(ref_);
  if (newref) {
    relativistic_ = newref->rel();
    paired_ = relativistic_ ? !geom->magnetism() : true;
  } else {
    relativistic_ = false;
    paired_ = true;
  }
  paired_ = idata_->get<bool>("paired", paired_);
  if (paired_ && !relativistic_)
    throw runtime_error("Individual spin-orbitals can only be visualized for 4-component methods.");

  // Determine which MOs to get
  const shared_ptr<const PTree> iorb = idata_->get_child_optional("orbitals");

  if (iorb) {
    // Subtracting one so that orbitals are input in 1-based format but are stored in C format (0-based)
    for (auto& i : *iorb)
      orbitals_.push_back(lexical_cast<int>(i->data()) - 1);
  } else {
    const int mult = paired_ ? 1 : 2;
    for (int i=0; i!=mult*geom_->nbasis(); ++i)
      orbitals_.push_back(i);
  }
  norb_ = orbitals_.size();


  // Determine gridpoints where density will be computed
  const bool angstrom = idata->get<bool>("angstrom", false);
  ngrid_dim_ = idata->get_array<size_t,3>("ngrid", {{61, 61, 61}});
  ngrid_ = ngrid_dim_[0]*ngrid_dim_[1]*ngrid_dim_[2];
  array<double,3> start_pos = idata->get_array<double,3>("start_pos", {{-99.0, -99.0, -99.0}});
  inc_size_ = idata->get_array<double,3>("inc_size", {{0.25, 0.25, 0.25}});

  // By default, assign positions covering the molecule + 4.0 Bohr on each side
  if (start_pos == array<double,3>({{-99.0, -99.0, -99.0}})) {
    array<double,3> max_pos = geom_->atoms(0)->position();
    array<double,3> min_pos = geom_->atoms(0)->position();
    for (auto& i : geom_->atoms()) {
      for (int j=0; j!=3; ++j) {
        if (i->position(j) > max_pos[j])
          max_pos[j] = i->position(j);
        if (i->position(j) < min_pos[j])
          min_pos[j] = i->position(j);
      }
    }
    for (int i=0; i!=3; ++i) {
      start_pos[i] = min_pos[i] - 4.0;
      inc_size_[i] = (max_pos[i] - min_pos[i] + 8.0) / (ngrid_dim_[i] - 1);
    }
  }

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
    // GIAO or 4-component wavefunction
    const double scale = relativistic_ ? 1.0 : 2.0;
    const int ncol = (relativistic_ && paired_) ? 2 : 1;

    // TODO Fix RelReference so this is not needed
    // Ugly patch to compensate for RelReference::nclosed() giving number of spin-orbitals after Dirac but spatial orbitals after RHF_London or ZCASSCF
    const int patch = (newref->nact() == 0) ? 1 : 2;

    for (int i=0; i!=norb_; ++i) {
      density_.push_back(newref->relcoeff()->form_density_rhf(ncol, ncol*orbitals_[i], scale));
    }
    density_.push_back(newref->relcoeff()->form_density_rhf(patch*newref->nclosed(), 0, scale));

  } else {
    // Conventional non-relativistic wavefunction
    // TODO This could be optimized, since we are storing real matrices as complex.  I guess it would require templating...
    for (int i=0; i!=norb_; ++i) {
      density_.push_back(make_shared<const ZMatrix>(*ref_->coeff()->form_density_rhf(1, orbitals_[i]), 1.0));
    }
    density_.push_back(make_shared<const ZMatrix>(*ref_->coeff()->form_density_rhf(ref_->nclosed(), 0), 1.0));
  }

  const string mtype = relativistic_ ? "relativistic" : "non-relativistic";
  const string stype = paired_ ? "spatial orbital" : "spin-orbital";
  cout << "Printing " << norb_+1 << " " << mtype << " " << stype << " densities at " << ngrid_ << " gridpoint" << ((ngrid_ > 1) ? "s" : "") << ". " << endl;
  if (relativistic_) cout << "Note that orbital printing ignores the small components of relativistic MOs."  << endl;

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

  cout << "Orbital printout computation finished; generating output files." << endl;

  print();

}


void MOPrint::computepoint(const size_t pos) {
  shared_ptr<ZMatrix> ao_density;
  shared_ptr<ZMatrix> input_ovlp;
  array<double,3> tmp;

  if (!geom_->london()) {

    // Standard basis
    if (pos != ngrid_) {
      // for each point in space
      tmp = {{ coords_[3*pos], coords_[3*pos+1], coords_[3*pos+2] }};
      auto ovlp = make_shared<Overlap_Point>(geom_, tmp);
      input_ovlp = make_shared<ZMatrix>(*ovlp->compute(), 1.0);
    } else {
      // total integrated overlap is also stored
      auto ovlp = make_shared<Overlap>(geom_);
      input_ovlp = make_shared<ZMatrix>(*ovlp, 1.0);
    }
  } else {

    // GIAO version
    if (pos != ngrid_) {
      // for each point in space
      tmp = {{ coords_[3*pos], coords_[3*pos+1], coords_[3*pos+2] }};
      auto ovlp = make_shared<Overlap_Point_London>(geom_, tmp);
      input_ovlp = ovlp->compute();
    } else {
      // total integrated overlap is also stored
      input_ovlp = make_shared<ZOverlap>(geom_);
    }
  }

  // First build the overlap matrix in AO basis
  if (!relativistic_) {
    ao_density = input_ovlp;
  } else {
    const int n = geom_->nbasis();
    ao_density = make_shared<ZMatrix>(4*n, 4*n);

    // Assumes RKB (or RMB) basis
    // small component neglected, although we could compute and add in the appropriate blocks
    ao_density->add_block(1.0, 0*n, 0*n, n, n, *input_ovlp);
    ao_density->add_block(1.0, 1*n, 1*n, n, n, *input_ovlp);
  }

  // Now compute total MO density using AO contributions
  for (int i=0; i!=norb_+1; ++i) {
    const complex<double> out = density_[i]->dot_product(*ao_density);
    assert(std::abs(std::imag(out)) < 1.0e-8);
    points_[(norb_+1)*pos+i] = std::real(out);
  }

}


void MOPrint::print() const {
  vector<double> density_sum = {};
  density_sum.resize(norb_+1);
  const bool cube_format = idata_->get<bool>("cube", true);

  if (cube_format) {
    for (int i=0; i<=norb_; ++i) {
      const string title = (i == norb_) ? "density" : "mo_" + to_string(orbitals_[i]+1);
      Muffle muffle(title + ".cub");
      cout << "BAGEL generated cube file." << endl;
      if (i == norb_)
        cout << "Full electronic density" << endl;
      else
        cout << "Molecular orbital " << orbitals_[i]+1 << endl;

      // Number of atoms, and origin of the volumetric cata
      cout << fixed << setprecision(6) << setw(5) << geom_->natom() << setw(12) << coords_[0] << setw(12) << coords_[1] << setw(12) << coords_[2] << endl;

      // Three axis vectors, and the number of gridpoints along each
      // TODO This could be generalized to allow non-perpendicular axis vectors, but probably is not necessary
      cout << setw(5) << ngrid_dim_[0] << setw(12) << inc_size_[0] << setw(12) << 0.0          << setw(12) << 0.0 << endl;
      cout << setw(5) << ngrid_dim_[1] << setw(12) << 0.0          << setw(12) << inc_size_[1] << setw(12) << 0.0 << endl;
      cout << setw(5) << ngrid_dim_[2] << setw(12) << 0.0          << setw(12) << 0.0          << setw(12) << inc_size_[2] << endl;

      // Atomic coordinates
      for (int j=0; j!=geom_->natom(); ++j) {
        cout << setw(5) << geom_->atoms(j)->atom_number() << setw(12) << 0.0;
        for (int k=0; k!=3; ++k)
          cout << setw(12) << geom_->atoms(j)->position(k);
        //cout << setw(12) << geom_->atoms(j)->position(0) << setw(12) << geom_->atoms(j)->position(1) << setw(12) << geom_->atoms(j)->position(2) << endl;
        cout << endl;
      }

      // And now the actual data
      int cnt = 0;
      cout << scientific << setprecision(5);
      for (int j=0; j!=ngrid_; ++j) {
        cout << setw(13) << points_[(norb_+1)*j + i];
        density_sum[i] += points_[(norb_+1)*j + i];
        if (cnt++ % 6 == 5) cout << endl;
        if (j % ngrid_dim_[2] == ngrid_dim_[2]-1) {
          cout << endl;
          cnt = 0;
        }
      }
    }
  } else {
    cout << fixed << setprecision(10);
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
  }


  cout << fixed;
  const double scale = inc_size_[0] * inc_size_[1] * inc_size_[2];
  for (int j=0; j!=norb_; ++j)
    cout << "Sum of all gridpoints for orbital " << orbitals_[j]+1 << " = " << density_sum[j]*scale << ".  Integrated orbital density = " << points_[(norb_+1)*ngrid_+j] << "." << endl;
  cout << "Sum of all gridpoints for total density = " << density_sum.back()*scale << ".  Total integrated density = " << points_.back() << "." << endl;
  cout << fixed << setprecision(5);

}

