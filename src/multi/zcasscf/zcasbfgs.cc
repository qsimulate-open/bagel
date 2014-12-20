//
// BAGEL - Parallel electron correlation program.
// Filename: zcasbfgs.cc
// Copyright (C) 2014 Jefferson E. Bates
//
// Author: Toru Shiozaki <jefferson.bates@northwestern.edu>
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

#include <src/multi/zcasscf/zcasbfgs.h>
#include <src/multi/zcasscf/zqvec.h>

using namespace std;
using namespace bagel;


complex<double> ZCASBFGS::find_level_shift(shared_ptr<const ZRotFile> rotmat) const {
  /* returns the smallest Hessian value that is not below -mc^2 to be used as a level shift
     This particular choice of level shift parameter ensures that the initial diagonal guess has Np negative values
     where Np is the number of positronic orbitals */
  double csq = c__*c__;
  complex<double> l0 = rotmat->data(0);

  for (int j = 1; j != rotmat->size(); ++j) {
    if (l0.real() > rotmat->data(j).real() && csq + rotmat->data(j).real() > 0)
      l0 = rotmat->data(j);
  }
  complex<double> level_shift;
  if (real(l0) < 0) {
    double scale = idata_->get<double>("scalefac", 1.05);
    level_shift = l0 * scale;
    cout << " " << endl;
    cout << setprecision(8) << "Smallest Hessian Element (excluding positrons) = " << l0 << endl;
    cout << setprecision(2) << "Scaling Factor                                 = " << scale << endl;
    cout << setprecision(8) << "Level Shift                                    = " << level_shift << endl << endl;
  } else {
    cout << " level shift is not negative " << endl;
    level_shift = complex<double> (0.0,0.0);
  }

  return level_shift;
}


tuple<shared_ptr<ZRotFile>, shared_ptr<ZRotFile>> ZCASBFGS::optimize_subspace_rotations(vector<double> energy, shared_ptr<const ZRotFile> grad, shared_ptr<const ZRotFile> rot, shared_ptr<SRBFGS<ZRotFile>> srbfgs, bool optimize_electrons) {
  // function to optimize only a single subspace of orbitals neglecting the coupling to the other
  const int nvirtnr = nvirt_ - nneg_/2;

  // copy the inputs
  shared_ptr<ZRotFile> newgrad;
  shared_ptr<ZRotFile> newrot = rot->copy();
  if (optimize_electrons) {
    newgrad = copy_electronic_rotations(grad);
  } else {
    newgrad = copy_positronic_rotations(grad);
  }

  const bool tight = idata_->get<bool>("tight", false);
  const int limmem = idata_->get<int>("limited_memory", 0);
  const int hebden = idata_->get<bool>("hebden", false);
  auto reset = srbfgs->check_step(energy, newgrad, newrot, tight, limmem);
  if (reset) {
    cout << " STEP DOES NOT MEET PROPER CRITERIA " << endl;
    // TODO : implement step rejection or another clever alternative ; if step rejection is used energy array will need to be output since it will require modification
  }

  shared_ptr<ZRotFile> a;
  if (optimize_electrons) {
    if (!hebden)
      a = srbfgs->more_sorensen_extrapolate(newgrad, newrot);
    else
      a = srbfgs->extrapolate(newgrad, newrot);
  } else {
    // positronic optimization results in a negative level shift so use Hebden method
    a = srbfgs->extrapolate(newgrad, newrot);
  }

  if (optimize_electrons) {
    kramers_adapt(a, nclosed_, nact_, nvirtnr);
  } else {
    kramers_adapt(a, nclosed_, nact_, nneg_/2);
  }
  cout << setprecision(6) << " Subspace gradient rms  = " << newgrad->rms() << endl;

  return make_tuple(a, newgrad);
}


shared_ptr<ZMatrix> ZCASBFGS::compute_unitary_rotation(vector<double>& subspace_energy, shared_ptr<SRBFGS<ZRotFile>> subspace_bfgs, shared_ptr<ZMatrix> displacement_history, const int nvirt_subspace, shared_ptr<const ZMatrix> cfockao, shared_ptr<ZRotFile>& grad, const bool optimize_electrons) {
    // grad is modified by optimize_subspace_rotations and hence has to be passed by reference
    // be sure nneg_/2 is being passed for nvirt_subspace when e-p rotations are active

    // get energy
    if (nact_) {
      // use state averaged energy to update trust radius
      assert(fci_->energy().size() > 0);
      double sa_en = 0.0;
      for (auto& i : fci_->energy())
        sa_en += i;
      sa_en /= double((fci_->energy()).size());
      subspace_energy.push_back(sa_en);
      if (energy_.size() > 0) prev_energy_ = energy_;
      energy_ = fci_->energy();
    } else {
      assert(nstate_ == 1 && energy_.size() == 1);
      const int subspace_size = subspace_energy.size(); // location fixed by initial size of array since c++ is 0 based indexing
      subspace_energy.resize(subspace_size+1); // increase size of array by 1
      subspace_energy[subspace_size] = geom_->nuclear_repulsion();
      auto mo = make_shared<ZMatrix>(*coeff_ % (*cfockao+*hcore_) * *coeff_);
      for (int i = 0; i != nclosed_*2; ++i)
        subspace_energy[subspace_size] += 0.5*mo->element(i,i).real();
      energy_[0] = subspace_energy[subspace_size];
    }

    // compute rotation via extrapolation
    shared_ptr<ZRotFile> xlog;
    shared_ptr<ZRotFile> subspace_rot;
    Timer more_sorensen_timer(0);
    cout << " " << endl;
    cout << " -------  Step Restricted BFGS Extrapolation  ------- " << endl;
    if (optimize_electrons)
      cout << " --- Optimizing electrons --- " << endl;
    else
      cout << " --- Optimizing positrons --- " << endl;
    // grad is altered during optimization of subspace rotations
    xlog    = make_shared<ZRotFile>(displacement_history->log(4), nclosed_*2, nact_*2, nvirt_subspace*2);
    tie(subspace_rot, grad) = optimize_subspace_rotations(subspace_energy, grad, xlog, subspace_bfgs, optimize_electrons);
    cout << " ---------------------------------------------------- " << endl << endl;
    more_sorensen_timer.tick_print("More-Sorensen/Hebden extrapolation");
    kramers_adapt(subspace_rot, nclosed_, nact_, nvirt_subspace);

    // Rotate orbitals
    shared_ptr<ZMatrix> amat = subspace_rot->unpack<ZMatrix>();

    // multiply -1 from the formula taken care of in extrap. multiply -i to make amat hermite (will be compensated)
    *amat *= 1.0 * complex<double>(0.0, -1.0);

    // restore the matrix from RotFile
    VectorB teig(amat->ndim());
    amat->diagonalize(teig);
    auto amat_sav = amat->copy();
    for (int i = 0; i != amat->ndim(); ++i) {
      complex<double> ex = exp(complex<double>(0.0, teig(i)));
      for_each(amat->element_ptr(0,i), amat->element_ptr(0,i+1), [&ex](complex<double>& a) { a *= ex; });
    }
    auto expa = make_shared<ZMatrix>(*amat ^ *amat_sav);
    // enforce time-reversal symmetry and unitarity
    kramers_adapt(expa, nvirt_subspace);
    expa->purify_unitary();

    return expa;
}
