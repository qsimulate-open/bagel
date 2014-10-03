//
// BAGEL - Parallel electron correlation program.
// Filename: pscf.cc
// Copyright (C) 2014 Toru Shiozaki
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


#include <iomanip>
#include <algorithm>
#include <src/util/timer.h>
#include <src/periodic/pscf.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(PSCF)

PSCF::PSCF(const shared_ptr<const PTree> idata, const shared_ptr<const Geometry> geom, const shared_ptr<const Reference> re)
  : PSCF_base(idata, geom, re), dodf_(idata->get<bool>("df",true)) {
  cout << "  *** Periodic Hartree--Fock ***" << endl << endl;
  if (dodf_)
    throw runtime_error("Periodic code does not work with density fitting yet!");

  if (nocc_ != noccB_)
    throw runtime_error("PSCF only works for closed shell systems.");

  cout << indent << "=== V(unit cell) in direct space ===" << endl << indent << endl;
  cout << indent << fixed << setprecision(10) << setw(15) << lattice_->volume() << endl << endl;

  lattice_->print_primitive_vectors();
  lattice_->print_lattice_coordinates();
  lattice_->print_primitive_kvectors();

}

void PSCF::compute() {

  Timer pscftime;

  shared_ptr<const PData> koverlap = overlap_->ft(lattice_->lattice_vectors(), lattice_->lattice_kvectors());
  shared_ptr<const PData> ktildex  = tildex_->ft(lattice_->lattice_vectors(), lattice_->lattice_kvectors());
  shared_ptr<PData> kcoeff;

  const int nkblock = lattice_->num_lattice_kvectors();
  const int blocksize = overlap_->blocksize();

  if (coeff_ == nullptr) {
    shared_ptr<const PData> kfock = hcore_->ft(lattice_->lattice_vectors(), lattice_->lattice_kvectors());
    auto intermediate = make_shared<PData>(blocksize, nkblock);
    for (int i = 0; i != nkblock; ++i) {
      const ZMatrix kblock = *((*ktildex)(i)) % *((*kfock)(i)) * *((*ktildex)(i));
      (*intermediate)[i] = make_shared<ZMatrix>(kblock);
    }
    kcoeff = make_shared<PData>(blocksize, nkblock);
    for (int i = 0; i != nkblock; ++i) {
      const ZMatrix kblock = *((*ktildex)(i)) * *((*intermediate)(i));
      (*kcoeff)[i] = make_shared<ZMatrix>(kblock);
    }
  } else {
    throw runtime_error("Working on it...");
  }

  shared_ptr<const PData> aodensity = kcoeff->form_density_rhf(nocc_);

  cout << indent << "=== Nuclear Repulsion ===" << endl << indent << endl;
  cout << indent << fixed << setprecision(10) << setw(15) << geom_->nuclear_repulsion() << endl << endl;

  cout << indent << "=== Lattice Nuclear Repulsion ===" << endl << indent << endl;
  cout << indent << fixed << setprecision(10) << setw(15) << lattice_->nuclear_repulsion() << endl << endl;

  cout << indent << "    * DIIS with orbital gradients will be used." << endl << endl;
  pscftime.tick_print("PSCF startup");
  cout << endl;

  cout << indent << "=== PSCF iteration (" + geom_->basisfile() + ") ===" << endl << indent << endl;

  for (int iter = 0; iter !=  max_iter_; ++iter) {

  }


}
