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
#include <src/math/diis.h>
#include <src/periodic/pscf.h>
#include <src/periodic/poverlap.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(PSCF)

PSCF::PSCF(const shared_ptr<const PTree> idata, const shared_ptr<const Geometry> geom, const shared_ptr<const Reference> re)
  : PSCF_base(idata, geom, re), dodf_(idata->get<bool>("df",true)) {
  cout << "  *** Periodic Hartree--Fock ***" << endl << endl;
  if (!dodf_)
    throw runtime_error("Periodic SCF only works with density fitting!");

  if (nocc_ != noccB_)
    throw runtime_error("PSCF only works for closed shell systems.");

  cout << indent << "=== V(unit cell) in direct space ===" << endl << indent << endl;
  cout << indent << fixed << setprecision(10) << setw(15) << lattice_->volume() << endl << endl;

  lattice_->print_primitive_vectors();
  lattice_->print_lattice_coordinates();
  lattice_->print_primitive_kvectors();

//lattice_->print_lattice_vectors();

/********************** TEST FOURIER TRANSFORM ************************/
#if 0
  POverlap overlap(lattice_);
//overlap.print("Overlap", 100);

  shared_ptr<const PData> ft_overlap
  = overlap.ft(lattice_->lattice_vectors(), lattice_->lattice_kvectors());
  ft_overlap->print("FT-Overlap", 100);
  cout << "****** Check Hermitianity ******" << endl;
  for (int i = 0; i != ft_overlap->nblock(); ++i)
    cout << setprecision(15) << (*(*ft_overlap)(i) - *((*ft_overlap)(i)->transpose_conjg())).norm()/(*ft_overlap)(i)->size() << endl;

  shared_ptr<const PData> ift_overlap = ft_overlap->ift(lattice_->lattice_vectors(), lattice_->lattice_kvectors());
//ift_overlap->print("IFT-Overlap", 100);

  shared_ptr<const PData> ft2_overlap = ift_overlap->ft(lattice_->lattice_vectors(), lattice_->lattice_kvectors());
  ft2_overlap->print("I(IFT)-Overlap", 100);
#endif
/**********************************************************************/

}

void PSCF::compute() {

  Timer pscftime;

  shared_ptr<const PData> khcore   = hcore_->ft(lattice_->lattice_vectors(), lattice_->lattice_kvectors());

  const int nkblock = lattice_->num_lattice_kvectors();
  const int blocksize = hcore_->blocksize();
  auto kcoeff = make_shared<PCoeff>(blocksize, nkblock);
  shared_ptr<const PData> coeff;

  if (coeff_ == nullptr) {
    shared_ptr<const PData> kfock = khcore;
    auto intermediate = make_shared<PData>(blocksize, nkblock);
    for (int i = 0; i != nkblock; ++i) {
      const ZMatrix kblock = *((*ktildex_)(i)) % *((*kfock)(i)) * *((*ktildex_)(i));
      (*intermediate)[i] = make_shared<ZMatrix>(kblock);
      (*kcoeff)[i] = make_shared<ZMatrix>(*((*ktildex_)(i)) * *((*intermediate)(i)));
      coeff = kcoeff->ift(lattice_->lattice_vectors(), lattice_->lattice_kvectors());
      coeff_ = make_shared<const PCoeff>(*coeff);
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
  DIIS<ZMatrix, ZMatrix> diis(diis_size_);

  const int gamma = lattice_->gamma_point();

  for (int iter = 0; iter !=  max_iter_; ++iter) {
    auto c = make_shared<PCoeff>(*coeff);
    shared_ptr<const PData> pdensity = c->form_density_rhf(nocc_);
    auto fock = make_shared<const PFock>(lattice_, hcore_, pdensity);
    shared_ptr<const PData> kfock = fock->ft(lattice_->lattice_vectors(), lattice_->lattice_kvectors());
    vector<complex<double>> energy(nkblock);
    double max_imag = 0.0;
    auto fock0 = make_shared<ZMatrix>(*((*kfock)(gamma)));
    auto error_vector = make_shared<ZMatrix>(*fock0 * *((*aodensity)(gamma)) * *((*koverlap_)(gamma))
                                          - *((*koverlap_)(gamma)) * *((*aodensity)(gamma)) * *fock0);
    const double error = error_vector->rms();

    cout << indent << setw(5) << iter;
    for (int i = 0; i != nkblock; ++i) {
      energy[i] = 0.5 * ((*((*khcore)(i)) + *((*kfock)(i)) * *((*aodensity)(i))).trace()) + lattice_->nuclear_repulsion();
      assert(energy[i].imag() < 1e-8);
      if (abs(energy[i].imag()) > max_imag) max_imag = abs(energy[i].imag());
    }
    cout << setw(17) << setprecision(8) << error << setw(15) << setprecision(2) << pscftime.tick();
    if (abs(max_imag) > 1e-12) {
      cout << "  *** Warning *** Im(E) = " << setw(15) << fixed << setprecision(12) << max_imag << endl;
    } else {
      cout << endl;
    }

    if (error < thresh_scf_) {
      cout << indent << endl << indent << "  * PSCF iteration converged." << endl << endl;
      cout << indent << endl << indent << "    Eigenvalues with nbasis = " << blocksize << endl << endl;
      for (int i = 0; i != nkblock; ++i) {
        for (int n = 0; n != blocksize; ++n) cout << setprecision(9) << (*eig_[i])(n) << "    ";
        cout << endl;
      }
      break;
    } else if (iter == max_iter_-1) {
      cout << indent << endl << indent << "  * Max iteration reached in SOSCF." << endl << endl;
      break;
    }

    auto intermediate = make_shared<PData>(blocksize, nkblock);
    for (int i = 0; i != nkblock; ++i) {
      if (iter >= diis_start_) {
        fock0 = diis.extrapolate({(*kfock)(i), error_vector});
      } else {
        fock0 = make_shared<ZMatrix>(*((*kfock)(i)));
      }
      ZMatrix kblock = *((*ktildex_)(i)) % *fock0 * *((*ktildex_)(i));
      //cout << i << "   " << setprecision(15) << (kblock - *(kblock.transpose_conjg())).norm()/kblock.size() << endl;
      kblock.fill_upper_conjg(); // very dangerous!
      (*intermediate)[i] = make_shared<ZMatrix>(kblock);
      (*intermediate)[i]->diagonalize(*eig_[i]);
      (*kcoeff)[i] = make_shared<ZMatrix>(*((*ktildex_)(i)) * *((*intermediate)(i)));
    }

    aodensity = kcoeff->form_density_rhf(nocc_);
    coeff = kcoeff->ift(lattice_->lattice_vectors(), lattice_->lattice_kvectors());
  }

  coeff_ = make_shared<const PCoeff>(*coeff);

}
