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
//#include <src/util/math/diis.h>
#include <src/periodic/pscf.h>
#include <src/periodic/poverlap.h>
//#include <src/periodic/pfmm.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(PSCF)

PSCF::PSCF(const shared_ptr<const PTree> idata, const shared_ptr<const Geometry> geom, const shared_ptr<const Reference> re)
  : PSCF_base(idata, geom, re) {
  cout << "  *** Periodic Hartree--Fock ***" << endl << endl;
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
#if 0
  cout << "DO FMM = " << dofmm_ << endl;

  vector<array<double, 3>> primvecs(3);
  primvecs[0] = {{1.0, 0.0, 0.0}};
  primvecs[1] = {{0.0, 1.0, 0.0}};
  primvecs[2] = {{0.0, 0.0, 1.0}};

  vector<shared_ptr<const Atom>> atoms = geom->atoms();
  shared_ptr<const Geometry> newgeom = geom->periodic(atoms);
  auto scell = make_shared<const SimulationCell>(newgeom, primvecs);
  cout << "*** Simulation Cell ***" << endl;
  scell->print();
  PFMM test(scell, true, fmm_lmax_, fmm_ws_, fmm_extent_);
  for (int l = 0; l < test.max_rank(); ++l) {
    for (int m = 0; m <= l; ++m) { // Mlm = -Ml-m
      const int imul = l * l + m + l;
      const double mlm = (test.mlm(imul)).real();
//      if (abs(mlm) > 1e-8)
//        cout << "l = " << l << "  m = " << m << "  mlm = " << setw(20) << setprecision(14) << mlm << endl;
      if (l % 2 == 0 && m % 4 == 0)
        cout << l << "   " << m << "  " << setw(20) << setprecision(14) << mlm << endl;
    }
  }
#endif


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
      (*intermediate)[i]->diagonalize(*eig_[i]);
      (*kcoeff)[i] = make_shared<ZMatrix>(*((*ktildex_)(i)) * *((*intermediate)(i)));
    }
    coeff = kcoeff->ift(lattice_->lattice_vectors(), lattice_->lattice_kvectors());
    coeff_ = make_shared<const PCoeff>(*coeff);
  } else {
    throw runtime_error("Working on it...");
  }

  shared_ptr<const PData> kdensity = kcoeff->form_density_rhf(nocc_);

  const int gamma = lattice_->gamma_point();
  auto olddensity = make_shared<const ZMatrix>(blocksize, blocksize);

  cout << indent << "=== Nuclear Repulsion ===" << endl << indent << endl;
  cout << indent << fixed << setprecision(10) << setw(15) << geom_->nuclear_repulsion() << endl << endl;

  cout << indent << "=== Lattice Nuclear Repulsion ===" << endl << indent << endl;
  cout << indent << fixed << setprecision(10) << setw(15) << lattice_->nuclear_repulsion() << endl << endl;

  cout << indent << "    * DIIS with orbital gradients will be used." << endl << endl;
  pscftime.tick_print("PSCF startup");
  cout << endl;

  cout << indent << "=== PSCF iteration (" + geom_->basisfile() + ") ===" << endl << indent << endl;
  //DIIS<ZMatrix, ZMatrix> diis(diis_size_);

  for (int iter = 0; iter !=  max_iter_; ++iter) {
    auto c = make_shared<PCoeff>(*coeff);
    shared_ptr<const PData> pdensity = c->form_density_rhf(nocc_);
    shared_ptr<const PFock> fock;
    if (!dofmm_) {
      fock = make_shared<const PFock>(lattice_, hcore_, pdensity);
    } else {
      fock = make_shared<const PFock>(lattice_, hcore_, pdensity, fmm_);
    }
    shared_ptr<const PData> kfock = fock->ft(lattice_->lattice_vectors(), lattice_->lattice_kvectors());
    auto fock0 = make_shared<ZMatrix>(*((*kfock)(gamma)));
    double error = 0;
    for (int i =0; i != nkblock; ++i) {
      auto error_vector = make_shared<ZMatrix>(*((*kfock)(i)) * *((*kdensity)(i)) * *((*koverlap_)(i))
                                            - *((*koverlap_)(i)) * *((*kdensity)(i)) * *((*kfock)(i)));
      error += error_vector->rms();
    }
    //auto diis_vector = make_shared<const ZMatrix>(*((*kdensity)(gamma)) - *olddensity);
    //auto diis_vector = make_shared<ZMatrix>(*fock0 * *((*kdensity)(gamma)) * *((*koverlap_)(gamma))
    //                                      - *((*koverlap_)(gamma)) * *((*kdensity)(gamma)) * *fock0);

    complex<double> energy;
    double charge = 0.0;
    for (int i = 0; i != lattice_->num_lattice_vectors(); ++i) {
      energy += 0.5 * ((*((*fock)(i)) + *((*hcore_)(i))) * *((*pdensity)(i))).trace();
      //assert(energy.imag() < 1e-8);
      if (energy.imag() >= 1e-8)
        cout << "*** Warning: energy.imag() >= 1e-8 " << setprecision(9) << energy.imag() << endl;
      for (int j = 0; j != blocksize; ++j)
        for (int k = 0; k != blocksize; ++k)
          charge += ((*overlap_)(lattice_->num_lattice_vectors() - i - 1)->element(j, k) * (*pdensity)(i)->element(j, k)).real();
    }
    cout << "*** sr = " << setprecision(1) << charge << endl;
    for (int i = 0; i != nkblock; ++i) {
      double charge = 0.0;
      for (int j = 0; j != blocksize; ++j)
        for (int k = 0; k != blocksize; ++k)
          charge += ((*koverlap_)(i)->element(j, k) * (*kdensity)(i)->element(j, k)).real();
      cout << "SP = " << setprecision(1) << charge << endl;
    }

    cout << "   #ele = " << lattice_->nele();
    //cout << "SP = " << setprecision(1) << charge << "       #ele = " << lattice_->nele();
    energy_ = energy.real() + lattice_->nuclear_repulsion() + fock->correction();
    cout << indent << setw(5) << iter << setw(30) << fixed << setprecision(8) << energy_ << "   "
                                      << setw(17) << error << setw(15) << setprecision(2) << pscftime.tick();
    if (abs(energy.imag()) > 1e-12) {
      cout << "  *** Warning *** Im(E) = " << setw(15) << fixed << setprecision(12) << energy.imag() << endl;
    } else {
      cout << endl;
    }

    if (error < thresh_scf_) {
      cout << indent << endl << indent << "  * PSCF iteration converged." << endl << endl;
      cout << indent << endl << indent << "    Eigenvalues with nbasis = " << blocksize << endl << endl;
#if 1
      for (int i = 0; i != nkblock; ++i) {
        for (int n = 0; n != blocksize; ++n) cout << setprecision(9) << (*eig_[i])(n) << "    ";
        cout << endl;
      }
#endif
      break;
    } else if (iter == max_iter_-1) {
      cout << indent << endl << indent << "  * Max iteration reached in PSCF." << endl << endl;
      break;
    }

    auto intermediate = make_shared<PData>(blocksize, nkblock);
//    if (iter >= diis_start_)
//      fock0 = diis.extrapolate({(*kfock)(gamma), diis_vector});

    for (int i = 0; i != nkblock; ++i) {
      //if (iter < diis_start_) *fock0 = *(*kfock)(i);
      *fock0 = *(*kfock)(i);
      ZMatrix kblock = *((*ktildex_)(i)) % *fock0 * *((*ktildex_)(i));
      //cout << i << "   " << setprecision(15) << (kblock - *(kblock.transpose_conjg())).norm()/kblock.size() << endl;
      (*intermediate)[i] = make_shared<ZMatrix>(kblock);
      (*intermediate)[i]->diagonalize(*eig_[i]);
      (*kcoeff)[i] = make_shared<ZMatrix>(*((*ktildex_)(i)) * *((*intermediate)(i)));
    }

    olddensity = make_shared<const ZMatrix>(*(*kdensity)(gamma));
    kdensity = kcoeff->form_density_rhf(nocc_);
    coeff = kcoeff->ift(lattice_->lattice_vectors(), lattice_->lattice_kvectors());
  }

  coeff_ = make_shared<const PCoeff>(*coeff);

}
