//
// BAGEL - Parallel electron correlation program.
// Filename: scf_london.cc
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

#include <src/london/scf_london.h>
#include <src/prop/multipole.h>
#include <src/scf/atomicdensities.h>

using namespace bagel;
using namespace std;

BOOST_CLASS_EXPORT_IMPLEMENT(SCF_London)


SCF_London::SCF_London(const shared_ptr<const PTree> idata, const shared_ptr<const Geometry_London> geom, const shared_ptr<const Reference> re)
 : SCF_base_London(idata, geom, re, !idata->get<bool>("df",true)), dodf_(idata->get<bool>("df",true)), restarted_(false) {

  cout << indent << "*** RHF ***" << endl << endl;
  if (nocc_ != noccB_) throw runtime_error("Closed shell SCF was called with nact != 0");

  // For the moment, I can't be bothered to test the level shifting apparatus for UHF and ROHF cases.
  // In the future, this should probably be moved to SCF_base and designed to work properly there
  lshift_ = idata->get<double>("levelshift", 0.0);
  if (lshift_ != 0.0) {
    cout << "  level shift : " << setprecision(3) << lshift_ << endl << endl;
    levelshift_ = make_shared<ShiftVirtual<DistZMatrix>>(nocc_, lshift_);
  }
}


void SCF_London::compute() {
  Timer scftime;

  shared_ptr<const ZMatrix> previous_fock = hcore_;
  shared_ptr<const ZMatrix> aodensity_;

  shared_ptr<const DistZMatrix> tildex = tildex_->distmatrix();
  shared_ptr<const DistZMatrix> tildextc = tildex_->transpose_conjg()->distmatrix();
  shared_ptr<const DistZMatrix> hcore = hcore_->distmatrix();
  shared_ptr<const DistZMatrix> overlap = overlap_->distmatrix();
  shared_ptr<const DistZMatrix> coeff;
  shared_ptr<const DistZMatrix> aodensity;

  if (!restarted_) {
    if (coeff_ == nullptr) {
      shared_ptr<const DistZMatrix> fock = hcore;
      if (dodf_ && cgeom_->spherical()) {
        auto aden = make_shared<const AtomicDensities>(cgeom_);
        auto zaden = std::make_shared<ZMatrix>(*aden, 1.0);
        auto focka = make_shared<const Fock_London<1>>(cgeom_, hcore_, zaden, schwarz_);
        //focka->print("FockA", 20);
        fock = focka->distmatrix();
      }
      // TODO Do you need the transpose conjugate of *tildex here?
      DistZMatrix intermediate = *tildex % *fock * *tildex;
      intermediate.diagonalize(eig());
      coeff = make_shared<const DistZMatrix>(*tildex * intermediate);
      //coeff->print("coeff, as first assigned", 20);
    } else {
      shared_ptr<const ZMatrix> focka;
      if (!dodf_) {
        //throw runtime_error("Only worrying about density-fitted HF for now");
        shared_ptr<const ZMatrix> halfaodensity = coeff->form_density_rhf(nocc_);
        aodensity_ = make_shared<ZMatrix>(*halfaodensity * 2.0);
        focka = make_shared<const Fock_London<0>>(cgeom_, hcore_, aodensity_, schwarz_);
        } else {
        focka = make_shared<const Fock_London<1>>(cgeom_, hcore_, nullptr, coeff_->slice(0, nocc_), do_grad_, true/*rhf*/);
        //focka->print("FockA", 20);
      }
      DistZMatrix intermediate = *tildex % *focka->distmatrix() * *tildex;
      intermediate.diagonalize(eig());
      coeff = make_shared<const DistZMatrix>(*tildex * intermediate);
    }
    coeff_ = make_shared<const ZCoeff>(*coeff->matrix());

    cout << indent << "=== Nuclear Repulsion ===" << endl << indent << endl;
    cout << indent << fixed << setprecision(10) << setw(15) << cgeom_->nuclear_repulsion() << endl << endl;
    cout << indent << "    * DIIS with orbital gradients will be used." << endl << endl;
    scftime.tick_print("SCF startup");
    cout << endl;
    cout << indent << "=== RHF iteration (" + cgeom_->basisfile() + ") ===" << endl << indent << endl;

    diis_ = make_shared<DIIS<DistZMatrix,ZMatrix>>(diis_size_);
  } else {
    coeff = coeff_->distmatrix();
  }

  if (!dodf_) {
    //throw runtime_error("Only worrying about density-fitted HF for now");
    shared_ptr<const ZMatrix> halfaodensity = coeff->form_density_rhf(nocc_);
    aodensity_ = make_shared<ZMatrix>(*halfaodensity * 2.0);
    aodensity = aodensity_->distmatrix();
  } else {
    //coeff->print("coeff", 20);
    shared_ptr<const ZMatrix> halfaodensity = coeff->form_density_rhf(nocc_);
    aodensity = make_shared<DistZMatrix>(*halfaodensity * 2.0);
  }

  // starting SCF iteration
  shared_ptr<const ZMatrix> densitychange = aodensity_;

  for (int iter = 0; iter != max_iter_; ++iter) {
    //cout << endl << endl << endl << "STARTING ITERATION " << iter << endl << endl << endl;
    Timer pdebug(1);

    if (restart_) {
      stringstream ss; ss << "scf_" << iter;
      OArchive archive(ss.str());
      archive << static_cast<Method*>(this);
    }

    if (!dodf_) {
      //throw runtime_error("Only worrying about density-fitted HF for now");
      previous_fock = make_shared<Fock_London<0>>(cgeom_, previous_fock, densitychange, schwarz_);
      mpi__->broadcast(const_pointer_cast<ZMatrix>(previous_fock)->data(), previous_fock->size(), 0);
    } else {
      //coeff_->slice(0, nocc_)->print("slice of coeff_ being used",20);
      previous_fock = make_shared<Fock_London<1>>(cgeom_, hcore_, nullptr, coeff_->slice(0, nocc_), do_grad_, true/*rhf*/);
    }
    shared_ptr<const DistZMatrix> fock = previous_fock->distmatrix();

    // TODO Is this the best place to convert to real?
    const complex<double> zenergy  = 0.5*aodensity->dot_product(*hcore+*fock) + cgeom_->nuclear_repulsion();
    //cout << "Energy of iteration " << iter << " = " << zenergy << endl;
    energy_  = real(zenergy);
    //aodensity->print("aodensity", 20);
    //hcore->print("hcore", 20);
    //fock->print("fock", 20);
    if (abs(imag(zenergy))>1.0e-8) throw logic_error("Energy should be real!");

    pdebug.tick_print("Fock build");

    auto error_vector = make_shared<const DistZMatrix>(*fock**aodensity**overlap - *overlap**aodensity**fock);
    const double error = error_vector->rms();
    //error_vector->print("Error Vector", 20);
    //cout << endl << endl << "ERROR = " << error_vector->rms() << endl << endl;

    cout << indent << setw(5) << iter << setw(20) << fixed << setprecision(8) << energy_ << "   "
                                      << setw(17) << error << setw(15) << setprecision(2) << scftime.tick() << endl;

    if (error < thresh_scf_) {
      cout << indent << endl << indent << "  * SCF iteration converged." << endl << endl;
      if (do_grad_) half_ = dynamic_pointer_cast<const Fock_London<1>>(previous_fock)->half();
      break;
    } else if (iter == max_iter_-1) {
      cout << indent << endl << indent << "  * Max iteration reached in SCF." << endl << endl;
      break;
    }

    if (diis_ || iter >= diis_start_) {
      fock = diis_->extrapolate(make_pair(fock, error_vector));
      pdebug.tick_print("DIIS");
    }

    DistZMatrix intermediate(*coeff % *fock * *coeff);
    //fock->print("Fock Matrix by extrapolation from diis_", 20);
    //intermediate.print("intermediate using coeff & extrapolated Fock", 20);

    if (levelshift_)
      levelshift_->shift(intermediate);

    intermediate.diagonalize(eig());

    //for (int i=0; i!=eig_.size(); i++) cout << "eigenvalue " << i << " = " << eig_[i] << endl;

    pdebug.tick_print("Diag");

    //intermediate.print("intermediate after diagonalization", 20);
    coeff = make_shared<const DistZMatrix>(*coeff * intermediate);
    coeff_ = make_shared<const ZCoeff>(*coeff->matrix());


    if (!dodf_) {
      //throw runtime_error("Only worrying about density-fitted HF for now");
      shared_ptr<const ZMatrix> halfaodensity = coeff->form_density_rhf(nocc_);
      shared_ptr<const ZMatrix> new_density = make_shared<ZMatrix>(*halfaodensity * 2.0);
      densitychange = make_shared<ZMatrix>(*new_density - *aodensity_);
      aodensity_ = new_density;
      aodensity = aodensity_->distmatrix();
    } else {
      //aodensity = coeff->form_density_rhf(nocc_);
      //coeff->print("coeff from this iteration", 20);
      shared_ptr<const DistZMatrix> halfaodensity = coeff->form_density_rhf(nocc_);
      aodensity = make_shared<DistZMatrix>(*halfaodensity * 2.0);
      //aodensity->print("the new aodensity", 20);
    }
    pdebug.tick_print("Post process");
  }


  //TODO If we're going to compute dipole moments in SCF_London, this is the place to do it
#if 0
  if (!cgeom_->external()) {
    if (dodf_) aodensity_ = aodensity->matrix();
    Multipole mu(cgeom_, aodensity_, multipole_print_);
    mu.compute();
  }
#endif
}


shared_ptr<const Reference> SCF_London::conv_to_ref() const {
  cout << endl << "CAUTION:  Reference class has not been properly set up for London orbital basis." << endl; // TODO
  //auto out = make_shared<Reference>(geom_, coeff(), nocc(), 0, coeff_->mdim()-nocc(), energy());
  //out->set_eig(eig_);
  //return out;
  return nullptr;
}
