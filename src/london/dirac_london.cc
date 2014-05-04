//
// BAGEL - Parallel electron correlation program.
// Filename: dirac_london.cc
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

#include <src/util/constants.h>
#include <src/london/dirac_london.h>
#include <src/london/dfock_london.h>
#include <src/math/zmatrix.h>
#include <src/math/matrix.h>
#include <src/math/diis.h>
#include <src/rel/relreference.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(Dirac_London)

Dirac_London::Dirac_London(const shared_ptr<const PTree> idata, const shared_ptr<const Geometry_London> cgeom,
             const shared_ptr<const Reference> re) : Method(idata, cgeom, re) {
  gaunt_ = idata->get<bool>("gaunt", false);
  breit_ = idata->get<bool>("breit", gaunt_);
  robust_ = idata->get<bool>("robust", false);

  // when computing gradient, we store half-transform integrals
  do_grad_ = idata->get<bool>("gradient", false);

  cgeom_ = cgeom->relativistic(gaunt_);
  common_init(idata);
}


void Dirac_London::common_init(const shared_ptr<const PTree> idata) {
  cout << "  *** Dirac HF ***" << endl << endl;

  // reading input keywords
  max_iter_ = idata->get<int>("maxiter", 100);
  max_iter_ = idata->get<int>("maxiter_scf", max_iter_);
  diis_start_ = idata->get<int>("diis_start", 1);
  thresh_scf_ = idata->get<double>("thresh", 1.0e-8);
  thresh_scf_ = idata->get<double>("thresh_scf", thresh_scf_);
  thresh_overlap_ = idata_->get<double>("thresh_overlap", 1.0e-8);
  ncharge_ = idata->get<int>("charge", 0);
  nele_ = cgeom_->nele()-ncharge_;

  hcore_ = make_shared<const RelHcore>(cgeom_);
  overlap_ = make_shared<const RelOverlap>(cgeom_);
  s12_ = overlap_->tildex(thresh_overlap_);

  nneg_ = s12_->mdim()/2;
  assert(s12_->mdim() % 2 == 0);

  if (breit_ && !gaunt_) throw runtime_error("Breit cannot be turned on if Gaunt is off");
}


void Dirac_London::compute() {
  Timer scftime;
  string indent = "  ";

  shared_ptr<const DistZMatrix> hcore = hcore_->distmatrix();
  shared_ptr<const DistZMatrix> distovl = overlap_->distmatrix();
  shared_ptr<const DistZMatrix> s12 = s12_->distmatrix();
  eig_ = unique_ptr<double[]>(new double[hcore->ndim()]);

  // making initial guess
  shared_ptr<const DistZMatrix> coeff = initial_guess(s12, hcore);
  shared_ptr<const DistZMatrix> aodensity = coeff->form_density_rhf(nele_, nneg_);

  cout << indent << "=== Nuclear Repulsion ===" << endl << indent << endl;
  cout << indent << fixed << setprecision(10) << setw(15) << cgeom_->nuclear_repulsion() << endl << endl;
  cout << indent << "    * DIIS with orbital gradients will be used." << endl << endl;
  scftime.tick_print("SCF startup");
  cout << endl;
  cout << indent << "=== Dirac RHF iteration (" + cgeom_->basisfile() + ", RKB) ===" << endl << indent << endl;

  DIIS<DistZMatrix, ZMatrix> diis(5);

  for (int iter = 0; iter != max_iter_; ++iter) {
    Timer ptime(1);

    auto fock = make_shared<DFock_London>(cgeom_, hcore_, coeff->matrix()->slice(nneg_, nele_+nneg_), gaunt_, breit_, do_grad_, robust_);

// TODO I have a feeling that the code should not need this, but sometimes there are slight errors. still looking on it.
#if 0
    fock->hermite();
#endif
    // distribute
    shared_ptr<const DistZMatrix> distfock = fock->distmatrix();

    // compute energy here
    const complex<double> prod = aodensity->dot_product(*hcore+*distfock); // identical to Tr(D^+ F)
    if (fabs(prod.imag()) > 1.0e-12) {
      stringstream ss; ss << "imaginary part of energy is nonzero!! Perhaps Fock is not Hermite for some reasons " << setprecision(10) << prod.imag();
//    throw runtime_error(ss.str());
      cout << ss.str() << endl;
    }
    energy_ = 0.5*prod.real() + cgeom_->nuclear_repulsion();

    auto error_vector = make_shared<const DistZMatrix>(*distfock**aodensity**distovl - *distovl**aodensity**distfock);
    const double error = error_vector->rms();

    ptime.tick_print("Fock build");
    cout << indent << setw(5) << iter << setw(20) << fixed << setprecision(8) << energy_
         << "   " << setw(17) << error << setw(15) << setprecision(2) << scftime.tick() << endl;

    if (error < thresh_scf_ && iter > 0) {
      cout << indent << endl << indent << "  * SCF iteration converged." << endl << endl;
      // when computing gradient, we store half-transform integrals to avoid recomputation
      if (do_grad_) {
        throw std::logic_error("Gradient integrals have not been implemented with London orbitals yet.");
        half_ = fock->half();
      }
      break;
    } else if (iter == max_iter_-1) {
      cout << indent << endl << indent << "  * Max iteration reached in SCF." << endl << endl;
      throw runtime_error("Max iteration reached in Dirac--Fock SCF");
    }

    if (iter >= diis_start_) {
      distfock = diis.extrapolate(make_pair(distfock, error_vector));
      ptime.tick_print("DIIS");
    }

    DistZMatrix intermediate(*coeff % *distfock * *coeff);
    intermediate.diagonalize(eig_.get());
    coeff = make_shared<DistZMatrix>(*coeff * intermediate);

    aodensity = coeff->form_density_rhf(nele_, nneg_);

  }

  coeff_ = coeff->matrix();
}


//Print non dirac sea eigenvalues
void Dirac_London::print_eig() const {
  const int n = cgeom_->nbasis();
  for (int i = 0*n; i != 4*n; ++i) cout << setprecision(10) << setw(15) << eig_[i] <<  endl;
}


shared_ptr<const Reference> Dirac_London::conv_to_ref() const {
  cout << endl << "CAUTION:  Reference class has not been properly set up for London orbital basis." << endl << endl; // TODO
  return nullptr;
  /*
  // we store only positive state coefficients
  const size_t npos = coeff_->mdim() - nneg_;
  auto out =  make_shared<RelReference>(cgeom_, coeff_, energy_, nneg_, nele_, npos-nele_, gaunt_, breit_);
  vector<double> eig(eig_.get()+nneg_, eig_.get()+nneg_+npos);
  vector<double> eigm(eig_.get(), eig_.get()+nneg_);
  eig.insert(eig.end(), eigm.begin(), eigm.end());
  out->set_eig(eig);
  return out;
  */
}


shared_ptr<const DistZMatrix> Dirac_London::initial_guess(const shared_ptr<const DistZMatrix> s12, const shared_ptr<const DistZMatrix> hcore) const {
  const int n = cgeom_->nbasis();
  unique_ptr<double[]> eig(new double[hcore->ndim()]);

  shared_ptr<const DistZMatrix> coeff;
  if (!ref_) {
    DistZMatrix interm = *s12 % *hcore * *s12;
    interm.diagonalize(eig.get());
    coeff = make_shared<const DistZMatrix>(*s12 * interm);
  } else if (dynamic_pointer_cast<const RelReference>(ref_)) {
    throw logic_error("Not worrying about Reference for gauge-invariant Dirac-Fock yet");
    /*
    auto relref = dynamic_pointer_cast<const RelReference>(ref_);
    shared_ptr<ZMatrix> fock = make_shared<DFock_London>(cgeom_, hcore_, relref->relcoeff()->slice(0, nele_), gaunt_, breit_, *//*store_half*//*false, robust_);
    DistZMatrix interm = *s12 % *fock->distmatrix() * *s12;
    interm.diagonalize(eig.get());
    coeff = make_shared<const DistZMatrix>(*s12 * interm);
    */
  } else if (ref_->coeff()->ndim() == n) {
    throw logic_error("Not worrying about Reference for gauge-invariant Dirac-Fock yet");
    // non-relativistic reference.
    /*
    const int nocc = ref_->nocc();
    shared_ptr<ZMatrix> fock;
    if (nocc*2 == nele_) {
      auto ocoeff = make_shared<ZMatrix>(n*4, 2*nocc);
      ocoeff->add_real_block(1.0, 0,    0, n, nocc, ref_->coeff()->data());
      ocoeff->add_real_block(1.0, n, nocc, n, nocc, ref_->coeff()->data());
      fock = make_shared<DFock_London>(cgeom_, hcore_, ocoeff, gaunt_, breit_, *//*store_half*//*false, robust_);
    } else {
      const int nocca = ref_->noccA();
      const int noccb = ref_->noccB();
      assert(nocca+noccb == nele_);
      auto ocoeff = make_shared<ZMatrix>(n*4, nocca+noccb);
     ocoeff->add_real_block(1.0, 0,     0, n, nocca, ref_->coeffA()->data());
      ocoeff->add_real_block(1.0, n, nocca, n, noccb, ref_->coeffB()->data());
      fock = make_shared<DFock_London>(cgeom_, hcore_, ocoeff, gaunt_, breit_, *//*store_half*//*false, robust_);
    }
    DistZMatrix interm = *s12 % *fock->distmatrix() * *s12;
    interm.diagonalize(eig.get());
    coeff = make_shared<const DistZMatrix>(*s12 * interm);
    */
  } else {
    assert(ref_->coeff()->ndim() == n*4);
    throw logic_error("not yet implemented");
  }
  return coeff;
}
