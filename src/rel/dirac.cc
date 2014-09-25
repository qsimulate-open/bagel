//
// BAGEL - Parallel electron correlation program.
// Filename: dirac.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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
#include <src/rel/dirac.h>
#include <src/rel/dfock.h>
#include <src/rel/relhcore.h>
#include <src/london/relhcore_london.h>
#include <src/rel/reloverlap.h>
#include <src/london/reloverlap_london.h>
#include <src/math/zmatrix.h>
#include <src/math/matrix.h>
#include <src/math/diis.h>
#include <src/rel/relreference.h>
#include <src/prop/momentum_london.h>
#include <src/prop/momentum_point.h>

using namespace std;
using namespace bagel;

Dirac::Dirac(const shared_ptr<const PTree> idata, const shared_ptr<const Geometry> geom,
             const shared_ptr<const Reference> re) : Method(idata, geom, re) {
  gaunt_ = idata->get<bool>("gaunt", false);
  breit_ = idata->get<bool>("breit", gaunt_);
  robust_ = idata->get<bool>("robust", false);

  // when computing gradient, we store half-transform integrals
  do_grad_ = idata->get<bool>("gradient", false);
  if (do_grad_ && geom_->magnetism()) throw runtime_error("Gradient integrals have not been implemented for a GIAO basis.");

  geom_ = geom->relativistic(gaunt_);
  common_init(idata);
}


void Dirac::common_init(const shared_ptr<const PTree> idata) {
  Timer init_time;
  cout << "  *** Dirac HF ***" << endl << endl;

  // reading input keywords
  max_iter_ = idata->get<int>("maxiter", 100);
  max_iter_ = idata->get<int>("maxiter_scf", max_iter_);
  diis_start_ = idata->get<int>("diis_start", 1);
  thresh_scf_ = idata->get<double>("thresh", 1.0e-8);
  thresh_scf_ = idata->get<double>("thresh_scf", thresh_scf_);
  thresh_overlap_ = idata_->get<double>("thresh_overlap", 1.0e-8);
  ncharge_ = idata->get<int>("charge", 0);
  nele_ = geom_->nele()-ncharge_;

  if (!geom_->magnetism()) {
    hcore_ = make_shared<const RelHcore>(geom_);
    overlap_ = make_shared<const RelOverlap>(geom_);
  } else {
    hcore_ = make_shared<const RelHcore_London>(geom_);
    overlap_ = make_shared<const RelOverlap_London>(geom_);
  }
  init_time.tick_print("1-electron integrals");

  s12_ = overlap_->tildex(thresh_overlap_);

  nneg_ = s12_->mdim()/2;
  assert(s12_->mdim() % 2 == 0);

  if (breit_ && !gaunt_) throw runtime_error("Breit cannot be turned on if Gaunt is off");
}


void Dirac::compute() {
  Timer scftime;
  string indent = "  ";

  shared_ptr<const DistZMatrix> hcore = hcore_->distmatrix();
  shared_ptr<const DistZMatrix> distovl = overlap_->distmatrix();
  shared_ptr<const DistZMatrix> s12 = s12_->distmatrix();
  eig_ = VectorB(hcore->ndim());

  // making initial guess
  shared_ptr<const DistZMatrix> coeff = initial_guess(s12, hcore);
  shared_ptr<const DistZMatrix> aodensity = coeff->form_density_rhf(nele_, nneg_);

  cout << indent << "=== Nuclear Repulsion ===" << endl << indent << endl;
  cout << indent << fixed << setprecision(10) << setw(15) << geom_->nuclear_repulsion() << endl << endl;
  cout << indent << "    * DIIS with orbital gradients will be used." << endl << endl;
  scftime.tick_print("SCF startup");
  cout << endl;
  cout << indent << "=== Dirac RHF iteration (" + geom_->basisfile() + ", " << (geom_->magnetism() ? "RMB" : "RKB") << ") ===" << endl << indent << endl;

  DIIS<DistZMatrix, ZMatrix> diis(5);

  for (int iter = 0; iter != max_iter_; ++iter) {
    Timer ptime(1);

    auto fock = make_shared<DFock>(geom_, hcore_, coeff->matrix()->slice_copy(nneg_, nele_+nneg_), gaunt_, breit_, do_grad_, robust_);

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
    energy_ = 0.5*prod.real() + geom_->nuclear_repulsion();

    auto error_vector = make_shared<const DistZMatrix>(*distfock**aodensity**distovl - *distovl**aodensity**distfock);
    const double error = error_vector->rms();

    ptime.tick_print("Fock build");
    cout << indent << setw(5) << iter << setw(20) << fixed << setprecision(12) << energy_
         << "   " << setw(17) << error << setw(15) << setprecision(2) << scftime.tick() << endl;

    if (error < thresh_scf_ && iter > 0) {
      cout << indent << endl << indent << "  * SCF iteration converged." << endl << endl;
      // when computing gradient, we store half-transform integrals to avoid recomputation
      if (do_grad_) half_ = fock->half();
      break;
    } else if (iter == max_iter_-1) {
      cout << indent << endl << indent << "  * Max iteration reached in SCF." << endl << endl;
      throw runtime_error("Max iteration reached in Dirac--Fock SCF");
    }

    if (iter >= diis_start_) {
      distfock = diis.extrapolate({distfock, error_vector});
      ptime.tick_print("DIIS");
    }

    DistZMatrix intermediate(*coeff % *distfock * *coeff);
    intermediate.diagonalize(eig_);
    coeff = make_shared<DistZMatrix>(*coeff * intermediate);

    aodensity = coeff->form_density_rhf(nele_, nneg_);

  }

  coeff_ = coeff->matrix();


  // Trying to compute charge current
  if (geom_->magnetism()) {
#if 1
    // Integrated over all space
    const int n = geom_->nbasis();
    auto mom = make_shared<Momentum_London>(geom_);
    array<shared_ptr<ZMatrix>,3> pi = mom->compute();
    auto current_x = make_shared<ZMatrix>(4*n, 4*n);
    auto current_y = make_shared<ZMatrix>(4*n, 4*n);
    auto current_z = make_shared<ZMatrix>(4*n, 4*n);

    const complex<double> re(-0.5,  0.0);
    const complex<double> im( 0.0, -0.5);

    current_x->add_block( re, 0*n, 2*n, n, n, *pi[0]);
    current_x->add_block( im, 0*n, 2*n, n, n, *pi[1]);
    current_x->add_block(-re, 0*n, 3*n, n, n, *pi[2]);
    current_x->add_block( re, 1*n, 2*n, n, n, *pi[2]);
    current_x->add_block( re, 1*n, 3*n, n, n, *pi[0]);
    current_x->add_block(-im, 1*n, 3*n, n, n, *pi[1]);

    current_x->add_block( re, 2*n, 0*n, n, n, *pi[0]);
    current_x->add_block(-im, 2*n, 0*n, n, n, *pi[1]);
    current_x->add_block( re, 2*n, 1*n, n, n, *pi[2]);
    current_x->add_block(-re, 3*n, 0*n, n, n, *pi[2]);
    current_x->add_block( re, 3*n, 1*n, n, n, *pi[0]);
    current_x->add_block( im, 3*n, 1*n, n, n, *pi[1]);

    current_y->add_block(-im, 0*n, 2*n, n, n, *pi[0]);
    current_y->add_block( re, 0*n, 2*n, n, n, *pi[1]);
    current_y->add_block( im, 0*n, 3*n, n, n, *pi[2]);
    current_y->add_block( im, 1*n, 2*n, n, n, *pi[2]);
    current_y->add_block( im, 1*n, 3*n, n, n, *pi[0]);
    current_y->add_block( re, 1*n, 3*n, n, n, *pi[1]);

    current_y->add_block( im, 2*n, 0*n, n, n, *pi[0]);
    current_y->add_block( re, 2*n, 0*n, n, n, *pi[1]);
    current_y->add_block(-im, 2*n, 1*n, n, n, *pi[2]);
    current_y->add_block(-im, 3*n, 0*n, n, n, *pi[2]);
    current_y->add_block(-im, 3*n, 1*n, n, n, *pi[0]);
    current_y->add_block( re, 3*n, 1*n, n, n, *pi[1]);

    current_z->add_block( re, 0*n, 2*n, n, n, *pi[2]);
    current_z->add_block( re, 0*n, 3*n, n, n, *pi[0]);
    current_z->add_block(-im, 0*n, 3*n, n, n, *pi[1]);
    current_z->add_block(-re, 1*n, 2*n, n, n, *pi[0]);
    current_z->add_block(-im, 1*n, 2*n, n, n, *pi[1]);
    current_z->add_block( re, 1*n, 3*n, n, n, *pi[2]);

    current_z->add_block( re, 2*n, 0*n, n, n, *pi[2]);
    current_z->add_block(-re, 2*n, 1*n, n, n, *pi[0]);
    current_z->add_block( im, 2*n, 1*n, n, n, *pi[1]);
    current_z->add_block( re, 3*n, 0*n, n, n, *pi[0]);
    current_z->add_block( im, 3*n, 0*n, n, n, *pi[1]);
    current_z->add_block( re, 3*n, 1*n, n, n, *pi[2]);

    overlap_->print("overlap matrix in AO basis", 40);
    current_x->print("integrated x-current in AO basis", 40);
    current_y->print("integrated y-current in AO basis", 40);
    current_z->print("integrated z-current in AO basis", 40);

    //shared_ptr<const ZMatrix> ocoeff1 = coeff_->slice_copy(nneg_, nneg_+nele_);
    //shared_ptr<const ZMatrix> ocoeff2 = coeff_->slice_copy(0, nele_);
    //shared_ptr<const ZMatrix> ocoeff = ocoeff1->merge(ocoeff2);
    //shared_ptr<const ZMatrix> ocoeff = coeff_->slice_copy(nneg_, 2*nneg_);
    shared_ptr<const ZMatrix> ocoeff = coeff_->slice_copy(nneg_, nneg_+nele_);
    ocoeff->print("MO Coefficient matrix", 40);

    auto ovlmo = make_shared<ZMatrix>(*ocoeff % *overlap_ * *ocoeff);
    auto x_current = make_shared<ZMatrix>(*ocoeff % *current_x * *ocoeff);
    auto y_current = make_shared<ZMatrix>(*ocoeff % *current_y * *ocoeff);
    auto z_current = make_shared<ZMatrix>(*ocoeff % *current_z * *ocoeff);

    ovlmo->print("overlap matrix in MO basis", 40);
    x_current->print("integrated x-current in MO basis", 40);
    y_current->print("integrated y-current in MO basis", 40);
    z_current->print("integrated z-current in MO basis", 40);

    const complex<double> ovltot = ovlmo->trace();
    const complex<double> xtot = x_current->trace();
    const complex<double> ytot = y_current->trace();
    const complex<double> ztot = z_current->trace();

    cout << "total self-overlap = " << ovltot << endl;
    cout << "total x-current = " << xtot << endl;
    cout << "total y-current = " << ytot << endl;
    cout << "total z-current = " << ztot << endl;

    const complex<double> ovlprod = aodensity->dot_product(*overlap_);
    const complex<double> xprod = aodensity->dot_product(*current_x);
    const complex<double> yprod = aodensity->dot_product(*current_y);
    const complex<double> zprod = aodensity->dot_product(*current_z);

    cout << endl;
    cout << "recalculated self-overlap = " << ovlprod << endl;
    cout << "recalculated x-current = " << xprod << endl;
    cout << "recalculated y-current = " << yprod << endl;
    cout << "recalculated z-current = " << zprod << endl;
#endif

#if 1
    // at a specific point
    const int n = geom_->nbasis();
    const array<double,3> loc = idata_->get_array<double,3>("point", {{0.0, 0.0, 0.0}});
    auto mom = make_shared<Momentum_Point>(geom_, loc);
    array<shared_ptr<ZMatrix>,3> pi = mom->compute();
    auto current_x = make_shared<ZMatrix>(4*n, 4*n);
    auto current_y = make_shared<ZMatrix>(4*n, 4*n);
    auto current_z = make_shared<ZMatrix>(4*n, 4*n);

    const complex<double> re(-0.5,  0.0);
    const complex<double> im( 0.0, -0.5);

    current_x->add_block( re, 0*n, 2*n, n, n, *pi[0]);
    current_x->add_block( im, 0*n, 2*n, n, n, *pi[1]);
    current_x->add_block(-re, 0*n, 3*n, n, n, *pi[2]);
    current_x->add_block( re, 1*n, 2*n, n, n, *pi[2]);
    current_x->add_block( re, 1*n, 3*n, n, n, *pi[0]);
    current_x->add_block(-im, 1*n, 3*n, n, n, *pi[1]);

    current_x->add_block( re, 2*n, 0*n, n, n, *pi[0]);
    current_x->add_block(-im, 2*n, 0*n, n, n, *pi[1]);
    current_x->add_block( re, 2*n, 1*n, n, n, *pi[2]);
    current_x->add_block(-re, 3*n, 0*n, n, n, *pi[2]);
    current_x->add_block( re, 3*n, 1*n, n, n, *pi[0]);
    current_x->add_block( im, 3*n, 1*n, n, n, *pi[1]);

    current_y->add_block(-im, 0*n, 2*n, n, n, *pi[0]);
    current_y->add_block( re, 0*n, 2*n, n, n, *pi[1]);
    current_y->add_block( im, 0*n, 3*n, n, n, *pi[2]);
    current_y->add_block( im, 1*n, 2*n, n, n, *pi[2]);
    current_y->add_block( im, 1*n, 3*n, n, n, *pi[0]);
    current_y->add_block( re, 1*n, 3*n, n, n, *pi[1]);

    current_y->add_block( im, 2*n, 0*n, n, n, *pi[0]);
    current_y->add_block( re, 2*n, 0*n, n, n, *pi[1]);
    current_y->add_block(-im, 2*n, 1*n, n, n, *pi[2]);
    current_y->add_block(-im, 3*n, 0*n, n, n, *pi[2]);
    current_y->add_block(-im, 3*n, 1*n, n, n, *pi[0]);
    current_y->add_block( re, 3*n, 1*n, n, n, *pi[1]);

    current_z->add_block( re, 0*n, 2*n, n, n, *pi[2]);
    current_z->add_block( re, 0*n, 3*n, n, n, *pi[0]);
    current_z->add_block(-im, 0*n, 3*n, n, n, *pi[1]);
    current_z->add_block(-re, 1*n, 2*n, n, n, *pi[0]);
    current_z->add_block(-im, 1*n, 2*n, n, n, *pi[1]);
    current_z->add_block( re, 1*n, 3*n, n, n, *pi[2]);

    current_z->add_block( re, 2*n, 0*n, n, n, *pi[2]);
    current_z->add_block(-re, 2*n, 1*n, n, n, *pi[0]);
    current_z->add_block( im, 2*n, 1*n, n, n, *pi[1]);
    current_z->add_block( re, 3*n, 0*n, n, n, *pi[0]);
    current_z->add_block( im, 3*n, 0*n, n, n, *pi[1]);
    current_z->add_block( re, 3*n, 1*n, n, n, *pi[2]);

    overlap_->print("overlap matrix in AO basis", 40);
    current_x->print("point x-current in AO basis", 40);
    current_y->print("point y-current in AO basis", 40);
    current_z->print("point z-current in AO basis", 40);

    //shared_ptr<const ZMatrix> ocoeff1 = coeff_->slice_copy(nneg_, nneg_+nele_);
    //shared_ptr<const ZMatrix> ocoeff2 = coeff_->slice_copy(0, nele_);
    //shared_ptr<const ZMatrix> ocoeff = ocoeff1->merge(ocoeff2);
    //shared_ptr<const ZMatrix> ocoeff = coeff_->slice_copy(nneg_, 2*nneg_);
    shared_ptr<const ZMatrix> ocoeff = coeff_->slice_copy(nneg_, nneg_+nele_);
    ocoeff->print("MO Coefficient matrix", 40);

    auto ovlmo = make_shared<ZMatrix>(*ocoeff % *overlap_ * *ocoeff);
    auto x_current = make_shared<ZMatrix>(*ocoeff % *current_x * *ocoeff);
    auto y_current = make_shared<ZMatrix>(*ocoeff % *current_y * *ocoeff);
    auto z_current = make_shared<ZMatrix>(*ocoeff % *current_z * *ocoeff);

    ovlmo->print("overlap matrix in MO basis", 40);
    x_current->print("point x-current in MO basis", 40);
    y_current->print("point y-current in MO basis", 40);
    z_current->print("point z-current in MO basis", 40);

    const complex<double> ovltot = ovlmo->trace();
    const complex<double> xtot = x_current->trace();
    const complex<double> ytot = y_current->trace();
    const complex<double> ztot = z_current->trace();

    cout << "total self-overlap = " << ovltot << endl;
    cout << "total x-current = " << xtot << endl;
    cout << "total y-current = " << ytot << endl;
    cout << "total z-current = " << ztot << endl;

    const complex<double> ovlprod = aodensity->dot_product(*overlap_);
    const complex<double> xprod = aodensity->dot_product(*current_x);
    const complex<double> yprod = aodensity->dot_product(*current_y);
    const complex<double> zprod = aodensity->dot_product(*current_z);

    cout << endl;
    cout << "recalculated self-overlap = " << ovlprod << endl;
    cout << "recalculated x-current = " << xprod << endl;
    cout << "recalculated y-current = " << yprod << endl;
    cout << "recalculated z-current = " << zprod << endl;
#endif


#if 1
  // To print out many points

  vector<array<double,3>> coords;
  const int n = geom_->nbasis();

  for (int i=0; i<=20; ++i) {
    for (int j=0; j<=20; ++j) {
      array<double,3> current = {{ -5.0+0.5*i, -5.0+0.5*j, 1.00  }};
      coords.push_back(current);
    }
  }

  cout << "       x                y                z             Re(pi_x)         Im(pi_x)         Re(pi_y)         Im(pi_y)         Re(pi_z)         Im(pi_z)" << endl;

  for (int i=0; i!=coords.size(); ++i) {

    auto mom = make_shared<Momentum_Point>(geom_, coords[i]);
    array<shared_ptr<ZMatrix>,3> pi = mom->compute();
    auto current_x = make_shared<ZMatrix>(4*n, 4*n);
    auto current_y = make_shared<ZMatrix>(4*n, 4*n);
    auto current_z = make_shared<ZMatrix>(4*n, 4*n);

    const complex<double> re(-0.5,  0.0);
    const complex<double> im( 0.0, -0.5);

    current_x->add_block( re, 0*n, 2*n, n, n, *pi[0]);
    current_x->add_block( im, 0*n, 2*n, n, n, *pi[1]);
    current_x->add_block(-re, 0*n, 3*n, n, n, *pi[2]);
    current_x->add_block( re, 1*n, 2*n, n, n, *pi[2]);
    current_x->add_block( re, 1*n, 3*n, n, n, *pi[0]);
    current_x->add_block(-im, 1*n, 3*n, n, n, *pi[1]);

    current_x->add_block( re, 2*n, 0*n, n, n, *pi[0]);
    current_x->add_block(-im, 2*n, 0*n, n, n, *pi[1]);
    current_x->add_block( re, 2*n, 1*n, n, n, *pi[2]);
    current_x->add_block(-re, 3*n, 0*n, n, n, *pi[2]);
    current_x->add_block( re, 3*n, 1*n, n, n, *pi[0]);
    current_x->add_block( im, 3*n, 1*n, n, n, *pi[1]);

    current_y->add_block(-im, 0*n, 2*n, n, n, *pi[0]);
    current_y->add_block( re, 0*n, 2*n, n, n, *pi[1]);
    current_y->add_block( im, 0*n, 3*n, n, n, *pi[2]);
    current_y->add_block( im, 1*n, 2*n, n, n, *pi[2]);
    current_y->add_block( im, 1*n, 3*n, n, n, *pi[0]);
    current_y->add_block( re, 1*n, 3*n, n, n, *pi[1]);

    current_y->add_block( im, 2*n, 0*n, n, n, *pi[0]);
    current_y->add_block( re, 2*n, 0*n, n, n, *pi[1]);
    current_y->add_block(-im, 2*n, 1*n, n, n, *pi[2]);
    current_y->add_block(-im, 3*n, 0*n, n, n, *pi[2]);
    current_y->add_block(-im, 3*n, 1*n, n, n, *pi[0]);
    current_y->add_block( re, 3*n, 1*n, n, n, *pi[1]);

    current_z->add_block( re, 0*n, 2*n, n, n, *pi[2]);
    current_z->add_block( re, 0*n, 3*n, n, n, *pi[0]);
    current_z->add_block(-im, 0*n, 3*n, n, n, *pi[1]);
    current_z->add_block(-re, 1*n, 2*n, n, n, *pi[0]);
    current_z->add_block(-im, 1*n, 2*n, n, n, *pi[1]);
    current_z->add_block( re, 1*n, 3*n, n, n, *pi[2]);

    current_z->add_block( re, 2*n, 0*n, n, n, *pi[2]);
    current_z->add_block(-re, 2*n, 1*n, n, n, *pi[0]);
    current_z->add_block( im, 2*n, 1*n, n, n, *pi[1]);
    current_z->add_block( re, 3*n, 0*n, n, n, *pi[0]);
    current_z->add_block( im, 3*n, 0*n, n, n, *pi[1]);
    current_z->add_block( re, 3*n, 1*n, n, n, *pi[2]);

    cout << fixed << setprecision(12);
    const complex<double> xprod = aodensity->dot_product(*current_x);
    const complex<double> yprod = aodensity->dot_product(*current_y);
    const complex<double> zprod = aodensity->dot_product(*current_z);

    cout << ((coords[i][0] < 0) ? "" : " ") << coords[i][0] << "  "
         << ((coords[i][1] < 0) ? "" : " ") << coords[i][1] << "  "
         << ((coords[i][2] < 0) ? "" : " ") << coords[i][2] << "  "
         << ((xprod.real() < 0) ? "" : " ") << xprod.real() << "  "
         << ((xprod.imag() < 0) ? "" : " ") << xprod.imag() << "  "
         << ((yprod.real() < 0) ? "" : " ") << yprod.real() << "  "
         << ((yprod.imag() < 0) ? "" : " ") << yprod.imag() << "  "
         << ((zprod.real() < 0) ? "" : " ") << zprod.real() << "  "
         << ((zprod.imag() < 0) ? "" : " ") << zprod.imag() << "  "
         << endl;

  }

#endif


    //print_eig();
  }
}


//Print non dirac sea eigenvalues
void Dirac::print_eig() const {
  const int n = geom_->nbasis();
  for (int i = 0*n; i != 4*n; ++i) cout << setprecision(10) << setw(15) << eig_[i] <<  endl;
}


shared_ptr<const Reference> Dirac::conv_to_ref() const {
  // we store only positive state coefficients
  const size_t npos = coeff_->mdim() - nneg_;
  auto out = make_shared<RelReference>(geom_, coeff_, energy_, nneg_, nele_, npos-nele_, gaunt_, breit_);
  vector<double> eigp(eig_.begin()+nneg_, eig_.end());
  vector<double> eigm(eig_.begin(), eig_.begin()+nneg_);
  VectorB eig(eig_.size());
  copy(eigp.begin(), eigp.end(), eig.begin());
  copy(eigm.begin(), eigm.end(), eig.begin()+eigp.size());
  out->set_eig(eig);
  return out;
}


shared_ptr<const DistZMatrix> Dirac::initial_guess(const shared_ptr<const DistZMatrix> s12, const shared_ptr<const DistZMatrix> hcore) const {
  const int n = geom_->nbasis();
  VectorB eig(hcore->ndim());

  shared_ptr<const DistZMatrix> coeff;
  if (!ref_) {
    // No reference; starting from hcore
    DistZMatrix interm = *s12 % *hcore * *s12;
    interm.diagonalize(eig);
    coeff = make_shared<const DistZMatrix>(*s12 * interm);

  } else if (dynamic_pointer_cast<const RelReference>(ref_)) {
    auto relref = dynamic_pointer_cast<const RelReference>(ref_);

    if (relref->rel()) {
      // Relativistic (4-component) reference
      shared_ptr<ZMatrix> fock = make_shared<DFock>(geom_, hcore_, relref->relcoeff()->slice_copy(0, nele_), gaunt_, breit_, /*store_half*/false, robust_);
      DistZMatrix interm = *s12 % *fock->distmatrix() * *s12;
      interm.diagonalize(eig);
      coeff = make_shared<const DistZMatrix>(*s12 * interm);
    } else {
      // Non-relativistic, GIAO-based reference
      assert(geom_->magnetism());
      const int nocc = ref_->nocc();
      shared_ptr<ZMatrix> fock;
      assert(nocc*2 == nele_);
      auto ocoeff = make_shared<ZMatrix>(n*4, 2*nocc);
      ocoeff->add_block(1.0, 0,    0, n, nocc, relref->relcoeff()->slice(0,nocc));
      ocoeff->add_block(1.0, n, nocc, n, nocc, relref->relcoeff()->slice(0,nocc));
      fock = make_shared<DFock>(geom_, hcore_, ocoeff, gaunt_, breit_, /*store_half*/false, robust_);
      DistZMatrix interm = *s12 % *fock->distmatrix() * *s12;
      interm.diagonalize(eig);
      coeff = make_shared<const DistZMatrix>(*s12 * interm);
    }
  } else if (ref_->coeff()->ndim() == n) {
    // Non-relativistic, real reference
    assert(!geom_->magnetism());
    const int nocc = ref_->nocc();
    shared_ptr<ZMatrix> fock;
    if (nocc*2 == nele_) {
      // RHF
      auto ocoeff = make_shared<ZMatrix>(n*4, 2*nocc);
      ocoeff->add_real_block(1.0, 0,    0, n, nocc, ref_->coeff()->slice(0,nocc));
      ocoeff->add_real_block(1.0, n, nocc, n, nocc, ref_->coeff()->slice(0,nocc));
      fock = make_shared<DFock>(geom_, hcore_, ocoeff, gaunt_, breit_, /*store_half*/false, robust_);
    } else if (ref_->noccB() != 0) {
      // UHF & ROHF
      const int nocca = ref_->noccA();
      const int noccb = ref_->noccB();
      assert(nocca+noccb == nele_);
      auto ocoeff = make_shared<ZMatrix>(n*4, nocca+noccb);
      ocoeff->add_real_block(1.0, 0,     0, n, nocca, ref_->coeffA()->slice(0,nocca));
      ocoeff->add_real_block(1.0, n, nocca, n, noccb, ref_->coeffB()->slice(nocca,nocca+noccb));
      fock = make_shared<DFock>(geom_, hcore_, ocoeff, gaunt_, breit_, /*store_half*/false, robust_);
    } else {
      // CASSCF
      auto ocoeff = make_shared<ZMatrix>(n*4, 2*nele_);
      ocoeff->add_real_block(1.0, 0,     0, n, nele_, ref_->coeff()->slice(0,nele_));
      ocoeff->add_real_block(1.0, n, nele_, n, nele_, ref_->coeff()->slice(0,nele_));
      fock = make_shared<DFock>(geom_, hcore_, ocoeff, gaunt_, breit_, /*store_half*/false, robust_);
    }
    DistZMatrix interm = *s12 % *fock->distmatrix() * *s12;
    interm.diagonalize(eig);
    coeff = make_shared<const DistZMatrix>(*s12 * interm);
  } else {
    assert(ref_->coeff()->ndim() == n*4);
    throw logic_error("Invalid Reference provided for Dirac.  (Initial guess not implemented.)");
  }
  return coeff;
}
