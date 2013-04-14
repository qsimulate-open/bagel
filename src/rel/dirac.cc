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
// the Free Software Foundation; either version 2, or (at your option)
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
#include <src/util/zmatrix.h>
#include <src/util/matrix.h>
#include <src/util/diis.h>
#include <src/rel/kramers.h>

using namespace std;
using namespace bagel;

Dirac::Dirac(const multimap<string, string>& idata, const shared_ptr<const Geometry> geom,
             const shared_ptr<const Reference> re) : ref_(re) {
  gaunt_ = read_input<bool>(idata, "gaunt", true);
  breit_ = read_input<bool>(idata, "breit", gaunt_);
  geom_ = geom->relativistic(gaunt_);
  common_init(idata);
}


Dirac::Dirac(const multimap<string, string>& idata, const shared_ptr<const Geometry> geom,
             const shared_ptr<const RelReference> re) : relref_(re) {
  gaunt_ = read_input<bool>(idata, "gaunt", true);
  breit_ = read_input<bool>(idata, "breit", gaunt_);
  geom_ = geom->relativistic(gaunt_);
  common_init(idata);
}


void Dirac::common_init(const multimap<string, string>& idata) {
  cout << "  *** Dirac HF ***" << endl << endl;

  hcore_ = shared_ptr<const RelHcore>(new RelHcore(geom_));
  overlap_ = shared_ptr<const RelOverlap>(new RelOverlap(geom_, false));
  s12_ = shared_ptr<const RelOverlap>(new RelOverlap(geom_, true));
  // reading input keywords
  max_iter_ = read_input<int>(idata, "maxiter", 100);
  max_iter_ = read_input<int>(idata, "maxiter_scf", max_iter_);
  diis_start_ = read_input<int>(idata, "diis_start", 1);
  thresh_scf_ = read_input<double>(idata, "thresh", 1.0e-8);
  thresh_scf_ = read_input<double>(idata, "thresh_scf", thresh_scf_);
  ncharge_ = read_input<int>(idata, "charge", 0);
  nele_ = geom_->nele()-ncharge_;
  nneg_ = geom_->nbasis()*2;

  if (breit_ && !gaunt_) throw runtime_error("Breit cannot be turned on if Gaunt is off");
}


void Dirac::compute() {
  Timer scftime;
  string indent = "  ";
  const int n = geom_->nbasis();

  shared_ptr<const DistZMatrix> hcore = hcore_->distmatrix();
  shared_ptr<const DistZMatrix> distovl = overlap_->distmatrix();
  shared_ptr<const DistZMatrix> s12 = s12_->distmatrix();
  unique_ptr<double[]> eig(new double[hcore->ndim()]);

  // making initial guess
  shared_ptr<const DistZMatrix> coeff = initial_guess(s12, hcore);
  shared_ptr<const DistZMatrix> aodensity = coeff->form_density_rhf(nele_, nneg_);

  cout << indent << "=== Nuclear Repulsion ===" << endl << indent << endl;
  cout << indent << fixed << setprecision(10) << setw(15) << geom_->nuclear_repulsion() << endl << endl;
  cout << indent << "    * DIIS with orbital gradients will be used." << endl << endl;
  scftime.tick_print("SCF startup");
  cout << endl;
  cout << indent << "=== Dirac RHF iteration (" + geom_->basisfile() + ", RKB) ===" << endl << indent << endl;

  DIIS<DistZMatrix, ZMatrix> diis(5);

  for (int iter = 0; iter != max_iter_; ++iter) {
    Timer ptime(1);

    // fock construction
    shared_ptr<ZMatrix> fock(new DFock(geom_, hcore_, coeff->matrix()->slice(nneg_, nele_+nneg_), gaunt_, breit_));
// TODO I have a feeling that the code should not need this, but sometimes there are slight errors. still looking on it.
#if 0
    fock->hermite();
#endif
    // distribute
    shared_ptr<const DistZMatrix> distfock = fock->distmatrix();

    // compute energy here
    const complex<double> prod = aodensity->zdotc(*hcore+*distfock);
    if (fabs(prod.imag()) > 1.0e-12) {
      stringstream ss; ss << "imaginary part of energy is nonzero!! Perhaps Fock is not Hermite for some reasons " << setprecision(10) << prod.imag();
//    throw runtime_error(ss.str());
      cout << ss.str() << endl;
    }
    energy_ = 0.5*prod.real() + geom_->nuclear_repulsion();

    shared_ptr<const DistZMatrix> error_vector(new DistZMatrix(*distfock**aodensity**distovl - *distovl**aodensity**distfock));
    const double error = error_vector->rms();

    ptime.tick_print("Fock build");
    cout << indent << setw(5) << iter << setw(20) << fixed << setprecision(8) << energy_
         << "   " << setw(17) << error << setw(15) << setprecision(2) << scftime.tick() << endl;

    if (error < thresh_scf_) {
      cout << indent << endl << indent << "  * SCF iteration converged." << endl << endl;
#if 1
      unique_ptr<double[]> eig2(new double[hcore->ndim()]);
      shared_ptr<Kramers> kramers(new Kramers(n));
      DistZMatrix int2(*coeff % *kramers * *coeff);
      int2.diagonalize(eig2.get());
      print_eig(eig2);
      DistZMatrix int3(int2 % *distfock * int2);
      int3.diagonalize(eig2.get());
      print_eig(eig2);
#endif
      break;
    } else if (iter == max_iter_-1) {
      cout << indent << endl << indent << "  * Max iteration reached in SCF." << endl << endl;
      break;
    }

    if (iter >= diis_start_) {
      distfock = diis.extrapolate(make_pair(distfock, error_vector));
      ptime.tick_print("DIIS");
    }

    DistZMatrix intermediate(*coeff % *distfock * *coeff);
    intermediate.diagonalize(eig.get());
    coeff = shared_ptr<DistZMatrix>(new DistZMatrix(*coeff * intermediate));

    aodensity = coeff->form_density_rhf(nele_, nneg_);

  }

  coeff_ = coeff->matrix();
//print_eig(eig);
}


//Print non dirac sea eigenvalues
void Dirac::print_eig(const unique_ptr<double[]>& eig) {
  const int n = geom_->nbasis();
  for (int i = 2*n; i != 4*n; ++(++i)) cout << setprecision(10) << setw(15) << eig[i] <<  "    " << eig[i+1] << endl;
}


shared_ptr<RelReference> Dirac::conv_to_ref() const {
  // we store only positive state coefficients
  const int npos = geom_->nbasis()*2;
  return shared_ptr<RelReference>(new RelReference(geom_, coeff_->slice(nneg_, nneg_+npos), energy_));
}


shared_ptr<const DistZMatrix> Dirac::initial_guess(const shared_ptr<const DistZMatrix> s12, const shared_ptr<const DistZMatrix> hcore) const {
  const int n = geom_->nbasis();
  unique_ptr<double[]> eig(new double[hcore->ndim()]);

  shared_ptr<const DistZMatrix> coeff;
  if (!ref_ && !relref_) {
    DistZMatrix interm = *s12 % *hcore * *s12;
    interm.diagonalize(eig.get());
    coeff = shared_ptr<const DistZMatrix>(new DistZMatrix(*s12 * interm));
  } else if (relref_) {
    shared_ptr<ZMatrix> fock(new DFock(geom_, hcore_, relref_->coeff()->slice(0, nele_), gaunt_, breit_));
    DistZMatrix interm = *s12 % *fock->distmatrix() * *s12;
    interm.diagonalize(eig.get());
    coeff = shared_ptr<const DistZMatrix>(new DistZMatrix(*s12 * interm));
  } else if (ref_->coeff()->ndim() == n) {
    // non-relativistic reference.
    const int nocc = ref_->nocc();
    shared_ptr<ZMatrix> fock;
    if (nocc*2 == nele_) {
      shared_ptr<ZMatrix> ocoeff(new ZMatrix(n*4, 2*nocc));
      ocoeff->add_real_block(complex<double>(1.0,0.0), 0,    0, n, nocc, ref_->coeff()->data());
      ocoeff->add_real_block(complex<double>(1.0,0.0), n, nocc, n, nocc, ref_->coeff()->data());
      fock = shared_ptr<ZMatrix>(new DFock(geom_, hcore_, ocoeff, gaunt_, breit_));
    } else {
      const int nocca = ref_->noccA();
      const int noccb = ref_->noccB();
      assert(nocca+noccb == nele_);
      shared_ptr<ZMatrix> ocoeff(new ZMatrix(n*4, nocca+noccb));
      ocoeff->add_real_block(complex<double>(1.0,0.0), 0,     0, n, nocca, ref_->coeffA()->data());
      ocoeff->add_real_block(complex<double>(1.0,0.0), n, nocca, n, noccb, ref_->coeffB()->data());
      fock = shared_ptr<ZMatrix>(new DFock(geom_, hcore_, ocoeff, gaunt_, breit_));
    }
    DistZMatrix interm = *s12 % *fock->distmatrix() * *s12;
    interm.diagonalize(eig.get());
    coeff = shared_ptr<const DistZMatrix>(new DistZMatrix(*s12 * interm));
  } else {
    assert(ref_->coeff()->ndim() == n*4);
    throw logic_error("not yet implemented");
  }
  return coeff;
}
