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
#include <src/rel/relcoeff.h>
#include <src/util/zmatrix.h>
#include <src/util/matrix.h>
#include <src/util/diis.h>

using namespace std;
using namespace bagel;

Dirac::Dirac(const multimap<string, string>& idata_, const shared_ptr<const Geometry> geom,
             const shared_ptr<const Reference> re) : geom_(geom->relativistic()) {
  hcore_ = shared_ptr<const RelHcore>(new RelHcore(geom_));
  overlap_ = shared_ptr<const RelOverlap>(new RelOverlap(geom_, false));
  s12_ = shared_ptr<const RelOverlap>(new RelOverlap(geom_, true));
  // reading input keywords
  max_iter_ = read_input<int>(idata_, "maxiter", 100);
  max_iter_ = read_input<int>(idata_, "maxiter_scf", max_iter_);
  diis_start_ = read_input<int>(idata_, "diis_start", 1);
  thresh_scf_ = read_input<double>(idata_, "thresh", 1.0e-8);
  thresh_scf_ = read_input<double>(idata_, "thresh_scf", thresh_scf_);
}

void Dirac::compute() {
  Timer scftime;
  string indent = "  ";
  const int n = geom_->nbasis();

  shared_ptr<const DistZMatrix> hcore = hcore_->distmatrix();
  shared_ptr<const DistZMatrix> distovl = overlap_->distmatrix();
  shared_ptr<const DistZMatrix> s12 = s12_->distmatrix();

  // distributed hcore and overlap
  DistZMatrix interm = *s12 % *hcore * *s12;
  unique_ptr<double[]> eig(new double[hcore->ndim()]);
  interm.diagonalize(eig.get()); 

  // make a relativistic version of Coeff class (c.f. coeff.h in src/scf)
  // only implement form_density_rhf..
  const int nele = geom_->nele();
  const int nrows = geom_->nbasis();
  const int column = 2 * geom_->nbasis();
  const int nneg = 2 * geom_->nbasis(); 
  
  // coefficient matrix
  shared_ptr<const DistZMatrix> coeff(new DistZMatrix(*s12 * interm));
  shared_ptr<const DistZMatrix> aodensity = coeff->form_density_rhf(nele, nneg);

  cout << indent << "=== Nuclear Repulsion ===" << endl << indent << endl;
  cout << indent << fixed << setprecision(10) << setw(15) << geom_->nuclear_repulsion() << endl << endl;
  cout << indent << "    * DIIS with orbital gradients will be used." << endl << endl;
  scftime.tick_print("SCF startup");
  cout << endl;
  cout << indent << "=== Dirac RHF iteration (" + geom_->basisfile() + ", RKB) ===" << endl << indent << endl;

  DIIS<DistZMatrix, ZMatrix> diis(5);

  for (int iter = 0; iter != max_iter_; ++iter) {
    Timer ptime(1);

    // TODO fock construction here. Fock construction requires a local copy of coeff
#if 1
    shared_ptr<const ZMatrix> fock(new DFock(geom_, hcore_, coeff->matrix()->slice(nneg, nele+nneg)));
#else
    shared_ptr<const ZMatrix> fock = hcore->matrix();
#endif
    // distribute
    shared_ptr<const DistZMatrix> distfock = fock->distmatrix();

    // compute energy here
    const complex<double> prod = aodensity->zdotc(*hcore+*distfock);
    if (fabs(prod.imag()) > 1.0e-12) {
      stringstream ss; ss << "imaginary part of energy is nonzero!! Perhaps Fock is not Hermite for some reasons " << prod.imag();
      throw runtime_error(ss.str());
    }
    energy_ = 0.5*prod.real() + geom_->nuclear_repulsion();

    shared_ptr<const DistZMatrix> error_vector(new DistZMatrix(*distfock**aodensity**distovl - *distovl**aodensity**distfock));
    const double error = error_vector->rms();

    ptime.tick_print("Fock build");
    cout << indent << setw(5) << iter << setw(20) << fixed << setprecision(8) << energy_
         << "   " << setw(17) << error << setw(15) << setprecision(2) << scftime.tick() << endl;

    if (error < thresh_scf_) {
      cout << indent << endl << indent << "  * SCF iteration converged." << endl << endl;
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
    coeff = shared_ptr<RelCoeff>(new RelCoeff(*coeff * intermediate));

    aodensity = coeff->form_density_rhf(nele, nneg); 

  }

  print_eig(eig);
}


//Print non dirac sea eigenvalues
void Dirac::print_eig(const unique_ptr<double[]>& eig) {
  const int n = geom_->nbasis();
  for (int i = 2*n; i != 4*n; ++i) cout << setprecision(10) << setw(15) << eig[i] << endl;
}


shared_ptr<Reference> Dirac::conv_to_ref() const {
  assert(false);
  return shared_ptr<Reference>();
}

