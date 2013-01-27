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
#include <src/util/zmatrix.h>
#include <src/util/matrix.h>
#include <src/util/diis.h>

using namespace std;
using namespace bagel;

void Dirac::compute() {
  Timer scftime;
  string indent = "  ";

  shared_ptr<ZMatrix> hcore = hcore_construct();
  shared_ptr<ZMatrix> s12 = s12_construct();
  const int n = geom_->nbasis();
  shared_ptr<ZMatrix> overlap(new ZMatrix(4*n, 4*n));
  {
    // Maybe you would want to make a DOverlap class.
    shared_ptr<Matrix> ovl(new Overlap(*overlap_));
    shared_ptr<Matrix> k12(new Matrix(*kinetic_ * (0.5/(c__*c__))));
    complex<double> one(1.0);
    overlap->copy_real_block(one, 0, 0, n, n, ovl);
    overlap->copy_real_block(one, n, n, n, n, ovl);
    overlap->copy_real_block(one, 2*n, 2*n, n, n, k12); 
    overlap->copy_real_block(one, 3*n, 3*n, n, n, k12); 
  }

  ZMatrix interm = *s12 % *hcore * *s12;
  unique_ptr<double[]> eig(new double[hcore->ndim()]);
  interm.diagonalize(eig.get()); 

  // make a relativistic version of Coeff class (c.f. coeff.h in src/scf)
  // only implement form_density_rhf..
  const int nele = geom_->nele();
  const int nrows = geom_->nbasis();
  const int column = 2 * geom_->nbasis();
  const int nneg = 2 * geom_->nbasis(); 
  
  // coefficient matrix
  shared_ptr<const ZMatrix> coeff(new ZMatrix(*s12 * interm));
  shared_ptr<const ZMatrix> ocoeff = coeff->slice(nneg, nneg+nele);
  shared_ptr<const ZMatrix> aodensity(new ZMatrix(*ocoeff ^ *ocoeff));

  cout << indent << "=== Nuclear Repulsion ===" << endl << indent << endl;
  cout << indent << fixed << setprecision(10) << setw(15) << geom_->nuclear_repulsion() << endl << endl;
  cout << indent << "    * DIIS with orbital gradients will be used." << endl << endl;
  scftime.tick_print("SCF startup");
  cout << endl;
  cout << indent << "=== Dirac RHF iteration (" + geom_->basisfile() + ", RKB) ===" << endl << indent << endl;

  DIIS<ZMatrix> diis(5);
  // TODO these are going to be read from input >>>>>>>>>>>>>
  const int max_iter_ = 100;
  const int diis_start_ = 0;
  const double thresh_scf_ = 1.0e-8;
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  for (int iter = 0; iter != max_iter_; ++iter) {
    Timer ptime(1);

    // TODO fock construction here
#if 0
    shared_ptr<ZMatrix> fock(new DFock(...));
#else
    shared_ptr<ZMatrix> fock = hcore;
#endif

    // compute energy here
    const complex<double> prod = aodensity->zdotc(*hcore+*fock);
    if (fabs(prod.imag()) > 1.0e-12) {
      stringstream ss; ss << "imaginary part of energy is nonzero!! Perhaps Fock is not Hermite for some reasons " << prod.imag();
      throw runtime_error(ss.str());
    }
    energy_ = 0.5*aodensity->zdotc(*hcore+*fock).real() + geom_->nuclear_repulsion();

    shared_ptr<const ZMatrix> error_vector(new ZMatrix(*fock**aodensity**overlap - *overlap**aodensity**fock));
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

#if 0
    // TODO use DIIS
    if (iter >= diis_start_) {
      fock = diis.extrapolate(make_pair(fock, error_vector));
      pdebug.tick_print("DIIS");
    }
#endif

    ZMatrix intermediate(*coeff % *fock * *coeff);
    intermediate.diagonalize(eig.get());
    coeff = shared_ptr<const ZMatrix>(new ZMatrix(*coeff * intermediate));
    ocoeff = coeff->slice(nneg, nneg+nele);
    aodensity = shared_ptr<const ZMatrix>(new ZMatrix(*ocoeff ^ *ocoeff));

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


shared_ptr<ZMatrix> Dirac::hcore_construct() {
  const int n = geom_->nbasis();

  shared_ptr<ZMatrix> out(new ZMatrix(4*n, 4*n));
  shared_ptr<ZMatrix> znai(new ZMatrix(2*n, 2*n));
  shared_ptr<ZMatrix> zkinetic(new ZMatrix(2*n, 2*n));

  array<shared_ptr<ZMatrix>,4> zsmallnai;
  for (auto& i : zsmallnai)
    i = znai->clone(); 

  const complex<double> coeff1 (1.0, 0.0);
  const complex<double> coeffi (0.0, 1.0);

  znai->copy_real_block(coeff1, 0, 0, n, n, nai_);
  znai->copy_real_block(coeff1, n, n, n, n, nai_);
  zkinetic->copy_real_block(coeff1, 0, 0, n, n, kinetic_);
  zkinetic->copy_real_block(coeff1, n, n, n, n, kinetic_);

  zsmallnai[0]->copy_real_block(coeff1, 0, 0, n, n, (*smallnai_)[0]);
  zsmallnai[0]->copy_real_block(coeff1, n, n, n, n, (*smallnai_)[0]);
  zsmallnai[1]->copy_real_block(coeffi, 0, 0, n, n, (*smallnai_)[1]);
  zsmallnai[1]->copy_real_block(-coeffi, n, n, n, n, (*smallnai_)[1]);
  zsmallnai[2]->copy_real_block(coeffi, 0, n, n, n, (*smallnai_)[2]);
  zsmallnai[2]->copy_real_block(coeffi, n, 0, n, n, (*smallnai_)[2]);
  zsmallnai[3]->copy_real_block(coeff1, 0, n, n, n, (*smallnai_)[3]);
  zsmallnai[3]->copy_real_block(-coeff1, n, 0, n, n, (*smallnai_)[3]);

  shared_ptr<ZMatrix> smallnai(new ZMatrix(*zsmallnai[0] + *zsmallnai[1] + *zsmallnai[2] + *zsmallnai[3]));

  // RKB hcore: T is off diagonal block matrices, V is first main diagonal, and 1/4m^2c^2W-T is second main diagonal
  const complex<double> w(0.25/(c__*c__), 0.0);
  out->zero();
  out->copy_block(0, 0, 2*n, 2*n, znai);
  out->copy_block(0, 2*n, 2*n, 2*n, zkinetic);
  out->copy_block(2*n, 0, 2*n, 2*n, zkinetic);
  out->copy_block(2*n, 2*n, 2*n, 2*n, shared_ptr<ZMatrix>(new ZMatrix(*smallnai * w - *zkinetic)));

  return out;
}


shared_ptr<ZMatrix> Dirac::s12_construct() {
  const int n = geom_->nbasis();
  const complex<double> coeff1 (1.0, 0.0);

  shared_ptr<ZMatrix> out(new ZMatrix(4*n, 4*n));
  shared_ptr<Matrix> ovl(new Overlap(*overlap_));
  shared_ptr<Matrix> k12(new Matrix(*kinetic_ * (0.5/(c__*c__))));
  ovl->inverse_half();
  k12->inverse_half();

  out->copy_real_block(coeff1, 0, 0, n, n, ovl);
  out->copy_real_block(coeff1, n, n, n, n, ovl);
  out->copy_real_block(coeff1, 2*n, 2*n, n, n, k12); 
  out->copy_real_block(coeff1, 3*n, 3*n, n, n, k12); 

  return out;
}
