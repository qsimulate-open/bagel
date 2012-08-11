//
// Newint - Parallel electron correlation program.
// Filename: dimer.cc
// Copyright (C) 2012 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//

#include <tuple>

#include <src/scf/geometry.h>
#include <src/scf/coeff.h>
#include <src/scf/matrix1e.h>
#include <src/scf/overlap.h>
#include <src/scf/fock.h>
#include <src/dimer/dimer.h>

using namespace std;

/************************************************************************************
*  Dimer::Dimer(shared_ptr<Geometry> A, vector<double> displacement)                *
*                                                                                   *
************************************************************************************/
Dimer::Dimer(shared_ptr<const Geometry> A, array<double,3> displacement) : nbasis_(2*A->nbasis()) {
   /************************************************************
   *  Set up variables that will contain the organized info    *
   ************************************************************/
   shared_ptr<const Geometry> geomB(new const Geometry((*A), displacement));

   geompair_ = make_pair(A, geomB);

   cout << " ===== Constructing Dimer geometry ===== " << endl;
   construct_geometry();
}

Dimer::Dimer(shared_ptr<const Reference> A, array<double,3> displacement) : nbasis_(2*A->geom()->nbasis())
{
   /************************************************************
   *  Set up variables that will contain the organized info    *
   ************************************************************/
   coeffs_.push_back(A->coeff());

   shared_ptr<const Geometry> geomA = A->geom();
   shared_ptr<const Geometry> geomB(new const Geometry((*geomA), displacement));

   geompair_ = make_pair(geomA, geomB);

   cout << " ===== Constructing Dimer geometry ===== " << endl;
   construct_geometry();

   cout << " ===== Constructing Dimer reference ===== " << endl;
   construct_coeff();

   int nclo = 2*A->nclosed();
   int nact = 2*A->nact();
   int nvirt = 2*A->nvirt();

   superreference_ = shared_ptr<Reference>(new Reference(supergeometry_, supercoeff_, nclo, nact, nvirt ));
}

void Dimer::construct_geometry() {
   vector<shared_ptr<const Geometry> > geo_vec;
   geo_vec.push_back(geompair_.first);
   geo_vec.push_back(geompair_.second);

   supergeometry_ = shared_ptr<Geometry>(new Geometry(geo_vec));
}

/* Note to self: assuming everything is a homodimer so the MOs get striped */
void Dimer::construct_coeff() {
   assert( coeffs_.size() <= 2 );

   if( coeffs_.empty() ) {
      throw runtime_error("Attempting to concatenate coefficients on a Dimer with no defined coefficients");
   }

   supercoeff_ = shared_ptr<Coeff>(new Coeff(supergeometry_));

   double *Adata = coeffs_.front()->data();
   double *Bdata = coeffs_.back()->data();
   double *Sdata = supercoeff_->data();

   int num_basis = coeffs_.front()->nbasis();

   for(int ii = 0; ii != num_basis; ++ii) {
      Sdata = copy(Adata, Adata + num_basis, Sdata);
      /* Fill twice, one for trailing zeros of MO on A, then for leading zeros on MO of B */
      fill(Sdata, Sdata + 2*num_basis, 0.0); Sdata += 2*num_basis; 
      Sdata = copy(Bdata, Bdata + num_basis, Sdata);
      Adata += num_basis; Bdata += num_basis;
   }
}

shared_ptr<Coeff> Dimer::overlap() {
/* What I need to do is use a supergeo to make a big overlap matrix
   and then transform it with supercoeffs into the MO basis */

   if(coeffs_.size() == 0) {
      throw runtime_error("Attempting to compute overlaps in a Dimer without a coefficient matrix");
   }

   /* super overlap in AO basis */
   Overlap ovlp(supergeometry_);

   /* transform to MO basis with supercoeff */
   shared_ptr<Coeff> novlp(new Coeff( (*supercoeff_) % ovlp * (*supercoeff_) ));

   return novlp;
}

void Dimer::orthonormalize() {
   shared_ptr<Coeff> ovlp = overlap();
   Matrix1e S = *ovlp;

   unique_ptr<double[]> eig(new double[nbasis_]);
   S.diagonalize(eig.get());

   Matrix1e S_1_2(S);
   double *S12_data = S_1_2.data();
   for( int ii = 0; ii != nbasis_; ++ii) {
      double scale = 1.0/sqrt(eig[ii]);
      for( int jj = 0; jj != nbasis_; ++jj) {
         *S12_data++ *= scale;
      }
   }

   S_1_2 = S_1_2 * *(S.transpose());

   supercoeff_ = shared_ptr<Coeff>(new Coeff(*supercoeff_ * S_1_2));
   int nclo = superreference_->nclosed();
   int nact =  superreference_->nact();
   int nvirt = superreference_->nvirt();

   superreference_ = shared_ptr<Reference>(new Reference(supergeometry_, supercoeff_, nclo, nact, nvirt ));
}

double Dimer::energy() {
   shared_ptr<Matrix1e> ao_density = supercoeff_->form_density_rhf(superreference_->nclosed());
   shared_ptr<Fock<1> > hcore(new Fock<1>(supergeometry_));
   shared_ptr<Fock<1> > fock(new Fock<1>(supergeometry_, hcore, ao_density, supergeometry_->schwarz()));

   Matrix1e hcore_fock = (*hcore + *fock);
   double energy = ao_density->ddot(*hcore_fock.transpose());
   energy = 0.5*energy + supergeometry_->nuclear_repulsion();

   return energy;
}

