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
#include <src/dimer/dimer.h>

using namespace std;

/************************************************************************************
*  Dimer::Dimer(shared_ptr<Geometry> A, vector<double> displacement)                *
*                                                                                   *
************************************************************************************/
Dimer::Dimer(shared_ptr<const Geometry> A, tuple<double,double,double> displacement) {
   /************************************************************
   *  Set up variables that will contain the organized info    *
   ************************************************************/
   shared_ptr<const Geometry> geomB(new const Geometry((*A), displacement));

   geompair_ = make_pair(A, geomB);
}

Dimer::Dimer(shared_ptr<const Reference> A, tuple<double,double,double> displacement)
{
   /************************************************************
   *  Set up variables that will contain the organized info    *
   ************************************************************/
   coeffs_.push_back(A->coeff());

   shared_ptr<const Geometry> geomA = A->geom();
   shared_ptr<const Geometry> geomB(new const Geometry((*geomA), displacement));

   geompair_ = make_pair(geomA, geomB);
}

shared_ptr<const Geometry> Dimer::geometry() {
   vector<shared_ptr<const Geometry> > geo_vec;
   geo_vec.push_back(geompair_.first);
   geo_vec.push_back(geompair_.second);

   shared_ptr<const Geometry> geo_out(new const Geometry(geo_vec));

   return geo_out;
}

shared_ptr<const Coeff> Dimer::coefficients() {
   shared_ptr<const Geometry> cgeo = geometry();

   return coefficients(cgeo);
}

shared_ptr<const Coeff> Dimer::coefficients(shared_ptr<const Geometry> cgeo) {
   assert( coeffs_.size() <= 2 );

   vector<shared_ptr<const Coeff> > out_cvec;

   if( coeffs_.empty() ) {
      throw runtime_error("Attempting to concatenate coefficients on a Dimer with no defined coefficients");
   }
   else if (coeffs_.size() == 1) {
      out_cvec.push_back(coeffs_.front());
      out_cvec.push_back(coeffs_.front());
   }
   else {
      out_cvec = coeffs_;
   }

   shared_ptr<const Coeff> out(new const Coeff(out_cvec));

   return out;
}

shared_ptr<const Coeff> Dimer::overlap(int nocc1, int nocc2) {
/* What I need to do is use a supergeo to make a big overlap matrix
   and then transform it with supercoeffs into the MO basis */
   if( nocc1 == -1 ) { nocc1 = geompair_.first->nbasis(); }
   if( nocc2 == -1 ) { nocc2 = geompair_.second->nbasis(); }

   if(coeffs_.size() == 0) {
      throw runtime_error("Attempting to compute overlaps in a Dimer without a coefficient matrix");
   }

   vector<shared_ptr<const Coeff> > coeff_vec;
   shared_ptr<Matrix1e> leftslice_m = coeffs_.front()->slice(0,nocc1);
   shared_ptr<const Coeff> leftslice_c(new const Coeff(*leftslice_m));
   
   coeff_vec.push_back(leftslice_c);

   shared_ptr<Matrix1e> rightslice_m = coeffs_.back()->slice(0,nocc2);
   shared_ptr<const Coeff> rightslice_c(new const Coeff(*rightslice_m));
   coeff_vec.push_back(rightslice_c);

   /* super geo */
   shared_ptr<const Geometry> supergeo = geometry();

   /* super coeff */
   shared_ptr<const Coeff> supercoeff(new const Coeff(coeff_vec));

   /* super overlap in AO basis */
   Overlap ovlp(supergeo);

   /* transform to MO basis with supercoeff */
   shared_ptr<const Coeff> novlp(new const Coeff( (*supercoeff) % ovlp * (*supercoeff) ));

   return novlp;
}
