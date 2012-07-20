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
Dimer::Dimer(shared_ptr<Geometry> A, tuple<double,double,double> displacement) {
   /************************************************************
   *  Set up variables that will contain the organized info    *
   ************************************************************/
   shared_ptr<Geometry> geomB(new Geometry((*A), displacement));

   geompair_ = make_pair(A, geomB);
}

Dimer::Dimer(shared_ptr<const Reference> A, tuple<double,double,double> displacement)
{
   /************************************************************
   *  Set up variables that will contain the organized info    *
   ************************************************************/
   coeffs_.push_back(A->coeff());

   shared_ptr<Geometry> geomA(new Geometry(*(A->geom()), make_tuple(0.0,0.0,0.0))) ;
   shared_ptr<Geometry> geomB(new Geometry((*geomA), displacement));

   geompair_ = make_pair(geomA, geomB);
}

shared_ptr<Geometry> Dimer::geometry() {
   vector<shared_ptr<Geometry> > geo_vec;
   geo_vec.push_back(geompair_.first);
   geo_vec.push_back(geompair_.second);

   shared_ptr<Geometry> geo_out(new Geometry(geo_vec));

   return geo_out;
}

shared_ptr<const Geometry> Dimer::const_geometry() {
   vector<shared_ptr<Geometry> > geo_vec;
   geo_vec.push_back(geompair_.first);
   geo_vec.push_back(geompair_.second);

   shared_ptr<const Geometry> geo_out(new const Geometry(geo_vec));

   return geo_out;
}

shared_ptr<Coeff> Dimer::coefficients() {
   shared_ptr<const Geometry> cgeo = const_geometry();

   return coefficients(cgeo);
}

shared_ptr<Coeff> Dimer::coefficients(shared_ptr<const Geometry> cgeo) {
   assert( coeffs_.size() >= 0 && coeffs_.size() <= 2 );

   shared_ptr<Coeff> supercoeff(new Coeff(cgeo));
   int nbasis = coeffs_.front()->nbasis();

   double* catdata =supercoeff->data();

   if( coeffs_.empty() ) {
      throw runtime_error("Attempting to concatenate coefficients on a Dimer with no defined coefficients");
   }
   else if (coeffs_.size() == 1) {
      double* cdata0 = coeffs_.front()->data();

      for(int i = 0; i < nbasis; ++i) {
         for(int j = 0; j < nbasis; ++j, ++catdata, ++cdata0) {
            *catdata = *cdata0;
         }
         for(int j = 0; j < nbasis; ++j, ++catdata) {
            *catdata = 0.0;
         }
      }
      cdata0 = coeffs_.front()->data();
      for(int i = 0; i < nbasis; ++i) {
         for(int j = 0; j< nbasis; ++j, ++catdata) {
            *catdata = 0.0;
         }
         for(int j = 0; j< nbasis; ++j, ++catdata, ++cdata0) {
            *catdata = *cdata0;
         }
      }
   }
   else {
      double* cdata0 = coeffs_[0]->data();
      double* cdata1 = coeffs_[1]->data();

      for(int i = 0; i < nbasis; ++i) {
         for(int j = 0; j < nbasis; ++j, ++catdata, ++cdata0) {
            *catdata = *cdata0;
         }
         for(int j = 0; j < nbasis; ++j, ++catdata) {
            *catdata = 0.0;
         }
      }
      for(int i = 0; i< nbasis; ++i) {
         for(int j = 0; j< nbasis; ++j, ++catdata) {
            *catdata = 0.0;
         }
         for(int j = 0; j< nbasis; ++j, ++catdata, ++cdata1) {
            *catdata = *cdata1;
         }
      }

   }

   return supercoeff;
}

Matrix1e Dimer::overlap() {
/* What I need to do is use a supergeo to make a big overlap matrix
   and then transform it with supercoeffs into the MO basis */

   /* super geo */
   shared_ptr<const Geometry> supergeo = const_geometry();

   /* super coeff */
   shared_ptr<Coeff> supercoeff = coefficients(supergeo);

   /* super overlap in AO basis */
   Overlap ovlp(supergeo);

   /* transform to MO basis with supercoeff */
   Matrix1e novlp = (*supercoeff) % ovlp * (*supercoeff);

   return novlp;
}
