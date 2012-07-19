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

Dimer::Dimer(shared_ptr<Reference> A, tuple<double,double,double> displacement) :
   coeff_(A->coeff())
   {
   /************************************************************
   *  Set up variables that will contain the organized info    *
   ************************************************************/

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
