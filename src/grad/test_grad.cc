//
// Newint - Parallel electron correlation program.
// Filename: test_grad.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
// Maintainer: Shiozaki group
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

#include <src/grad/gradeval_hf.h>
#include <src/wfn/reference.h>
#include <iostream>


using namespace std;

void test_grad(shared_ptr<Reference> ref) {

  const size_t start = ::clock();
  cout << "  testing grad.." << endl;

  // target quantity here ... ==========

  GradEval_HF g1(ref);
  vector<double> grad = g1.compute();


  cout << endl << "  * Nuclear energy gradient" << endl << endl;
  for (int i = 0; i != ref->geom()->natom(); ++i) {
    cout << "    o Atom " << setw(3) << i << endl;
    cout << "        x  " << setprecision(10) << setw(20) << fixed << grad[3*i+0] << endl;
    cout << "        y  " << setprecision(10) << setw(20) << fixed << grad[3*i+1] << endl;
    cout << "        z  " << setprecision(10) << setw(20) << fixed << grad[3*i+2] << endl;
  }

  cout << setw(50) << left << "  * Gradient computed with " << setprecision(3) << right <<
          setw(10) << (::clock() - start)/static_cast<double>(CLOCKS_PER_SEC) << endl << endl;

}


