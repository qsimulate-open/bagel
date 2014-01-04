//
// BAGEL - Parallel electron correlation program.
// Filename: bagel_interface.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Ryan Reynolds <rreynoldschem@u.northwestern.edu>
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

// This file defines a function that can be called from Bagel, using its inputs to compute the ERI using the functional architecture from my test code.
// It loops over the various contracted function and angular momentum combinations to compute an entire cartesian batch in series, for easy comparison to Bagel's output.

#include <cassert>
#include <complex>
#include <array>
#include <memory>
#include <vector>
#include <utility>
#include "ericompute.h"
#include "src/molecule/shell.h"

using namespace std;

namespace ryan {

// Returns a vector identifying all the possible Cartesian functions with total angular momentum L
vector<vector<int>> assign_angular (const int L) {
  vector<vector<int>> out = {};
  int x, y, z;
  for (int i=0; i<=L; i++) {
    z = i;
    for (int j=0; j<=(L-i); j++) {
      y = j;
      x = L - i - j;
      vector<int> next = {x, y, z};
      out.push_back (next);
    }
  }
#if 0
  cout << endl << "Assigning angular momentum: Total = " << L << endl;
  for (int i=0; i!=out.size(); i++) {
    cout << out[i][0] << " " << out[i][1] << " " << out[i][2] << endl;
  }
#endif
  return out;
}


// Returns a vector identifying all possible Spherical functions with total angular momentum L
// The elements of this vector identify linear combinations of the Cartesian functions produced by assign_angular (const int L)
vector<vector<double>> spherical_combinations (const int L, const bool convert) {
  vector<vector<double>> out = {};
  if (convert) {
    out.resize (2*L+1);
    if (L == 0) {
      const double one = 1.0;
      //        000      <-- Cartesian orbital indices
      out[0] = {one}; // s
    }
    if (L == 1) {
      const double one = 1.0;
      //        100  010  001      <-- Cartesian orbital indices
      out[0] = {one, 0.0, 0.0}; // p_x
      out[1] = {0.0, one, 0.0}; // p_y
      out[2] = {0.0, 0.0, one}; // p_z
    }
    if (L == 2) {
      const double one = 1.0;
      const double c0 = (sqrt(3.0)/2.0);
      const double c2 = 0.5;
      const double c1 = sqrt(3.0);
      //        200  110  020  101  011  002      <-- Cartesian orbital indices
      out[0] = { c0, 0.0, -c0, 0.0, 0.0, 0.0}; // d_x2-y2
      out[1] = {0.0,  c1, 0.0, 0.0, 0.0, 0.0}; // d_xy
      out[2] = {0.0, 0.0, 0.0,  c1, 0.0, 0.0}; // d_xz
      out[3] = {0.0, 0.0, 0.0, 0.0,  c1, 0.0}; // d_yz
      out[4] = {-c2, 0.0, -c2, 0.0, 0.0, one}; // d_z2
    }
    if (L == 3) {
      const double one = 1.0;
      const double c3 = 1.0 * sqrt(15);
      const double c6 = (1.5/sqrt(5.0)) * sqrt(5.0);
      const double c4 = sqrt(1.2) * sqrt(5.0);
      const double c5 = (sqrt(6.0)/4.0);
      const double c7 = (sqrt(1.2)/4.0) * sqrt(5.0);
      const double c2 = (sqrt(3.0)/2.0) * sqrt(5.0);
      const double c0 = (sqrt(10.0)/4.0);
      const double c1 = (1.5/sqrt(2.0)) * sqrt(5.0);
      //        300  210  120  030  201  111  021  102  012  003      <-- Cartesian orbital indices
      out[0] = { c0, 0.0, -c1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // f_x3-3xy2
      out[1] = {0.0,  c1, 0.0, -c0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // f_3x2y-y3
      out[2] = {0.0, 0.0, 0.0, 0.0,  c2, 0.0, -c2, 0.0, 0.0, 0.0}; // f_zx2-zy2
      out[3] = {0.0, 0.0, 0.0, 0.0, 0.0,  c3, 0.0, 0.0, 0.0, 0.0}; // f_xyz
      out[4] = {-c5, 0.0, -c7, 0.0, 0.0, 0.0, 0.0,  c4, 0.0, 0.0}; // f_xz2
      out[5] = {0.0, -c7, 0.0, -c5, 0.0, 0.0, 0.0, 0.0,  c4, 0.0}; // f_yz2
      out[6] = {0.0, 0.0, 0.0, 0.0, -c6, 0.0, -c6, 0.0, 0.0, one}; // f_z3
    }
    if (L > 3) throw runtime_error ("The test code can only do Cartesian-to-spherical conversion for f-type orbitals and below.");
  } else {
    // If not doing the spherical conversion, just assign all coefficients as one or zero
    const int size = (((L+1)*(L+2))/2);
    out.resize (size);
    for (int i=0; i!=size; i++) {
      out[i].resize(size);
      for (int j=0; j!=size; j++) {
        if (i==j) out[i][j] = 1.0;
        else out[i][j] = 0.0;
      }
    }
  }
  return out;
}


vector<pair<vector<int>,complex<double>>> get_comparison_ERI (const array<shared_ptr<const bagel::Shell>,4>& basisinfo) {

  const bool normalize_basis = 0;
  const bool scale_input = 0;
  const bool orthogonalize = 0;

  const vector<double> field = {  0.0032,  0.0006, -0.0051};

  // Declare vectors to be used below
  vector<double> positions = {};
  vector<int> full_angular = {};
  vector<int> angular = {};
  vector<int> ncontracted = {};
  vector<bool> spherical = {};
  vector<vector<int>> nprimitive = {};
  vector<vector<double>> exponents = {};
  vector<vector<double>> contraction_coefficients = {};
  nprimitive.resize(4);
  exponents.resize(4);
  contraction_coefficients.resize(4);

  // Pull input data from the four bage::Shells
  for (int i=0; i!=4; i++) {
    spherical.push_back (basisinfo[i]->spherical());
//    spherical.push_back (false);


    full_angular.push_back (basisinfo[i]->angular_number());
    for (int k=0; k!=3; k++) {
      positions.push_back (basisinfo[i]->position(k));
    }
    const int nexponents = basisinfo[i]->num_primitive();
    int counter = 0;
    ncontracted.push_back (basisinfo[i]->num_contracted());
    for (int j=0; j!=nexponents; j++) {
      exponents[i].push_back (basisinfo[i]->exponents(j));
    }
    for (int j=0; j!=ncontracted[i]; j++) {
      const int start = basisinfo[i]->contraction_ranges(j).first;
      const int finish = basisinfo[i]->contraction_ranges(j).second;
      const int nfunctions = finish - start;
      nprimitive[i].push_back (nfunctions);
      counter = counter + nfunctions;
      for (int k=start; k!=finish; k++) {
        contraction_coefficients[i].push_back (basisinfo[i]->contractions(j)[k]);
      }
    }
    assert (counter == nexponents);
  }
  assert (spherical[0]==spherical[1]);
  assert (spherical[2]==spherical[3]);

  vector<int> basis_size = {};
  vector<vector<vector<int>>> cartesian_all = {};
  vector<vector<vector<double>>> spherical_all = {};
  for (int i=0; i!=4; i++) {
    basis_size.push_back (((full_angular[i]+1)*(full_angular[i]+2))/2);
    cartesian_all.push_back (assign_angular (full_angular[i]));
    spherical_all.push_back (spherical_combinations(full_angular[i],spherical[i]));
  }

#if 0
  for (int i=0; i!=4; i++) {
    cout << "MO coefficients for center " << i+1 << ":" << endl;
    for (int j=0; j!=spherical_all[i].size(); j++) {
      for (int k=0; k!=spherical_all[i][j].size(); k++) {
        cout << spherical_all[i][j][k] << "  ";
        cout << endl;
      }
    }
  }
#endif

  const int ang0 = spherical_all[0].size();
  const int ang1 = spherical_all[1].size();
  const int ang2 = spherical_all[2].size();
  const int ang3 = spherical_all[3].size();

  // Check that the right values have been assigned for the number of angular momentum possibilities
#if 0
  assert (spherical[0] || ang0 == basis_size[0]);
  assert (spherical[1] || ang1 == basis_size[1]);
  assert (spherical[2] || ang2 == basis_size[2]);
  assert (spherical[3] || ang3 == basis_size[3]);
  assert (!spherical[0] || ang0 == (2*full_angular[0] + 1) );
  assert (!spherical[1] || ang1 == (2*full_angular[1] + 1) );
  assert (!spherical[2] || ang2 == (2*full_angular[2] + 1) );
  assert (!spherical[3] || ang3 == (2*full_angular[3] + 1) );
#endif

  const int fnc0 = ncontracted[0];
  const int fnc1 = ncontracted[1];
  const int fnc2 = ncontracted[2];
  const int fnc3 = ncontracted[3];

  const int total = ang0*ang1*ang2*ang3*fnc0*fnc1*fnc2*fnc3;
  vector<pair<vector<int>,complex<double>>> out = {};

  const int nbasis_contracted = basis_size[0] + basis_size[1] + basis_size[2] + basis_size[3];

  // Iterate over all possible combinations of contracted basis functions and angular momentum for a given shell
  for (int p=0; p!=fnc3; p++) {
    for (int l=0; l!=ang3; l++) {
      for (int o=0; o!=fnc2; o++) {
        for (int k=0; k!=ang2; k++) {
          for (int n=0; n!=fnc1; n++) {
            for (int j=0; j!=ang1; j++) {
              for (int m=0; m!=fnc0; m++) {
                for (int i=0; i!=ang0; i++) {

                  vector<int> indices = {i,j,k,l,m,n,o,p};
                  const int nbasis = nprimitive[0][m] + nprimitive[1][n] + nprimitive[2][o] + nprimitive[3][p];
                  vector<double> positions_now = {};
                  vector<int> angular_now = {};
                  vector<double> exponents_now = {};
                  vector<double> contractions_now = {};
                  vector<int> nprimitive_now = {};
                  vector<vector<complex<double>>> orbitals = {};
                  orbitals.resize(4);

                  for (int q=0; q!=4; q++) {
                    for (int x=0; x!=basis_size[q]; x++) {
                      const int s = indices[q];
                      const int t = indices[4+q];
                      nprimitive_now.push_back(nprimitive[q][t]);
                      for (int r=0; r!=nprimitive[q][t]; r++) {
                        for (int u=0; u!=3; u++) {
                          positions_now.push_back (positions[3*q+u]);
                          angular_now.push_back (cartesian_all[q][x][u]);
                        }
                      }
                      int position = 0;
                      for (int v=0; v!=t; v++) position += nprimitive[q][v];
                      for (int w=0; w!=nprimitive[q][t]; w++) {
                        exponents_now.push_back (exponents[q][position+w]);
                        contractions_now.push_back (contraction_coefficients[q][position+w]);
                      }
                      for (int y=0; y!=4; y++) {
                        if (y==q) orbitals[y].push_back(spherical_all[q][s][x]);
                        else orbitals[y].push_back(0.0);
                      }
                    }
                  }
                  for (int q=0; q!=4; q++) {
                    assert (orbitals[q].size()==nbasis_contracted);
                    //cout << "Mol. Orbital " << q+1 << ":  ";
                    //for (int z=0; z!=orbitals[q].size(); z++) {
                    //  cout << orbitals[q][z] << "  ";
                    //}
                    //cout << endl;
                  }

                  // Compute the ERI for this particular term!
                  complex<double> eri = compute_eri (nbasis_contracted, normalize_basis, scale_input, orthogonalize, field,
                              positions_now, angular_now, exponents_now, contractions_now, nprimitive_now, orbitals[0], orbitals[1], orbitals[2], orbitals[3]);
                  pair<vector<int>,complex<double>> result (indices, eri);
                  out.push_back(result);

                }
              }
            }
          }
        }
      }
    }
  }

  assert (total == out.size());

  return out;

}

}
