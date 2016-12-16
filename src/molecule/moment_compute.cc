//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: moment_compute.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <src/molecule/moment_compute.h>
#include <src/molecule/carsph_shell.h>
#include <src/util/math/matop.h>

using namespace std;
using namespace bagel;


// for each primitive basis function
array<shared_ptr<const Matrix>,6> MomentCompute::mblock(const Shell& shell, const double exponent) {

  const int angular_number = shell.angular_number();

  const int norig = (angular_number+1) * (angular_number+2) / 2;
  const int ninc = norig + angular_number + 2;
  const int ndec = norig - angular_number - 1;

  array<shared_ptr<Matrix>,6> mcart;
  for (int i = 0; i != 3; ++i) {
    mcart[i]   = make_shared<Matrix>(ninc, norig, true);
    mcart[3+i] = shell.aux_decrement() ? make_shared<Matrix>(ndec, norig, true) : nullptr;
  }

  for (int z = 0, column = 0; z <= angular_number; ++z) {
    for (int y = 0; y <= angular_number-z; ++y, ++column) {
      const int x = angular_number - y - z;

      // three components of the angular momentum
      const array<int,3> index = {{x, y, z}};

      const double talph = 2.0 * exponent;

      // k tells us which dimension of the momentum operator we're using
      for (int k = 0; k != 3; ++k) {

        // -i a_x phi^(x-1)
        if (index[k] != 0) {
          array<int,3> newindex = index;
          --newindex[k];
          const size_t row = newindex[2]*(angular_number) - newindex[2]*(newindex[2]-1)/2 + newindex[1];
          mcart[3+k]->element(row, column) = -index[k];
        }

        // +i 2alpha phi^(x+1)
        {
          array<int,3> newindex = index;
          ++newindex[k];
          const size_t row = newindex[2]*(angular_number+2) - newindex[2]*(newindex[2]-1)/2 + newindex[1];
          mcart[k]->element(row, column) = talph;
        }
      }
    }
  }

  array<shared_ptr<const Matrix>,6> out;

  // convert this block from cartesian to shell.spherical
  if (shell.spherical()) {
    shared_ptr<Matrix> carsphmatrix = carsph_matrix(angular_number);
    for (int i = 0; i != 6; ++i)
      out[i] = mcart[i] ? make_shared<Matrix>(*mcart[i] * *carsphmatrix) : nullptr;
  } else {
    copy(mcart.begin(), mcart.end(), out.begin());
  }

  return out;
}


// for each primitive basis function
array<shared_ptr<const ZMatrix>,9> MomentCompute::mblock(const Shell& shell, const double exponent, const array<double,3> magnetic_field, const bool london) {

  const int angular_number = shell.angular_number();

  const int norig = (angular_number+1) * (angular_number+2) / 2;
  const int ninc = norig + angular_number + 2;
  const int ndec = norig - angular_number - 1;
  const complex<double> imag(0.0, 1.0);

  array<shared_ptr<ZMatrix>,9> mcart;
  for (int i = 0; i != 3; ++i) {
    mcart[i]   = make_shared<ZMatrix>(ninc, norig, true);
    mcart[3+i] = shell.aux_decrement() ? make_shared<ZMatrix>(ndec, norig, true) : nullptr;
    mcart[6+i] = shell.aux_same() ? make_shared<ZMatrix>(norig, norig, true) : nullptr;
  }

  for (int z = 0, column = 0; z <= angular_number; ++z) {
    for (int y = 0; y <= angular_number-z; ++y, ++column) {
      const int x = angular_number - y - z;

      // three components of the angular momentum
      const array<int,3> index = {{x, y, z}};
      const array<int,3> fwd  = {{1, 2, 0}};
      const array<int,3> back = {{2, 0, 1}};

      const complex<double> tialph = imag * 2.0 * exponent;
      const array<complex<double>,3> halfb = {{0.5*magnetic_field[0], 0.5*magnetic_field[1], 0.5*magnetic_field[2]}};

      // k tells us which dimension of the momentum operator we're using
      for (int k = 0; k != 3; ++k) {

        // -i a_x phi^(x-1)
        if (index[k] != 0) {
          array<int,3> newindex = index;
          --newindex[k];
          const size_t row = newindex[2]*(angular_number) - newindex[2]*(newindex[2]-1)/2 + newindex[1];
          mcart[3+k]->element(row, column) = -static_cast<double>(index[k])*imag;
        }

        // +i 2alpha phi^(x+1)
        {
          array<int,3> newindex = index;
          ++newindex[k];
          const size_t row = newindex[2]*(angular_number+2) - newindex[2]*(newindex[2]-1)/2 + newindex[1];
          mcart[k]->element(row, column) = tialph;
        }

        // + 1/2 B_y phi^(z+1)
        {
          array<int,3> newindex = index;
          ++newindex[back[k]];
          const size_t row = newindex[2]*(angular_number+2) - newindex[2]*(newindex[2]-1)/2 + newindex[1];
          mcart[k]->element(row, column) = halfb[fwd[k]];
        }

        // - 1/2 B_z phi^(y+1)
        {
          array<int,3> newindex = index;
          ++newindex[fwd[k]];
          const size_t row = newindex[2]*(angular_number+2) - newindex[2]*(newindex[2]-1)/2 + newindex[1];
          mcart[k]->element(row, column) = -halfb[back[k]];
        }

        // + 1/2 (B_y R_z - B_z R_y) phi
        if (!london) {
          array<int,3> newindex = index;
          const size_t row = newindex[2]*(angular_number+1) - newindex[2]*(newindex[2]-1)/2 + newindex[1];
          mcart[6+k]->element(row, column) = halfb[fwd[k]]*shell.position(back[k]) - halfb[back[k]]*shell.position(fwd[k]);
        }

      }
    }
  }

  array<shared_ptr<const ZMatrix>,9> out;
  for (int i = 0; i != 9; ++i)
    if (mcart[i])
      mcart[i]->scale(-imag);

  // convert this block from cartesian to shell.spherical
  if (shell.spherical()) {
    auto carsphmatrix = make_shared<ZMatrix>(*carsph_matrix(angular_number), 1.0);
    for (int i = 0; i != 9; ++i)
      out[i] = mcart[i] ? make_shared<ZMatrix>(*mcart[i] * *carsphmatrix) : nullptr;
  } else {
    copy(mcart.begin(), mcart.end(), out.begin());
  }

  return out;
}


array<shared_ptr<const Matrix>,3> MomentCompute::call(const Shell& shell) {

  const int angular_number = shell.angular_number();
  const int ncart = (angular_number+1) * (angular_number+2) / 2;
  const int ninc = ncart + (angular_number + 2);
  const int ndec = ncart - (angular_number + 1);
  const int n = ninc + ndec;
  const int norig = shell.spherical() ? 2*angular_number+1 : ncart;

  // build the momentum transformation matrix for primitive functions
  // each exponent gets 2 blocks, one for L+1 & one for L-1
  array<shared_ptr<Matrix>,3> tmp;
  for (int i = 0; i != 3; ++i)
    tmp[i] = make_shared<Matrix>(n*shell.num_primitive(), norig*shell.num_primitive(), true);

  for (int j = 0; j != shell.num_primitive(); ++j) {
    array<shared_ptr<const Matrix>,6> thisblock = mblock(shell, shell.exponents(j));
    for (int i = 0; i != 3; ++i) {
      tmp[i]->copy_block(j*ninc, j*norig, ninc, norig, thisblock[i]);
      if (shell.aux_decrement())
        tmp[i]->copy_block(ninc*shell.num_primitive()+j*ndec, j*norig, ndec, norig, thisblock[3+i]);
    }
  }

  // build the contraction matrix
  Matrix contract(norig*shell.num_primitive(), norig*shell.num_contracted(), true);
  for (int j = 0; j != shell.num_contracted(); ++j)
    for (int k = shell.contraction_ranges(j).first; k != shell.contraction_ranges(j).second; ++k)
      for (int l = 0; l != norig; ++l)
        contract(k*norig+l, j*norig+l) = shell.contractions(j)[k];

  // contract
  array<shared_ptr<const Matrix>,3> out;
  for (int i = 0; i != 3; ++i)
    out[i] = make_shared<Matrix>(*tmp[i] * contract);
  return out;
}


array<shared_ptr<const ZMatrix>,3> MomentCompute::call(const Shell& shell, const array<double,3> magnetic_field, const bool london) {

  const int angular_number = shell.angular_number();
  const int ncart = (angular_number+1) * (angular_number+2) / 2;
  const int ninc = ncart + (angular_number + 2);
  const int nsame = london ? 0 : ncart;
  const int ndec = ncart - (angular_number + 1);
  const int norig = shell.spherical() ? 2*angular_number+1 : ncart;
  const int n = ninc + ndec + nsame;

  // build the momentum transformation matrix for primitive functions
  // each exponent gets 2-3 blocks, for L+1, L-1, and if needed, L+0 (in that order)
  array<shared_ptr<ZMatrix>,3> tmp;
  for (int i = 0; i != 3; ++i)
    tmp[i] = make_shared<ZMatrix>(n*shell.num_primitive(), norig*shell.num_primitive(), true);

  for (int j = 0; j != shell.num_primitive(); ++j) {
    array<shared_ptr<const ZMatrix>,9> thisblock = mblock(shell, shell.exponents(j), magnetic_field, london);
    for (int i = 0; i != 3; ++i) {
      tmp[i]->copy_block(j*ninc, j*norig, ninc, norig, thisblock[i]);
      if (shell.aux_decrement())
        tmp[i]->copy_block(ninc*shell.num_primitive()+j*ndec, j*norig, ndec, norig, thisblock[3+i]);
      if (shell.aux_same())
        tmp[i]->copy_block((ninc+ndec)*shell.num_primitive()+j*nsame, j*norig, nsame, norig, thisblock[6+i]);
    }
  }

  // build the contraction matrix
  ZMatrix contract(norig*shell.num_primitive(), norig*shell.num_contracted(), true);
  for (int j = 0; j != shell.num_contracted(); ++j)
    for (int k = shell.contraction_ranges(j).first; k != shell.contraction_ranges(j).second; ++k)
      for (int l = 0; l != norig; ++l)
        contract(k*norig+l, j*norig+l) = shell.contractions(j)[k];

  // contract
  array<shared_ptr<const ZMatrix>,3> out;
  for (int i = 0; i != 3; ++i)
    out[i] = make_shared<ZMatrix>(*tmp[i] * contract);
  return out;
}

