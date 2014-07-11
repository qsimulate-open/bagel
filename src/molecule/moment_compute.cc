//
// BAGEL - Parallel electron correlation program.
// Filename: moment_compute.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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

#include <src/molecule/moment_compute.h>
#include <src/molecule/carsph_shell.h>
#include <src/math/matop.h>

using namespace std;
using namespace bagel;


// for each primitive basis function
array<shared_ptr<const Matrix>,6> MomentCompute::mblock(const Shell& shell, const double exponent) {

  const int angular_number = shell.angular_number();

  const int norig = (angular_number+1) * (angular_number+2) / 2;
  const int ninc = norig + angular_number + 2;
  const int ndec = norig - angular_number - 1;

  assert(ninc == (shell.aux_increment()->angular_number()+1) * (shell.aux_increment()->angular_number()+2) / 2);
  if (shell.aux_decrement())
    assert(ndec == (shell.aux_decrement()->angular_number()+1) * (shell.aux_decrement()->angular_number()+2) / 2);
  else
    assert(ndec == 0);

  array<shared_ptr<Matrix>,6> mcart;
  for (int i = 0; i != 3; ++i)
    mcart[i] = make_shared<Matrix>(ninc, norig, true);

  for (int i = 3; i != 6; ++i)
    mcart[i] = shell.aux_decrement() ? make_shared<Matrix>(ndec, norig, true) : nullptr;

  for (int i = 0, column = 0; i <= angular_number; ++i) {
    const int z = i;
    for (int j = 0; j <= angular_number-i; ++j, ++column) {
      const int y = j;
      const int x = angular_number - i - j;

      // three components of the angular momentum
      array<int,3> index = {{x, y, z}};

      assert(column == index[2]*(angular_number+1) - index[2]*(index[2]-1)/2 + index[1]);
      const double talph = 2.0 * exponent;

      // k tells us which dimension of the momentum operator we're using
      for (int k = 0; k != 3; ++k) {

        // -i a_x phi^(x-1)
        if (index[k] != 0) {
          array<int,3> newindex = index;
          --newindex[k];
          assert(newindex[0] + newindex[1] + newindex[2] == angular_number - 1);
          const size_t row = newindex[2]*(angular_number) - newindex[2]*(newindex[2]-1)/2 + newindex[1];
          mcart[3+k]->element(row, column) = -index[k];
        }

        // +i 2alpha phi^(x+1)
        {
          array<int,3> newindex = index;
          ++newindex[k];
          assert(newindex[0] + newindex[1] + newindex[2] == angular_number + 1);
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
  const complex<double> imag (0.0, 1.0);

  assert(ninc == (shell.aux_increment()->angular_number()+1) * (shell.aux_increment()->angular_number()+2) / 2);
  if (shell.aux_decrement())
    assert(ndec == (shell.aux_decrement()->angular_number()+1) * (shell.aux_decrement()->angular_number()+2) / 2);
  else
    assert (ndec == 0);
  if (shell.aux_same())
    assert(norig == (shell.aux_same()->angular_number()+1) * (shell.aux_same()->angular_number()+2) / 2);
  else
    assert(london);

  array<shared_ptr<ZMatrix>,9> mcart;
  for (int i = 0; i != 3; ++i)
    mcart[i] = make_shared<ZMatrix>(ninc, norig, true);

  for (int i = 3; i != 6; ++i)
    mcart[i] = shell.aux_decrement() ? make_shared<ZMatrix>(ndec, norig, true) : nullptr;

  for (int i = 6; i != 9; ++i)
    mcart[i] = shell.aux_same() ? make_shared<ZMatrix>(norig, norig, true) : nullptr;

  for (int i = 0, column = 0; i <= angular_number; ++i) {
    const int z = i;
    for (int j = 0; j <= angular_number-i; ++j, ++column) {
      const int y = j;
      const int x = angular_number - i - j;

      // three components of the angular momentum
      array<int,3> index = {{x, y, z}};
      const array<int,3> fwd  = {{1, 2, 0}};
      const array<int,3> back = {{2, 0, 1}};

      assert(column == index[2]*(angular_number+1) - index[2]*(index[2]-1)/2 + index[1]);
      const complex<double> tialph = imag * 2.0 * exponent;
      const array<complex<double>,3> halfb = {{0.5*magnetic_field[0], 0.5*magnetic_field[1], 0.5*magnetic_field[2]}};

      // k tells us which dimension of the momentum operator we're using
      for (int k = 0; k != 3; ++k) {

        // -i a_x phi^(x-1)
        if (index[k] != 0) {
          array<int,3> newindex = index;
          --newindex[k];
          assert(newindex[0] + newindex[1] + newindex[2] == angular_number - 1);
          const size_t row = newindex[2]*(angular_number) - newindex[2]*(newindex[2]-1)/2 + newindex[1];
          mcart[3+k]->element(row, column) = -static_cast<double>(index[k])*imag;
        }

        // +i 2alpha phi^(x+1)
        {
          array<int,3> newindex = index;
          ++newindex[k];
          assert(newindex[0] + newindex[1] + newindex[2] == angular_number + 1);
          const size_t row = newindex[2]*(angular_number+2) - newindex[2]*(newindex[2]-1)/2 + newindex[1];
          mcart[k]->element(row, column) = tialph;
        }

        // + 1/2 B_y phi^(z+1)
        {
          array<int,3> newindex = index;
          ++newindex[back[k]];
          assert(newindex[0] + newindex[1] + newindex[2] == angular_number + 1);
          const size_t row = newindex[2]*(angular_number+2) - newindex[2]*(newindex[2]-1)/2 + newindex[1];
          mcart[k]->element(row, column) = halfb[fwd[k]];
        }

        // - 1/2 B_z phi^(y+1)
        {
          array<int,3> newindex = index;
          ++newindex[fwd[k]];
          assert(newindex[0] + newindex[1] + newindex[2] == angular_number + 1);
          const size_t row = newindex[2]*(angular_number+2) - newindex[2]*(newindex[2]-1)/2 + newindex[1];
          mcart[k]->element(row, column) = -halfb[back[k]];
        }

        // + 1/2 (B_y R_z - B_z R_y) phi
        if (!london) {
          array<int,3> newindex = index;
          assert(newindex[0] + newindex[1] + newindex[2] == angular_number);
          const size_t row = newindex[2]*(angular_number+1) - newindex[2]*(newindex[2]-1)/2 + newindex[1];
          mcart[6+k]->element(row, column) = halfb[fwd[k]]*shell.position(back[k]) - halfb[back[k]]*shell.position(fwd[k]);
        }

      }
    }
  }

  array<shared_ptr<const ZMatrix>,9> out;
  for (int i = 0; i != 9; ++i)
    if (mcart[i])
      mcart[i]->scale(complex<double>(0.0, -1.0));

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

  assert(ninc == (shell.aux_increment()->angular_number()+1) * (shell.aux_increment()->angular_number()+2)/2);
  if (shell.aux_decrement())
    assert(ndec == (shell.aux_decrement()->angular_number()+1) * (shell.aux_decrement()->angular_number()+2)/2);
  else
    assert(ndec == 0);
  assert(shell.aux_increment()->num_primitive() == shell.num_primitive());
  if (shell.aux_decrement())
    assert(shell.aux_decrement()->num_primitive() == shell.num_primitive());
  assert(!shell.aux_same());

  // build the momentum transformation matrix for primitive functions
  // each exponent gets 2 blocks, one for L+1 & one for L-1
  array<shared_ptr<Matrix>,3> tmp;
  for (int i = 0; i != 3; ++i)
    tmp[i] = make_shared<Matrix>(n*shell.num_primitive(), norig*shell.num_primitive(), true);

  for (int j = 0; j != shell.num_primitive(); ++j) {
    array<shared_ptr<const Matrix>,6> thisblock = mblock(shell, shell.exponents(j));
    for (int i = 0; i != 3; ++i)
      tmp[i]->copy_block(j*ninc, j*norig, ninc, norig, thisblock[i]);
    if (shell.aux_decrement())
      for (int i = 0; i != 3; ++i)
        tmp[i]->copy_block(ninc*shell.num_primitive()+j*ndec, j*norig, ndec, norig, thisblock[3+i]);
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

  assert(ninc == (shell.aux_increment()->angular_number()+1) * (shell.aux_increment()->angular_number()+2)/2);
  if (shell.aux_decrement())
    assert(ndec == (shell.aux_decrement()->angular_number()+1) * (shell.aux_decrement()->angular_number()+2)/2);
  else
    assert(ndec == 0);
  if (shell.aux_same())
    assert(nsame == (shell.aux_same()->angular_number()+1) * (shell.aux_same()->angular_number()+2)/2);
  else
    assert(nsame == 0);
  assert(shell.aux_increment()->num_primitive() == shell.num_primitive());
  if (shell.aux_decrement())
    assert(shell.aux_decrement()->num_primitive() == shell.num_primitive());
  if (!london)
    assert(shell.aux_same()->num_primitive() == shell.num_primitive());
  else
    assert(!shell.aux_same());

  // build the momentum transformation matrix for primitive functions
  // each exponent gets 2-3 blocks, for L+1, L-1, and if needed, L+0 (in that order)
  array<shared_ptr<ZMatrix>,3> tmp;
  for (int i = 0; i != 3; ++i)
    tmp[i] = make_shared<ZMatrix>(n*shell.num_primitive(), norig*shell.num_primitive(), true);
  for (int j = 0; j != shell.num_primitive(); ++j) {
    array<shared_ptr<const ZMatrix>,9> thisblock = mblock(shell, shell.exponents(j), magnetic_field, london);
    for (int i = 0; i != 3; ++i)
      tmp[i]->copy_block(j*ninc, j*norig, ninc, norig, thisblock[i]);
    if (shell.aux_decrement())
      for (int i = 0; i != 3; ++i)
        tmp[i]->copy_block(ninc*shell.num_primitive()+j*ndec, j*norig, ndec, norig, thisblock[3+i]);
    if (shell.aux_same())
      for (int i = 0; i != 3; ++i)
        tmp[i]->copy_block((ninc+ndec)*shell.num_primitive()+j*nsame, j*norig, nsame, norig, thisblock[6+i]);
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
    out[i] = make_shared<const ZMatrix>(*tmp[i] * contract);
  return out;
}

