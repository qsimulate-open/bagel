//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: carsph_shell.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@northwestern.edu>
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


#ifndef __SRC_MOLECULE_CARSPH_SHELL_H
#define __SRC_MOLECULE_CARSPH_SHELL_H

#include <memory>
#include <src/util/math/matrix.h>

namespace bagel {

static std::shared_ptr<Matrix> carsph_matrix (const int i) {
  const int m = (i+1)*(i+2)/2;
  const int n = 2*i+1;
  std::shared_ptr<Matrix> out = std::make_shared<Matrix> (m, n, true);

  constexpr double d0 = 0.8660254037844386;
  constexpr double d1 = 1.7320508075688772;
  constexpr double d2 = 0.5;
  constexpr double f0 = 0.79056941504209488;
  constexpr double f1 = 2.3717082451262845;
  constexpr double f2 = 1.9364916731037085;
  constexpr double f3 = 3.872983346207417;
  constexpr double f4 = 2.4494897427831779;
  constexpr double f5 = 0.61237243569579447;
  constexpr double f6 = 1.5;
  constexpr double g0 = 0.73950997288745202;
  constexpr double g1 = 4.4370598373247123;
  constexpr double g2 = 2.9580398915498081;
  constexpr double g3 = 2.0916500663351889;
  constexpr double g4 = 6.2749501990055663;
  constexpr double g5 = 3.3541019662496847;
  constexpr double g6 = 0.55901699437494745;
  constexpr double g7 = 6.7082039324993694;
  constexpr double g8 = 1.1180339887498949;
  constexpr double g9 = 3.1622776601683795;
  constexpr double g10 = 2.3717082451262845;
  constexpr double g11 = 3;
  constexpr double g12 = 0.375;
  constexpr double g13 = 0.75;
  constexpr double h0 = 0.70156076002011403;
  constexpr double h1 = 7.0156076002011405;
  constexpr double h2 = 3.5078038001005702;
  constexpr double h3 = 2.2185299186623562;
  constexpr double h4 = 13.311179511974137;
  constexpr double h5 = 8.8741196746494246;
  constexpr double h6 = 4.1833001326703778;
  constexpr double h7 = 12.549900398011133;
  constexpr double h8 = 0.52291251658379723;
  constexpr double h9 = 1.5687375497513916;
  constexpr double h10 = 5.123475382979799;
  constexpr double h11 = 2.5617376914898995;
  constexpr double h12 = 10.246950765959598;
  constexpr double h13 = 3.872983346207417;
  constexpr double h14 = 5.8094750193111251;
  constexpr double h15 = 0.48412291827592713;
  constexpr double h16 = 0.96824583655185426;
  constexpr double h17 = 5;
  constexpr double h18 = 1.875;
  constexpr double h19 = 3.75;
  constexpr double i0 = 0.67169328938139616;
  constexpr double i1 = 10.075399340720942;
  constexpr double i2 = 4.0301597362883772;
  constexpr double i3 = 13.433865787627923;
  constexpr double i4 = 2.3268138086232857;
  constexpr double i5 = 23.268138086232856;
  constexpr double i6 = 11.634069043116428;
  constexpr double i7 = 4.9607837082461073;
  constexpr double i8 = 29.764702249476645;
  constexpr double i9 = 0.49607837082461076;
  constexpr double i10 = 2.9764702249476644;
  constexpr double i11 = 19.843134832984429;
  constexpr double i12 = 1.984313483298443;
  constexpr double i13 = 7.245688373094719;
  constexpr double i14 = 21.737065119284157;
  constexpr double i15 = 2.7171331399105196;
  constexpr double i16 = 8.1513994197315593;
  constexpr double i17 = 0.45285552331841994;
  constexpr double i18 = 0.90571104663683988;
  constexpr double i19 = 14.491376746189438;
  constexpr double i20 = 1.8114220932736798;
  constexpr double i21 = 4.5825756949558398;
  constexpr double i22 = 11.456439237389599;
  constexpr double i23 = 2.8641098093473998;
  constexpr double i24 = 5.7282196186947996;
  constexpr double i25 = 7.5;
  constexpr double i26 = 5.625;
  constexpr double i27 = 11.25;
  constexpr double i28 = 0.3125;
  constexpr double i29 = 0.9375;

  static constexpr std::array<double,1> css = { 1.0 };
  static constexpr std::array<double,9> csp = { 1.0, 0.0, 0.0,
                                           0.0, 1.0, 0.0,
                                           0.0, 0.0, 1.0 };
  static constexpr std::array<double,30> csd = {  d0, 0.0, -d0, 0.0, 0.0, 0.0,
                                            0.0,  d1, 0.0, 0.0, 0.0, 0.0,
                                            0.0, 0.0, 0.0,  d1, 0.0, 0.0,
                                            0.0, 0.0, 0.0, 0.0,  d1, 0.0,
                                            -d2, 0.0, -d2, 0.0, 0.0, 1.0 };
  static constexpr std::array<double,170> csf = { f0, 0.0, -f1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                            0.0,  f1, 0.0, -f0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                            0.0, 0.0, 0.0, 0.0,  f2, 0.0, -f2, 0.0, 0.0, 0.0,
                                            0.0, 0.0, 0.0, 0.0, 0.0,  f3, 0.0, 0.0, 0.0, 0.0,
                                            -f5, 0.0, -f5, 0.0, 0.0, 0.0, 0.0,  f4, 0.0, 0.0,
                                            0.0, -f5, 0.0, -f5, 0.0, 0.0, 0.0, 0.0,  f4, 0.0,
                                            0.0, 0.0, 0.0, 0.0, -f6, 0.0, -f6, 0.0, 0.0, 1.0 };
  static constexpr std::array<double,135> csg = {  g0, 0.0, -g1, 0.0,  g0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                             0.0,  g2, 0.0, -g2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                             0.0, 0.0, 0.0, 0.0, 0.0,  g3, 0.0, -g4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  g4, 0.0, -g3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                             -g6, 0.0, 0.0, 0.0,  g6, 0.0, 0.0, 0.0, 0.0,  g5, 0.0, -g5, 0.0, 0.0, 0.0,
                                             0.0, -g8, 0.0, -g8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  g7, 0.0, 0.0, 0.0, 0.0,
                                             0.0, 0.0, 0.0, 0.0, 0.0, -g10, 0.0, -g10, 0.0, 0.0, 0.0, 0.0,  g9, 0.0, 0.0,
                                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -g10, 0.0, -g10, 0.0, 0.0, 0.0, 0.0,  g9, 0.0,
                                             g12, 0.0, g13, 0.0, g12, 0.0, 0.0, 0.0, 0.0, -g11, 0.0, -g11, 0.0, 0.0, 1.0 };
  static constexpr std::array<double,231> csh = {  h0, 0.0, -h1, 0.0,  h2, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0,
                                             0.0,  h2, 0.0, -h1, 0.0,  h0,  0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0,
                                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   h3, 0.0, -h4, 0.0,  h3,  0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0,
                                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0,  h5, 0.0, -h5, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0,
                                             -h8, 0.0, (h9-h8), 0.0,  h9, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0,   h6, 0.0, -h7, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0,
                                             0.0, -h9, 0.0, (h8-h9), 0.0,  h8,  0.0, 0.0, 0.0, 0.0, 0.0,  0.0,  h7, 0.0, -h6, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0,
                                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  -h11, 0.0, 0.0, 0.0, h11,  0.0, 0.0, 0.0, 0.0, h10,  0.0, -h10, 0.0, 0.0, 0.0,
                                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0, -h10, 0.0, -h10, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0,  h12, 0.0, 0.0, 0.0, 0.0,
                                             h15, 0.0, h16, 0.0, h15, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0,  -h14, 0.0, -h14, 0.0, 0.0,  0.0, 0.0, h13, 0.0, 0.0,
                                             0.0, h15, 0.0, h16, 0.0, h15,  0.0, 0.0, 0.0, 0.0, 0.0,  0.0, -h14, 0.0, -h14, 0.0,  0.0, 0.0, 0.0, h13, 0.0,
                                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  h18, 0.0, h19, 0.0, h18,  0.0, 0.0, 0.0, 0.0, -h17,  0.0, -h17, 0.0, 0.0, 1.0 };
  static constexpr std::array<double,364> csi = {  i0, 0.0, -i1, 0.0,  i1, 0.0,  -i0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                             0.0,  i2, 0.0, -i3, 0.0,  i2,  0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0,  i4, 0.0, -i5, 0.0,   i6, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0,  i6, 0.0, -i5,  0.0,  i4, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                             -i9, 0.0, (i10-i9), 0.0, (i10-i9), 0.0,  -i9, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0,  i7, 0.0, -i8,  0.0,  i7, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                             0.0, -i12, 0.0, 0.0, 0.0, i12,  0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, i11, 0.0,  -i11, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0, -i15, 0.0, (i16-i15), 0.0,  i16, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, i13, 0.0, -i14,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, -i16, 0.0, (i15-i16),  0.0, i15, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, i14, 0.0,  -i13, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                             i17, 0.0, (i18-i17), 0.0, (i17-i18), 0.0,  -i17, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, -i13, 0.0, 0.0,  0.0, i13, 0.0, 0.0, 0.0,  0.0, i13, 0.0, -i13, 0.0, 0.0, 0.0,
                                             0.0, i18, 0.0, i20, 0.0, i18,  0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, -i19, 0.0,  -i19, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, i19, 0.0, 0.0, 0.0, 0.0,
                                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0, i23, 0.0, i24, 0.0,  i23, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, -i22, 0.0, -i22,  0.0, 0.0, 0.0, 0.0, i21, 0.0, 0.0,
                                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, i23, 0.0, i24,  0.0, i23, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, -i22, 0.0,  -i22, 0.0, 0.0, 0.0, 0.0, i21, 0.0,
                                             -i28, 0.0, -i29, 0.0, -i29, 0.0,  -i28, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, i26, 0.0, i27,  0.0, i26, 0.0, 0.0, 0.0,  0.0, -i25, 0.0, -i25, 0.0, 0.0, 1.0 };
  switch (i) {
    case 0:
      assert(m*n == 1);
      out->copy_block(0, 0, m, n, css.data());
      break;
    case 1:
      assert(m*n == 9);
      out->copy_block(0, 0, m, n, csp.data());
      break;
    case 2:
      assert(m*n == 30);
      out->copy_block(0, 0, m, n, csd.data());
      break;
    case 3:
      assert(m*n == 70);
      out->copy_block(0, 0, m, n, csf.data());
      break;
    case 4:
      assert(m*n == 135);
      out->copy_block(0, 0, m, n, csg.data());
      break;
    case 5:
      assert(m*n == 231);
      out->copy_block(0, 0, m, n, csh.data());
      break;
    case 6:
      assert(m*n == 364);
#ifndef COMPILE_J_ORB
      throw std::runtime_error("Relativistic calculations with i-type orbital basis functions require j-type integrals for the small component.  Recompile with -DCOMPILE_J_ORB to use this feature.");
#endif
      out->copy_block(0, 0, m, n, csi.data());
      break;
    case 7:
      assert(m*n == 540);
      throw std::runtime_error("Relativistic calculations cannot use j-type orbital basis functions.  (k-type would be needed for the small component.)");
      break;
    default:
      throw std::runtime_error("Angular momentum index not recognized");
  }

  return out;
}

}

#endif
