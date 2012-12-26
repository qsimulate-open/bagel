//
// BAGEL - Parallel electron correlation program.
// Filename: vrr.cc
// Copyright (C) 2009 Toru Shiozaki
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

#include <src/rysint/eribatch.h>
#include <src/rysint/_vrr_drv.h>

using namespace std;
using namespace bagel;


void ERIBatch::perform_VRR() {
#ifndef LIBINT_INTERFACE
  const int acsize = asize_ * csize_;
  const int a = basisinfo_[0]->angular_number();
  const int b = basisinfo_[1]->angular_number();
  const int c = basisinfo_[2]->angular_number();
  const int d = basisinfo_[3]->angular_number();
  const int isize = (amax_+1) * (cmax_+1); 
  double* const workx = stack_->get(isize*rank_*3);
  double* const worky = workx + isize*rank_;
  double* const workz = worky + isize*rank_;
  const int hashkey = (a << 24) + (b << 16) + (c << 8) + d;
  switch (hashkey) {
  case 0 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_0_0_0_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 256 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_0_0_1_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 257 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_0_0_1_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 512 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_0_0_2_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 513 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_0_0_2_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 514 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_0_0_2_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 768 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_0_0_3_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 769 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_0_0_3_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 770 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_0_0_3_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 771 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_0_0_3_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 1024 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_0_0_4_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 1025 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_0_0_4_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 1026 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_0_0_4_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 1027 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_0_0_4_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 1028 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_0_0_4_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 1280 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_0_0_5_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 1281 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_0_0_5_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 1282 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_0_0_5_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 1283 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_0_0_5_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 1284 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_0_0_5_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 1285 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_0_0_5_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 1536 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_0_0_6_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 1537 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_0_0_6_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 1538 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_0_0_6_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 1539 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_0_0_6_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 1540 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_0_0_6_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 1541 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_0_0_6_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 1542 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_0_0_6_6(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16777216 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_0_0_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16777472 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_0_1_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16777473 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_0_1_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16777728 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_0_2_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16777729 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_0_2_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16777730 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_0_2_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16777984 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_0_3_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16777985 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_0_3_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16777986 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_0_3_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16777987 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_0_3_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16778240 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_0_4_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16778241 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_0_4_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16778242 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_0_4_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16778243 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_0_4_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16778244 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_0_4_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16778496 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_0_5_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16778497 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_0_5_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16778498 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_0_5_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16778499 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_0_5_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16778500 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_0_5_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16778501 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_0_5_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16778752 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_0_6_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16778753 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_0_6_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16778754 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_0_6_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16778755 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_0_6_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16778756 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_0_6_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16778757 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_0_6_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16778758 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_0_6_6(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16842752 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_1_0_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16843008 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_1_1_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16843009 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_1_1_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16843264 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_1_2_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16843265 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_1_2_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16843266 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_1_2_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16843520 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_1_3_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16843521 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_1_3_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16843522 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_1_3_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16843523 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_1_3_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16843776 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_1_4_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16843777 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_1_4_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16843778 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_1_4_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16843779 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_1_4_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16843780 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_1_4_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16844032 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_1_5_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16844033 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_1_5_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16844034 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_1_5_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16844035 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_1_5_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16844036 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_1_5_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16844037 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_1_5_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16844288 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_1_6_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16844289 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_1_6_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16844290 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_1_6_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16844291 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_1_6_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16844292 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_1_6_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16844293 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_1_6_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 16844294 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_1_1_6_6(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33554432 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_0_0_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33554688 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_0_1_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33554689 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_0_1_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33554944 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_0_2_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33554945 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_0_2_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33554946 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_0_2_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33555200 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_0_3_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33555201 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_0_3_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33555202 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_0_3_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33555203 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_0_3_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33555456 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_0_4_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33555457 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_0_4_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33555458 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_0_4_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33555459 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_0_4_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33555460 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_0_4_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33555712 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_0_5_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33555713 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_0_5_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33555714 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_0_5_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33555715 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_0_5_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33555716 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_0_5_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33555717 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_0_5_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33555968 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_0_6_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33555969 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_0_6_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33555970 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_0_6_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33555971 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_0_6_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33555972 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_0_6_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33555973 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_0_6_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33555974 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_0_6_6(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33619968 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_1_0_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33620224 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_1_1_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33620225 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_1_1_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33620480 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_1_2_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33620481 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_1_2_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33620482 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_1_2_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33620736 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_1_3_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33620737 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_1_3_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33620738 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_1_3_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33620739 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_1_3_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33620992 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_1_4_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33620993 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_1_4_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33620994 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_1_4_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33620995 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_1_4_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33620996 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_1_4_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33621248 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_1_5_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33621249 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_1_5_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33621250 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_1_5_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33621251 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_1_5_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33621252 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_1_5_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33621253 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_1_5_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33621504 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_1_6_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33621505 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_1_6_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33621506 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_1_6_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33621507 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_1_6_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33621508 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_1_6_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33621509 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_1_6_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33621510 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_1_6_6(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33685504 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_2_0_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33685760 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_2_1_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33685761 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_2_1_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33686016 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_2_2_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33686017 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_2_2_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33686018 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_2_2_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33686272 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_2_3_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33686273 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_2_3_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33686274 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_2_3_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33686275 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_2_3_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33686528 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_2_4_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33686529 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_2_4_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33686530 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_2_4_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33686531 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_2_4_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33686532 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_2_4_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33686784 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_2_5_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33686785 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_2_5_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33686786 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_2_5_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33686787 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_2_5_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33686788 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_2_5_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33686789 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_2_5_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33687040 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_2_6_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33687041 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_2_6_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33687042 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_2_6_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33687043 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_2_6_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33687044 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_2_6_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33687045 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_2_6_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 33687046 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_2_2_6_6(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50331648 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_0_0_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50331904 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_0_1_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50331905 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_0_1_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50332160 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_0_2_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50332161 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_0_2_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50332162 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_0_2_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50332416 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_0_3_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50332417 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_0_3_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50332418 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_0_3_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50332419 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_0_3_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50332672 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_0_4_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50332673 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_0_4_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50332674 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_0_4_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50332675 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_0_4_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50332676 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_0_4_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50332928 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_0_5_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50332929 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_0_5_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50332930 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_0_5_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50332931 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_0_5_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50332932 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_0_5_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50332933 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_0_5_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50333184 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_0_6_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50333185 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_0_6_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50333186 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_0_6_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50333187 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_0_6_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50333188 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_0_6_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50333189 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_0_6_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50333190 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_0_6_6(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50397184 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_1_0_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50397440 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_1_1_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50397441 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_1_1_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50397696 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_1_2_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50397697 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_1_2_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50397698 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_1_2_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50397952 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_1_3_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50397953 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_1_3_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50397954 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_1_3_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50397955 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_1_3_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50398208 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_1_4_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50398209 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_1_4_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50398210 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_1_4_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50398211 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_1_4_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50398212 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_1_4_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50398464 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_1_5_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50398465 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_1_5_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50398466 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_1_5_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50398467 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_1_5_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50398468 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_1_5_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50398469 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_1_5_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50398720 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_1_6_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50398721 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_1_6_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50398722 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_1_6_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50398723 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_1_6_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50398724 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_1_6_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50398725 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_1_6_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50398726 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_1_6_6(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50462720 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_2_0_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50462976 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_2_1_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50462977 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_2_1_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50463232 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_2_2_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50463233 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_2_2_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50463234 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_2_2_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50463488 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_2_3_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50463489 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_2_3_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50463490 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_2_3_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50463491 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_2_3_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50463744 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_2_4_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50463745 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_2_4_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50463746 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_2_4_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50463747 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_2_4_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50463748 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_2_4_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50464000 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_2_5_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50464001 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_2_5_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50464002 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_2_5_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50464003 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_2_5_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50464004 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_2_5_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50464005 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_2_5_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50464256 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_2_6_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50464257 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_2_6_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50464258 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_2_6_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50464259 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_2_6_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50464260 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_2_6_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50464261 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_2_6_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50464262 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_2_6_6(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50528256 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_3_0_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50528512 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_3_1_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50528513 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_3_1_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50528768 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_3_2_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50528769 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_3_2_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50528770 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_3_2_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50529024 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_3_3_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50529025 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_3_3_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50529026 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_3_3_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50529027 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_3_3_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50529280 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_3_4_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50529281 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_3_4_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50529282 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_3_4_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50529283 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_3_4_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50529284 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_3_4_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50529536 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_3_5_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50529537 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_3_5_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50529538 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_3_5_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50529539 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_3_5_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50529540 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_3_5_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50529541 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_3_5_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50529792 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_3_6_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50529793 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_3_6_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50529794 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_3_6_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50529795 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_3_6_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50529796 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_3_6_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50529797 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_3_6_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 50529798 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_3_3_6_6(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67108864 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_0_0_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67109120 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_0_1_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67109121 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_0_1_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67109376 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_0_2_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67109377 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_0_2_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67109378 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_0_2_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67109632 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_0_3_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67109633 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_0_3_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67109634 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_0_3_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67109635 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_0_3_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67109888 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_0_4_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67109889 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_0_4_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67109890 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_0_4_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67109891 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_0_4_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67109892 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_0_4_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67110144 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_0_5_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67110145 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_0_5_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67110146 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_0_5_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67110147 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_0_5_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67110148 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_0_5_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67110149 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_0_5_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67110400 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_0_6_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67110401 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_0_6_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67110402 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_0_6_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67110403 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_0_6_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67110404 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_0_6_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67110405 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_0_6_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67110406 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_0_6_6(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67174400 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_1_0_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67174656 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_1_1_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67174657 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_1_1_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67174912 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_1_2_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67174913 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_1_2_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67174914 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_1_2_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67175168 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_1_3_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67175169 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_1_3_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67175170 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_1_3_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67175171 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_1_3_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67175424 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_1_4_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67175425 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_1_4_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67175426 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_1_4_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67175427 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_1_4_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67175428 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_1_4_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67175680 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_1_5_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67175681 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_1_5_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67175682 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_1_5_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67175683 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_1_5_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67175684 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_1_5_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67175685 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_1_5_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67175936 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_1_6_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67175937 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_1_6_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67175938 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_1_6_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67175939 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_1_6_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67175940 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_1_6_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67175941 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_1_6_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67175942 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_1_6_6(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67239936 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_2_0_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67240192 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_2_1_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67240193 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_2_1_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67240448 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_2_2_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67240449 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_2_2_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67240450 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_2_2_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67240704 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_2_3_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67240705 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_2_3_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67240706 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_2_3_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67240707 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_2_3_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67240960 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_2_4_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67240961 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_2_4_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67240962 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_2_4_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67240963 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_2_4_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67240964 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_2_4_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67241216 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_2_5_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67241217 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_2_5_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67241218 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_2_5_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67241219 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_2_5_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67241220 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_2_5_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67241221 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_2_5_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67241472 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_2_6_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67241473 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_2_6_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67241474 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_2_6_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67241475 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_2_6_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67241476 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_2_6_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67241477 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_2_6_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67241478 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_2_6_6(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67305472 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_3_0_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67305728 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_3_1_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67305729 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_3_1_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67305984 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_3_2_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67305985 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_3_2_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67305986 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_3_2_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67306240 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_3_3_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67306241 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_3_3_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67306242 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_3_3_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67306243 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_3_3_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67306496 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_3_4_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67306497 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_3_4_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67306498 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_3_4_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67306499 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_3_4_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67306500 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_3_4_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67306752 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_3_5_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67306753 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_3_5_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67306754 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_3_5_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67306755 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_3_5_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67306756 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_3_5_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67306757 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_3_5_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67307008 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_3_6_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67307009 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_3_6_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67307010 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_3_6_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67307011 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_3_6_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67307012 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_3_6_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67307013 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_3_6_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67307014 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_3_6_6(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67371008 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_4_0_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67371264 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_4_1_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67371265 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_4_1_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67371520 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_4_2_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67371521 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_4_2_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67371522 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_4_2_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67371776 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_4_3_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67371777 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_4_3_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67371778 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_4_3_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67371779 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_4_3_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67372032 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_4_4_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67372033 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_4_4_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67372034 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_4_4_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67372035 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_4_4_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67372036 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_4_4_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67372288 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_4_5_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67372289 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_4_5_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67372290 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_4_5_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67372291 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_4_5_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67372292 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_4_5_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67372293 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_4_5_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67372544 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_4_6_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67372545 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_4_6_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67372546 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_4_6_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67372547 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_4_6_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67372548 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_4_6_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67372549 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_4_6_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 67372550 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_4_4_6_6(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83886080 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_0_0_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83886336 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_0_1_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83886337 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_0_1_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83886592 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_0_2_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83886593 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_0_2_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83886594 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_0_2_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83886848 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_0_3_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83886849 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_0_3_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83886850 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_0_3_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83886851 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_0_3_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83887104 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_0_4_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83887105 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_0_4_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83887106 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_0_4_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83887107 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_0_4_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83887108 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_0_4_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83887360 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_0_5_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83887361 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_0_5_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83887362 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_0_5_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83887363 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_0_5_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83887364 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_0_5_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83887365 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_0_5_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83887616 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_0_6_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83887617 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_0_6_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83887618 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_0_6_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83887619 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_0_6_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83887620 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_0_6_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83887621 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_0_6_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83887622 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_0_6_6(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83951616 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_1_0_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83951872 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_1_1_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83951873 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_1_1_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83952128 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_1_2_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83952129 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_1_2_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83952130 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_1_2_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83952384 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_1_3_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83952385 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_1_3_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83952386 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_1_3_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83952387 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_1_3_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83952640 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_1_4_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83952641 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_1_4_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83952642 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_1_4_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83952643 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_1_4_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83952644 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_1_4_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83952896 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_1_5_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83952897 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_1_5_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83952898 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_1_5_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83952899 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_1_5_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83952900 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_1_5_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83952901 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_1_5_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83953152 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_1_6_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83953153 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_1_6_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83953154 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_1_6_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83953155 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_1_6_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83953156 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_1_6_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83953157 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_1_6_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 83953158 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_1_6_6(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84017152 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_2_0_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84017408 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_2_1_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84017409 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_2_1_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84017664 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_2_2_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84017665 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_2_2_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84017666 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_2_2_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84017920 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_2_3_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84017921 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_2_3_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84017922 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_2_3_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84017923 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_2_3_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84018176 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_2_4_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84018177 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_2_4_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84018178 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_2_4_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84018179 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_2_4_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84018180 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_2_4_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84018432 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_2_5_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84018433 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_2_5_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84018434 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_2_5_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84018435 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_2_5_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84018436 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_2_5_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84018437 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_2_5_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84018688 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_2_6_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84018689 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_2_6_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84018690 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_2_6_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84018691 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_2_6_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84018692 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_2_6_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84018693 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_2_6_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84018694 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_2_6_6(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84082688 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_3_0_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84082944 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_3_1_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84082945 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_3_1_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84083200 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_3_2_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84083201 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_3_2_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84083202 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_3_2_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84083456 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_3_3_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84083457 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_3_3_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84083458 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_3_3_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84083459 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_3_3_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84083712 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_3_4_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84083713 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_3_4_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84083714 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_3_4_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84083715 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_3_4_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84083716 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_3_4_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84083968 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_3_5_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84083969 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_3_5_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84083970 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_3_5_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84083971 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_3_5_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84083972 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_3_5_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84083973 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_3_5_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84084224 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_3_6_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84084225 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_3_6_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84084226 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_3_6_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84084227 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_3_6_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84084228 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_3_6_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84084229 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_3_6_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84084230 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_3_6_6(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84148224 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_4_0_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84148480 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_4_1_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84148481 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_4_1_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84148736 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_4_2_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84148737 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_4_2_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84148738 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_4_2_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84148992 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_4_3_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84148993 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_4_3_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84148994 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_4_3_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84148995 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_4_3_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84149248 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_4_4_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84149249 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_4_4_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84149250 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_4_4_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84149251 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_4_4_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84149252 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_4_4_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84149504 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_4_5_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84149505 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_4_5_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84149506 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_4_5_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84149507 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_4_5_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84149508 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_4_5_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84149509 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_4_5_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84149760 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_4_6_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84149761 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_4_6_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84149762 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_4_6_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84149763 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_4_6_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84149764 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_4_6_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84149765 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_4_6_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84149766 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_4_6_6(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84213760 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_5_0_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84214016 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_5_1_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84214017 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_5_1_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84214272 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_5_2_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84214273 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_5_2_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84214274 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_5_2_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84214528 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_5_3_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84214529 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_5_3_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84214530 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_5_3_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84214531 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_5_3_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84214784 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_5_4_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84214785 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_5_4_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84214786 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_5_4_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84214787 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_5_4_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84214788 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_5_4_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84215040 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_5_5_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84215041 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_5_5_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84215042 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_5_5_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84215043 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_5_5_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84215044 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_5_5_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84215045 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_5_5_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84215296 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_5_6_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84215297 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_5_6_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84215298 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_5_6_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84215299 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_5_6_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84215300 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_5_6_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84215301 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_5_6_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 84215302 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_5_5_6_6(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100663296 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_0_0_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100663552 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_0_1_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100663553 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_0_1_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100663808 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_0_2_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100663809 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_0_2_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100663810 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_0_2_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100664064 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_0_3_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100664065 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_0_3_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100664066 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_0_3_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100664067 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_0_3_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100664320 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_0_4_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100664321 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_0_4_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100664322 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_0_4_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100664323 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_0_4_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100664324 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_0_4_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100664576 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_0_5_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100664577 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_0_5_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100664578 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_0_5_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100664579 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_0_5_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100664580 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_0_5_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100664581 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_0_5_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100664832 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_0_6_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100664833 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_0_6_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100664834 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_0_6_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100664835 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_0_6_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100664836 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_0_6_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100664837 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_0_6_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100664838 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_0_6_6(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100728832 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_1_0_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100729088 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_1_1_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100729089 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_1_1_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100729344 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_1_2_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100729345 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_1_2_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100729346 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_1_2_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100729600 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_1_3_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100729601 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_1_3_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100729602 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_1_3_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100729603 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_1_3_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100729856 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_1_4_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100729857 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_1_4_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100729858 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_1_4_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100729859 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_1_4_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100729860 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_1_4_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100730112 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_1_5_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100730113 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_1_5_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100730114 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_1_5_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100730115 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_1_5_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100730116 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_1_5_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100730117 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_1_5_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100730368 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_1_6_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100730369 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_1_6_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100730370 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_1_6_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100730371 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_1_6_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100730372 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_1_6_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100730373 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_1_6_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100730374 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_1_6_6(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100794368 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_2_0_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100794624 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_2_1_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100794625 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_2_1_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100794880 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_2_2_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100794881 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_2_2_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100794882 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_2_2_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100795136 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_2_3_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100795137 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_2_3_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100795138 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_2_3_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100795139 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_2_3_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100795392 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_2_4_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100795393 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_2_4_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100795394 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_2_4_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100795395 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_2_4_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100795396 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_2_4_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100795648 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_2_5_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100795649 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_2_5_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100795650 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_2_5_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100795651 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_2_5_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100795652 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_2_5_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100795653 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_2_5_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100795904 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_2_6_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100795905 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_2_6_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100795906 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_2_6_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100795907 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_2_6_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100795908 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_2_6_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100795909 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_2_6_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100795910 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_2_6_6(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100859904 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_3_0_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100860160 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_3_1_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100860161 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_3_1_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100860416 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_3_2_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100860417 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_3_2_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100860418 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_3_2_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100860672 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_3_3_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100860673 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_3_3_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100860674 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_3_3_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100860675 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_3_3_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100860928 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_3_4_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100860929 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_3_4_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100860930 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_3_4_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100860931 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_3_4_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100860932 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_3_4_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100861184 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_3_5_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100861185 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_3_5_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100861186 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_3_5_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100861187 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_3_5_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100861188 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_3_5_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100861189 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_3_5_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100861440 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_3_6_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100861441 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_3_6_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100861442 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_3_6_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100861443 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_3_6_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100861444 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_3_6_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100861445 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_3_6_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100861446 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_3_6_6(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100925440 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_4_0_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100925696 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_4_1_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100925697 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_4_1_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100925952 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_4_2_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100925953 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_4_2_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100925954 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_4_2_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100926208 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_4_3_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100926209 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_4_3_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100926210 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_4_3_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100926211 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_4_3_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100926464 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_4_4_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100926465 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_4_4_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100926466 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_4_4_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100926467 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_4_4_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100926468 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_4_4_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100926720 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_4_5_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100926721 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_4_5_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100926722 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_4_5_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100926723 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_4_5_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100926724 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_4_5_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100926725 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_4_5_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100926976 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_4_6_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100926977 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_4_6_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100926978 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_4_6_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100926979 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_4_6_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100926980 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_4_6_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100926981 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_4_6_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100926982 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_4_6_6(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100990976 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_5_0_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100991232 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_5_1_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100991233 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_5_1_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100991488 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_5_2_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100991489 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_5_2_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100991490 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_5_2_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100991744 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_5_3_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100991745 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_5_3_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100991746 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_5_3_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100991747 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_5_3_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100992000 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_5_4_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100992001 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_5_4_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100992002 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_5_4_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100992003 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_5_4_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100992004 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_5_4_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100992256 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_5_5_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100992257 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_5_5_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100992258 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_5_5_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100992259 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_5_5_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100992260 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_5_5_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100992261 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_5_5_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100992512 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_5_6_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100992513 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_5_6_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100992514 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_5_6_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100992515 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_5_6_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100992516 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_5_6_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100992517 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_5_6_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 100992518 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_5_6_6(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 101056512 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_6_0_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 101056768 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_6_1_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 101056769 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_6_1_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 101057024 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_6_2_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 101057025 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_6_2_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 101057026 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_6_2_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 101057280 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_6_3_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 101057281 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_6_3_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 101057282 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_6_3_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 101057283 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_6_3_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 101057536 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_6_4_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 101057537 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_6_4_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 101057538 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_6_4_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 101057539 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_6_4_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 101057540 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_6_4_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 101057792 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_6_5_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 101057793 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_6_5_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 101057794 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_6_5_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 101057795 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_6_5_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 101057796 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_6_5_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 101057797 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_6_5_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 101058048 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_6_6_0(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 101058049 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_6_6_1(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 101058050 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_6_6_2(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 101058051 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_6_6_3(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 101058052 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_6_6_4(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 101058053 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_6_6_5(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  case 101058054 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      VRR_Driver::vrr_driver_6_6_6_6(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz);
    } break;
  }
  stack_->release(rank_*isize*3, workx);

#endif
}