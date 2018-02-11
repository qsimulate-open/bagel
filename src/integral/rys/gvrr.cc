//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gvrr.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <src/integral/rys/gradbatch.h>
#include <src/integral/rys/_gvrr_drv.h>
#include <src/util/math/comb.h>
#include <src/util/f77.h>

using namespace std;
using namespace bagel;

static const Comb comb;


void GradBatch::perform_VRR() {
#ifndef LIBINT_INTERFACE
  const int a = basisinfo_[0]->angular_number();
  const int b = basisinfo_[1]->angular_number();
  const int c = basisinfo_[2]->angular_number();
  const int d = basisinfo_[3]->angular_number();
  const int acsize = (a+1)*(a+2)*(b+1)*(b+2)*(c+1)*(c+2)*(d+1)*(d+2)/16;

  const int isize = (amax_ + 1) * (cmax_ + 1);
  double* const workx = stack_->get(isize*rank_*3);
  double* const worky = workx + isize*rank_;
  double* const workz = worky + isize*rank_;

  const int a2 = a+2;
  const int b2 = b+2;
  const int c2 = c+2;
  const int d2 = d+2;

  double* const transx = stack_->get((amax_+1)*a2*b2);
  double* const transy = stack_->get((amax_+1)*a2*b2);
  double* const transz = stack_->get((amax_+1)*a2*b2);
  double* const trans2x = stack_->get((cmax_+1)*c2*d2);
  double* const trans2y = stack_->get((cmax_+1)*c2*d2);
  double* const trans2z = stack_->get((cmax_+1)*c2*d2);
  fill_n(transx,  (amax_+1)*a2*b2, 0.0);
  fill_n(transy,  (amax_+1)*a2*b2, 0.0);
  fill_n(transz,  (amax_+1)*a2*b2, 0.0);
  fill_n(trans2x, (cmax_+1)*c2*d2, 0.0);
  fill_n(trans2y, (cmax_+1)*c2*d2, 0.0);
  fill_n(trans2z, (cmax_+1)*c2*d2, 0.0);
  // for usual integrals
  for (int ib = 0, k = 0; ib <= b+1; ++ib) {
    for (int ia = 0; ia <= a+1; ++ia, ++k) {
      if (ia == a+1 && ib == b+1) continue;
      for (int i = ia; i <= ia+ib; ++i) {
        transx[i + (amax_+1)*k] = comb(ib, ia+ib-i) * pow(AB_[0], ia+ib-i);
        transy[i + (amax_+1)*k] = comb(ib, ia+ib-i) * pow(AB_[1], ia+ib-i);
        transz[i + (amax_+1)*k] = comb(ib, ia+ib-i) * pow(AB_[2], ia+ib-i);
      }
    }
  }
  for (int id = 0, k = 0; id <= d+1; ++id) {
    for (int ic = 0; ic <= c+1; ++ic, ++k) {
      if (ic == c+1 && id == d+1) continue;
      for (int i = ic; i <= ic+id; ++i) {
        trans2x[i + (cmax_+1)*k] = comb(id, ic+id-i) * pow(CD_[0], ic+id-i);
        trans2y[i + (cmax_+1)*k] = comb(id, ic+id-i) * pow(CD_[1], ic+id-i);
        trans2z[i + (cmax_+1)*k] = comb(id, ic+id-i) * pow(CD_[2], ic+id-i);
      }
    }
  }
  double* const intermediate = stack_->get(b2*a2*(cmax_+1)*rank_);
  double* const final_x  = stack_->get(b2*a2*c2*d2*rank_);
  double* const final_y  = stack_->get(b2*a2*c2*d2*rank_);
  double* const final_z  = stack_->get(b2*a2*c2*d2*rank_);
  double* const final_xa = stack_->get(b2*a2*c2*d2*rank_);
  double* const final_xb = stack_->get(b2*a2*c2*d2*rank_);
  double* const final_xc = stack_->get(b2*a2*c2*d2*rank_);
  double* const final_ya = stack_->get(b2*a2*c2*d2*rank_);
  double* const final_yb = stack_->get(b2*a2*c2*d2*rank_);
  double* const final_yc = stack_->get(b2*a2*c2*d2*rank_);
  double* const final_za = stack_->get(b2*a2*c2*d2*rank_);
  double* const final_zb = stack_->get(b2*a2*c2*d2*rank_);
  double* const final_zc = stack_->get(b2*a2*c2*d2*rank_);
  const array<bool,4> dummy{{basisinfo_[0]->dummy(), basisinfo_[1]->dummy(), basisinfo_[2]->dummy(), basisinfo_[3]->dummy()}};
  const int hashkey = (a << 24) + (b << 16) + (c << 8) + d;
  switch (hashkey) {
  case 0 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,0,0,1>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 256 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,1,0,2>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 257 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,1,1,2>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 512 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,2,0,2>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 513 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,2,1,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 514 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,2,2,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 768 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,3,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 769 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,3,1,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 770 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,3,2,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 771 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,3,3,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 1024 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,4,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 1025 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,4,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 1026 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,4,2,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 1027 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,4,3,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 1028 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,4,4,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 1280 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,5,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 1281 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,5,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 1282 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,5,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 1283 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,5,3,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 1284 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,5,4,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 1285 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,5,5,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 1536 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,6,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 1537 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,6,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 1538 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,6,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 1539 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,6,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 1540 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,6,4,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 1541 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,6,5,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 1542 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,6,6,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#ifdef COMPILE_J_ORB
  case 1792 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,7,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 1793 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,7,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 1794 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,7,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 1795 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,7,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 1796 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,7,4,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 1797 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,7,5,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 1798 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,7,6,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 1799 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,7,7,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
  case 16777216 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,0,0,2>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16777472 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,1,0,2>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16777473 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,1,1,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16777728 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,2,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16777729 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,2,1,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16777730 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,2,2,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16777984 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,3,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16777985 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,3,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16777986 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,3,2,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16777987 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,3,3,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16778240 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,4,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16778241 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,4,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16778242 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,4,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16778243 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,4,3,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16778244 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,4,4,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16778496 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,5,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16778497 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,5,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16778498 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,5,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16778499 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,5,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16778500 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,5,4,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16778501 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,5,5,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16778752 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,6,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16778753 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,6,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16778754 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,6,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16778755 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,6,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16778756 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,6,4,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16778757 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,6,5,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16778758 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,6,6,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#ifdef COMPILE_J_ORB
  case 16779008 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,7,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 16779009 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,7,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 16779010 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,7,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 16779011 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,7,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 16779012 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,7,4,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 16779013 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,7,5,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 16779014 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,7,6,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 16779015 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,7,7,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
  case 16842752 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,0,0,2>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16843008 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,1,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16843009 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,1,1,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16843264 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,2,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16843265 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,2,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16843266 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,2,2,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16843520 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,3,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16843521 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,3,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16843522 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,3,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16843523 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,3,3,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16843776 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,4,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16843777 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,4,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16843778 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,4,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16843779 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,4,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16843780 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,4,4,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16844032 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,5,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16844033 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,5,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16844034 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,5,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16844035 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,5,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16844036 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,5,4,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16844037 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,5,5,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16844288 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,6,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16844289 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,6,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16844290 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,6,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16844291 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,6,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16844292 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,6,4,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16844293 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,6,5,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 16844294 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,6,6,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#ifdef COMPILE_J_ORB
  case 16844544 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,7,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 16844545 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,7,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 16844546 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,7,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 16844547 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,7,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 16844548 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,7,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 16844549 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,7,5,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 16844550 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,7,6,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 16844551 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,7,7,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
  case 33554432 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,0,0,2>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33554688 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,1,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33554689 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,1,1,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33554944 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,2,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33554945 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,2,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33554946 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,2,2,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33555200 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,3,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33555201 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,3,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33555202 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,3,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33555203 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,3,3,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33555456 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,4,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33555457 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,4,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33555458 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,4,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33555459 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,4,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33555460 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,4,4,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33555712 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,5,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33555713 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,5,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33555714 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,5,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33555715 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,5,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33555716 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,5,4,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33555717 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,5,5,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33555968 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,6,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33555969 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,6,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33555970 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,6,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33555971 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,6,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33555972 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,6,4,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33555973 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,6,5,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33555974 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,6,6,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#ifdef COMPILE_J_ORB
  case 33556224 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,7,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 33556225 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,7,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 33556226 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,7,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 33556227 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,7,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 33556228 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,7,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 33556229 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,7,5,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 33556230 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,7,6,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 33556231 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,7,7,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
  case 33619968 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,0,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33620224 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,1,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33620225 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,1,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33620480 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,2,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33620481 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,2,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33620482 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,2,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33620736 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,3,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33620737 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,3,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33620738 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,3,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33620739 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,3,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33620992 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,4,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33620993 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,4,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33620994 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,4,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33620995 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,4,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33620996 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,4,4,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33621248 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,5,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33621249 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,5,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33621250 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,5,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33621251 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,5,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33621252 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,5,4,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33621253 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,5,5,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33621504 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,6,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33621505 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,6,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33621506 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,6,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33621507 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,6,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33621508 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,6,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33621509 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,6,5,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33621510 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,6,6,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#ifdef COMPILE_J_ORB
  case 33621760 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,7,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 33621761 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,7,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 33621762 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,7,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 33621763 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,7,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 33621764 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,7,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 33621765 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,7,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 33621766 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,7,6,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 33621767 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,7,7,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
  case 33685504 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,0,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33685760 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,1,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33685761 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,1,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33686016 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,2,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33686017 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,2,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33686018 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,2,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33686272 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,3,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33686273 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,3,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33686274 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,3,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33686275 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,3,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33686528 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,4,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33686529 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,4,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33686530 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,4,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33686531 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,4,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33686532 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,4,4,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33686784 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,5,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33686785 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,5,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33686786 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,5,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33686787 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,5,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33686788 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,5,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33686789 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,5,5,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33687040 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,6,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33687041 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,6,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33687042 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,6,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33687043 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,6,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33687044 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,6,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33687045 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,6,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 33687046 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,6,6,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#ifdef COMPILE_J_ORB
  case 33687296 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,7,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 33687297 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,7,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 33687298 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,7,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 33687299 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,7,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 33687300 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,7,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 33687301 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,7,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 33687302 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,7,6,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 33687303 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,7,7,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
  case 50331648 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,0,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50331904 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,1,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50331905 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,1,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50332160 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,2,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50332161 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,2,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50332162 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,2,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50332416 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,3,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50332417 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,3,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50332418 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,3,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50332419 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,3,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50332672 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,4,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50332673 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,4,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50332674 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,4,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50332675 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,4,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50332676 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,4,4,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50332928 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,5,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50332929 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,5,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50332930 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,5,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50332931 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,5,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50332932 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,5,4,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50332933 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,5,5,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50333184 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,6,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50333185 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,6,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50333186 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,6,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50333187 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,6,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50333188 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,6,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50333189 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,6,5,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50333190 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,6,6,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#ifdef COMPILE_J_ORB
  case 50333440 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,7,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 50333441 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,7,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 50333442 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,7,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 50333443 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,7,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 50333444 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,7,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 50333445 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,7,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 50333446 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,7,6,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 50333447 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,7,7,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
  case 50397184 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,0,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50397440 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,1,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50397441 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,1,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50397696 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,2,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50397697 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,2,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50397698 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,2,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50397952 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,3,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50397953 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,3,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50397954 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,3,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50397955 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,3,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50398208 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,4,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50398209 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,4,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50398210 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,4,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50398211 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,4,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50398212 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,4,4,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50398464 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,5,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50398465 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,5,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50398466 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,5,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50398467 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,5,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50398468 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,5,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50398469 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,5,5,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50398720 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,6,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50398721 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,6,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50398722 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,6,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50398723 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,6,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50398724 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,6,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50398725 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,6,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50398726 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,6,6,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#ifdef COMPILE_J_ORB
  case 50398976 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,7,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 50398977 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,7,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 50398978 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,7,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 50398979 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,7,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 50398980 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,7,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 50398981 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,7,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 50398982 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,7,6,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 50398983 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,7,7,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
  case 50462720 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,0,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50462976 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,1,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50462977 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,1,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50463232 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,2,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50463233 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,2,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50463234 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,2,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50463488 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,3,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50463489 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,3,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50463490 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,3,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50463491 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,3,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50463744 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,4,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50463745 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,4,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50463746 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,4,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50463747 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,4,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50463748 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,4,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50464000 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,5,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50464001 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,5,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50464002 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,5,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50464003 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,5,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50464004 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,5,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50464005 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,5,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50464256 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,6,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50464257 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,6,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50464258 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,6,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50464259 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,6,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50464260 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,6,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50464261 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,6,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50464262 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,6,6,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#ifdef COMPILE_J_ORB
  case 50464512 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,7,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 50464513 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,7,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 50464514 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,7,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 50464515 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,7,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 50464516 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,7,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 50464517 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,7,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 50464518 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,7,6,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 50464519 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,7,7,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
  case 50528256 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,0,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50528512 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,1,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50528513 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,1,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50528768 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,2,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50528769 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,2,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50528770 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,2,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50529024 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,3,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50529025 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,3,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50529026 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,3,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50529027 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,3,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50529280 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,4,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50529281 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,4,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50529282 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,4,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50529283 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,4,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50529284 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,4,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50529536 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,5,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50529537 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,5,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50529538 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,5,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50529539 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,5,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50529540 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,5,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50529541 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,5,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50529792 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,6,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50529793 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,6,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50529794 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,6,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50529795 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,6,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50529796 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,6,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50529797 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,6,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 50529798 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,6,6,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#ifdef COMPILE_J_ORB
  case 50530048 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,7,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 50530049 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,7,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 50530050 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,7,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 50530051 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,7,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 50530052 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,7,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 50530053 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,7,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 50530054 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,7,6,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 50530055 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,7,7,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
  case 67108864 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,0,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67109120 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,1,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67109121 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,1,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67109376 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,2,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67109377 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,2,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67109378 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,2,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67109632 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,3,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67109633 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,3,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67109634 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,3,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67109635 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,3,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67109888 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,4,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67109889 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,4,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67109890 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,4,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67109891 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,4,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67109892 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,4,4,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67110144 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,5,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67110145 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,5,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67110146 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,5,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67110147 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,5,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67110148 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,5,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67110149 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,5,5,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67110400 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,6,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67110401 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,6,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67110402 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,6,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67110403 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,6,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67110404 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,6,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67110405 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,6,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67110406 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,6,6,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#ifdef COMPILE_J_ORB
  case 67110656 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,7,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67110657 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,7,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67110658 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,7,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67110659 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,7,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67110660 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,7,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67110661 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,7,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67110662 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,7,6,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67110663 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,7,7,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
  case 67174400 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,0,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67174656 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,1,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67174657 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,1,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67174912 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,2,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67174913 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,2,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67174914 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,2,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67175168 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,3,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67175169 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,3,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67175170 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,3,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67175171 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,3,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67175424 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,4,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67175425 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,4,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67175426 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,4,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67175427 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,4,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67175428 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,4,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67175680 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,5,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67175681 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,5,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67175682 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,5,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67175683 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,5,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67175684 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,5,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67175685 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,5,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67175936 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,6,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67175937 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,6,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67175938 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,6,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67175939 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,6,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67175940 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,6,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67175941 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,6,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67175942 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,6,6,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#ifdef COMPILE_J_ORB
  case 67176192 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,7,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67176193 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,7,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67176194 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,7,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67176195 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,7,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67176196 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,7,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67176197 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,7,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67176198 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,7,6,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67176199 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,7,7,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
  case 67239936 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,0,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67240192 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,1,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67240193 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,1,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67240448 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,2,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67240449 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,2,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67240450 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,2,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67240704 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,3,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67240705 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,3,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67240706 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,3,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67240707 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,3,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67240960 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,4,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67240961 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,4,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67240962 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,4,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67240963 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,4,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67240964 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,4,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67241216 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,5,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67241217 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,5,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67241218 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,5,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67241219 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,5,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67241220 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,5,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67241221 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,5,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67241472 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,6,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67241473 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,6,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67241474 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,6,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67241475 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,6,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67241476 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,6,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67241477 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,6,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67241478 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,6,6,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#ifdef COMPILE_J_ORB
  case 67241728 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,7,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67241729 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,7,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67241730 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,7,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67241731 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,7,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67241732 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,7,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67241733 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,7,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67241734 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,7,6,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67241735 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,7,7,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
  case 67305472 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,0,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67305728 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,1,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67305729 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,1,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67305984 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,2,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67305985 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,2,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67305986 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,2,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67306240 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,3,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67306241 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,3,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67306242 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,3,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67306243 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,3,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67306496 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,4,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67306497 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,4,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67306498 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,4,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67306499 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,4,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67306500 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,4,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67306752 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,5,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67306753 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,5,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67306754 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,5,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67306755 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,5,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67306756 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,5,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67306757 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,5,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67307008 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,6,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67307009 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,6,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67307010 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,6,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67307011 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,6,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67307012 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,6,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67307013 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,6,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67307014 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,6,6,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#ifdef COMPILE_J_ORB
  case 67307264 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,7,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67307265 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,7,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67307266 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,7,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67307267 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,7,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67307268 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,7,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67307269 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,7,5,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67307270 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,7,6,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67307271 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,7,7,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
  case 67371008 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,0,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67371264 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,1,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67371265 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,1,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67371520 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,2,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67371521 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,2,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67371522 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,2,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67371776 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,3,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67371777 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,3,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67371778 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,3,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67371779 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,3,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67372032 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,4,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67372033 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,4,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67372034 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,4,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67372035 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,4,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67372036 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,4,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67372288 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,5,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67372289 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,5,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67372290 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,5,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67372291 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,5,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67372292 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,5,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67372293 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,5,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67372544 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,6,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67372545 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,6,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67372546 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,6,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67372547 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,6,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67372548 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,6,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67372549 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,6,5,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 67372550 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,6,6,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#ifdef COMPILE_J_ORB
  case 67372800 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,7,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67372801 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,7,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67372802 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,7,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67372803 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,7,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67372804 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,7,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67372805 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,7,5,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67372806 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,7,6,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 67372807 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,7,7,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
  case 83886080 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,0,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83886336 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,1,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83886337 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,1,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83886592 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,2,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83886593 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,2,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83886594 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,2,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83886848 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,3,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83886849 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,3,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83886850 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,3,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83886851 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,3,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83887104 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,4,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83887105 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,4,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83887106 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,4,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83887107 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,4,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83887108 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,4,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83887360 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,5,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83887361 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,5,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83887362 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,5,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83887363 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,5,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83887364 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,5,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83887365 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,5,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83887616 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,6,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83887617 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,6,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83887618 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,6,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83887619 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,6,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83887620 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,6,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83887621 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,6,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83887622 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,6,6,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#ifdef COMPILE_J_ORB
  case 83887872 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,7,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 83887873 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,7,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 83887874 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,7,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 83887875 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,7,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 83887876 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,7,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 83887877 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,7,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 83887878 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,7,6,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 83887879 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,7,7,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
  case 83951616 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,0,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83951872 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,1,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83951873 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,1,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83952128 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,2,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83952129 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,2,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83952130 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,2,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83952384 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,3,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83952385 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,3,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83952386 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,3,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83952387 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,3,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83952640 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,4,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83952641 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,4,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83952642 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,4,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83952643 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,4,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83952644 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,4,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83952896 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,5,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83952897 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,5,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83952898 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,5,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83952899 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,5,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83952900 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,5,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83952901 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,5,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83953152 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,6,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83953153 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,6,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83953154 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,6,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83953155 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,6,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83953156 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,6,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83953157 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,6,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 83953158 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,6,6,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#ifdef COMPILE_J_ORB
  case 83953408 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,7,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 83953409 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,7,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 83953410 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,7,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 83953411 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,7,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 83953412 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,7,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 83953413 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,7,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 83953414 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,7,6,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 83953415 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,7,7,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
  case 84017152 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,0,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84017408 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,1,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84017409 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,1,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84017664 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,2,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84017665 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,2,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84017666 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,2,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84017920 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,3,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84017921 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,3,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84017922 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,3,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84017923 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,3,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84018176 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,4,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84018177 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,4,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84018178 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,4,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84018179 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,4,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84018180 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,4,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84018432 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,5,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84018433 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,5,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84018434 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,5,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84018435 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,5,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84018436 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,5,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84018437 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,5,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84018688 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,6,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84018689 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,6,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84018690 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,6,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84018691 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,6,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84018692 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,6,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84018693 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,6,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84018694 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,6,6,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#ifdef COMPILE_J_ORB
  case 84018944 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,7,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 84018945 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,7,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 84018946 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,7,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 84018947 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,7,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 84018948 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,7,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 84018949 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,7,5,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 84018950 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,7,6,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 84018951 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,7,7,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
  case 84082688 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,0,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84082944 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,1,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84082945 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,1,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84083200 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,2,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84083201 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,2,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84083202 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,2,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84083456 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,3,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84083457 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,3,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84083458 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,3,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84083459 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,3,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84083712 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,4,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84083713 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,4,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84083714 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,4,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84083715 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,4,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84083716 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,4,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84083968 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,5,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84083969 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,5,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84083970 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,5,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84083971 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,5,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84083972 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,5,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84083973 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,5,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84084224 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,6,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84084225 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,6,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84084226 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,6,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84084227 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,6,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84084228 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,6,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84084229 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,6,5,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84084230 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,6,6,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#ifdef COMPILE_J_ORB
  case 84084480 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,7,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 84084481 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,7,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 84084482 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,7,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 84084483 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,7,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 84084484 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,7,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 84084485 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,7,5,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 84084486 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,7,6,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 84084487 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,7,7,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
  case 84148224 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,0,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84148480 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,1,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84148481 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,1,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84148736 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,2,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84148737 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,2,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84148738 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,2,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84148992 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,3,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84148993 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,3,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84148994 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,3,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84148995 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,3,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84149248 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,4,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84149249 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,4,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84149250 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,4,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84149251 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,4,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84149252 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,4,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84149504 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,5,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84149505 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,5,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84149506 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,5,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84149507 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,5,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84149508 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,5,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84149509 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,5,5,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84149760 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,6,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84149761 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,6,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84149762 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,6,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84149763 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,6,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84149764 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,6,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84149765 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,6,5,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84149766 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,6,6,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#ifdef COMPILE_J_ORB
  case 84150016 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,7,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 84150017 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,7,1,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 84150018 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,7,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 84150019 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,7,3,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 84150020 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,7,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 84150021 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,7,5,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 84150022 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,7,6,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 84150023 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,7,7,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
  case 84213760 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,0,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84214016 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,1,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84214017 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,1,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84214272 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,2,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84214273 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,2,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84214274 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,2,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84214528 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,3,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84214529 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,3,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84214530 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,3,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84214531 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,3,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84214784 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,4,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84214785 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,4,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84214786 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,4,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84214787 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,4,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84214788 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,4,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84215040 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,5,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84215041 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,5,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84215042 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,5,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84215043 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,5,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84215044 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,5,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84215045 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,5,5,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84215296 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,6,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84215297 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,6,1,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84215298 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,6,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84215299 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,6,3,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84215300 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,6,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84215301 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,6,5,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 84215302 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,6,6,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#ifdef COMPILE_J_ORB
  case 84215552 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,7,0,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 84215553 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,7,1,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 84215554 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,7,2,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 84215555 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,7,3,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 84215556 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,7,4,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 84215557 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,7,5,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 84215558 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,7,6,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 84215559 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,7,7,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
  case 100663296 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,0,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100663552 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,1,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100663553 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,1,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100663808 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,2,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100663809 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,2,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100663810 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,2,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100664064 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,3,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100664065 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,3,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100664066 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,3,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100664067 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,3,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100664320 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,4,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100664321 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,4,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100664322 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,4,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100664323 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,4,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100664324 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,4,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100664576 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,5,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100664577 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,5,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100664578 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,5,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100664579 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,5,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100664580 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,5,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100664581 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,5,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100664832 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,6,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100664833 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,6,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100664834 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,6,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100664835 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,6,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100664836 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,6,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100664837 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,6,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100664838 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,6,6,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#ifdef COMPILE_J_ORB
  case 100665088 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,7,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100665089 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,7,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100665090 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,7,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100665091 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,7,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100665092 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,7,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100665093 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,7,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100665094 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,7,6,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100665095 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,7,7,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
  case 100728832 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,0,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100729088 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,1,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100729089 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,1,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100729344 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,2,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100729345 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,2,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100729346 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,2,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100729600 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,3,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100729601 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,3,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100729602 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,3,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100729603 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,3,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100729856 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,4,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100729857 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,4,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100729858 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,4,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100729859 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,4,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100729860 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,4,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100730112 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,5,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100730113 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,5,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100730114 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,5,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100730115 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,5,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100730116 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,5,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100730117 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,5,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100730368 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,6,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100730369 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,6,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100730370 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,6,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100730371 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,6,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100730372 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,6,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100730373 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,6,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100730374 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,6,6,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#ifdef COMPILE_J_ORB
  case 100730624 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,7,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100730625 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,7,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100730626 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,7,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100730627 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,7,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100730628 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,7,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100730629 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,7,5,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100730630 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,7,6,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100730631 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,7,7,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
  case 100794368 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,0,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100794624 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,1,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100794625 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,1,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100794880 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,2,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100794881 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,2,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100794882 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,2,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100795136 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,3,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100795137 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,3,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100795138 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,3,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100795139 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,3,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100795392 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,4,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100795393 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,4,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100795394 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,4,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100795395 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,4,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100795396 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,4,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100795648 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,5,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100795649 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,5,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100795650 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,5,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100795651 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,5,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100795652 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,5,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100795653 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,5,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100795904 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,6,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100795905 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,6,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100795906 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,6,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100795907 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,6,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100795908 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,6,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100795909 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,6,5,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100795910 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,6,6,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#ifdef COMPILE_J_ORB
  case 100796160 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,7,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100796161 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,7,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100796162 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,7,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100796163 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,7,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100796164 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,7,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100796165 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,7,5,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100796166 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,7,6,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100796167 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,7,7,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
  case 100859904 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,0,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100860160 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,1,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100860161 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,1,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100860416 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,2,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100860417 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,2,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100860418 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,2,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100860672 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,3,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100860673 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,3,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100860674 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,3,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100860675 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,3,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100860928 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,4,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100860929 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,4,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100860930 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,4,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100860931 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,4,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100860932 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,4,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100861184 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,5,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100861185 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,5,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100861186 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,5,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100861187 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,5,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100861188 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,5,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100861189 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,5,5,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100861440 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,6,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100861441 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,6,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100861442 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,6,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100861443 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,6,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100861444 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,6,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100861445 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,6,5,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100861446 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,6,6,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#ifdef COMPILE_J_ORB
  case 100861696 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,7,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100861697 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,7,1,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100861698 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,7,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100861699 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,7,3,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100861700 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,7,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100861701 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,7,5,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100861702 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,7,6,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100861703 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,7,7,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
  case 100925440 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,0,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100925696 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,1,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100925697 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,1,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100925952 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,2,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100925953 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,2,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100925954 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,2,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100926208 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,3,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100926209 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,3,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100926210 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,3,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100926211 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,3,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100926464 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,4,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100926465 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,4,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100926466 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,4,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100926467 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,4,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100926468 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,4,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100926720 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,5,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100926721 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,5,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100926722 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,5,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100926723 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,5,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100926724 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,5,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100926725 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,5,5,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100926976 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,6,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100926977 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,6,1,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100926978 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,6,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100926979 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,6,3,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100926980 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,6,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100926981 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,6,5,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100926982 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,6,6,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#ifdef COMPILE_J_ORB
  case 100927232 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,7,0,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100927233 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,7,1,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100927234 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,7,2,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100927235 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,7,3,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100927236 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,7,4,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100927237 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,7,5,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100927238 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,7,6,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100927239 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,7,7,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
  case 100990976 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,0,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100991232 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,1,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100991233 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,1,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100991488 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,2,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100991489 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,2,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100991490 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,2,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100991744 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,3,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100991745 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,3,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100991746 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,3,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100991747 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,3,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100992000 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,4,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100992001 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,4,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100992002 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,4,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100992003 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,4,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100992004 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,4,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100992256 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,5,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100992257 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,5,1,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100992258 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,5,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100992259 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,5,3,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100992260 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,5,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100992261 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,5,5,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100992512 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,6,0,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100992513 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,6,1,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100992514 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,6,2,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100992515 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,6,3,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100992516 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,6,4,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100992517 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,6,5,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 100992518 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,6,6,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#ifdef COMPILE_J_ORB
  case 100992768 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,7,0,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100992769 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,7,1,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100992770 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,7,2,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100992771 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,7,3,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100992772 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,7,4,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100992773 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,7,5,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100992774 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,7,6,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 100992775 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,7,7,14>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
  case 101056512 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,0,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 101056768 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,1,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 101056769 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,1,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 101057024 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,2,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 101057025 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,2,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 101057026 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,2,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 101057280 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,3,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 101057281 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,3,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 101057282 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,3,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 101057283 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,3,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 101057536 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,4,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 101057537 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,4,1,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 101057538 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,4,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 101057539 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,4,3,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 101057540 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,4,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 101057792 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,5,0,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 101057793 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,5,1,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 101057794 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,5,2,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 101057795 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,5,3,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 101057796 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,5,4,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 101057797 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,5,5,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 101058048 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,6,0,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 101058049 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,6,1,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 101058050 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,6,2,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 101058051 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,6,3,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 101058052 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,6,4,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 101058053 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,6,5,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
  case 101058054 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,6,6,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#ifdef COMPILE_J_ORB
  case 101058304 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,7,0,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 101058305 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,7,1,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 101058306 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,7,2,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 101058307 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,7,3,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 101058308 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,7,4,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 101058309 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,7,5,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 101058310 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,7,6,14>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 101058311 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,7,7,14>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117440512 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,0,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117440768 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,1,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117440769 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,1,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117441024 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,2,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117441025 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,2,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117441026 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,2,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117441280 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,3,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117441281 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,3,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117441282 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,3,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117441283 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,3,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117441536 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,4,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117441537 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,4,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117441538 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,4,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117441539 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,4,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117441540 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,4,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117441792 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,5,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117441793 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,5,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117441794 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,5,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117441795 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,5,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117441796 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,5,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117441797 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,5,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117442048 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,6,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117442049 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,6,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117442050 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,6,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117442051 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,6,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117442052 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,6,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117442053 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,6,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117442054 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,6,6,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117442304 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,7,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117442305 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,7,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117442306 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,7,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117442307 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,7,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117442308 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,7,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117442309 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,7,5,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117442310 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,7,6,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117442311 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,0,7,7,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117506048 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,0,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117506304 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,1,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117506305 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,1,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117506560 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,2,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117506561 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,2,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117506562 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,2,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117506816 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,3,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117506817 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,3,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117506818 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,3,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117506819 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,3,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117507072 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,4,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117507073 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,4,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117507074 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,4,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117507075 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,4,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117507076 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,4,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117507328 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,5,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117507329 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,5,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117507330 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,5,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117507331 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,5,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117507332 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,5,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117507333 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,5,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117507584 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,6,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117507585 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,6,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117507586 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,6,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117507587 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,6,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117507588 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,6,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117507589 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,6,5,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117507590 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,6,6,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117507840 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,7,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117507841 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,7,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117507842 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,7,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117507843 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,7,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117507844 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,7,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117507845 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,7,5,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117507846 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,7,6,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117507847 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,1,7,7,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117571584 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,0,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117571840 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,1,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117571841 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,1,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117572096 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,2,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117572097 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,2,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117572098 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,2,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117572352 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,3,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117572353 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,3,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117572354 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,3,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117572355 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,3,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117572608 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,4,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117572609 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,4,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117572610 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,4,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117572611 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,4,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117572612 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,4,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117572864 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,5,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117572865 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,5,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117572866 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,5,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117572867 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,5,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117572868 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,5,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117572869 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,5,5,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117573120 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,6,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117573121 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,6,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117573122 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,6,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117573123 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,6,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117573124 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,6,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117573125 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,6,5,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117573126 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,6,6,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117573376 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,7,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117573377 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,7,1,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117573378 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,7,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117573379 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,7,3,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117573380 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,7,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117573381 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,7,5,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117573382 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,7,6,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117573383 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,2,7,7,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117637120 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,0,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117637376 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,1,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117637377 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,1,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117637632 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,2,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117637633 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,2,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117637634 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,2,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117637888 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,3,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117637889 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,3,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117637890 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,3,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117637891 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,3,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117638144 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,4,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117638145 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,4,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117638146 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,4,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117638147 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,4,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117638148 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,4,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117638400 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,5,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117638401 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,5,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117638402 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,5,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117638403 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,5,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117638404 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,5,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117638405 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,5,5,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117638656 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,6,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117638657 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,6,1,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117638658 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,6,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117638659 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,6,3,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117638660 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,6,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117638661 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,6,5,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117638662 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,6,6,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117638912 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,7,0,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117638913 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,7,1,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117638914 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,7,2,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117638915 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,7,3,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117638916 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,7,4,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117638917 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,7,5,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117638918 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,7,6,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117638919 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,3,7,7,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117702656 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,0,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117702912 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,1,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117702913 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,1,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117703168 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,2,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117703169 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,2,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117703170 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,2,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117703424 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,3,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117703425 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,3,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117703426 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,3,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117703427 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,3,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117703680 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,4,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117703681 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,4,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117703682 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,4,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117703683 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,4,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117703684 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,4,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117703936 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,5,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117703937 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,5,1,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117703938 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,5,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117703939 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,5,3,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117703940 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,5,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117703941 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,5,5,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117704192 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,6,0,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117704193 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,6,1,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117704194 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,6,2,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117704195 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,6,3,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117704196 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,6,4,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117704197 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,6,5,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117704198 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,6,6,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117704448 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,7,0,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117704449 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,7,1,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117704450 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,7,2,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117704451 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,7,3,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117704452 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,7,4,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117704453 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,7,5,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117704454 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,7,6,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117704455 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,4,7,7,14>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117768192 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,0,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117768448 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,1,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117768449 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,1,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117768704 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,2,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117768705 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,2,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117768706 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,2,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117768960 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,3,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117768961 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,3,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117768962 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,3,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117768963 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,3,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117769216 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,4,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117769217 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,4,1,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117769218 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,4,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117769219 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,4,3,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117769220 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,4,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117769472 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,5,0,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117769473 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,5,1,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117769474 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,5,2,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117769475 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,5,3,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117769476 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,5,4,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117769477 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,5,5,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117769728 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,6,0,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117769729 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,6,1,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117769730 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,6,2,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117769731 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,6,3,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117769732 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,6,4,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117769733 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,6,5,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117769734 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,6,6,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117769984 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,7,0,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117769985 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,7,1,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117769986 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,7,2,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117769987 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,7,3,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117769988 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,7,4,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117769989 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,7,5,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117769990 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,7,6,14>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117769991 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,5,7,7,14>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117833728 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,0,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117833984 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,1,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117833985 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,1,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117834240 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,2,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117834241 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,2,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117834242 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,2,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117834496 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,3,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117834497 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,3,1,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117834498 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,3,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117834499 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,3,3,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117834752 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,4,0,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117834753 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,4,1,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117834754 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,4,2,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117834755 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,4,3,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117834756 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,4,4,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117835008 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,5,0,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117835009 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,5,1,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117835010 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,5,2,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117835011 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,5,3,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117835012 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,5,4,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117835013 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,5,5,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117835264 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,6,0,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117835265 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,6,1,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117835266 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,6,2,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117835267 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,6,3,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117835268 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,6,4,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117835269 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,6,5,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117835270 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,6,6,14>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117835520 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,7,0,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117835521 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,7,1,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117835522 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,7,2,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117835523 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,7,3,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117835524 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,7,4,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117835525 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,7,5,14>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117835526 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,7,6,14>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117835527 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,6,7,7,15>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117899264 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,0,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117899520 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,1,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117899521 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,1,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117899776 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,2,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117899777 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,2,1,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117899778 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,2,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117900032 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,3,0,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117900033 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,3,1,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117900034 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,3,2,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117900035 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,3,3,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117900288 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,4,0,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117900289 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,4,1,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117900290 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,4,2,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117900291 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,4,3,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117900292 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,4,4,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117900544 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,5,0,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117900545 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,5,1,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117900546 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,5,2,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117900547 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,5,3,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117900548 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,5,4,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117900549 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,5,5,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117900800 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,6,0,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117900801 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,6,1,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117900802 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,6,2,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117900803 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,6,3,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117900804 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,6,4,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117900805 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,6,5,14>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117900806 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,6,6,14>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117901056 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,7,0,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117901057 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,7,1,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117901058 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,7,2,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117901059 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,7,3,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117901060 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,7,4,14>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117901061 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,7,5,14>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117901062 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,7,6,15>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
#ifdef COMPILE_J_ORB
  case 117901063 :
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<7,7,7,7,15>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz, dummy);
    } break;
#endif
  default :
    assert(false);   // hashkey not found
  }
  stack_->release(b2*a2*c2*d2*rank_, final_zc);
  stack_->release(b2*a2*c2*d2*rank_, final_zb);
  stack_->release(b2*a2*c2*d2*rank_, final_za);
  stack_->release(b2*a2*c2*d2*rank_, final_yc);
  stack_->release(b2*a2*c2*d2*rank_, final_yb);
  stack_->release(b2*a2*c2*d2*rank_, final_ya);
  stack_->release(b2*a2*c2*d2*rank_, final_xc);
  stack_->release(b2*a2*c2*d2*rank_, final_xb);
  stack_->release(b2*a2*c2*d2*rank_, final_xa);
  stack_->release(b2*a2*c2*d2*rank_, final_z);
  stack_->release(b2*a2*c2*d2*rank_, final_y);
  stack_->release(b2*a2*c2*d2*rank_, final_x);

  stack_->release(b2*a2*(cmax_+1)*rank_, intermediate);

  stack_->release((cmax_+1)*c2*d2, trans2z);
  stack_->release((cmax_+1)*c2*d2, trans2y);
  stack_->release((cmax_+1)*c2*d2, trans2x);

  stack_->release((amax_+1)*a2*b2, transz);
  stack_->release((amax_+1)*a2*b2, transy);
  stack_->release((amax_+1)*a2*b2, transx);
  stack_->release(rank_*isize*3, workx);

#endif
}
