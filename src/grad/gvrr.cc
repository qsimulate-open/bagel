//
// BAGEL - Parallel electron correlation program.
// Filename: gvrr.cc
// Copyright (C) 2012 Toru Shiozaki
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

#include <src/grad/gradbatch.h>
#include <src/grad/_gvrr_drv.h>
#include <src/util/comb.h>
#include <src/util/f77.h>

using namespace std;
using namespace bagel;

static const Comb comb;


void GradBatch::perform_VRR() {
#ifndef LIBINT_INTERFACE
  const int acsize = size_block_ / primsize_;
  const int a = basisinfo_[0]->angular_number();
  const int b = basisinfo_[1]->angular_number();
  const int c = basisinfo_[2]->angular_number();
  const int d = basisinfo_[3]->angular_number();

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
  fill(transx,  transx +(amax_+1)*a2*b2, 0.0);
  fill(transy,  transy +(amax_+1)*a2*b2, 0.0);
  fill(transz,  transz +(amax_+1)*a2*b2, 0.0);
  fill(trans2x, trans2x+(cmax_+1)*c2*d2, 0.0);
  fill(trans2y, trans2y+(cmax_+1)*c2*d2, 0.0);
  fill(trans2z, trans2z+(cmax_+1)*c2*d2, 0.0);
  // for usual integrals
  for (int ib = 0, k = 0; ib <= b+1; ++ib) {
    for (int ia = 0; ia <= a+1; ++ia, ++k) {
      if (ia == a+1 && ib == b+1) continue;
      for (int i = ia; i <= ia+ib; ++i) {
        transx[i + (amax_+1)*k] = comb.c(ib, ia+ib-i) * pow(AB_[0], ia+ib-i);
        transy[i + (amax_+1)*k] = comb.c(ib, ia+ib-i) * pow(AB_[1], ia+ib-i);
        transz[i + (amax_+1)*k] = comb.c(ib, ia+ib-i) * pow(AB_[2], ia+ib-i);
      }   
    }   
  }
  for (int id = 0, k = 0; id <= d+1; ++id) {
    for (int ic = 0; ic <= c+1; ++ic, ++k) {
      if (ic == c+1 && id == d+1) continue;
      for (int i = ic; i <= ic+id; ++i) {
        trans2x[i + (cmax_+1)*k] = comb.c(id, ic+id-i) * pow(CD_[0], ic+id-i);
        trans2y[i + (cmax_+1)*k] = comb.c(id, ic+id-i) * pow(CD_[1], ic+id-i);
        trans2z[i + (cmax_+1)*k] = comb.c(id, ic+id-i) * pow(CD_[2], ic+id-i);
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
  if (a == 0 && b == 0 && c == 0 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,0,0,1>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 0 && b == 0 && c == 1 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,1,0,2>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 0 && b == 0 && c == 1 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,1,1,2>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 0 && b == 0 && c == 2 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,2,0,2>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 0 && b == 0 && c == 2 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,2,1,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 0 && b == 0 && c == 2 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,2,2,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 0 && b == 0 && c == 3 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,3,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 0 && b == 0 && c == 3 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,3,1,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 0 && b == 0 && c == 3 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,3,2,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 0 && b == 0 && c == 3 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,3,3,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 0 && b == 0 && c == 4 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,4,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 0 && b == 0 && c == 4 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,4,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 0 && b == 0 && c == 4 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,4,2,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 0 && b == 0 && c == 4 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,4,3,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 0 && b == 0 && c == 4 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,4,4,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 0 && b == 0 && c == 5 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,5,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 0 && b == 0 && c == 5 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,5,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 0 && b == 0 && c == 5 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,5,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 0 && b == 0 && c == 5 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,5,3,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 0 && b == 0 && c == 5 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,5,4,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 0 && b == 0 && c == 5 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,5,5,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 0 && b == 0 && c == 6 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,6,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 0 && b == 0 && c == 6 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,6,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 0 && b == 0 && c == 6 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,6,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 0 && b == 0 && c == 6 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,6,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 0 && b == 0 && c == 6 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,6,4,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 0 && b == 0 && c == 6 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,6,5,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 0 && b == 0 && c == 6 && d == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<0,0,6,6,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 0 && c == 0 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,0,0,2>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 0 && c == 1 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,1,0,2>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 0 && c == 1 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,1,1,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 0 && c == 2 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,2,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 0 && c == 2 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,2,1,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 0 && c == 2 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,2,2,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 0 && c == 3 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,3,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 0 && c == 3 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,3,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 0 && c == 3 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,3,2,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 0 && c == 3 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,3,3,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 0 && c == 4 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,4,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 0 && c == 4 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,4,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 0 && c == 4 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,4,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 0 && c == 4 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,4,3,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 0 && c == 4 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,4,4,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 0 && c == 5 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,5,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 0 && c == 5 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,5,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 0 && c == 5 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,5,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 0 && c == 5 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,5,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 0 && c == 5 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,5,4,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 0 && c == 5 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,5,5,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 0 && c == 6 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,6,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 0 && c == 6 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,6,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 0 && c == 6 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,6,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 0 && c == 6 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,6,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 0 && c == 6 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,6,4,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 0 && c == 6 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,6,5,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 0 && c == 6 && d == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,0,6,6,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 1 && c == 0 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,0,0,2>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 1 && c == 1 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,1,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 1 && c == 1 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,1,1,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 1 && c == 2 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,2,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 1 && c == 2 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,2,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 1 && c == 2 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,2,2,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 1 && c == 3 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,3,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 1 && c == 3 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,3,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 1 && c == 3 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,3,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 1 && c == 3 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,3,3,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 1 && c == 4 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,4,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 1 && c == 4 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,4,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 1 && c == 4 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,4,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 1 && c == 4 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,4,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 1 && c == 4 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,4,4,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 1 && c == 5 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,5,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 1 && c == 5 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,5,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 1 && c == 5 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,5,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 1 && c == 5 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,5,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 1 && c == 5 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,5,4,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 1 && c == 5 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,5,5,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 1 && c == 6 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,6,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 1 && c == 6 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,6,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 1 && c == 6 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,6,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 1 && c == 6 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,6,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 1 && c == 6 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,6,4,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 1 && c == 6 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,6,5,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 1 && b == 1 && c == 6 && d == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<1,1,6,6,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 0 && c == 0 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,0,0,2>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 0 && c == 1 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,1,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 0 && c == 1 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,1,1,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 0 && c == 2 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,2,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 0 && c == 2 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,2,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 0 && c == 2 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,2,2,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 0 && c == 3 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,3,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 0 && c == 3 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,3,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 0 && c == 3 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,3,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 0 && c == 3 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,3,3,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 0 && c == 4 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,4,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 0 && c == 4 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,4,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 0 && c == 4 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,4,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 0 && c == 4 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,4,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 0 && c == 4 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,4,4,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 0 && c == 5 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,5,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 0 && c == 5 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,5,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 0 && c == 5 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,5,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 0 && c == 5 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,5,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 0 && c == 5 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,5,4,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 0 && c == 5 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,5,5,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 0 && c == 6 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,6,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 0 && c == 6 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,6,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 0 && c == 6 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,6,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 0 && c == 6 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,6,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 0 && c == 6 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,6,4,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 0 && c == 6 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,6,5,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 0 && c == 6 && d == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,0,6,6,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 1 && c == 0 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,0,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 1 && c == 1 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,1,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 1 && c == 1 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,1,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 1 && c == 2 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,2,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 1 && c == 2 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,2,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 1 && c == 2 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,2,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 1 && c == 3 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,3,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 1 && c == 3 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,3,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 1 && c == 3 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,3,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 1 && c == 3 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,3,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 1 && c == 4 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,4,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 1 && c == 4 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,4,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 1 && c == 4 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,4,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 1 && c == 4 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,4,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 1 && c == 4 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,4,4,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 1 && c == 5 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,5,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 1 && c == 5 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,5,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 1 && c == 5 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,5,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 1 && c == 5 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,5,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 1 && c == 5 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,5,4,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 1 && c == 5 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,5,5,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 1 && c == 6 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,6,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 1 && c == 6 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,6,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 1 && c == 6 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,6,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 1 && c == 6 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,6,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 1 && c == 6 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,6,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 1 && c == 6 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,6,5,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 1 && c == 6 && d == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,1,6,6,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 2 && c == 0 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,0,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 2 && c == 1 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,1,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 2 && c == 1 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,1,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 2 && c == 2 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,2,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 2 && c == 2 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,2,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 2 && c == 2 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,2,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 2 && c == 3 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,3,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 2 && c == 3 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,3,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 2 && c == 3 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,3,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 2 && c == 3 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,3,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 2 && c == 4 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,4,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 2 && c == 4 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,4,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 2 && c == 4 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,4,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 2 && c == 4 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,4,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 2 && c == 4 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,4,4,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 2 && c == 5 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,5,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 2 && c == 5 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,5,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 2 && c == 5 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,5,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 2 && c == 5 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,5,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 2 && c == 5 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,5,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 2 && c == 5 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,5,5,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 2 && c == 6 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,6,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 2 && c == 6 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,6,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 2 && c == 6 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,6,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 2 && c == 6 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,6,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 2 && c == 6 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,6,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 2 && c == 6 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,6,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 2 && b == 2 && c == 6 && d == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<2,2,6,6,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 0 && c == 0 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,0,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 0 && c == 1 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,1,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 0 && c == 1 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,1,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 0 && c == 2 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,2,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 0 && c == 2 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,2,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 0 && c == 2 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,2,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 0 && c == 3 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,3,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 0 && c == 3 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,3,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 0 && c == 3 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,3,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 0 && c == 3 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,3,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 0 && c == 4 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,4,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 0 && c == 4 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,4,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 0 && c == 4 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,4,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 0 && c == 4 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,4,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 0 && c == 4 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,4,4,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 0 && c == 5 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,5,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 0 && c == 5 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,5,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 0 && c == 5 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,5,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 0 && c == 5 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,5,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 0 && c == 5 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,5,4,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 0 && c == 5 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,5,5,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 0 && c == 6 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,6,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 0 && c == 6 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,6,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 0 && c == 6 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,6,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 0 && c == 6 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,6,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 0 && c == 6 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,6,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 0 && c == 6 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,6,5,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 0 && c == 6 && d == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,0,6,6,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 1 && c == 0 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,0,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 1 && c == 1 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,1,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 1 && c == 1 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,1,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 1 && c == 2 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,2,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 1 && c == 2 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,2,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 1 && c == 2 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,2,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 1 && c == 3 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,3,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 1 && c == 3 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,3,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 1 && c == 3 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,3,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 1 && c == 3 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,3,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 1 && c == 4 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,4,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 1 && c == 4 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,4,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 1 && c == 4 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,4,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 1 && c == 4 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,4,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 1 && c == 4 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,4,4,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 1 && c == 5 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,5,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 1 && c == 5 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,5,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 1 && c == 5 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,5,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 1 && c == 5 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,5,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 1 && c == 5 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,5,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 1 && c == 5 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,5,5,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 1 && c == 6 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,6,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 1 && c == 6 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,6,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 1 && c == 6 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,6,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 1 && c == 6 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,6,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 1 && c == 6 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,6,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 1 && c == 6 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,6,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 1 && c == 6 && d == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,1,6,6,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 2 && c == 0 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,0,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 2 && c == 1 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,1,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 2 && c == 1 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,1,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 2 && c == 2 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,2,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 2 && c == 2 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,2,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 2 && c == 2 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,2,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 2 && c == 3 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,3,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 2 && c == 3 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,3,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 2 && c == 3 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,3,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 2 && c == 3 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,3,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 2 && c == 4 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,4,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 2 && c == 4 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,4,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 2 && c == 4 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,4,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 2 && c == 4 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,4,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 2 && c == 4 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,4,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 2 && c == 5 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,5,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 2 && c == 5 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,5,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 2 && c == 5 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,5,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 2 && c == 5 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,5,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 2 && c == 5 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,5,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 2 && c == 5 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,5,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 2 && c == 6 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,6,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 2 && c == 6 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,6,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 2 && c == 6 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,6,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 2 && c == 6 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,6,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 2 && c == 6 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,6,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 2 && c == 6 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,6,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 2 && c == 6 && d == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,2,6,6,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 3 && c == 0 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,0,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 3 && c == 1 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,1,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 3 && c == 1 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,1,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 3 && c == 2 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,2,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 3 && c == 2 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,2,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 3 && c == 2 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,2,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 3 && c == 3 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,3,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 3 && c == 3 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,3,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 3 && c == 3 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,3,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 3 && c == 3 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,3,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 3 && c == 4 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,4,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 3 && c == 4 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,4,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 3 && c == 4 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,4,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 3 && c == 4 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,4,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 3 && c == 4 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,4,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 3 && c == 5 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,5,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 3 && c == 5 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,5,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 3 && c == 5 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,5,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 3 && c == 5 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,5,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 3 && c == 5 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,5,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 3 && c == 5 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,5,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 3 && c == 6 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,6,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 3 && c == 6 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,6,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 3 && c == 6 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,6,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 3 && c == 6 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,6,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 3 && c == 6 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,6,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 3 && c == 6 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,6,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 3 && b == 3 && c == 6 && d == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<3,3,6,6,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 0 && c == 0 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,0,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 0 && c == 1 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,1,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 0 && c == 1 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,1,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 0 && c == 2 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,2,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 0 && c == 2 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,2,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 0 && c == 2 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,2,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 0 && c == 3 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,3,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 0 && c == 3 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,3,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 0 && c == 3 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,3,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 0 && c == 3 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,3,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 0 && c == 4 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,4,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 0 && c == 4 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,4,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 0 && c == 4 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,4,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 0 && c == 4 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,4,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 0 && c == 4 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,4,4,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 0 && c == 5 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,5,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 0 && c == 5 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,5,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 0 && c == 5 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,5,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 0 && c == 5 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,5,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 0 && c == 5 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,5,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 0 && c == 5 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,5,5,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 0 && c == 6 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,6,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 0 && c == 6 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,6,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 0 && c == 6 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,6,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 0 && c == 6 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,6,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 0 && c == 6 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,6,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 0 && c == 6 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,6,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 0 && c == 6 && d == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,0,6,6,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 1 && c == 0 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,0,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 1 && c == 1 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,1,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 1 && c == 1 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,1,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 1 && c == 2 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,2,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 1 && c == 2 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,2,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 1 && c == 2 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,2,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 1 && c == 3 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,3,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 1 && c == 3 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,3,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 1 && c == 3 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,3,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 1 && c == 3 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,3,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 1 && c == 4 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,4,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 1 && c == 4 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,4,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 1 && c == 4 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,4,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 1 && c == 4 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,4,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 1 && c == 4 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,4,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 1 && c == 5 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,5,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 1 && c == 5 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,5,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 1 && c == 5 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,5,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 1 && c == 5 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,5,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 1 && c == 5 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,5,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 1 && c == 5 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,5,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 1 && c == 6 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,6,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 1 && c == 6 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,6,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 1 && c == 6 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,6,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 1 && c == 6 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,6,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 1 && c == 6 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,6,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 1 && c == 6 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,6,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 1 && c == 6 && d == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,1,6,6,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 2 && c == 0 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,0,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 2 && c == 1 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,1,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 2 && c == 1 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,1,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 2 && c == 2 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,2,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 2 && c == 2 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,2,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 2 && c == 2 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,2,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 2 && c == 3 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,3,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 2 && c == 3 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,3,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 2 && c == 3 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,3,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 2 && c == 3 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,3,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 2 && c == 4 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,4,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 2 && c == 4 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,4,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 2 && c == 4 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,4,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 2 && c == 4 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,4,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 2 && c == 4 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,4,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 2 && c == 5 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,5,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 2 && c == 5 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,5,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 2 && c == 5 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,5,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 2 && c == 5 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,5,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 2 && c == 5 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,5,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 2 && c == 5 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,5,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 2 && c == 6 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,6,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 2 && c == 6 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,6,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 2 && c == 6 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,6,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 2 && c == 6 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,6,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 2 && c == 6 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,6,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 2 && c == 6 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,6,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 2 && c == 6 && d == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,2,6,6,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 3 && c == 0 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,0,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 3 && c == 1 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,1,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 3 && c == 1 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,1,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 3 && c == 2 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,2,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 3 && c == 2 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,2,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 3 && c == 2 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,2,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 3 && c == 3 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,3,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 3 && c == 3 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,3,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 3 && c == 3 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,3,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 3 && c == 3 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,3,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 3 && c == 4 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,4,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 3 && c == 4 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,4,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 3 && c == 4 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,4,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 3 && c == 4 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,4,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 3 && c == 4 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,4,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 3 && c == 5 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,5,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 3 && c == 5 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,5,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 3 && c == 5 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,5,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 3 && c == 5 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,5,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 3 && c == 5 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,5,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 3 && c == 5 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,5,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 3 && c == 6 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,6,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 3 && c == 6 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,6,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 3 && c == 6 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,6,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 3 && c == 6 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,6,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 3 && c == 6 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,6,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 3 && c == 6 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,6,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 3 && c == 6 && d == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,3,6,6,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 4 && c == 0 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,0,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 4 && c == 1 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,1,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 4 && c == 1 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,1,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 4 && c == 2 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,2,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 4 && c == 2 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,2,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 4 && c == 2 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,2,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 4 && c == 3 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,3,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 4 && c == 3 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,3,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 4 && c == 3 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,3,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 4 && c == 3 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,3,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 4 && c == 4 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,4,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 4 && c == 4 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,4,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 4 && c == 4 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,4,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 4 && c == 4 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,4,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 4 && c == 4 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,4,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 4 && c == 5 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,5,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 4 && c == 5 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,5,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 4 && c == 5 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,5,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 4 && c == 5 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,5,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 4 && c == 5 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,5,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 4 && c == 5 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,5,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 4 && c == 6 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,6,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 4 && c == 6 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,6,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 4 && c == 6 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,6,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 4 && c == 6 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,6,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 4 && c == 6 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,6,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 4 && c == 6 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,6,5,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 4 && b == 4 && c == 6 && d == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<4,4,6,6,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 0 && c == 0 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,0,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 0 && c == 1 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,1,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 0 && c == 1 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,1,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 0 && c == 2 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,2,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 0 && c == 2 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,2,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 0 && c == 2 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,2,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 0 && c == 3 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,3,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 0 && c == 3 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,3,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 0 && c == 3 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,3,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 0 && c == 3 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,3,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 0 && c == 4 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,4,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 0 && c == 4 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,4,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 0 && c == 4 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,4,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 0 && c == 4 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,4,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 0 && c == 4 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,4,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 0 && c == 5 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,5,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 0 && c == 5 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,5,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 0 && c == 5 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,5,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 0 && c == 5 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,5,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 0 && c == 5 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,5,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 0 && c == 5 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,5,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 0 && c == 6 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,6,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 0 && c == 6 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,6,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 0 && c == 6 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,6,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 0 && c == 6 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,6,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 0 && c == 6 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,6,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 0 && c == 6 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,6,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 0 && c == 6 && d == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,0,6,6,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 1 && c == 0 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,0,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 1 && c == 1 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,1,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 1 && c == 1 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,1,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 1 && c == 2 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,2,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 1 && c == 2 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,2,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 1 && c == 2 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,2,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 1 && c == 3 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,3,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 1 && c == 3 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,3,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 1 && c == 3 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,3,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 1 && c == 3 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,3,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 1 && c == 4 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,4,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 1 && c == 4 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,4,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 1 && c == 4 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,4,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 1 && c == 4 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,4,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 1 && c == 4 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,4,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 1 && c == 5 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,5,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 1 && c == 5 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,5,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 1 && c == 5 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,5,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 1 && c == 5 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,5,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 1 && c == 5 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,5,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 1 && c == 5 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,5,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 1 && c == 6 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,6,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 1 && c == 6 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,6,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 1 && c == 6 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,6,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 1 && c == 6 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,6,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 1 && c == 6 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,6,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 1 && c == 6 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,6,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 1 && c == 6 && d == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,1,6,6,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 2 && c == 0 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,0,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 2 && c == 1 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,1,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 2 && c == 1 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,1,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 2 && c == 2 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,2,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 2 && c == 2 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,2,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 2 && c == 2 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,2,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 2 && c == 3 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,3,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 2 && c == 3 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,3,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 2 && c == 3 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,3,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 2 && c == 3 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,3,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 2 && c == 4 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,4,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 2 && c == 4 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,4,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 2 && c == 4 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,4,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 2 && c == 4 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,4,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 2 && c == 4 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,4,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 2 && c == 5 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,5,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 2 && c == 5 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,5,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 2 && c == 5 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,5,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 2 && c == 5 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,5,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 2 && c == 5 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,5,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 2 && c == 5 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,5,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 2 && c == 6 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,6,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 2 && c == 6 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,6,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 2 && c == 6 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,6,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 2 && c == 6 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,6,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 2 && c == 6 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,6,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 2 && c == 6 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,6,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 2 && c == 6 && d == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,2,6,6,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 3 && c == 0 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,0,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 3 && c == 1 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,1,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 3 && c == 1 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,1,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 3 && c == 2 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,2,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 3 && c == 2 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,2,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 3 && c == 2 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,2,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 3 && c == 3 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,3,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 3 && c == 3 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,3,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 3 && c == 3 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,3,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 3 && c == 3 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,3,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 3 && c == 4 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,4,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 3 && c == 4 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,4,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 3 && c == 4 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,4,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 3 && c == 4 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,4,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 3 && c == 4 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,4,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 3 && c == 5 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,5,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 3 && c == 5 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,5,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 3 && c == 5 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,5,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 3 && c == 5 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,5,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 3 && c == 5 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,5,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 3 && c == 5 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,5,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 3 && c == 6 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,6,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 3 && c == 6 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,6,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 3 && c == 6 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,6,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 3 && c == 6 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,6,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 3 && c == 6 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,6,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 3 && c == 6 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,6,5,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 3 && c == 6 && d == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,3,6,6,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 4 && c == 0 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,0,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 4 && c == 1 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,1,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 4 && c == 1 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,1,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 4 && c == 2 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,2,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 4 && c == 2 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,2,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 4 && c == 2 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,2,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 4 && c == 3 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,3,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 4 && c == 3 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,3,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 4 && c == 3 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,3,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 4 && c == 3 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,3,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 4 && c == 4 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,4,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 4 && c == 4 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,4,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 4 && c == 4 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,4,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 4 && c == 4 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,4,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 4 && c == 4 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,4,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 4 && c == 5 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,5,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 4 && c == 5 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,5,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 4 && c == 5 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,5,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 4 && c == 5 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,5,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 4 && c == 5 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,5,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 4 && c == 5 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,5,5,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 4 && c == 6 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,6,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 4 && c == 6 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,6,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 4 && c == 6 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,6,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 4 && c == 6 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,6,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 4 && c == 6 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,6,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 4 && c == 6 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,6,5,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 4 && c == 6 && d == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,4,6,6,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 5 && c == 0 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,0,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 5 && c == 1 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,1,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 5 && c == 1 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,1,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 5 && c == 2 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,2,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 5 && c == 2 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,2,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 5 && c == 2 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,2,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 5 && c == 3 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,3,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 5 && c == 3 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,3,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 5 && c == 3 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,3,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 5 && c == 3 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,3,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 5 && c == 4 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,4,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 5 && c == 4 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,4,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 5 && c == 4 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,4,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 5 && c == 4 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,4,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 5 && c == 4 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,4,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 5 && c == 5 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,5,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 5 && c == 5 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,5,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 5 && c == 5 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,5,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 5 && c == 5 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,5,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 5 && c == 5 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,5,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 5 && c == 5 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,5,5,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 5 && c == 6 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,6,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 5 && c == 6 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,6,1,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 5 && c == 6 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,6,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 5 && c == 6 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,6,3,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 5 && c == 6 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,6,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 5 && c == 6 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,6,5,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 5 && b == 5 && c == 6 && d == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<5,5,6,6,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 0 && c == 0 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,0,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 0 && c == 1 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,1,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 0 && c == 1 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,1,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 0 && c == 2 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,2,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 0 && c == 2 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,2,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 0 && c == 2 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,2,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 0 && c == 3 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,3,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 0 && c == 3 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,3,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 0 && c == 3 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,3,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 0 && c == 3 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,3,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 0 && c == 4 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,4,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 0 && c == 4 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,4,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 0 && c == 4 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,4,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 0 && c == 4 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,4,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 0 && c == 4 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,4,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 0 && c == 5 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,5,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 0 && c == 5 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,5,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 0 && c == 5 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,5,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 0 && c == 5 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,5,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 0 && c == 5 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,5,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 0 && c == 5 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,5,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 0 && c == 6 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,6,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 0 && c == 6 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,6,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 0 && c == 6 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,6,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 0 && c == 6 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,6,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 0 && c == 6 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,6,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 0 && c == 6 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,6,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 0 && c == 6 && d == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,0,6,6,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 1 && c == 0 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,0,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 1 && c == 1 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,1,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 1 && c == 1 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,1,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 1 && c == 2 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,2,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 1 && c == 2 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,2,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 1 && c == 2 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,2,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 1 && c == 3 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,3,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 1 && c == 3 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,3,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 1 && c == 3 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,3,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 1 && c == 3 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,3,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 1 && c == 4 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,4,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 1 && c == 4 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,4,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 1 && c == 4 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,4,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 1 && c == 4 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,4,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 1 && c == 4 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,4,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 1 && c == 5 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,5,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 1 && c == 5 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,5,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 1 && c == 5 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,5,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 1 && c == 5 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,5,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 1 && c == 5 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,5,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 1 && c == 5 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,5,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 1 && c == 6 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,6,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 1 && c == 6 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,6,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 1 && c == 6 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,6,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 1 && c == 6 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,6,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 1 && c == 6 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,6,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 1 && c == 6 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,6,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 1 && c == 6 && d == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,1,6,6,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 2 && c == 0 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,0,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 2 && c == 1 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,1,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 2 && c == 1 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,1,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 2 && c == 2 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,2,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 2 && c == 2 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,2,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 2 && c == 2 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,2,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 2 && c == 3 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,3,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 2 && c == 3 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,3,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 2 && c == 3 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,3,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 2 && c == 3 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,3,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 2 && c == 4 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,4,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 2 && c == 4 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,4,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 2 && c == 4 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,4,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 2 && c == 4 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,4,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 2 && c == 4 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,4,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 2 && c == 5 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,5,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 2 && c == 5 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,5,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 2 && c == 5 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,5,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 2 && c == 5 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,5,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 2 && c == 5 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,5,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 2 && c == 5 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,5,5,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 2 && c == 6 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,6,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 2 && c == 6 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,6,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 2 && c == 6 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,6,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 2 && c == 6 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,6,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 2 && c == 6 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,6,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 2 && c == 6 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,6,5,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 2 && c == 6 && d == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,2,6,6,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 3 && c == 0 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,0,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 3 && c == 1 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,1,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 3 && c == 1 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,1,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 3 && c == 2 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,2,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 3 && c == 2 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,2,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 3 && c == 2 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,2,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 3 && c == 3 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,3,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 3 && c == 3 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,3,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 3 && c == 3 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,3,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 3 && c == 3 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,3,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 3 && c == 4 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,4,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 3 && c == 4 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,4,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 3 && c == 4 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,4,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 3 && c == 4 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,4,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 3 && c == 4 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,4,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 3 && c == 5 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,5,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 3 && c == 5 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,5,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 3 && c == 5 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,5,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 3 && c == 5 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,5,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 3 && c == 5 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,5,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 3 && c == 5 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,5,5,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 3 && c == 6 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,6,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 3 && c == 6 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,6,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 3 && c == 6 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,6,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 3 && c == 6 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,6,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 3 && c == 6 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,6,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 3 && c == 6 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,6,5,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 3 && c == 6 && d == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,3,6,6,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 4 && c == 0 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,0,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 4 && c == 1 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,1,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 4 && c == 1 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,1,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 4 && c == 2 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,2,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 4 && c == 2 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,2,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 4 && c == 2 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,2,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 4 && c == 3 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,3,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 4 && c == 3 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,3,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 4 && c == 3 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,3,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 4 && c == 3 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,3,3,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 4 && c == 4 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,4,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 4 && c == 4 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,4,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 4 && c == 4 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,4,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 4 && c == 4 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,4,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 4 && c == 4 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,4,4,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 4 && c == 5 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,5,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 4 && c == 5 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,5,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 4 && c == 5 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,5,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 4 && c == 5 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,5,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 4 && c == 5 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,5,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 4 && c == 5 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,5,5,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 4 && c == 6 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,6,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 4 && c == 6 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,6,1,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 4 && c == 6 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,6,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 4 && c == 6 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,6,3,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 4 && c == 6 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,6,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 4 && c == 6 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,6,5,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 4 && c == 6 && d == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,4,6,6,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 5 && c == 0 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,0,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 5 && c == 1 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,1,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 5 && c == 1 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,1,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 5 && c == 2 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,2,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 5 && c == 2 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,2,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 5 && c == 2 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,2,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 5 && c == 3 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,3,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 5 && c == 3 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,3,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 5 && c == 3 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,3,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 5 && c == 3 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,3,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 5 && c == 4 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,4,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 5 && c == 4 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,4,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 5 && c == 4 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,4,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 5 && c == 4 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,4,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 5 && c == 4 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,4,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 5 && c == 5 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,5,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 5 && c == 5 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,5,1,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 5 && c == 5 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,5,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 5 && c == 5 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,5,3,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 5 && c == 5 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,5,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 5 && c == 5 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,5,5,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 5 && c == 6 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,6,0,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 5 && c == 6 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,6,1,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 5 && c == 6 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,6,2,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 5 && c == 6 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,6,3,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 5 && c == 6 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,6,4,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 5 && c == 6 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,6,5,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 5 && c == 6 && d == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,5,6,6,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 6 && c == 0 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,0,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 6 && c == 1 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,1,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 6 && c == 1 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,1,1,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 6 && c == 2 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,2,0,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 6 && c == 2 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,2,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 6 && c == 2 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,2,2,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 6 && c == 3 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,3,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 6 && c == 3 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,3,1,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 6 && c == 3 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,3,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 6 && c == 3 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,3,3,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 6 && c == 4 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,4,0,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 6 && c == 4 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,4,1,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 6 && c == 4 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,4,2,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 6 && c == 4 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,4,3,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 6 && c == 4 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,4,4,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 6 && c == 5 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,5,0,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 6 && c == 5 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,5,1,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 6 && c == 5 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,5,2,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 6 && c == 5 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,5,3,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 6 && c == 5 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,5,4,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 6 && c == 5 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,5,5,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 6 && c == 6 && d == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,6,0,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 6 && c == 6 && d == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,6,1,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 6 && c == 6 && d == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,6,2,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 6 && c == 6 && d == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,6,3,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 6 && c == 6 && d == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,6,4,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 6 && c == 6 && d == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,6,5,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
  } else if (a == 6 && b == 6 && c == 6 && d == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      gvrr_driver<6,6,6,6,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], size_block_,
                    exponents_.get()+ii*4, transx, transy, transz, trans2x, trans2y, trans2z, intermediate,
                    final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc, final_za, final_zb, final_zc, workx, worky, workz);
    }
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
