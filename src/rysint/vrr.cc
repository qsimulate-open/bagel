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
#include <src/rysint/int2d.h>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <src/parallel/resources.h>
#include <src/rysint/_vrr_drv.h>
#include <src/util/f77.h>

using namespace std;
using namespace bagel;


void ERIBatch::perform_VRR() {
  const int acsize = asize_ * csize_;
  if (amax_ == 0 && cmax_ == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<0,0,1>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 0 && cmax_ == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<0,1,1>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 0 && cmax_ == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<0,2,2>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 0 && cmax_ == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<0,3,2>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 0 && cmax_ == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<0,4,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 0 && cmax_ == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<0,5,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 0 && cmax_ == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<0,6,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 0 && cmax_ == 7) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<0,7,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 0 && cmax_ == 8) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<0,8,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 0 && cmax_ == 9) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<0,9,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 0 && cmax_ == 10) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<0,10,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 0 && cmax_ == 11) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<0,11,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 0 && cmax_ == 12) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<0,12,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 1 && cmax_ == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<1,0,1>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 1 && cmax_ == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<1,1,2>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 1 && cmax_ == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<1,2,2>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 1 && cmax_ == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<1,3,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 1 && cmax_ == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<1,4,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 1 && cmax_ == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<1,5,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 1 && cmax_ == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<1,6,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 1 && cmax_ == 7) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<1,7,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 1 && cmax_ == 8) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<1,8,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 1 && cmax_ == 9) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<1,9,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 1 && cmax_ == 10) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<1,10,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 1 && cmax_ == 11) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<1,11,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 1 && cmax_ == 12) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<1,12,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 2 && cmax_ == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<2,0,2>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 2 && cmax_ == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<2,1,2>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 2 && cmax_ == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<2,2,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 2 && cmax_ == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<2,3,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 2 && cmax_ == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<2,4,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 2 && cmax_ == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<2,5,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 2 && cmax_ == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<2,6,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 2 && cmax_ == 7) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<2,7,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 2 && cmax_ == 8) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<2,8,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 2 && cmax_ == 9) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<2,9,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 2 && cmax_ == 10) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<2,10,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 2 && cmax_ == 11) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<2,11,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 2 && cmax_ == 12) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<2,12,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 3 && cmax_ == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<3,0,2>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 3 && cmax_ == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<3,1,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 3 && cmax_ == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<3,2,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 3 && cmax_ == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<3,3,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 3 && cmax_ == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<3,4,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 3 && cmax_ == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<3,5,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 3 && cmax_ == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<3,6,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 3 && cmax_ == 7) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<3,7,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 3 && cmax_ == 8) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<3,8,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 3 && cmax_ == 9) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<3,9,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 3 && cmax_ == 10) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<3,10,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 3 && cmax_ == 11) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<3,11,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 3 && cmax_ == 12) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<3,12,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 4 && cmax_ == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<4,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 4 && cmax_ == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<4,1,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 4 && cmax_ == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<4,2,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 4 && cmax_ == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<4,3,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 4 && cmax_ == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<4,4,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 4 && cmax_ == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<4,5,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 4 && cmax_ == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<4,6,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 4 && cmax_ == 7) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<4,7,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 4 && cmax_ == 8) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<4,8,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 4 && cmax_ == 9) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<4,9,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 4 && cmax_ == 10) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<4,10,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 4 && cmax_ == 11) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<4,11,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 4 && cmax_ == 12) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<4,12,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 5 && cmax_ == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<5,0,3>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 5 && cmax_ == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<5,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 5 && cmax_ == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<5,2,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 5 && cmax_ == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<5,3,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 5 && cmax_ == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<5,4,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 5 && cmax_ == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<5,5,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 5 && cmax_ == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<5,6,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 5 && cmax_ == 7) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<5,7,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 5 && cmax_ == 8) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<5,8,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 5 && cmax_ == 9) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<5,9,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 5 && cmax_ == 10) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<5,10,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 5 && cmax_ == 11) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<5,11,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 5 && cmax_ == 12) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<5,12,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 6 && cmax_ == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<6,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 6 && cmax_ == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<6,1,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 6 && cmax_ == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<6,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 6 && cmax_ == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<6,3,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 6 && cmax_ == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<6,4,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 6 && cmax_ == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<6,5,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 6 && cmax_ == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<6,6,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 6 && cmax_ == 7) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<6,7,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 6 && cmax_ == 8) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<6,8,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 6 && cmax_ == 9) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<6,9,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 6 && cmax_ == 10) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<6,10,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 6 && cmax_ == 11) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<6,11,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 6 && cmax_ == 12) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<6,12,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 7 && cmax_ == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<7,0,4>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 7 && cmax_ == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<7,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 7 && cmax_ == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<7,2,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 7 && cmax_ == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<7,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 7 && cmax_ == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<7,4,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 7 && cmax_ == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<7,5,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 7 && cmax_ == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<7,6,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 7 && cmax_ == 7) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<7,7,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 7 && cmax_ == 8) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<7,8,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 7 && cmax_ == 9) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<7,9,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 7 && cmax_ == 10) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<7,10,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 7 && cmax_ == 11) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<7,11,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 7 && cmax_ == 12) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<7,12,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 8 && cmax_ == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<8,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 8 && cmax_ == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<8,1,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 8 && cmax_ == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<8,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 8 && cmax_ == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<8,3,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 8 && cmax_ == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<8,4,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 8 && cmax_ == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<8,5,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 8 && cmax_ == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<8,6,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 8 && cmax_ == 7) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<8,7,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 8 && cmax_ == 8) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<8,8,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 8 && cmax_ == 9) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<8,9,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 8 && cmax_ == 10) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<8,10,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 8 && cmax_ == 11) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<8,11,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 8 && cmax_ == 12) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<8,12,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 9 && cmax_ == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<9,0,5>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 9 && cmax_ == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<9,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 9 && cmax_ == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<9,2,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 9 && cmax_ == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<9,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 9 && cmax_ == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<9,4,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 9 && cmax_ == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<9,5,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 9 && cmax_ == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<9,6,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 9 && cmax_ == 7) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<9,7,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 9 && cmax_ == 8) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<9,8,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 9 && cmax_ == 9) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<9,9,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 9 && cmax_ == 10) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<9,10,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 9 && cmax_ == 11) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<9,11,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 9 && cmax_ == 12) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<9,12,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 10 && cmax_ == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<10,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 10 && cmax_ == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<10,1,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 10 && cmax_ == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<10,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 10 && cmax_ == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<10,3,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 10 && cmax_ == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<10,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 10 && cmax_ == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<10,5,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 10 && cmax_ == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<10,6,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 10 && cmax_ == 7) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<10,7,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 10 && cmax_ == 8) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<10,8,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 10 && cmax_ == 9) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<10,9,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 10 && cmax_ == 10) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<10,10,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 10 && cmax_ == 11) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<10,11,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 10 && cmax_ == 12) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<10,12,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 11 && cmax_ == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<11,0,6>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 11 && cmax_ == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<11,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 11 && cmax_ == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<11,2,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 11 && cmax_ == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<11,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 11 && cmax_ == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<11,4,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 11 && cmax_ == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<11,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 11 && cmax_ == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<11,6,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 11 && cmax_ == 7) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<11,7,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 11 && cmax_ == 8) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<11,8,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 11 && cmax_ == 9) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<11,9,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 11 && cmax_ == 10) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<11,10,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 11 && cmax_ == 11) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<11,11,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 11 && cmax_ == 12) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<11,12,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 12 && cmax_ == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<12,0,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 12 && cmax_ == 1) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<12,1,7>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 12 && cmax_ == 2) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<12,2,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 12 && cmax_ == 3) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<12,3,8>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 12 && cmax_ == 4) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<12,4,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 12 && cmax_ == 5) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<12,5,9>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 12 && cmax_ == 6) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<12,6,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 12 && cmax_ == 7) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<12,7,10>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 12 && cmax_ == 8) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<12,8,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 12 && cmax_ == 9) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<12,9,11>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 12 && cmax_ == 10) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<12,10,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 12 && cmax_ == 11) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<12,11,12>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  } else if (amax_ == 12 && cmax_ == 12) {
    for (int j = 0; j != screening_size_; ++j) {
      int ii = screening_[j];
      vrr_driver<12,12,13>(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);
    }
  }

}
