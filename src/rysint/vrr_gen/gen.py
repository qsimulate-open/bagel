#!/usr/bin/python

import math

ss = "\
//\n\
// BAGEL - Parallel electron correlation program.\n\
// Filename: vrr.cc\n\
// Copyright (C) 2009 Toru Shiozaki\n\
//\n\
// Author: Toru Shiozaki <shiozaki@northwestern.edu>\n\
// Maintainer: Shiozaki group\n\
//\n\
// This file is part of the BAGEL package.\n\
//\n\
// The BAGEL package is free software; you can redistribute it and/or modify\n\
// it under the terms of the GNU Library General Public License as published by\n\
// the Free Software Foundation; either version 2, or (at your option)\n\
// any later version.\n\
//\n\
// The BAGEL package is distributed in the hope that it will be useful,\n\
// but WITHOUT ANY WARRANTY; without even the implied warranty of\n\
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n\
// GNU Library General Public License for more details.\n\
//\n\
// You should have received a copy of the GNU Library General Public License\n\
// along with the BAGEL package; see COPYING.  If not, write to\n\
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.\n\
//\n\
\n\
#include <src/rysint/eribatch.h>\n\
#include <src/rysint/int2d.h>\n\
#include <cmath>\n\
#include <algorithm>\n\
#include <cstring>\n\
#include <src/parallel/resources.h>\n\
#include <src/rysint/_vrr_drv.h>\n\
#include <src/util/f77.h>\n\
\n\
using namespace std;\n\
using namespace bagel;\n\
\n\
\n\
void ERIBatch::perform_VRR() {\n\
  const int acsize = asize_ * csize_;\n"

for a in range(0,13):
  for c in range(0,13):
    rank = int(math.ceil((a+c+1)*0.5-0.001))

    if a == 0 and c == 0: 
      ss += "\
  if (amax_ == " + str(a) + " && cmax_ == " + str(c) + ") {\n\
    for (int j = 0; j != screening_size_; ++j) {\n\
      int ii = screening_[j];\n\
      vrr_driver<" + str(a) + "," + str(c) + "," +  str(rank) + ">(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],\n\
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),\n\
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);\n\
    }\n"
    else:
      ss += "\
  } else if (amax_ == " + str(a) + " && cmax_ == " + str(c) + ") {\n\
    for (int j = 0; j != screening_size_; ++j) {\n\
      int ii = screening_[j];\n\
      vrr_driver<" + str(a) + "," + str(c) + "," +  str(rank) + ">(data_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii],\n\
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),\n\
                    p_+ii*3, q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, amin_, cmin_, asize_);\n\
    }\n"

ss += "\
  }\n\
\n\
}"


print ss;
