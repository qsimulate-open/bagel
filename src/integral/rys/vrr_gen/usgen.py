#!/usr/bin/python

import math

filename = "usvrr.cc"

ss = "\
//\n\
// BAGEL - Parallel electron correlation program.\n\
// Filename: " + filename + "\n\
// Copyright (C) 2009 Toru Shiozaki\n\
//\n\
// Author: Toru Shiozaki <shiozaki@northwestern.edu>\n\
// Maintainer: Shiozaki group\n\
//\n\
// This file is part of the BAGEL package.\n\
//\n\
// The BAGEL package is free software; you can redistribute it and/or modify\n\
// it under the terms of the GNU Library General Public License as published by\n\
// the Free Software Foundation; either version 3, or (at your option)\n\
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
#include <src/integral/rys/slaterbatch.h>\n\
#include <src/integral/rys/_svrr_drv.h>\n\
\n\
using namespace std;\n\
using namespace bagel;\n\
\n\
#ifdef HAVE_LIBSLATER\n\
\n\
void SlaterBatch::perform_USVRR() {\n\
  const int acsize = asize_ * csize_;\n\
  const int a = basisinfo_[0]->angular_number();\n\
  const int b = basisinfo_[1]->angular_number();\n\
  const int c = basisinfo_[2]->angular_number();\n\
  const int d = basisinfo_[3]->angular_number();\n\
  const int isize = (amax_+1) * (cmax_+1);\n\
  double* const workx = stack_->get(isize*rank_*4);\n\
  double* const worky = workx + isize*rank_;\n\
  double* const workz = worky + isize*rank_;\n\
  double* const workx2 = workz + isize*rank_;\n\
  const int hashkey = (a << 24) + (b << 16) + (c << 8) + d;\n"

for a in range(0,7):
 for b in range(0,7):
  if a < b: continue
  for c in range(0,7):
   for d in range(0,7):
    if c < d: continue
    rank = int(math.ceil((a+b+c+d+2)*0.5-0.001))
    off = 1 << 8
    key = d+off*(c+off*(b+off*a))

    if a == 0 and c == 0:
     ss += "\
  switch (hashkey) {\n"
    ss += "\
  case " + str(key) + " :\n\
    for (int j = 0; j != screening_size_; ++j) {\n\
      int ii = screening_[j];\n\
      usvrr_driver<" + str(a) + "," + str(b) + "," + str(c) + "," +  str(d) + "," + str(rank) + ">(data_+ii*acsize, data2_+ii*acsize, roots_+ii*rank_, weights_+ii*rank_, coeff_[ii], coeffy_[ii],\n\
                    basisinfo_[0]->position(), basisinfo_[1]->position(), basisinfo_[2]->position(), basisinfo_[3]->position(),\n\
                    P_+ii*3, Q_+ii*3, xp_[ii], xq_[ii], amapping_, cmapping_, asize_, workx, worky, workz, workx2);\n\
    } break;\n"
ss += "\
  }\n\
  stack_->release(rank_*isize*4, workx);\n\
\n\
}\n\
\n\
#endif"


f = open(filename, "w")
f.write(ss)
