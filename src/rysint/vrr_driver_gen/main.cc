//
// Author: Toru Shiozaki
// Date  : April 2009
//

#include <string>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <fstream>
#include "../macros.h"

using namespace std;
using namespace boost;

int main() {

ofstream ofs;
ofs.open("vrr.cc");

string header = "\
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
#include <src/util/f77.h>\n\
\n\
using namespace std;\n\
using namespace bagel;\n\
\n\
";
ofs << header << endl;

for (int r = 4; r <= RYS_MAX; ++r) {
string out;
string rank = lexical_cast<string>(r);
out += "\
void ERIBatch::perform_VRR" + rank + "() {\n\
  const int isize = (amax_ + 1) * (cmax_ + 1);\n\
  const int worksize = " + rank + " * isize;\n\
  const int vrr_index = amax_ * ANG_VRR_END + cmax_;\n\
\n\
  double* const workx = stack_->get(worksize);\n\
  double* const worky = stack_->get(worksize);\n\
  double* const workz = stack_->get(worksize);\n\
  double iyiz[" + rank + "];\n\
\n\
  const int acsize = asize_ * csize_;\n\
  const double ax = basisinfo_[0]->position(0);\n\
  const double ay = basisinfo_[0]->position(1);\n\
  const double az = basisinfo_[0]->position(2);\n\
  const double bx = basisinfo_[1]->position(0);\n\
  const double by = basisinfo_[1]->position(1);\n\
  const double bz = basisinfo_[1]->position(2);\n\
  const double cx = basisinfo_[2]->position(0);\n\
  const double cy = basisinfo_[2]->position(1);\n\
  const double cz = basisinfo_[2]->position(2);\n\
  const double dx = basisinfo_[3]->position(0);\n\
  const double dy = basisinfo_[3]->position(1);\n\
  const double dz = basisinfo_[3]->position(2);\n\
  for (int j = 0; j != screening_size_; ++j) {\n\
    int ii = screening_[j];\n\
    int offset = ii * " + rank + ";\n\
    int data_offset_ii = ii * acsize;\n\
\n\
    double* current_data = &data_[data_offset_ii];\n\
\n\
    const int ii3 = 3 * ii;\n\
    const double cxp = xp_[ii];\n\
    const double cxq = xq_[ii];\n\
    const double oxp2 = 0.5 / cxp;\n\
    const double oxq2 = 0.5 / cxq;\n\
    const double opq = 1.0 / (cxp + cxq);\n\
    const array<double, 11> dparamx = {{p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};\n\
    Int2D cix(dparamx, &roots_[offset], " + rank + ", worksize, workx, vrr_->vrrfunc[vrr_index]);\n\
    cix.scale_data(&weights_[offset], coeff_[ii]);\n\
\n\
    const array<double, 11> dparamy = {{p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};\n\
    Int2D ciy(dparamy, &roots_[offset], " + rank + ", worksize, worky, vrr_->vrrfunc[vrr_index]);\n\
\n\
    const array<double, 11> dparamz = {{p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};\n\
    Int2D ciz(dparamz, &roots_[offset], " + rank + ", worksize, workz, vrr_->vrrfunc[vrr_index]);\n\
\n\
    for (int iz = 0; iz <= cmax_; ++iz) {\n\
      for (int iy = 0; iy <= cmax_ - iz; ++iy) {\n\
        const int iyz = cmax1_ * (iy + cmax1_ * iz);\n\
        for (int jz = 0; jz <= amax_; ++jz) {\n\
          const int offsetz = " + rank + " * (amax1_ * iz + jz);\n\
          for (int jy = 0; jy <= amax_ - jz; ++jy) {\n\
            const int offsety = " + rank + " * (amax1_ * iy + jy);\n\
            const int jyz = amax1_ * (jy + amax1_ * jz);\n\
            for (int i = 0; i != " + rank + "; ++i)\n\
              iyiz[i] = worky[offsety + i] * workz[offsetz + i];\n";
out +="\
            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) {\n\
              const int iposition = cmapping_[ix + iyz];\n\
              const int ipos_asize = iposition * asize_;\n\
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) {\n\
                const int offsetx = " + rank + " * (amax1_ * ix + jx);\n\
                const int jposition = amapping_[jx + jyz];\n\
                const int ijposition = jposition + ipos_asize;\n\
                current_data[ijposition] = ddot_(" + rank + ", iyiz, 1, workx+offsetx, 1);\n\
              }\n\
            }\n\
          }\n\
        }\n\
      }\n\
    }\n\
\n\
  }\n\
\n\
  stack_->release(worksize, workz);\n\
  stack_->release(worksize, worky);\n\
  stack_->release(worksize, workx);\n\
}\n\n";

ofs << out << endl;

}

ofs.close();

}
