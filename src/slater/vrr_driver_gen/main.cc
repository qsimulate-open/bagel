//
// BAGEL - Parallel electron correlation program.
// Filename: main.cc
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


#include <sstream>
#include <fstream>
using namespace std;

int main() {
  ofstream ofs;
  ofs.open("svrr.cc");
  string header = "\
//\n\
// BAGEL - Parallel electron correlation program.\n\
// Filename: svrr.cc\n\
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
#include <src/slater/slaterbatch.h>\n\
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

  stringstream ss;
  ss << header << endl;

  for (int rank = 3; rank != 14; ++rank) {
    ss << "void SlaterBatch::perform_SVRR" << rank << "() {" << endl;
    ss << "  const int isize = (amax_ + 1) * (cmax_ + 1);" << endl;
    ss << "  const int worksize = " << rank << " * isize;" << endl;
    ss << "  const int svrr_index = amax_ * ANG_VRR_END + cmax_;" << endl;
    ss << "" << endl;
    ss << "  const int acsize = asize_ * csize_;" << endl;
    ss << "" << endl;
    ss << "  double* workx = stack_->get(worksize);" << endl;
    ss << "  double* worky = stack_->get(worksize);" << endl;
    ss << "  double* workz = stack_->get(worksize);" << endl;
    ss << "  double iyiz[RYS_MAX];" << endl;
    ss << "" << endl;
    ss << "  double womt[RYS_MAX];" << endl;
    ss << "  double wt[RYS_MAX];" << endl;
    ss << "" << endl;
    ss << "  const double ax = basisinfo_[0]->position(0);" << endl;
    ss << "  const double ay = basisinfo_[0]->position(1);" << endl;
    ss << "  const double az = basisinfo_[0]->position(2);" << endl;
    ss << "  const double bx = basisinfo_[1]->position(0);" << endl;
    ss << "  const double by = basisinfo_[1]->position(1);" << endl;
    ss << "  const double bz = basisinfo_[1]->position(2);" << endl;
    ss << "  const double cx = basisinfo_[2]->position(0);" << endl;
    ss << "  const double cy = basisinfo_[2]->position(1);" << endl;
    ss << "  const double cz = basisinfo_[2]->position(2);" << endl;
    ss << "  const double dx = basisinfo_[3]->position(0);" << endl;
    ss << "  const double dy = basisinfo_[3]->position(1);" << endl;
    ss << "  const double dz = basisinfo_[3]->position(2);" << endl;
    ss << "" << endl;
    ss << "  // Perform VRR" << endl;
    ss << "  for (int j = 0; j != screening_size_; ++j) {" << endl;
    ss << "    const int ii = screening_[j];" << endl;
    ss << "    int offset = ii * " << rank << ";" << endl;
    ss << "    int data_offset_ii = ii * acsize;" << endl;
    ss << "" << endl;
    ss << "    double* current_data = data_ + data_offset_ii;" << endl;
    ss << "" << endl;
    ss << "    /// workx, worky, workz would be Ix, Iy, and Iz data" << endl;
    ss << "    const int ii3 = 3 * ii;" << endl;
    ss << "    const double cxp = xp_[ii];" << endl;
    ss << "    const double cxq = xq_[ii];" << endl;
    ss << "    const double oxp2 = 0.5 / cxp; " << endl;
    ss << "    const double oxq2 = 0.5 / cxq; " << endl;
    ss << "    const double opq = 1.0 / (cxp + cxq);" << endl;
    ss << "    const array<double, 11> dparamx = {{p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};" << endl;
    ss << "    Int2D<" << rank << "> cix(dparamx, roots_+offset, worksize, workx, vrr_->vrrfunc[svrr_index]);" << endl;
    ss << "" << endl;
    ss << "    for (int i = 0; i != " << rank << "; ++i) {" << endl;
    ss << "      wt[i] = weights_[offset + i] * roots_[offset + i];" << endl;
    ss << "      womt[i] = weights_[offset + i] - wt[i];" << endl;
    ss << "    }    " << endl;
    ss << "    cix.scale_data(womt, coeff_[ii]);" << endl;
    ss << "" << endl;
    ss << "    const array<double, 11> dparamy = {{p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};" << endl;
    ss << "    Int2D<" << rank << "> ciy(dparamy, roots_+offset, worksize, worky, vrr_->vrrfunc[svrr_index]);" << endl;
    ss << "" << endl;
    ss << "    const array<double, 11> dparamz = {{p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};" << endl;
    ss << "    Int2D<" << rank << "> ciz(dparamz, roots_+offset, worksize, workz, vrr_->vrrfunc[svrr_index]);" << endl;
    ss << "" << endl;
    ss << "    /// assembly process" << endl;
    ss << "" << endl;
    ss << "    for (int iz = 0; iz <= cmax_; ++iz) {" << endl;
    ss << "      for (int iy = 0; iy <= cmax_ - iz; ++iy) {" << endl;
    ss << "        const int iyz = cmax1_ * (iy + cmax1_ * iz); " << endl;
    ss << "        for (int jz = 0; jz <= amax_; ++jz) {" << endl;
    ss << "          const int offsetz = " << rank << " * ((amax_ + 1) * iz + jz); " << endl;
    ss << "          for (int jy = 0; jy <= amax_-jz; ++jy) {" << endl;
    ss << "            const int offsety = " << rank << " * ((amax_ + 1) * iy + jy); " << endl;
    ss << "            const int jyz = amax1_ * (jy + amax1_ * jz); " << endl;
    ss << "            for (int i = 0; i != " << rank << "; ++i) iyiz[i] = worky[offsety + i] * workz[offsetz + i];" << endl;
    ss << "" << endl;
    ss << "            for (int ix = max(0, cmin_-iy-iz); ix <= cmax_-iy-iz; ++ix) {" << endl;
    ss << "              const int iposition = cmapping_[ix + iyz];" << endl;
    ss << "              const int ipos_asize = iposition * asize_;" << endl;
    ss << "              for (int jx = max(0, amin_-jy-jz); jx <= amax_-jy-jz; ++jx) {" << endl;
    ss << "                const int jposition = amapping_[jx + jyz];" << endl;
    ss << "                const int ijposition = jposition + ipos_asize;" << endl;
    ss << "" << endl;
    ss << "                const int offsetx = " << rank << " * ((amax_ + 1) * ix + jx); " << endl;
    ss << "                current_data[ijposition] = 0.0; " << endl;
    ss << "                for (int i = 0; i != " << rank << "; ++i) {" << endl;
    ss << "                  current_data[ijposition] += iyiz[i] * workx[offsetx + i];" << endl;
    ss << "                }    " << endl;
    ss << "              }    " << endl;
    ss << "            }    " << endl;
    ss << "          }    " << endl;
    ss << "        }    " << endl;
    ss << "      }    " << endl;
    ss << "    }    " << endl;
    ss << "" << endl;
    ss << "  } // end of primsize loop" << endl;
    ss << "" << endl;
    ss << "  stack_->release(worksize, workz);" << endl;
    ss << "  stack_->release(worksize, worky);" << endl;
    ss << "  stack_->release(worksize, workx);" << endl;
    ss << "}" << endl;
    ss << "" << endl;
    ss << endl;
    ss << endl;
  }
  ofs << ss.str();
  ofs.close();
}
