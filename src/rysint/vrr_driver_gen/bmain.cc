//
// Author: Toru Shiozaki
// Date  : April 2009
//

#include <string>
#include <fstream>
#include <sstream>

using namespace std;

int main() {

  ofstream ofs;
  ofs.open("bvrr.cc");

  string header = "\
//\n\
// BAGEL - Parallel electron correlation program.\n\
// Filename: bvrr.cc\n\
// Copyright (C) 2012 Toru Shiozaki\n\
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
#include <src/rysint/breitbatch.h>\n\
#include <src/rysint/int2d.h>\n\
#include <src/parallel/resources.h>\n\
#include <src/util/f77.h>\n\
\n\
using namespace std;\n\
using namespace bagel;\n\
\n\
";
  ofs << header << endl;

  stringstream ss;

  for (int rank = 1; rank <= 13; ++rank) {
    ss << "void BreitBatch::perform_VRR" << rank << "() {" << endl;
    ss << "  const int amax2 = amax1_+1;" << endl;
    ss << "  const int cmax2 = cmax1_+1;" << endl;
    ss << "  const int isize = amax2 * cmax2;" << endl;
    ss << "  const int worksize = " << rank << " * isize;" << endl;
    ss << "  // CAUTION we need up to amax1_, cmax1_" << endl;
    ss << "  const int vrr_index = amax1_ * ANG_VRR_END + cmax1_;" << endl;
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
    ss << "  double* const workx = stack_->get(worksize*3);" << endl;
    ss << "  double* const worky = workx + worksize;" << endl;
    ss << "  double* const workz = worky + worksize;" << endl;
    ss << "  double* const worktx = stack_->get(worksize*3);" << endl;
    ss << "  double* const workty = worktx + worksize;" << endl;
    ss << "  double* const worktz = workty + worksize;" << endl;
    ss << "  double* const worksx = stack_->get(worksize*3);" << endl;
    ss << "  double* const worksy = worksx + worksize;" << endl;
    ss << "  double* const worksz = worksy + worksize;" << endl;
    ss << "" << endl;
    ss << "  const int acsize = size_block_ / primsize_;" << endl;
    ss << "" << endl;
    ss << "  for (int j = 0; j != screening_size_; ++j) {" << endl;
    ss << "    const int ii = screening_[j];" << endl;
    ss << "    const size_t offset = ii * " << rank << ";" << endl;
    ss << "    const size_t data_offset_ii = ii * acsize;" << endl;
    ss << "" << endl;
    ss << "    const int ii3 = 3 * ii; " << endl;
    ss << "    const double cxp = xp_[ii];" << endl;
    ss << "    const double cxq = xq_[ii];" << endl;
    ss << "    const double oxp2 = 0.5 / cxp;" << endl;
    ss << "    const double oxq2 = 0.5 / cxq;" << endl;
    ss << "    const double opq = 1.0 / (cxp + cxq);" << endl;
    ss << "" << endl;
    ss << "    const array<double, 11> dparamx = {{p_[ii3],   q_[ii3],   ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};" << endl;
    ss << "    Int2D<" << rank << "> cix(dparamx, &roots_[offset], worksize, workx, vrr_->vrrfunc[vrr_index]);" << endl;
    ss << "    cix.scale_data(&weights_[offset], coeff_[ii]);" << endl;
    ss << "" << endl;
    ss << "    const array<double, 11> dparamy = {{p_[ii3+1], q_[ii3+1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};" << endl;
    ss << "    Int2D<" << rank << "> ciy(dparamy, &roots_[offset], worksize, worky, vrr_->vrrfunc[vrr_index]);" << endl;
    ss << "" << endl;
    ss << "    const array<double, 11> dparamz = {{p_[ii3+2], q_[ii3+2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};" << endl;
    ss << "    Int2D<" << rank << "> ciz(dparamz, &roots_[offset], worksize, workz, vrr_->vrrfunc[vrr_index]);" << endl;
    ss << "" << endl;
    ss << "    const double pq[3] = {p_[ii3]-q_[ii3], p_[ii3+1]-q_[ii3+1], p_[ii3+2]-q_[ii3+2]};" << endl;
    ss << "" << endl;
    ss << "    // next compute \tidle{I}_x,y,z up to amax_, cmax_" << endl;
    ss << "    for (int ic = 0; ic <= cmax1_; ++ic)" << endl;
    ss << "      for (int ia = 0; ia <= amax1_; ++ia)" << endl;
    ss << "        for (int i = 0; i != " << rank << "; ++i) {" << endl;
    ss << "          worktx[i+" << rank << "*(ia+amax2*ic)] = pq[0]*workx[i+" << rank << "*(ia+amax2*ic)] + (ia==0 ? 0.0 : oxp2*workx[i+" << rank << "*(ia+amax2*ic)]) - (ic==0 ? 0.0 : oxq2*workx[i+" << rank << "*(ia+amax2*(ic-1))]); " << endl;
    ss << "          workty[i+" << rank << "*(ia+amax2*ic)] = pq[1]*worky[i+" << rank << "*(ia+amax2*ic)] + (ia==0 ? 0.0 : oxp2*worky[i+" << rank << "*(ia+amax2*ic)]) - (ic==0 ? 0.0 : oxq2*worky[i+" << rank << "*(ia+amax2*(ic-1))]); " << endl;
    ss << "          worktz[i+" << rank << "*(ia+amax2*ic)] = pq[2]*workz[i+" << rank << "*(ia+amax2*ic)] + (ia==0 ? 0.0 : oxp2*workz[i+" << rank << "*(ia+amax2*ic)]) - (ic==0 ? 0.0 : oxq2*workz[i+" << rank << "*(ia+amax2*(ic-1))]); " << endl;
    ss << "        }   " << endl;
    ss << "    // then compute \tilde{\tilde{I}}_x,y,z up to amax_-1, cmax_-1" << endl;
    ss << "    for (int ic = 0; ic != cmax1_; ++ic)" << endl;
    ss << "      for (int ia = 0; ia != amax1_; ++ia)" << endl;
    ss << "        for (int i = 0; i != " << rank << "; ++i) {" << endl;
    ss << "          worksx[i+" << rank << "*(ia+amax2*ic)] = worktx[i+" << rank << "*((ia+1)+amax2*ic)] - worktx[i+" << rank << "*(ia+amax2*(ic+1))] + (ax - cx)*worktx[i+" << rank << "*(ia+amax2*ic)];" << endl;
    ss << "          worksy[i+" << rank << "*(ia+amax2*ic)] = workty[i+" << rank << "*((ia+1)+amax2*ic)] - workty[i+" << rank << "*(ia+amax2*(ic+1))] + (ay - cy)*workty[i+" << rank << "*(ia+amax2*ic)];" << endl;
    ss << "          worksz[i+" << rank << "*(ia+amax2*ic)] = worktz[i+" << rank << "*((ia+1)+amax2*ic)] - worktz[i+" << rank << "*(ia+amax2*(ic+1))] + (az - cz)*worktz[i+" << rank << "*(ia+amax2*ic)];" << endl;
    ss << "        }   " << endl;
    ss << "" << endl;
    ss << "    double* const dataxx = &data_[data_offset_ii];" << endl;
    ss << "    double* const dataxy = dataxx + size_block_; " << endl;
    ss << "    double* const datayy = dataxy + size_block_; " << endl;
    ss << "    double* const dataxz = datayy + size_block_; " << endl;
    ss << "    double* const datayz = dataxz + size_block_; " << endl;
    ss << "    double* const datazz = datayz + size_block_; " << endl;
    ss << "" << endl;
    ss << "    double* const iyiz_nn = stack_->get(" << rank << "*6);" << endl;
    ss << "    double* const iyiz_tn = iyiz_nn + " << rank << ";" << endl;
    ss << "    double* const iyiz_nt = iyiz_tn + " << rank << ";" << endl;
    ss << "    double* const iyiz_tt = iyiz_nt + " << rank << ";" << endl;
    ss << "    double* const iyiz_sn = iyiz_tt + " << rank << ";" << endl;
    ss << "    double* const iyiz_ns = iyiz_sn + " << rank << ";" << endl;
    ss << "" << endl;
    ss << "" << endl;
    ss << "    // assemble up to amax_, cmax_" << endl;
    ss << "    for (int iz = 0; iz <= cmax_; ++iz) {" << endl;
    ss << "      for (int iy = 0; iy <= cmax_ - iz; ++iy) {" << endl;
    ss << "        for (int jz = 0; jz <= amax_; ++jz) {" << endl;
    ss << "          for (int jy = 0; jy <= amax_ - jz; ++jy) {" << endl;
    ss << "            const int offsetz = " << rank << " * (amax2 * iz + jz);" << endl;
    ss << "            const int offsety = " << rank << " * (amax2 * iy + jy);" << endl;
    ss << "" << endl;
    ss << "            const int iyz = cmax1_ * (iy + cmax1_ * iz);" << endl;
    ss << "            const int jyz = amax1_ * (jy + amax1_ * jz);" << endl;
    ss << "" << endl;
    ss << "            for (int i = 0; i != " << rank << "; ++i) {" << endl;
    ss << "              iyiz_nn[i] = worky [offsety + i] * workz [offsetz + i]; " << endl;
    ss << "              iyiz_tn[i] = workty[offsety + i] * workz [offsetz + i] * (1.0-roots_[offset+i]);" << endl;
    ss << "              iyiz_nt[i] = worky [offsety + i] * worktz[offsetz + i] * (1.0-roots_[offset+i]);" << endl;
    ss << "              iyiz_tt[i] = workty[offsety + i] * worktz[offsetz + i] * (1.0-roots_[offset+i]);" << endl;
    ss << "              iyiz_sn[i] = worksy[offsety + i] * workz [offsetz + i]; " << endl;
    ss << "              iyiz_ns[i] = worky [offsety + i] * worksz[offsetz + i]; " << endl;
    ss << "            }   " << endl;
    ss << "" << endl;
    ss << "            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) {" << endl;
    ss << "              const int iposition = cmapping_[ix + iyz];" << endl;
    ss << "              const int ipos_asize = iposition * asize_;" << endl;
    ss << "              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) {" << endl;
    ss << "                const int offsetx = " << rank << " * (amax2 * ix + jx);" << endl;
    ss << "                const int jposition = amapping_[jx + jyz];" << endl;
    ss << "                const int ijposition = jposition + ipos_asize;" << endl;
    ss << "" << endl;
    ss << "                dataxx[ijposition] = ddot_(" << rank << ", iyiz_nn, 1, worksx+offsetx, 1); " << endl;
    ss << "                dataxy[ijposition] = ddot_(" << rank << ", iyiz_tn, 1, worktx+offsetx, 1); " << endl;
    ss << "                dataxz[ijposition] = ddot_(" << rank << ", iyiz_nt, 1, worktx+offsetx, 1); " << endl;
    ss << "                datayy[ijposition] = ddot_(" << rank << ", iyiz_sn, 1, workx +offsetx, 1); " << endl;
    ss << "                datazz[ijposition] = ddot_(" << rank << ", iyiz_ns, 1, workx +offsetx, 1); " << endl;
    ss << "                datayz[ijposition] = ddot_(" << rank << ", iyiz_tt, 1, workx +offsetx, 1); " << endl;
    ss << "              }   " << endl;
    ss << "            }   " << endl;
    ss << "          }   " << endl;
    ss << "        }   " << endl;
    ss << "      }   " << endl;
    ss << "    }   " << endl;
    ss << "" << endl;
    ss << "    stack_->release(" << rank << "*6, iyiz_nn);" << endl;
    ss << "  }" << endl;
    ss << "" << endl;
    ss << "  stack_->release(worksize*3, worksx);" << endl;
    ss << "  stack_->release(worksize*3, worktx);" << endl;
    ss << "  stack_->release(worksize*3, workx);" << endl;
    ss << "" << endl;
    ss << "}" << endl;
    ss << endl;
    ss << endl;
  }

  ofs << ss.str();
  ofs.close();

}
