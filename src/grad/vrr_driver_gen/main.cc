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
  ofs.open("gvrr.cc");
  stringstream ss;

  ss << "//" << endl;
  ss << "// BAGEL - Parallel electron correlation program." << endl;
  ss << "// Filename: gvrr.cc" << endl;
  ss << "// Copyright (C) 2012 Toru Shiozaki" << endl;
  ss << "//" << endl;
  ss << "// Author: Toru Shiozaki <shiozaki@northwestern.edu>" << endl;
  ss << "// Maintainer: Shiozaki group" << endl;
  ss << "//" << endl;
  ss << "// This file is part of the BAGEL package." << endl;
  ss << "//" << endl;
  ss << "// The BAGEL package is free software; you can redistribute it and/or modify" << endl;
  ss << "// it under the terms of the GNU Library General Public License as published by" << endl;
  ss << "// the Free Software Foundation; either version 2, or (at your option)" << endl;
  ss << "// any later version." << endl;
  ss << "//" << endl;
  ss << "// The BAGEL package is distributed in the hope that it will be useful," << endl;
  ss << "// but WITHOUT ANY WARRANTY; without even the implied warranty of" << endl;
  ss << "// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << endl;
  ss << "// GNU Library General Public License for more details." << endl;
  ss << "//" << endl;
  ss << "// You should have received a copy of the GNU Library General Public License" << endl;
  ss << "// along with the BAGEL package; see COPYING.  If not, write to" << endl;
  ss << "// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA." << endl;
  ss << "//" << endl;
  ss << "" << endl;
  ss << "" << endl;
  ss << "#include <src/grad/gradbatch.h>" << endl;
  ss << "#include <src/rysint/int2d.h>" << endl;
  ss << "#include <src/util/f77.h>" << endl;
  ss << "#include <src/util/comb.h>" << endl;
  ss << "" << endl;
  ss << "using namespace std;" << endl;
  ss << "using namespace bagel;" << endl;
  ss << "" << endl;
  ss << "static const Comb comb;" << endl;
  ss << "" << endl;

  for (int rank = 1; rank != 14; ++rank) {
    ss << "void GradBatch::perform_VRR" << rank << "() {" << endl;
    ss << "  const int isize = (amax_ + 1) * (cmax_ + 1); " << endl;
    ss << "  const int worksize = " << rank << " * isize;" << endl;
    ss << "  const int vrr_index = amax_ * ANG_VRR_END + cmax_;" << endl;
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
    ss << "  // transformation matrix" << endl;
    ss << "  const int a = basisinfo_[0]->angular_number();" << endl;
    ss << "  const int b = basisinfo_[1]->angular_number();" << endl;
    ss << "  const int c = basisinfo_[2]->angular_number();" << endl;
    ss << "  const int d = basisinfo_[3]->angular_number();" << endl;
    ss << "  assert(a+b+1 == amax_ && c+d+1 == cmax_);" << endl;
    ss << "" << endl;
    ss << "  const int a2 = a+2;" << endl;
    ss << "  const int b2 = b+2;" << endl;
    ss << "  const int c2 = c+2;" << endl;
    ss << "  const int d2 = d+2;" << endl;
    ss << "" << endl;
    ss << "  double* const workx = stack_->get(worksize*3);" << endl;
    ss << "  double* const worky = workx + worksize;" << endl;
    ss << "  double* const workz = worky + worksize;" << endl;
    ss << "  double* const transx = stack_->get((amax_+1)*a2*b2);" << endl;
    ss << "  double* const transy = stack_->get((amax_+1)*a2*b2);" << endl;
    ss << "  double* const transz = stack_->get((amax_+1)*a2*b2);" << endl;
    ss << "  double* const trans2x = stack_->get((cmax_+1)*c2*d2);" << endl;
    ss << "  double* const trans2y = stack_->get((cmax_+1)*c2*d2);" << endl;
    ss << "  double* const trans2z = stack_->get((cmax_+1)*c2*d2);" << endl;
    ss << "  fill(transx,  transx +(amax_+1)*a2*b2, 0.0);" << endl;
    ss << "  fill(transy,  transy +(amax_+1)*a2*b2, 0.0);" << endl;
    ss << "  fill(transz,  transz +(amax_+1)*a2*b2, 0.0);" << endl;
    ss << "  fill(trans2x, trans2x+(cmax_+1)*c2*d2, 0.0);" << endl;
    ss << "  fill(trans2y, trans2y+(cmax_+1)*c2*d2, 0.0);" << endl;
    ss << "  fill(trans2z, trans2z+(cmax_+1)*c2*d2, 0.0);" << endl;
    ss << "  // for usual integrals" << endl;
    ss << "  for (int ib = 0, k = 0; ib <= b+1; ++ib) {" << endl;
    ss << "    for (int ia = 0; ia <= a+1; ++ia, ++k) {" << endl;
    ss << "      if (ia == a+1 && ib == b+1) continue;" << endl;
    ss << "      for (int i = ia; i <= ia+ib; ++i) {" << endl;
    ss << "        transx[i + (amax_+1)*k] = comb.c(ib, ia+ib-i) * pow(AB_[0], ia+ib-i);" << endl;
    ss << "        transy[i + (amax_+1)*k] = comb.c(ib, ia+ib-i) * pow(AB_[1], ia+ib-i);" << endl;
    ss << "        transz[i + (amax_+1)*k] = comb.c(ib, ia+ib-i) * pow(AB_[2], ia+ib-i);" << endl;
    ss << "      }   " << endl;
    ss << "    }   " << endl;
    ss << "  }" << endl;
    ss << "  for (int id = 0, k = 0; id <= d+1; ++id) {" << endl;
    ss << "    for (int ic = 0; ic <= c+1; ++ic, ++k) {" << endl;
    ss << "      if (ic == c+1 && id == d+1) continue;" << endl;
    ss << "      for (int i = ic; i <= ic+id; ++i) {" << endl;
    ss << "        trans2x[i + (cmax_+1)*k] = comb.c(id, ic+id-i) * pow(CD_[0], ic+id-i);" << endl;
    ss << "        trans2y[i + (cmax_+1)*k] = comb.c(id, ic+id-i) * pow(CD_[1], ic+id-i);" << endl;
    ss << "        trans2z[i + (cmax_+1)*k] = comb.c(id, ic+id-i) * pow(CD_[2], ic+id-i);" << endl;
    ss << "      }   " << endl;
    ss << "    }   " << endl;
    ss << "  }" << endl;
    ss << "  double* const intermediate = stack_->get(b2*a2*(cmax_+1)*" << rank << ");" << endl;
    ss << "  double* const final_x  = stack_->get(b2*a2*c2*d2*" << rank << ");" << endl;
    ss << "  double* const final_y  = stack_->get(b2*a2*c2*d2*" << rank << ");" << endl;
    ss << "  double* const final_z  = stack_->get(b2*a2*c2*d2*" << rank << ");" << endl;
    ss << "  double* const final_xa = stack_->get(b2*a2*c2*d2*" << rank << ");" << endl;
    ss << "  double* const final_xb = stack_->get(b2*a2*c2*d2*" << rank << ");" << endl;
    ss << "  double* const final_xc = stack_->get(b2*a2*c2*d2*" << rank << ");" << endl;
    ss << "  double* const final_ya = stack_->get(b2*a2*c2*d2*" << rank << ");" << endl;
    ss << "  double* const final_yb = stack_->get(b2*a2*c2*d2*" << rank << ");" << endl;
    ss << "  double* const final_yc = stack_->get(b2*a2*c2*d2*" << rank << ");" << endl;
    ss << "  double* const final_za = stack_->get(b2*a2*c2*d2*" << rank << ");" << endl;
    ss << "  double* const final_zb = stack_->get(b2*a2*c2*d2*" << rank << ");" << endl;
    ss << "  double* const final_zc = stack_->get(b2*a2*c2*d2*" << rank << ");" << endl;
    ss << "" << endl;
    ss << "  const int acsize = size_block_ / primsize_;" << endl;
    ss << "  assert(acsize == (a+1)*(b+1)*(c+1)*(d+1)*a2*b2*c2*d2/16 && size_block_*12 == size_alloc_);" << endl;
    ss << "" << endl;
    ss << "  for (int j = 0; j != screening_size_; ++j) {" << endl;
    ss << "    const int ii = screening_[j];" << endl;
    ss << "" << endl;
    ss << "    int offset = ii * " << rank << ";" << endl;
    ss << "    int data_offset_ii = ii * acsize;" << endl;
    ss << "    double* expo = exponents_.get() + ii*4;" << endl;
    ss << "" << endl;
    ss << "    const int ii3 = 3 * ii; " << endl;
    ss << "    const double cxp = xp_[ii];" << endl;
    ss << "    const double cxq = xq_[ii];" << endl;
    ss << "    const double oxp2 = 0.5 / cxp;" << endl;
    ss << "    const double oxq2 = 0.5 / cxq;" << endl;
    ss << "    const double opq = 1.0 / (cxp + cxq);" << endl;
    ss << "    const array<double, 11> dparamx = {{p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};" << endl;
    ss << "    Int2D<" << rank << "> cix(dparamx, &roots_[offset], worksize, workx, vrr_->vrrfunc[vrr_index]);" << endl;
    ss << "    cix.scale_data(&weights_[offset], coeff_[ii]);" << endl;
    ss << "" << endl;
    ss << "    // first (0:a+b, 0, 0:c+d, 0)-> (0:a+1, 0:b+1, 0:c+1, 0:d+1)" << endl;
    ss << "    for (int i = 0; i <= cmax_; ++i)" << endl;
    ss << "      dgemm_(\"N\", \"N\", " << rank << ", b2*a2, amax_+1, 1.0, workx+i*" << rank << "*(amax_+1), " << rank << ", transx, amax_+1, 0.0, intermediate+i*" << rank << "*b2*a2, " << rank << ");" << endl;
    ss << "    dgemm_(\"N\", \"N\", " << rank << "*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, " << rank << "*b2*a2, trans2x, cmax_+1, 0.0, final_x, " << rank << "*b2*a2);" << endl;
    ss << "" << endl;
    ss << "    const array<double, 11> dparamy = {{p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};" << endl;
    ss << "    Int2D<" << rank << "> ciy(dparamy, &roots_[offset], worksize, worky, vrr_->vrrfunc[vrr_index]);" << endl;
    ss << "" << endl;
    ss << "    for (int i = 0; i <= cmax_; ++i)" << endl;
    ss << "      dgemm_(\"N\", \"N\", " << rank << ", b2*a2, amax_+1, 1.0, worky+i*" << rank << "*(amax_+1), " << rank << ", transy, amax_+1, 0.0, intermediate+i*" << rank << "*b2*a2, " << rank << ");" << endl;
    ss << "    dgemm_(\"N\", \"N\", " << rank << "*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, " << rank << "*b2*a2, trans2y, cmax_+1, 0.0, final_y, " << rank << "*b2*a2);" << endl;
    ss << "" << endl;
    ss << "    const array<double, 11> dparamz = {{p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};" << endl;
    ss << "    Int2D<" << rank << "> ciz(dparamz, &roots_[offset], worksize, workz, vrr_->vrrfunc[vrr_index]);" << endl;
    ss << "" << endl;
    ss << "    for (int i = 0; i <= cmax_; ++i)" << endl;
    ss << "      dgemm_(\"N\", \"N\", " << rank << ", b2*a2, amax_+1, 1.0, workz+i*" << rank << "*(amax_+1), " << rank << ", transz, amax_+1, 0.0, intermediate+i*" << rank << "*b2*a2, " << rank << ");" << endl;
    ss << "    dgemm_(\"N\", \"N\", " << rank << "*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, " << rank << "*b2*a2, trans2z, cmax_+1, 0.0, final_z, " << rank << "*b2*a2);" << endl;
    ss << "" << endl;
    ss << "    for (int id = 0; id <= d; ++id) {" << endl;
    ss << "      for (int ic = 0; ic <= c; ++ic) {" << endl;
    ss << "        for (int ib = 0; ib <= b; ++ib) {" << endl;
    ss << "          for (int ia = 0; ia <= a; ++ia) {" << endl;
    ss << "            for (int r = 0; r != " << rank << "; ++r) {" << endl;
    ss << "                                                                                // v- this is a little dangerous, but perhaps the best" << endl;
    ss << "              final_xa[m<" << rank << ">(r,ia,ib,ic,id)] = 2.0*expo[0]*final_x[m<" << rank << ">(r,ia+1,ib,ic,id)] - ia*final_x[m<" << rank << ">(r,ia-1,ib,ic,id)];" << endl;
    ss << "              final_xb[m<" << rank << ">(r,ia,ib,ic,id)] = 2.0*expo[1]*final_x[m<" << rank << ">(r,ia,ib+1,ic,id)] - ib*final_x[m<" << rank << ">(r,ia,ib-1,ic,id)];" << endl;
    ss << "              final_xc[m<" << rank << ">(r,ia,ib,ic,id)] = 2.0*expo[2]*final_x[m<" << rank << ">(r,ia,ib,ic+1,id)] - ic*final_x[m<" << rank << ">(r,ia,ib,ic-1,id)];" << endl;
    ss << "" << endl;
    ss << "              final_ya[m<" << rank << ">(r,ia,ib,ic,id)] = 2.0*expo[0]*final_y[m<" << rank << ">(r,ia+1,ib,ic,id)] - ia*final_y[m<" << rank << ">(r,ia-1,ib,ic,id)];" << endl;
    ss << "              final_yb[m<" << rank << ">(r,ia,ib,ic,id)] = 2.0*expo[1]*final_y[m<" << rank << ">(r,ia,ib+1,ic,id)] - ib*final_y[m<" << rank << ">(r,ia,ib-1,ic,id)];" << endl;
    ss << "              final_yc[m<" << rank << ">(r,ia,ib,ic,id)] = 2.0*expo[2]*final_y[m<" << rank << ">(r,ia,ib,ic+1,id)] - ic*final_y[m<" << rank << ">(r,ia,ib,ic-1,id)];" << endl;
    ss << "" << endl;
    ss << "              final_za[m<" << rank << ">(r,ia,ib,ic,id)] = 2.0*expo[0]*final_z[m<" << rank << ">(r,ia+1,ib,ic,id)] - ia*final_z[m<" << rank << ">(r,ia-1,ib,ic,id)];" << endl;
    ss << "              final_zb[m<" << rank << ">(r,ia,ib,ic,id)] = 2.0*expo[1]*final_z[m<" << rank << ">(r,ia,ib+1,ic,id)] - ib*final_z[m<" << rank << ">(r,ia,ib-1,ic,id)];" << endl;
    ss << "              final_zc[m<" << rank << ">(r,ia,ib,ic,id)] = 2.0*expo[2]*final_z[m<" << rank << ">(r,ia,ib,ic+1,id)] - ic*final_z[m<" << rank << ">(r,ia,ib,ic-1,id)];" << endl;
    ss << "            }   " << endl;
    ss << "          }   " << endl;
    ss << "        }   " << endl;
    ss << "      }   " << endl;
    ss << "    }   " << endl;
    ss << "" << endl;
    ss << "    double* current_data0  = data_ + data_offset_ii;" << endl;
    ss << "    double* current_data1  = data_ + data_offset_ii + size_block_;" << endl;
    ss << "    double* current_data2  = data_ + data_offset_ii + size_block_* 2;" << endl;
    ss << "    double* current_data3  = data_ + data_offset_ii + size_block_* 3;" << endl;
    ss << "    double* current_data4  = data_ + data_offset_ii + size_block_* 4;" << endl;
    ss << "    double* current_data5  = data_ + data_offset_ii + size_block_* 5;" << endl;
    ss << "    double* current_data6  = data_ + data_offset_ii + size_block_* 6;" << endl;
    ss << "    double* current_data7  = data_ + data_offset_ii + size_block_* 7;" << endl;
    ss << "    double* current_data8  = data_ + data_offset_ii + size_block_* 8;" << endl;
    ss << "" << endl;
    ss << "    // CAUTION!" << endl;
    ss << "    // integrals in the 0(1(2(3(x2(x3(x0(x1))))))) order" << endl;
    ss << "    for (int icz = 0; icz <= c; ++icz) {" << endl;
    ss << "    for (int icy = 0; icy <= c - icz; ++icy) {" << endl;
    ss << "    const int icx = c - icz - icy;" << endl;
    ss << "" << endl;
    ss << "      for (int idz = 0; idz <= d; ++idz) {" << endl;
    ss << "      for (int idy = 0; idy <= d - idz; ++idy) {" << endl;
    ss << "      const int idx = d - idz - idy;" << endl;
    ss << "" << endl;
    ss << "        for (int iaz = 0; iaz <= a; ++iaz) {" << endl;
    ss << "        for (int iay = 0; iay <= a - iaz; ++iay) {" << endl;
    ss << "        const int iax = a - iaz - iay;" << endl;
    ss << "" << endl;
    ss << "          for (int ibz = 0; ibz <= b; ++ibz) {" << endl;
    ss << "          for (int iby = 0; iby <= b - ibz; ++iby) {" << endl;
    ss << "          const int ibx = b - ibz - iby;" << endl;
    ss << "" << endl;
    ss << "            for (int i = 0; i != " << rank << "; ++i) {" << endl;
    ss << "              *current_data0  += final_xa[m<" << rank << ">(i, iax, ibx, icx, idx)] * final_y [m<" << rank << ">(i, iay, iby, icy, idy)] * final_z [m<" << rank << ">(i, iaz, ibz, icz, idz)];" << endl;
    ss << "              *current_data1  += final_x [m<" << rank << ">(i, iax, ibx, icx, idx)] * final_ya[m<" << rank << ">(i, iay, iby, icy, idy)] * final_z [m<" << rank << ">(i, iaz, ibz, icz, idz)];" << endl;
    ss << "              *current_data2  += final_x [m<" << rank << ">(i, iax, ibx, icx, idx)] * final_y [m<" << rank << ">(i, iay, iby, icy, idy)] * final_za[m<" << rank << ">(i, iaz, ibz, icz, idz)];" << endl;
    ss << "              *current_data3  += final_xb[m<" << rank << ">(i, iax, ibx, icx, idx)] * final_y [m<" << rank << ">(i, iay, iby, icy, idy)] * final_z [m<" << rank << ">(i, iaz, ibz, icz, idz)];" << endl;
    ss << "              *current_data4  += final_x [m<" << rank << ">(i, iax, ibx, icx, idx)] * final_yb[m<" << rank << ">(i, iay, iby, icy, idy)] * final_z [m<" << rank << ">(i, iaz, ibz, icz, idz)];" << endl;
    ss << "              *current_data5  += final_x [m<" << rank << ">(i, iax, ibx, icx, idx)] * final_y [m<" << rank << ">(i, iay, iby, icy, idy)] * final_zb[m<" << rank << ">(i, iaz, ibz, icz, idz)];" << endl;
    ss << "              *current_data6  += final_xc[m<" << rank << ">(i, iax, ibx, icx, idx)] * final_y [m<" << rank << ">(i, iay, iby, icy, idy)] * final_z [m<" << rank << ">(i, iaz, ibz, icz, idz)];" << endl;
    ss << "              *current_data7  += final_x [m<" << rank << ">(i, iax, ibx, icx, idx)] * final_yc[m<" << rank << ">(i, iay, iby, icy, idy)] * final_z [m<" << rank << ">(i, iaz, ibz, icz, idz)];" << endl;
    ss << "              *current_data8  += final_x [m<" << rank << ">(i, iax, ibx, icx, idx)] * final_y [m<" << rank << ">(i, iay, iby, icy, idy)] * final_zc[m<" << rank << ">(i, iaz, ibz, icz, idz)];" << endl;
    ss << "            }   " << endl;
    ss << "            ++current_data0;" << endl;
    ss << "            ++current_data1;" << endl;
    ss << "            ++current_data2;" << endl;
    ss << "            ++current_data3;" << endl;
    ss << "            ++current_data4;" << endl;
    ss << "            ++current_data5;" << endl;
    ss << "            ++current_data6;" << endl;
    ss << "            ++current_data7;" << endl;
    ss << "            ++current_data8;" << endl;
    ss << "" << endl;
    ss << "          }}  " << endl;
    ss << "        }}  " << endl;
    ss << "      }}  " << endl;
    ss << "    }}  " << endl;
    ss << "" << endl;
    ss << "  }" << endl;
    ss << "" << endl;
    ss << "  stack_->release(b2*a2*c2*d2*" << rank << ", final_zc);" << endl;
    ss << "  stack_->release(b2*a2*c2*d2*" << rank << ", final_zb);" << endl;
    ss << "  stack_->release(b2*a2*c2*d2*" << rank << ", final_za);" << endl;
    ss << "  stack_->release(b2*a2*c2*d2*" << rank << ", final_yc);" << endl;
    ss << "  stack_->release(b2*a2*c2*d2*" << rank << ", final_yb);" << endl;
    ss << "  stack_->release(b2*a2*c2*d2*" << rank << ", final_ya);" << endl;
    ss << "  stack_->release(b2*a2*c2*d2*" << rank << ", final_xc);" << endl;
    ss << "  stack_->release(b2*a2*c2*d2*" << rank << ", final_xb);" << endl;
    ss << "  stack_->release(b2*a2*c2*d2*" << rank << ", final_xa);" << endl;
    ss << "  stack_->release(b2*a2*c2*d2*" << rank << ", final_z);" << endl;
    ss << "  stack_->release(b2*a2*c2*d2*" << rank << ", final_y);" << endl;
    ss << "  stack_->release(b2*a2*c2*d2*" << rank << ", final_x);" << endl;
    ss << "" << endl;
    ss << "  stack_->release(b2*a2*(cmax_+1)*" << rank << ", intermediate);" << endl;
    ss << "" << endl;
    ss << "  stack_->release((cmax_+1)*c2*d2, trans2z);" << endl;
    ss << "  stack_->release((cmax_+1)*c2*d2, trans2y);" << endl;
    ss << "  stack_->release((cmax_+1)*c2*d2, trans2x);" << endl;
    ss << "" << endl;
    ss << "  stack_->release((amax_+1)*a2*b2, transz);" << endl;
    ss << "  stack_->release((amax_+1)*a2*b2, transy);" << endl;
    ss << "  stack_->release((amax_+1)*a2*b2, transx);" << endl;
    ss << "" << endl;
    ss << "  stack_->release(worksize*3, workx);" << endl;
    ss << "}" << endl;
    ss << endl;
  }

  ofs << ss.str();
  return 0;
}
