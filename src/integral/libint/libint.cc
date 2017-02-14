//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: libint.cc
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


#ifdef LIBINT_INTERFACE

#include <stddef.h>
#include <iostream>
#include <src/util/f77.h>
#include <src/util/math/algo.h>
#include <src/integral/carsphlist.h>
#include <src/integral/libint/libint.h>
#include <boys.h>

#if LIBINT2_CGSHELL_ORDERING != LIBINT2_CGSHELL_ORDERING_BAGEL
# error "incompatible Libint2 library, must support BAGEL cartesian Gaussian shell ordering"
#endif

using namespace std;
using namespace bagel;

const static CarSphList carsphlist;
const static libint2::FmEval_Chebyshev7<double> fmeval(18);

Libint::Libint(const array<shared_ptr<const Shell>,4>& shells, const double dum,  shared_ptr<StackMem> stack)
 : RysInt(shells, stack) {

  tenno_ = 0;
  size_allocated_ = 0LU;
  size_alloc_ = 0LU;

  array<int,4> order {{ 0,1,2,3 }};


  // first count the number of dummys
  int center = 4;
  for (auto& i : basisinfo_) if (i->dummy()) --center;

  if (center == 3 && !(basisinfo_[0]->dummy() || basisinfo_[1]->dummy()))
    throw logic_error("dummy shell in an illegal position: 3 index gradient");
  if (center == 2 && (!(basisinfo_[0]->dummy() || basisinfo_[1]->dummy()) || !(basisinfo_[2]->dummy() || basisinfo_[3]->dummy())))
    throw logic_error("dummy shell in an illegal position: 2 index gradient");
  if (center < 2) throw logic_error("there are only one or less non-dummy basis in GLibint::GLibint");


  if (basisinfo_[0]->angular_number() < basisinfo_[1]->angular_number() || basisinfo_[0]->dummy()) {
    swap(basisinfo_[0], basisinfo_[1]);
    swap(order[0], order[1]);
    swap01_ = true;
  } else {
    swap01_ = false;
  }
  // swap 23 indices when needed
  if (basisinfo_[2]->angular_number() < basisinfo_[3]->angular_number() || basisinfo_[2]->dummy()) {
    swap(basisinfo_[2], basisinfo_[3]);
    swap(order[2], order[3]);
    swap23_ = true;
  } else {
    swap23_ = false;
  }

  swap0123_ = false;
  if ((basisinfo_[0]->angular_number()+basisinfo_[1]->angular_number() > basisinfo_[2]->angular_number()+basisinfo_[3]->angular_number())
      && center != 3) {
      swap0123_ = true;
      tie(basisinfo_[0], basisinfo_[1], basisinfo_[2], basisinfo_[3], swap01_, swap23_)
        = make_tuple(basisinfo_[2], basisinfo_[3], basisinfo_[0], basisinfo_[1], swap23_, swap01_);
    swap(order[0], order[2]);
    swap(order[1], order[3]);
  }

  array<int,4> invmap;
  for (int i = 0; i != 4; ++i) {
    invmap[order[i]] = i;
  }

  // get contracted basis functions from each basisinfo_
  array<vector<pair<vector<double>, vector<double>>>,4> ce;
  for (int s = 0; s != 4; ++s) {
    auto range = basisinfo_[s]->contraction_ranges().begin();
    vector<pair<vector<double>, vector<double>>> seg;
    for (auto i = basisinfo_[s]->contractions().begin(); i != basisinfo_[s]->contractions().end(); ++i, ++range) {
      const int start = range->first;
      const int end   = range->second;
      vector<double> c(i->begin()+start, i->begin()+end);
      vector<double> e(basisinfo_[s]->exponents().begin()+start, basisinfo_[s]->exponents().begin()+end);
      seg.push_back(make_pair(c,e));
    }
    ce[s] = seg;
  }

  array<double,3> A = basisinfo_[0]->position();
  array<double,3> B = basisinfo_[1]->position();
  array<double,3> C = basisinfo_[2]->position();
  array<double,3> D = basisinfo_[3]->position();

  array<int,4> am = {{ basisinfo_[0]->angular_number(), basisinfo_[1]->angular_number(), basisinfo_[2]->angular_number(), basisinfo_[3]->angular_number() }};
  const unsigned int amtot = am[0] + am[1] + am[2] + am[3];

  array<unsigned int,4> sam, cam;
  for (int i = 0; i != 4; ++i) sam[i] = (am[i]+1)*(am[i]+2)/2;
  for (int i = 0; i != 4; ++i) cam[i] = 2*am[i]+1;
  if (!spherical1_) { cam[0] = sam[0]; cam[1] = sam[1]; }
  if (!spherical2_) { cam[2] = sam[2]; cam[3] = sam[3]; }

  array<size_t, 4> dim = {{cam[0]*ce[0].size(), cam[1]*ce[1].size(), cam[2]*ce[2].size(), cam[3]*ce[3].size()}};
  array<size_t, 4> base;

  size_alloc_ = dim[0]*dim[1]*dim[2]*dim[3];
  size_final_ = size_alloc_;
  data_ = stack_->get(size_alloc_);

  stack_save_ = data_;

  const double AB_x = A[0] - B[0];
  const double AB_y = A[1] - B[1];
  const double AB_z = A[2] - B[2];
  const double AB2 = AB_x * AB_x + AB_y * AB_y + AB_z * AB_z;
  const double CD_x = C[0] - D[0];
  const double CD_y = C[1] - D[1];
  const double CD_z = C[2] - D[2];
  const double CD2 = CD_x * CD_x + CD_y * CD_y + CD_z * CD_z;

  double* const F = stack_->get(LIBINT_MAX_AM*4 + 6);

  base[0] = 0;
  for (auto i0 = ce[0].begin(); i0 != ce[0].end(); ++i0, base[0] += cam[0]) {
    base[1] = 0;
    for (auto i1 = ce[1].begin(); i1 != ce[1].end(); ++i1, base[1] += cam[1]) {
      base[2] = 0;
      for (auto i2 = ce[2].begin(); i2 != ce[2].end(); ++i2, base[2] += cam[2]) {
        base[3] = 0;
        for (auto i3 = ce[3].begin(); i3 != ce[3].end(); ++i3, base[3] += cam[3]) {

          size_t p0123 = 0;
          for (auto j0 = i0->first.begin(), k0 = i0->second.begin(); j0 != i0->first.end(); ++j0, ++k0) {
            for (auto j1 = i1->first.begin(), k1 = i1->second.begin(); j1 != i1->first.end(); ++j1, ++k1) {
              const double alpha0 = *k0;
              const double alpha1 = *k1;

              const double gammap = alpha0 + alpha1;
              const double oogammap = 1.0 / gammap;
              const double rhop = alpha0 * alpha1 * oogammap;
              const double Px = (alpha0 * A[0] + alpha1 * B[0]) * oogammap;
              const double Py = (alpha0 * A[1] + alpha1 * B[1]) * oogammap;
              const double Pz = (alpha0 * A[2] + alpha1 * B[2]) * oogammap;
              const double PAx = Px - A[0];
              const double PAy = Py - A[1];
              const double PAz = Pz - A[2];
              const double PBx = Px - B[0];
              const double PBy = Py - B[1];
              const double PBz = Pz - B[2];

              const double K1 = exp(- rhop * AB2);

              for (auto j2 = i2->first.begin(), k2 = i2->second.begin(); j2 != i2->first.end(); ++j2, ++k2) {
                for (auto j3 = i3->first.begin(), k3 = i3->second.begin(); j3 != i3->first.end(); ++j3, ++k3, ++p0123) {

                  Libint_t* erieval = stack_->libint_t_ptr(p0123);
                  erieval->veclen = 1;

                  const double alpha2 = *k2;
                  const double alpha3 = *k3;

                  const double c0 = *j0;
                  const double c1 = *j1;
                  const double c2 = *j2;
                  const double c3 = *j3;

#if LIBINT2_DEFINED(eri,PA_x)
                  erieval->PA_x[0] = PAx;
#endif
#if LIBINT2_DEFINED(eri,PA_y)
                  erieval->PA_y[0] = PAy;
#endif
#if LIBINT2_DEFINED(eri,PA_z)
                  erieval->PA_z[0] = PAz;
#endif
#if LIBINT2_DEFINED(eri,PB_x)
                  erieval->PB_x[0] = PBx;
#endif
#if LIBINT2_DEFINED(eri,PB_y)
                  erieval->PB_y[0] = PBy;
#endif
#if LIBINT2_DEFINED(eri,PB_z)
                  erieval->PB_z[0] = PBz;
#endif

#if LIBINT2_DEFINED(eri,AB_x)
                  erieval->AB_x[0] = AB_x;
#endif
#if LIBINT2_DEFINED(eri,AB_y)
                  erieval->AB_y[0] = AB_y;
#endif
#if LIBINT2_DEFINED(eri,AB_z)
                  erieval->AB_z[0] = AB_z;
#endif
#if LIBINT2_DEFINED(eri,oo2z)
                  erieval->oo2z[0] = 0.5*oogammap;
#endif
                  const double gammaq = alpha2 + alpha3;
                  const double oogammaq = 1.0 / gammaq;
                  const double rhoq = alpha2 * alpha3 * oogammaq;
                  const double gammapq = gammap * gammaq / (gammap + gammaq);
                  const double gammap_o_gammapgammaq = gammapq * oogammaq;
                  const double gammaq_o_gammapgammaq = gammapq * oogammap;
                  const double Qx = (alpha2 * C[0] + alpha3 * D[0]) * oogammaq;
                  const double Qy = (alpha2 * C[1] + alpha3 * D[1]) * oogammaq;
                  const double Qz = (alpha2 * C[2] + alpha3 * D[2]) * oogammaq;
                  const double QCx = Qx - C[0];
                  const double QCy = Qy - C[1];
                  const double QCz = Qz - C[2];
                  const double QDx = Qx - D[0];
                  const double QDy = Qy - D[1];
                  const double QDz = Qz - D[2];

#if LIBINT2_DEFINED(eri,QC_x)
                  erieval->QC_x[0] = QCx;
#endif
#if LIBINT2_DEFINED(eri,QC_y)
                  erieval->QC_y[0] = QCy;
#endif
#if LIBINT2_DEFINED(eri,QC_z)
                  erieval->QC_z[0] = QCz;
#endif
#if LIBINT2_DEFINED(eri,QD_x)
                  erieval->QD_x[0] = QDx;
#endif
#if LIBINT2_DEFINED(eri,QD_y)
                  erieval->QD_y[0] = QDy;
#endif
#if LIBINT2_DEFINED(eri,QD_z)
                  erieval->QD_z[0] = QDz;
#endif

#if LIBINT2_DEFINED(eri,CD_x)
                  erieval->CD_x[0] = CD_x;
#endif
#if LIBINT2_DEFINED(eri,CD_y)
                  erieval->CD_y[0] = CD_y;
#endif
#if LIBINT2_DEFINED(eri,CD_z)
                  erieval->CD_z[0] = CD_z;
#endif
#if LIBINT2_DEFINED(eri,oo2e)
                  erieval->oo2e[0] = 0.5*oogammaq;
#endif

    // Prefactors for interelectron transfer relation
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_0_x)
                  erieval->TwoPRepITR_pfac0_0_0_x[0] = - (alpha1*AB_x + alpha3*CD_x)*oogammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_0_y)
                  erieval->TwoPRepITR_pfac0_0_0_y[0] = - (alpha1*AB_y + alpha3*CD_y)*oogammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_0_z)
                  erieval->TwoPRepITR_pfac0_0_0_z[0] = - (alpha1*AB_z + alpha3*CD_z)*oogammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_0_x)
                  erieval->TwoPRepITR_pfac0_1_0_x[0] = - (alpha1*AB_x + alpha3*CD_x)*oogammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_0_y)
                  erieval->TwoPRepITR_pfac0_1_0_y[0] = - (alpha1*AB_y + alpha3*CD_y)*oogammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_0_z)
                  erieval->TwoPRepITR_pfac0_1_0_z[0] = - (alpha1*AB_z + alpha3*CD_z)*oogammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_1_x)
                  erieval->TwoPRepITR_pfac0_0_1_x[0] = (alpha0*AB_x + alpha2*CD_x)*oogammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_1_y)
                  erieval->TwoPRepITR_pfac0_0_1_y[0] = (alpha0*AB_y + alpha2*CD_y)*oogammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_1_z)
                  erieval->TwoPRepITR_pfac0_0_1_z[0] = (alpha0*AB_z + alpha2*CD_z)*oogammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_1_x)
                  erieval->TwoPRepITR_pfac0_1_1_x[0] = (alpha0*AB_x + alpha2*CD_x)*oogammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_1_y)
                  erieval->TwoPRepITR_pfac0_1_1_y[0] = (alpha0*AB_y + alpha2*CD_y)*oogammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_1_z)
                  erieval->TwoPRepITR_pfac0_1_1_z[0] = (alpha0*AB_z + alpha2*CD_z)*oogammaq;
#endif
#if LIBINT2_DEFINED(eri,eoz)
                  erieval->eoz[0] = gammaq*oogammap;
#endif
#if LIBINT2_DEFINED(eri,zoe)
                  erieval->zoe[0] = gammap*oogammaq;
#endif

                  const double PQx = Px - Qx;
                  const double PQy = Py - Qy;
                  const double PQz = Pz - Qz;
                  const double PQ2 = PQx * PQx + PQy * PQy + PQz * PQz;
                  const double Wx = (gammap_o_gammapgammaq * Px + gammaq_o_gammapgammaq * Qx);
                  const double Wy = (gammap_o_gammapgammaq * Py + gammaq_o_gammapgammaq * Qy);
                  const double Wz = (gammap_o_gammapgammaq * Pz + gammaq_o_gammapgammaq * Qz);

#if LIBINT2_DEFINED(eri,WP_x)
                  erieval->WP_x[0] = Wx - Px;
#endif
#if LIBINT2_DEFINED(eri,WP_y)
                  erieval->WP_y[0] = Wy - Py;
#endif
#if LIBINT2_DEFINED(eri,WP_z)
                  erieval->WP_z[0] = Wz - Pz;
#endif
#if LIBINT2_DEFINED(eri,WQ_x)
                  erieval->WQ_x[0] = Wx - Qx;
#endif
#if LIBINT2_DEFINED(eri,WQ_y)
                  erieval->WQ_y[0] = Wy - Qy;
#endif
#if LIBINT2_DEFINED(eri,WQ_z)
                  erieval->WQ_z[0] = Wz - Qz;
#endif
#if LIBINT2_DEFINED(eri,oo2ze)
                  erieval->oo2ze[0] = 0.5/(gammap+gammaq);
#endif
#if LIBINT2_DEFINED(eri,roz)
                  erieval->roz[0] = gammapq*oogammap;
#endif
#if LIBINT2_DEFINED(eri,roe)
                  erieval->roe[0] = gammapq*oogammaq;
#endif

                  const double K2 = exp(- rhoq * CD2);
                  double pfac = 2 * pow(M_PI, 2.5) * K1 * K2 / (gammap * gammaq * sqrt(gammap + gammaq));
                  pfac *= c0 * c1 * c2 * c3;

                  fmeval.eval(F,PQ2*gammapq,amtot);

#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(0))
                  erieval->LIBINT_T_SS_EREP_SS(0)[0] = pfac*F[0];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(1))
                  erieval->LIBINT_T_SS_EREP_SS(1)[0] = pfac*F[1];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(2))
                  erieval->LIBINT_T_SS_EREP_SS(2)[0] = pfac*F[2];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(3))
                  erieval->LIBINT_T_SS_EREP_SS(3)[0] = pfac*F[3];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(4))
                  erieval->LIBINT_T_SS_EREP_SS(4)[0] = pfac*F[4];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(5))
                  erieval->LIBINT_T_SS_EREP_SS(5)[0] = pfac*F[5];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(6))
                  erieval->LIBINT_T_SS_EREP_SS(6)[0] = pfac*F[6];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(7))
                  erieval->LIBINT_T_SS_EREP_SS(7)[0] = pfac*F[7];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(8))
                  erieval->LIBINT_T_SS_EREP_SS(8)[0] = pfac*F[8];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(9))
                  erieval->LIBINT_T_SS_EREP_SS(9)[0] = pfac*F[9];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(10))
                  erieval->LIBINT_T_SS_EREP_SS(10)[0] = pfac*F[10];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(11))
                  erieval->LIBINT_T_SS_EREP_SS(11)[0] = pfac*F[11];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(12))
                  erieval->LIBINT_T_SS_EREP_SS(12)[0] = pfac*F[12];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(13))
                  erieval->LIBINT_T_SS_EREP_SS(13)[0] = pfac*F[13];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(14))
                  erieval->LIBINT_T_SS_EREP_SS(14)[0] = pfac*F[14];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(15))
                  erieval->LIBINT_T_SS_EREP_SS(15)[0] = pfac*F[15];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(16))
                  erieval->LIBINT_T_SS_EREP_SS(16)[0] = pfac*F[16];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(17))
                  erieval->LIBINT_T_SS_EREP_SS(17)[0] = pfac*F[17];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(18))
                  erieval->LIBINT_T_SS_EREP_SS(18)[0] = pfac*F[18];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(19))
                  erieval->LIBINT_T_SS_EREP_SS(19)[0] = pfac*F[19];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(20))
                  erieval->LIBINT_T_SS_EREP_SS(20)[0] = pfac*F[20];
#endif
                }
              }
            }
          }

          stack_->libint_t_ptr(0)->contrdepth = p0123;
          if (center == 4) {
            LIBINT2_PREFIXED_NAME(libint2_build_eri)[am[0]][am[1]][am[2]][am[3]](stack_->libint_t_ptr(0));
          } else if (center == 3) {
            LIBINT2_PREFIXED_NAME(libint2_build_3eri)[am[0]][am[2]][am[3]](stack_->libint_t_ptr(0));
          } else {
            LIBINT2_PREFIXED_NAME(libint2_build_2eri)[am[0]][am[2]](stack_->libint_t_ptr(0));
          }
          double* ints = stack_->libint_t_ptr(0)->targets[0];

          if (spherical1_ || spherical2_) {
            const size_t batchsize = sam[0]*sam[1]*cam[2]*cam[3];
            double* area = stack_->get(batchsize);
            const int carsphindex = am[2] * ANG_HRR_END + am[3];
            const int carsphindex2 = am[0] * ANG_HRR_END + am[1];
            const int m = cam[2]*cam[3];
            const int n = sam[0]*sam[1];
            const int nn = cam[0]*cam[1];
            if (spherical2_) {
              carsphlist.carsphfunc_call(carsphindex, n, ints, area);
            } else {
              copy_n(ints, m*n, area);
            }
            blas::transpose(area, m, n, ints);
            if (spherical1_) {
              carsphlist.carsphfunc_call(carsphindex2, m, ints, area);
              copy_n(area, nn*m, ints);
            }
            stack_->release(batchsize, area);
          }

          // ijkl runs over xyz components
          int ijkl = 0;
          array<int,4> id;
          if (!(spherical1_ || spherical2_)) {
            for (int i = 0; i != cam[0]; ++i) {
              id[0] = i;
              for (int j = 0; j != cam[1]; ++j) {
                id[1] = j;
                for (int k = 0; k != cam[2]; ++k) {
                  id[2] = k;
                  for (int l = 0; l != cam[3]; ++l, ++ijkl) {
                    id[3] = l;
                    const size_t ele =
                                 id[invmap[0]]+base[invmap[0]] + dim[invmap[0]]* (id[invmap[1]]+base[invmap[1]] + dim[invmap[1]]*
                                (id[invmap[2]]+base[invmap[2]] + dim[invmap[2]]* (id[invmap[3]]+base[invmap[3]])));
                    data_[ele] = ints[ijkl];
                  }
                }
              }
            }
          } else {
            for (int i = 0; i != cam[2]; ++i) {
              id[2] = i;
              for (int j = 0; j != cam[3]; ++j) {
                id[3] = j;
                for (int k = 0; k != cam[0]; ++k) {
                  id[0] = k;
                  for (int l = 0; l != cam[1]; ++l, ++ijkl) {
                    id[1] = l;
                    const size_t ele =
                                 id[invmap[0]]+base[invmap[0]] + dim[invmap[0]]* (id[invmap[1]]+base[invmap[1]] + dim[invmap[1]]*
                                (id[invmap[2]]+base[invmap[2]] + dim[invmap[2]]* (id[invmap[3]]+base[invmap[3]])));
                    data_[ele] = ints[ijkl];
                  }
                }
              }
            }
          }

        }
      }
    }
  }

  stack_->release(LIBINT_MAX_AM*4 + 6, F);
}


#endif
