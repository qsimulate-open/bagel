//
// Newint - Parallel electron correlation program.
// Filename: prim_op.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#ifndef __SRC_SMITH_PRIM_OP_H
#define __SRC_SMITH_PRIM_OP_H

#include <memory>
#include <vector>
#include <src/util/f77.h>

namespace SMITH {

static void sort_indices2(const std::unique_ptr<double[]>& unsorted, std::unique_ptr<double[]>& sorted,
                          const int b, const int a, // according to unsorted
                          const int i, const int j, const double afac = 0.0, const double factor = 1.0) {
  if (i==0) {
    if (afac == 1.0) {
      daxpy_(a*b, factor, unsorted, 1, sorted, 1);
    } else if (afac == 0.0) {
      dcopy_(a*b, unsorted, 1, sorted, 1);
    } else {
      dscal_(a*b, afac, sorted, 1);
      daxpy_(a*b, factor, unsorted, 1, sorted, 1);
    }
  } else { 
    int id[2];
    int jd[2] = {b, a};
    long iall=0;
    for(int j0=0;j0<(int)a;++j0){
      id[1]=j0;
      for(int j1=0;j1<(int)b;++j1,++iall){
        id[0]=j1;
        long ib=id[i]+jd[i]*id[j];
        sorted[ib]=afac*sorted[ib]+unsorted[iall]*factor;
      }
    }
  }
} 


// CAUTION :: I have changed the convention from that in mpqc.
static void sort_indices4(const std::unique_ptr<double[]>& unsorted, std::unique_ptr<double[]>& sorted,
                          const int d,const int c,const int b,const int a, // according to unsorted
                          const int i,const int j,const int k,const int l,
                          const double afac = 0.0, const double factor = 1.0) {
  int id[4];
  int jd[4] = {d, c, b, a};

  long iall=0;
  for(int j0=0;j0<a;++j0){
    id[3]=j0;
    for(int j1=0;j1<b;++j1){
      id[2]=j1;
      for(int j2=0;j2<c;++j2){
        id[1]=j2;
        for (int j3=0;j3<d;++j3,++iall){
          id[0]=j3;
          long ib=id[i]+jd[i]*(id[j]+jd[j]*(id[k]+jd[k]*id[l]));
          sorted[ib]=afac*sorted[ib]+unsorted[iall]*factor;
        }
      }
    } 
  } 
};


template <typename T>
static std::vector<T> vec(T i0) {
  std::vector<T> out(1,i0);
  return out;
};
template <typename T>
static std::vector<T> vec(T i0, T i1) {
  std::vector<T> out(1,i0); out.push_back(i1);
  return out;
};
template <typename T>
static std::vector<T> vec(T i0, T i1, T i2) {
  std::vector<T> out(1,i0); out.push_back(i1); out.push_back(i2);
  return out;
};
template <typename T>
static std::vector<T> vec(T i0, T i1, T i2, T i3) {
  std::vector<T> out(1,i0); out.push_back(i1); out.push_back(i2); out.push_back(i3);
  return out;
};
template <typename T>
static std::vector<T> vec(T i0, T i1, T i2, T i3, T i4) {
  std::vector<T> out(1,i0); out.push_back(i1); out.push_back(i2); out.push_back(i3); out.push_back(i4);
  return out;
};
template <typename T>
static std::vector<T> vec(T i0, T i1, T i2, T i3, T i4, T i5) {
  std::vector<T> out(1,i0); out.push_back(i1); out.push_back(i2); out.push_back(i3); out.push_back(i4); out.push_back(i5);
  return out;
};
template <typename T>
static std::vector<T> vec(T i0, T i1, T i2, T i3, T i4, T i5, T i6) {
  std::vector<T> out(1,i0); out.push_back(i1); out.push_back(i2); out.push_back(i3); out.push_back(i4); out.push_back(i5);
                            out.push_back(i6);
  return out;
};
template <typename T>
static std::vector<T> vec(T i0, T i1, T i2, T i3, T i4, T i5, T i6, T i7) {
  std::vector<T> out(1,i0); out.push_back(i1); out.push_back(i2); out.push_back(i3); out.push_back(i4); out.push_back(i5);
                            out.push_back(i6); out.push_back(i7);
  return out;
};


}

#endif
