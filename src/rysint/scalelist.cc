//
// Newint - Parallel electron correlation program.
// Filename: scalelist.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#include <src/rysint/scalelist.h>

ScaleList::ScaleList() {
  scalefunc[1] = &scale_data_1;
  scalefunc[2] = &scale_data_2;
  scalefunc[3] = &scale_data_3;
  scalefunc[4] = &scale_data_4;
  scalefunc[5] = &scale_data_5;
  scalefunc[6] = &scale_data_6;
  scalefunc[7] = &scale_data_7;
  scalefunc[8] = &scale_data_8;
  scalefunc[9] = &scale_data_9;
  scalefunc[10] = &scale_data_10;
  scalefunc[11] = &scale_data_11;
  scalefunc[12] = &scale_data_12;
  scalefunc[13] = &scale_data_13;
}


ScaleList::~ScaleList() {

}


void ScaleList::scale_data_1(double* out, const double* a, const double c, const double* in, const int datasize) {
  const int outerloop = datasize;
  double ca = c * a[0];
  for (int i = 0; i != outerloop; ++i, out += 1, in += 1) {
    *out = *in * ca;
  } 
}


void ScaleList::scale_data_2(double* out, const double* a, const double c, const double* in, const int datasize) {
  const int outerloop = datasize / 2;
  int offset = 0;
  double ca[2];
  ca[0] = c * a[0];
  ca[1] = c * a[1];
  for (int i = 0; i != outerloop; ++i, offset += 2) {
    out[offset    ] = in[offset    ] * ca[0];
    out[offset + 1] = in[offset + 1] * ca[1];
  } 
}


void ScaleList::scale_data_3(double* out, const double* a, const double c, const double* in, const int datasize) {
  const int outerloop = datasize / 3;
  int offset = 0;
  double ca[3];
  ca[0] = c * a[0];
  ca[1] = c * a[1];
  ca[2] = c * a[2];
  for (int i = 0; i != outerloop; ++i, offset += 3) {
    out[offset    ] = in[offset    ] * ca[0];
    out[offset + 1] = in[offset + 1] * ca[1];
    out[offset + 2] = in[offset + 2] * ca[2];
  } 
}


void ScaleList::scale_data_4(double* out, const double* a, const double c, const double* in, const int datasize) {
  const int outerloop = datasize / 4;
  int offset = 0;
  double ca[4];
  ca[0] = c * a[0];
  ca[1] = c * a[1];
  ca[2] = c * a[2];
  ca[3] = c * a[3];
  for (int i = 0; i != outerloop; ++i, offset += 4) {
    out[offset    ] = in[offset    ] * ca[0];
    out[offset + 1] = in[offset + 1] * ca[1];
    out[offset + 2] = in[offset + 2] * ca[2];
    out[offset + 3] = in[offset + 3] * ca[3];
  } 
}


void ScaleList::scale_data_5(double* out, const double* a, const double c, const double* in, const int datasize) {
  const int outerloop = datasize / 5;
  int offset = 0;
  double ca[5];
  ca[0] = c * a[0];
  ca[1] = c * a[1];
  ca[2] = c * a[2];
  ca[3] = c * a[3];
  ca[4] = c * a[4];
  for (int i = 0; i != outerloop; ++i, offset += 5) {
    out[offset    ] = in[offset    ] * ca[0];
    out[offset + 1] = in[offset + 1] * ca[1];
    out[offset + 2] = in[offset + 2] * ca[2];
    out[offset + 3] = in[offset + 3] * ca[3];
    out[offset + 4] = in[offset + 4] * ca[4];
  } 
}


void ScaleList::scale_data_6(double* out, const double* a, const double c, const double* in, const int datasize) {
  const int outerloop = datasize / 6;
  int offset = 0;
  double ca[6];
  ca[0] = c * a[0]; 
  ca[1] = c * a[1]; 
  ca[2] = c * a[2]; 
  ca[3] = c * a[3]; 
  ca[4] = c * a[4]; 
  ca[5] = c * a[5]; 
  for (int i = 0; i != outerloop; ++i, offset += 6) {
    out[offset    ] = in[offset    ] * ca[0];
    out[offset + 1] = in[offset + 1] * ca[1];
    out[offset + 2] = in[offset + 2] * ca[2];
    out[offset + 3] = in[offset + 3] * ca[3];
    out[offset + 4] = in[offset + 4] * ca[4];
    out[offset + 5] = in[offset + 5] * ca[5];
  } 
}


void ScaleList::scale_data_7(double* out, const double* a, const double c, const double* in, const int datasize) {
  const int outerloop = datasize / 7;
  int offset = 0;
  double ca[7];
  ca[0] = c * a[0];
  ca[1] = c * a[1];
  ca[2] = c * a[2];
  ca[3] = c * a[3];
  ca[4] = c * a[4];
  ca[5] = c * a[5];
  ca[6] = c * a[6];
  for (int i = 0; i != outerloop; ++i, offset += 7) {
    out[offset    ] = in[offset    ] * ca[0];
    out[offset + 1] = in[offset + 1] * ca[1];
    out[offset + 2] = in[offset + 2] * ca[2];
    out[offset + 3] = in[offset + 3] * ca[3];
    out[offset + 4] = in[offset + 4] * ca[4];
    out[offset + 5] = in[offset + 5] * ca[5];
    out[offset + 6] = in[offset + 6] * ca[6];
  } 
}


void ScaleList::scale_data_8(double* out, const double* a, const double c, const double* in, const int datasize) {
  const int outerloop = datasize / 8;
  int offset = 0;
  double ca[8];
  ca[0] = c * a[0];
  ca[1] = c * a[1];
  ca[2] = c * a[2];
  ca[3] = c * a[3];
  ca[4] = c * a[4];
  ca[5] = c * a[5];
  ca[6] = c * a[6];
  ca[7] = c * a[7];
  for (int i = 0; i != outerloop; ++i, offset += 8) {
    out[offset    ] = in[offset    ] * ca[0];
    out[offset + 1] = in[offset + 1] * ca[1];
    out[offset + 2] = in[offset + 2] * ca[2];
    out[offset + 3] = in[offset + 3] * ca[3];
    out[offset + 4] = in[offset + 4] * ca[4];
    out[offset + 5] = in[offset + 5] * ca[5];
    out[offset + 6] = in[offset + 6] * ca[6];
    out[offset + 7] = in[offset + 7] * ca[7];
  } 
}


void ScaleList::scale_data_9(double* out, const double* a, const double c, const double* in, const int datasize) {
  const int outerloop = datasize / 9;
  int offset = 0;
  double ca[9];
  ca[0] = c * a[0];
  ca[1] = c * a[1];
  ca[2] = c * a[2];
  ca[3] = c * a[3];
  ca[4] = c * a[4];
  ca[5] = c * a[5];
  ca[6] = c * a[6];
  ca[7] = c * a[7];
  ca[8] = c * a[8];
  for (int i = 0; i != outerloop; ++i, offset += 9) {
    out[offset    ] = in[offset    ] * ca[0];
    out[offset + 1] = in[offset + 1] * ca[1];
    out[offset + 2] = in[offset + 2] * ca[2];
    out[offset + 3] = in[offset + 3] * ca[3];
    out[offset + 4] = in[offset + 4] * ca[4];
    out[offset + 5] = in[offset + 5] * ca[5];
    out[offset + 6] = in[offset + 6] * ca[6];
    out[offset + 7] = in[offset + 7] * ca[7];
    out[offset + 8] = in[offset + 8] * ca[8];
  } 
}


void ScaleList::scale_data_10(double* out, const double* a, const double c, const double* in, const int datasize) {
  const int outerloop = datasize / 10;
  int offset = 0;
  double ca[10];
  ca[0] = c * a[0];
  ca[1] = c * a[1];
  ca[2] = c * a[2];
  ca[3] = c * a[3];
  ca[4] = c * a[4];
  ca[5] = c * a[5];
  ca[6] = c * a[6];
  ca[7] = c * a[7];
  ca[8] = c * a[8];
  ca[9] = c * a[9];
  for (int i = 0; i != outerloop; ++i, offset += 10) {
    out[offset    ] = in[offset    ] * ca[0];
    out[offset + 1] = in[offset + 1] * ca[1];
    out[offset + 2] = in[offset + 2] * ca[2];
    out[offset + 3] = in[offset + 3] * ca[3];
    out[offset + 4] = in[offset + 4] * ca[4];
    out[offset + 5] = in[offset + 5] * ca[5];
    out[offset + 6] = in[offset + 6] * ca[6];
    out[offset + 7] = in[offset + 7] * ca[7];
    out[offset + 8] = in[offset + 8] * ca[8];
    out[offset + 9] = in[offset + 9] * ca[9];
  } 
}


void ScaleList::scale_data_11(double* out, const double* a, const double c, const double* in, const int datasize) {
  const int outerloop = datasize / 11;
  int offset = 0;
  double ca[11];
  ca[0] = c * a[0];
  ca[1] = c * a[1];
  ca[2] = c * a[2];
  ca[3] = c * a[3];
  ca[4] = c * a[4];
  ca[5] = c * a[5];
  ca[6] = c * a[6];
  ca[7] = c * a[7];
  ca[8] = c * a[8];
  ca[9] = c * a[9];
  ca[10] = c * a[10];
  for (int i = 0; i != outerloop; ++i, offset += 11) {
    out[offset    ] = in[offset    ] * ca[0];
    out[offset + 1] = in[offset + 1] * ca[1];
    out[offset + 2] = in[offset + 2] * ca[2];
    out[offset + 3] = in[offset + 3] * ca[3];
    out[offset + 4] = in[offset + 4] * ca[4];
    out[offset + 5] = in[offset + 5] * ca[5];
    out[offset + 6] = in[offset + 6] * ca[6];
    out[offset + 7] = in[offset + 7] * ca[7];
    out[offset + 8] = in[offset + 8] * ca[8];
    out[offset + 9] = in[offset + 9] * ca[9];
    out[offset + 10] = in[offset + 10] * ca[10];
  } 
}


void ScaleList::scale_data_12(double* out, const double* a, const double c, const double* in, const int datasize) {
  const int outerloop = datasize / 12;
  int offset = 0;
  double ca[12];
  ca[0] = c * a[0];
  ca[1] = c * a[1];
  ca[2] = c * a[2];
  ca[3] = c * a[3];
  ca[4] = c * a[4];
  ca[5] = c * a[5];
  ca[6] = c * a[6];
  ca[7] = c * a[7];
  ca[8] = c * a[8];
  ca[9] = c * a[9];
  ca[10] = c * a[10];
  ca[11] = c * a[11];
  for (int i = 0; i != outerloop; ++i, offset += 12) {
    out[offset    ] = in[offset    ] * ca[0];
    out[offset + 1] = in[offset + 1] * ca[1];
    out[offset + 2] = in[offset + 2] * ca[2];
    out[offset + 3] = in[offset + 3] * ca[3];
    out[offset + 4] = in[offset + 4] * ca[4];
    out[offset + 5] = in[offset + 5] * ca[5];
    out[offset + 6] = in[offset + 6] * ca[6];
    out[offset + 7] = in[offset + 7] * ca[7];
    out[offset + 8] = in[offset + 8] * ca[8];
    out[offset + 9] = in[offset + 9] * ca[9];
    out[offset + 10] = in[offset + 10] * ca[10];
    out[offset + 11] = in[offset + 11] * ca[11];
  } 
}


void ScaleList::scale_data_13(double* out, const double* a, const double c, const double* in, const int datasize) {
  const int outerloop = datasize / 13;
  int offset = 0;
  double ca[13];
  ca[0] = c * a[0];
  ca[1] = c * a[1];
  ca[2] = c * a[2];
  ca[3] = c * a[3];
  ca[4] = c * a[4];
  ca[5] = c * a[5];
  ca[6] = c * a[6];
  ca[7] = c * a[7];
  ca[8] = c * a[8];
  ca[9] = c * a[9];
  ca[10] = c * a[10];
  ca[11] = c * a[11];
  ca[12] = c * a[12];
  for (int i = 0; i != outerloop; ++i, offset += 13) {
    out[offset    ] = in[offset    ] * ca[0];
    out[offset + 1] = in[offset + 1] * ca[1];
    out[offset + 2] = in[offset + 2] * ca[2];
    out[offset + 3] = in[offset + 3] * ca[3];
    out[offset + 4] = in[offset + 4] * ca[4];
    out[offset + 5] = in[offset + 5] * ca[5];
    out[offset + 6] = in[offset + 6] * ca[6];
    out[offset + 7] = in[offset + 7] * ca[7];
    out[offset + 8] = in[offset + 8] * ca[8];
    out[offset + 9] = in[offset + 9] * ca[9];
    out[offset + 10] = in[offset + 10] * ca[10];
    out[offset + 11] = in[offset + 11] * ca[11];
    out[offset + 12] = in[offset + 12] * ca[12];
  } 
}


