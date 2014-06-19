//
// BAGEL - Parallel electron correlation program.
// Filename: matop.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
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

#ifndef __SRC_MATH_MATOP_H
#define __SRC_MATH_MATOP_H

#include <src/math/matrix.h>
#include <src/math/zmatrix.h>
#include <src/math/matview.h>

namespace bagel {

// operator+= and -=
inline Matrix&  operator+=(Matrix& a,  const Matrix& b)  { assert(a.size() == b.size()); blas::ax_plus_y_n( 1.0, b.data(), a.size(), a.data()); return a; }
inline Matrix&  operator+=(Matrix& a,  const MatView& b) { assert(a.size() == b.size()); blas::ax_plus_y_n( 1.0, b.data(), a.size(), a.data()); return a; }
inline MatView& operator+=(MatView& a, const Matrix& b)  { assert(a.size() == b.size()); blas::ax_plus_y_n( 1.0, b.data(), a.size(), a.data()); return a; }
inline MatView& operator+=(MatView& a, const MatView& b) { assert(a.size() == b.size()); blas::ax_plus_y_n( 1.0, b.data(), a.size(), a.data()); return a; }
inline Matrix&  operator-=(Matrix& a,  const Matrix& b)  { assert(a.size() == b.size()); blas::ax_plus_y_n(-1.0, b.data(), a.size(), a.data()); return a; }
inline Matrix&  operator-=(Matrix& a,  const MatView& b) { assert(a.size() == b.size()); blas::ax_plus_y_n(-1.0, b.data(), a.size(), a.data()); return a; }
inline MatView& operator-=(MatView& a, const Matrix& b)  { assert(a.size() == b.size()); blas::ax_plus_y_n(-1.0, b.data(), a.size(), a.data()); return a; }
inline MatView& operator-=(MatView& a, const MatView& b) { assert(a.size() == b.size()); blas::ax_plus_y_n(-1.0, b.data(), a.size(), a.data()); return a; }

inline ZMatrix&  operator+=(ZMatrix& a,  const ZMatrix& b)  { assert(a.size() == b.size()); blas::ax_plus_y_n( 1.0, b.data(), a.size(), a.data()); return a; }
inline ZMatrix&  operator+=(ZMatrix& a,  const ZMatView& b) { assert(a.size() == b.size()); blas::ax_plus_y_n( 1.0, b.data(), a.size(), a.data()); return a; }
inline ZMatView& operator+=(ZMatView& a, const ZMatrix& b)  { assert(a.size() == b.size()); blas::ax_plus_y_n( 1.0, b.data(), a.size(), a.data()); return a; }
inline ZMatView& operator+=(ZMatView& a, const ZMatView& b) { assert(a.size() == b.size()); blas::ax_plus_y_n( 1.0, b.data(), a.size(), a.data()); return a; }
inline ZMatrix&  operator-=(ZMatrix& a,  const ZMatrix& b)  { assert(a.size() == b.size()); blas::ax_plus_y_n(-1.0, b.data(), a.size(), a.data()); return a; }
inline ZMatrix&  operator-=(ZMatrix& a,  const ZMatView& b) { assert(a.size() == b.size()); blas::ax_plus_y_n(-1.0, b.data(), a.size(), a.data()); return a; }
inline ZMatView& operator-=(ZMatView& a, const ZMatrix& b)  { assert(a.size() == b.size()); blas::ax_plus_y_n(-1.0, b.data(), a.size(), a.data()); return a; }
inline ZMatView& operator-=(ZMatView& a, const ZMatView& b) { assert(a.size() == b.size()); blas::ax_plus_y_n(-1.0, b.data(), a.size(), a.data()); return a; }

// operator+ and -
inline Matrix operator+(const Matrix& a,  const Matrix& b)  { Matrix out(a); out += b; return out; }
inline Matrix operator+(const Matrix& a,  const MatView& b) { Matrix out(a); out += b; return out; }
inline Matrix operator+(const MatView& a, const Matrix& b)  { Matrix out(a); out += b; return out; }
inline Matrix operator+(const MatView& a, const MatView& b) { Matrix out(a); out += b; return out; }
inline Matrix operator-(const Matrix& a,  const Matrix& b)  { Matrix out(a); out -= b; return out; }
inline Matrix operator-(const Matrix& a,  const MatView& b) { Matrix out(a); out -= b; return out; }
inline Matrix operator-(const MatView& a, const Matrix& b)  { Matrix out(a); out -= b; return out; }
inline Matrix operator-(const MatView& a, const MatView& b) { Matrix out(a); out -= b; return out; }

inline ZMatrix operator+(const ZMatrix& a,  const ZMatrix& b)  { ZMatrix out(a); out += b; return out; }
inline ZMatrix operator+(const ZMatrix& a,  const ZMatView& b) { ZMatrix out(a); out += b; return out; }
inline ZMatrix operator+(const ZMatView& a, const ZMatrix& b)  { ZMatrix out(a); out += b; return out; }
inline ZMatrix operator+(const ZMatView& a, const ZMatView& b) { ZMatrix out(a); out += b; return out; }
inline ZMatrix operator-(const ZMatrix& a,  const ZMatrix& b)  { ZMatrix out(a); out -= b; return out; }
inline ZMatrix operator-(const ZMatrix& a,  const ZMatView& b) { ZMatrix out(a); out -= b; return out; }
inline ZMatrix operator-(const ZMatView& a, const ZMatrix& b)  { ZMatrix out(a); out -= b; return out; }
inline ZMatrix operator-(const ZMatView& a, const ZMatView& b) { ZMatrix out(a); out -= b; return out; }
}

#endif
