//
// BAGEL - Parallel electron correlation program.
// Filename: f77_interface.cc
// Copyright (C) 2011 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
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

#include <src/util/f77.h>
#include <config.h>

#ifndef HAVE_MKL_H
extern "C" {
  void zgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k,
              const std::complex<double>* alpha, const std::complex<double>* a, const int* lda, const std::complex<double>* b, const int* ldb,
              const std::complex<double>* beta, std::complex<double>* c, const int* ldc);
}
void zgemm3m_(const char* transa, const char* transb, const int* m, const int* n, const int* k,
              const std::complex<double>* alpha, const std::complex<double>* a, const int* lda, const std::complex<double>* b, const int* ldb,
              const std::complex<double>* beta, std::complex<double>* c, const int* ldc) { zgemm_(transa, transb, m,n,k, alpha, a, lda, b, ldb, beta, c, ldc); }
#endif
