//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: f77_interface.cc
// Copyright (C) 2011 Toru Shiozaki
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

#include <src/util/f77.h>
#include <bagel_config.h>

#ifndef HAVE_ZGEMM3M
void zgemm3m_(const char* transa, const char* transb, const int* m, const int* n, const int* k,
              const std::complex<double>* alpha, const std::complex<double>* a, const int* lda, const std::complex<double>* b, const int* ldb,
              const std::complex<double>* beta, std::complex<double>* c, const int* ldc) { zgemm_(transa, transb, m,n,k, alpha, a, lda, b, ldb, beta, c, ldc); }
#endif
