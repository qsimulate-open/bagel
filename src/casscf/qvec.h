//
// Newint - Parallel electron correlation program.
// Filename: qvec.h
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


#ifndef __NEWINT_SRC_CASSCF_QVEC
#define __NEWINT_SRC_CASSCF_QVEC

#include <memory>
#include <src/df/df.h>
#include <src/scf/coeff.h>
#include <src/fci/fci.h> // 2RDM and half-transformed integrals
#include <src/casscf/rotfile.h>

class Qvec : public QFile {
  protected:

  public:
    Qvec(const int n, const int m, std::shared_ptr<DensityFit> df, std::shared_ptr<Coeff> c, const size_t nclosed,
         std::shared_ptr<FCI> fci);
    Qvec(const QFile& a) : QFile(a) {};
    ~Qvec() {};

};


#endif
