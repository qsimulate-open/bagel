//
// BAGEL - Parallel electron correlation program.
// Filename: slatermem.h
// Copyright (C) 2013 Toru Shiozaki
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

#ifndef __SRC_RYSINT_SLATERMEM_H
#define __SRC_RYSINT_SLATERMEM_H

#include <memory>
#include <vector>

namespace bagel {

class SlaterMem {
  private:
    std::vector<std::unique_ptr<double[]>> datax_; 
    std::vector<std::unique_ptr<double[]>> dataw_; 
    void fill1();
    void fill2();
    void fill3();
    void fill4();
    void fill5();
    void fill6();
    void fill7();
    void fill8();
    void fill9();
    void fill10();
    void fill11();
    void fill12();
    void fill13();

  public:
    SlaterMem();
    const double* x(const int i) const { return datax_[i-1].get(); }
    const double* w(const int i) const { return dataw_[i-1].get(); }
};

extern const SlaterMem slatermem__;

}

#endif
