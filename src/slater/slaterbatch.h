//
// Newint - Parallel electron correlation program.
// Filename: slaterbatch.h
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

#ifndef __src_slater_slaterbatch_h
#define __src_slater_slaterbatch_h 

#include <cassert>
#include <vector>
#include <src/rysint/int2d.h>
#include <src/rysint/rysint.h>
#include <src/rysint/macros.h>
#include <src/slater/svrrlist.h>
#include <memory>
#include <tuple>

class SlaterBatch : public RysInt {

  protected:
    bool yukawa_;
    double gamma_;

    /// buffer and intermediate storage
    double *bkup2_;
    std::vector<std::tuple<int, double, double> > indexpair23_;

    void perform_SVRR1();
    void perform_SVRR2();
    void perform_SVRR3();
    void perform_SVRR4();
    void perform_SVRR5();
    void perform_SVRR6();
    void perform_SVRR7();
    void perform_SVRR8();
    void perform_SVRR9();
    void perform_SVRR10();
    void perform_SVRR11();
    void perform_SVRR12();
    void perform_SVRR13();

    void perform_USVRR1();
    void perform_USVRR2();
    void perform_USVRR3();
    void perform_USVRR4();
    void perform_USVRR5();
    void perform_USVRR6();
    void perform_USVRR7();
    void perform_USVRR8();
    void perform_USVRR9();
    void perform_USVRR10();
    void perform_USVRR11();
    void perform_USVRR12();
    void perform_USVRR13();

    void perform_SVRR();
    void perform_USVRR();

    void root_weight(const int primsize_);
    void root1_direct();
    void root2_direct();

  public:
    
    SlaterBatch(const std::vector<std::shared_ptr<Shell> >, const double, const double gamma, const bool doyukawa = false);
    ~SlaterBatch();

    /// compute a batch of integrals
    void compute();

};

#endif

