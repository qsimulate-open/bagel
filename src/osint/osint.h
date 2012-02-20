//
// Newint - Parallel electron correlation program.
// Filename: osint.h
// Copyright (C) 2009 Toru Shiozaki
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


#ifndef __src_osint_osint_h
#define __src_osint_osint_h

#include <vector>
#include <src/scf/shell.h>
#include <src/rysint/sortlist.h>
#include <memory>

class OSInt {
  protected:
    bool spherical_;

    std::vector<std::shared_ptr<Shell> > basisinfo_;

    double* data_;
    std::vector<double> xp_, xa_, xb_, rho_, p_; 
    std::vector<double> coeffsx_, coeffsy_, coeffsz_;
    std::vector<double> coefftx_, coeffty_, coefftz_;
    double AB_[3];

    int ang0_, ang1_, cont0_, cont1_, prim0_, prim1_;

    int amax_, amax1_, amin_, asize_, asize_final_, asize_intermediate_;
    std::vector<int> amapping_;

    bool swap01_;

    SortList sort_;

    virtual void perform_VRR(double*) {};
    void perform_contraction(const int, const double*, const int, const int, double*, 
                             const std::vector<std::vector<double> >&, const std::vector<std::pair<int, int> >&, const int, 
                             const std::vector<std::vector<double> >&, const std::vector<std::pair<int, int> >&, const int);

  public:
    OSInt(const std::vector<std::shared_ptr<Shell> >&);
    ~OSInt();

    virtual void compute() {}; 

    const double* data() const { return data_; }; 

};

#endif

