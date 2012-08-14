//
// Newint - Parallel electron correlation program.
// Filename: symmat.h
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


#ifndef __src_scf_symmat_h
#define __src_scf_symmat_h

#include <src/scf/geometry.h>
#include <src/scf/matrix1e.h>
#include <src/scf/petite.h>
#include <src/scf/symrot.h>

class SymMat : public Matrix1e {
  protected:
    std::shared_ptr<SymRotAbel> symrot_; 
    std::shared_ptr<Petite> petite_; 
    void computebatch(const std::vector<std::shared_ptr<const Shell> >&, const int, const int, const int) override {};

  public:
    SymMat(const std::shared_ptr<const Geometry>, const int);
    ~SymMat(); 

};

#endif
