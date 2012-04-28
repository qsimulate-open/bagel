//
// Newint - Parallel electron correlation program.
// Filename: werner4.h
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


// Implements second-order CASSCF (two step)
// Uses Werner-Knowles algorithm, but construcging 4-index MO integrals.
// This is, therefere, only for reference

#ifndef __SRC_CASSCF_WERNER4_H
#define __SRC_CASSCF_WERNER4_H

#include <src/casscf/werner.h>
#include <src/casscf/jkop.h>

class Werner4 : public WernerKnowles {
  protected:
    void common_init() {
      std::cout << "    * Using the two-step 4-index Werner-Knowles algorithm" << std::endl << std::endl;
    };

  public:
    Werner4(const std::multimap<std::string, std::string> idat, const std::shared_ptr<Geometry> geom, std::shared_ptr<Reference> ref)
      : WernerKnowles(idat, geom, ref) { };
    ~Werner4() {};

    void compute();

}; 


#endif
