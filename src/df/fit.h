//
// BAGEL - Parallel electron correlation program.
// Filename: fit.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
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


#ifndef __SRC_DF_FIT_H
#define __SRC_DF_FIT_H

#include <memory>
#include <src/df/df.h>
#include <src/rysint/eribatch.h>
#include <src/slater/slaterbatch.h>
#include <src/grad/gradbatch.h>
#include <src/slater/slaterbatch.h>

namespace bagel {

class ERIFit : public DensityFit {
  protected:
    std::pair<const double*, std::shared_ptr<RysInt> > compute_batch(std::array<std::shared_ptr<const Shell>,4>& input) override;

  public:
    ERIFit(const int nbas, const int naux,
       const std::vector<std::shared_ptr<const Atom> >& atoms,  const std::vector<std::vector<int> >& offsets,
       const std::vector<std::shared_ptr<const Atom> >& aux_atoms,  const std::vector<std::vector<int> >& aux_offsets, const double thr,
       const bool inverse, const double gam = 0.0) // gam is dummy
     : DensityFit(nbas, naux) {
       common_init(atoms, offsets, atoms, offsets, aux_atoms, aux_offsets, thr, inverse);
    };
    ~ERIFit() {};
    
};



class YukawaFit : public DensityFit {
  protected:
    const double gamma_;
    std::pair<const double*, std::shared_ptr<RysInt> > compute_batch(std::array<std::shared_ptr<const Shell>,4>& input) override;

  public:
    YukawaFit(const int nbas, const int naux,
       const std::vector<std::shared_ptr<const Atom> >& atoms,  const std::vector<std::vector<int> >& offsets,
       const std::vector<std::shared_ptr<const Atom> >& aux_atoms,  const std::vector<std::vector<int> >& aux_offsets, const double thr,
       const bool inverse, const double gam)
     : DensityFit(nbas, naux), gamma_(gam) {
       common_init(atoms, offsets, atoms, offsets, aux_atoms, aux_offsets, thr, inverse);
    };
    ~YukawaFit() {};
    
};

class SlaterFit : public DensityFit {
  protected:
    const double gamma_;
    std::pair<const double*, std::shared_ptr<RysInt> > compute_batch(std::array<std::shared_ptr<const Shell>,4>& input) override;

  public:
    SlaterFit(const int nbas, const int naux,
       const std::vector<std::shared_ptr<const Atom> >& atoms,  const std::vector<std::vector<int> >& offsets,
       const std::vector<std::shared_ptr<const Atom> >& aux_atoms,  const std::vector<std::vector<int> >& aux_offsets, const double thr,
       const bool inverse, const double gam)
     : DensityFit(nbas, naux), gamma_(gam) {
       common_init(atoms, offsets, atoms, offsets, aux_atoms, aux_offsets, thr, inverse);
    };
    ~SlaterFit() {};
    
};

}

#endif

