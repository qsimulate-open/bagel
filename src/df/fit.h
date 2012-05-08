//
// Newint - Parallel electron correlation program.
// Filename: fit.h
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


#ifndef __SRC_DF_FIT_H
#define __SRC_DF_FIT_H

#include <src/df/df.h>
#include <src/grad/gradbatch.h>
#include <src/slater/slaterbatch.h>

class ERIFit : public DensityFit {
  protected:
    // TODO if I turn on primitive screening, it is broken.
    const double* compute_batch(std::vector<std::shared_ptr<Shell> >& input) {
      ERIBatch eribatch(input, 0.0);
      eribatch.compute();
      return eribatch.data();
    };
  public:
    ERIFit(const int nbas, const int naux,
       const std::vector<std::shared_ptr<Atom> >& atoms,  const std::vector<std::vector<int> >& offsets,
       const std::vector<std::shared_ptr<Atom> >& aux_atoms,  const std::vector<std::vector<int> >& aux_offsets, const double thr,
       const bool inverse, const double gam = 0.0) // gam is dummy
     : DensityFit(nbas, naux) {
       common_init(atoms, offsets, atoms, offsets, aux_atoms, aux_offsets, thr, inverse);

    };
    ~ERIFit() {};
    
};



class YukawaFit : public DensityFit {
  protected:
    const double gamma_;

    const double* compute_batch(std::vector<std::shared_ptr<Shell> >& input) {
      SlaterBatch slaterbatch(input, 0.0, gamma_, true); // TODO true meas it computes Yukawa and Slater together, but Slater is discarded
      slaterbatch.compute();
      return slaterbatch.data2();
    };
  public:
    YukawaFit(const int nbas, const int naux,
       const std::vector<std::shared_ptr<Atom> >& atoms,  const std::vector<std::vector<int> >& offsets,
       const std::vector<std::shared_ptr<Atom> >& aux_atoms,  const std::vector<std::vector<int> >& aux_offsets, const double thr,
       const bool inverse, const double gam)
     : DensityFit(nbas, naux), gamma_(gam) {
       common_init(atoms, offsets, atoms, offsets, aux_atoms, aux_offsets, thr, inverse);

    };
    ~YukawaFit() {};
    
};

class SlaterFit : public DensityFit {
  protected:
    const double gamma_;

    const double* compute_batch(std::vector<std::shared_ptr<Shell> >& input) {
      SlaterBatch slaterbatch(input, 0.0, gamma_, false);
      slaterbatch.compute();
      return slaterbatch.data();
    };
  public:
    SlaterFit(const int nbas, const int naux,
       const std::vector<std::shared_ptr<Atom> >& atoms,  const std::vector<std::vector<int> >& offsets,
       const std::vector<std::shared_ptr<Atom> >& aux_atoms,  const std::vector<std::vector<int> >& aux_offsets, const double thr,
       const bool inverse, const double gam)
     : DensityFit(nbas, naux), gamma_(gam) {
       common_init(atoms, offsets, atoms, offsets, aux_atoms, aux_offsets, thr, inverse);

    };
    ~SlaterFit() {};
    
};


class GradFit : public DensityFit {
  protected:
    const double* compute_batch(std::vector<std::shared_ptr<Shell> >& input) {
      GradBatch gradbatch(input, 0.0);
      gradbatch.compute();
      return gradbatch.data();
    };
  public:
    GradFit(const int nbas, const int naux,
       const std::vector<std::shared_ptr<Atom> >& atoms,  const std::vector<std::vector<int> >& offsets,
       const std::vector<std::shared_ptr<Atom> >& aux_atoms,  const std::vector<std::vector<int> >& aux_offsets, const double thr,
       const bool inverse, const double gam = 0) // gam is dummy
     : DensityFit(nbas, naux) {
       common_init(atoms, offsets, atoms, offsets, aux_atoms, aux_offsets, thr, inverse);

    };
    ~GradFit() {};
    
};

#endif

