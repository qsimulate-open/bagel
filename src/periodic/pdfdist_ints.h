//
// BAGEL - Parallel electron correlation program.
// Filename: pdfdist_ints.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
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

#ifndef __SRC_PERIODIC_PDFDIST_INTS_H
#define __SRC_PERIODIC_PDFDIST_INTS_H

#include <src/df/df.h>

namespace bagel {

/*periodic version of DFDist_ints*/
class PDFDist_ints : public DFDist {
  protected:
    void compute_3index(const std::vector<std::shared_ptr<const Shell>>& ashell,   /*aux   */
                        const std::vector<std::shared_ptr<const Shell>>& b0shell,  /*cell 0*/
                        const std::vector<std::shared_ptr<const Shell>>& bgshell); /*cell g*/

  public:
    PDFDist_ints(const int nbas, const int naux, const std::vector<std::shared_ptr<const Atom>>& atoms_c0,
                 std::vector<std::shared_ptr<const Atom>>& atoms_cg, const std::vector<std::shared_ptr<const Atom>>& aux_atoms,
                 const double thr, const bool inverse, const std::shared_ptr<Matrix> data2 = nullptr)
    : DFDist(nbas, naux) {

      // 3index integrals made in DFBlock.
      std::vector<std::shared_ptr<const Shell>> ashell, b0shell, bgshell;
      for (auto& i : aux_atoms)     ashell.insert(ashell.end(),  i->shells().begin(), i->shells().end());
      for (auto& i : atoms_c0)     b0shell.insert(b0shell.end(), i->shells().begin(), i->shells().end());
      for (auto& i : atoms_cg)     bgshell.insert(bgshell.end(), i->shells().begin(), i->shells().end());

      // distribute auxiliary shells to each nodes
      int astart;
      std::vector<std::shared_ptr<const Shell>> myashell;
      std::tie(astart, myashell) = get_ashell(ashell);

      std::shared_ptr<const StaticDist> adist_shell = make_table(astart);
      std::shared_ptr<const StaticDist> adist_averaged = std::make_shared<const StaticDist>(naux_, mpi__->size());

      // make empty dfblocks
      const size_t asize  = std::accumulate(myashell.begin(),myashell.end(),0, [](const int& i, const std::shared_ptr<const Shell>& o) { return i+o->nbasis(); });
      const size_t b0size = std::accumulate(b0shell.begin(), b0shell.end(), 0, [](const int& i, const std::shared_ptr<const Shell>& o) { return i+o->nbasis(); });
      const size_t bgsize = std::accumulate(bgshell.begin(), bgshell.end(), 0, [](const int& i, const std::shared_ptr<const Shell>& o) { return i+o->nbasis(); });
      block_.push_back(std::make_shared<DFBlock>(adist_shell, adist_averaged, asize, b0size, bgsize, astart, 0, 0));

      // 3-index integrals
      compute_3index(myashell, b0shell, bgshell);

      // 2-index integrals
      if (data2)
        data2_ = data2;
      else
        compute_2index(ashell, thr, inverse);
    }
};

}

#endif
