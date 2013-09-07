
// BAGEL - Parallel electron correlation program.
// Filename: relmofile.h
// Copyright (C) 2013 Michael Caldwell
//
// Author: Michael Caldwell  <caldwell@u.northwestern.edu>
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


#ifndef __BAGEL_ZFCI_RELMOFILE_H
#define __BAGEL_ZFCI_RELMOFILE_H


#include <src/zfci/zmofile_base.h>
#include <src/rel/relreference.h>

namespace bagel {

class RelMOFile : public ZMOFile_Base {

  protected:
    int norb_rel_;
    std::shared_ptr<ZMatrix> core_dfock_;
    // creates integral files and returns the core energy.
    double create_Jiiii(const int, const int) override;
    std::tuple<std::shared_ptr<const ZMatrix>, std::shared_ptr<const ZMatrix>> kramers_block(std::shared_ptr<const ZMatrix> buf1e, std::shared_ptr<const ZMatrix> buf2e);
    void compress(std::shared_ptr<const ZMatrix> buf1e, std::shared_ptr<const ZMatrix> buf2e) override;

    std::shared_ptr<const Geometry> relgeom_;
    std::shared_ptr<const RelReference> relref;
  public:
    RelMOFile(const std::shared_ptr<const Reference>, const std::string method = std::string("KH"));
    RelMOFile(const std::shared_ptr<const Reference>, const std::shared_ptr<const Coeff>, const std::string method = std::string("KH"));
    std::shared_ptr<const ZMatrix> core_fock() const { return core_dfock_; };
    std::complex<double>* core_dfock_ptr() { return core_dfock_->data(); };
    const std::complex<double>* core_dfock_ptr() const { return core_dfock_->data(); };
};


class RelJop : public RelMOFile {
  protected:
    std::tuple<std::shared_ptr<const ZMatrix>, double> compute_mo1e(const int, const int) override;
    std::shared_ptr<const ZMatrix> compute_mo2e(const int, const int) override;
  public:
    RelJop(const std::shared_ptr<const Reference> b, const int c, const int d, const std::string f = std::string("KH"))
      : RelMOFile(b, f) { core_energy_ = create_Jiiii(c, d); }
};


class RelHtilde : public ZHtilde_Base, public RelMOFile {
  protected:
    std::tuple<std::shared_ptr<const ZMatrix>, double> compute_mo1e(const int, const int) override { return std::make_tuple(h1_tmp_, 0.0); };
    std::shared_ptr<const ZMatrix> compute_mo2e(const int, const int) override { return h2_tmp_; };
  public:
    RelHtilde(const std::shared_ptr<const Reference> b, const int c, const int d, std::shared_ptr<const ZMatrix> h1, std::shared_ptr<const ZMatrix> h2)
      : ZHtilde_Base(h1, std::move(h2)), RelMOFile(b) {
      core_energy_ = create_Jiiii(c, d);
    }
};

}

#endif
