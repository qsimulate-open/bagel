//
// BAGEL - Parallel electron correlation program.
// Filename: scf_base_london.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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

#ifndef __scf_london_base_london_h
#define __scf_london_base_london_h

#include <src/molecule/zoverlap.h>
#include <src/molecule/zhcore.h>
#include <src/london/zcoeff.h>
#include <src/london/fock_london.h>
#include <src/wfn/method.h>

namespace bagel {

class SCF_base_London : public Method_London {
  protected:
    std::shared_ptr<const ZMatrix> tildex_;
    std::shared_ptr<const ZOverlap> overlap_;
    std::shared_ptr<const ZHcore> hcore_;
    std::shared_ptr<const ZCoeff> coeff_;

    int max_iter_;

    int diis_start_;
    int diis_size_;

    double thresh_overlap_;
    double thresh_scf_;
    int multipole_print_;

    std::vector<double> schwarz_;
    void init_schwarz();

    std::vector<double> eig_;
    double energy_;

    int nocc_;
    int noccB_;

    const std::string indent = "  ";

    // when gradient is requested, we store half-transformed integrals
    // TODO so far only implemented in closed-shell SCF
    bool do_grad_;
    std::shared_ptr<ComplexDFHalfDist> half_;

    bool restart_;

  private:
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Method_London);
      ar & tildex_ & overlap_ & hcore_ & coeff_ & max_iter_ & diis_start_ & diis_size_
         & thresh_overlap_ & thresh_scf_ & multipole_print_ & schwarz_ & eig_ & energy_
         & nocc_ & noccB_ & do_grad_ & restart_;
    }

  public:
    SCF_base_London() { }
    SCF_base_London(const std::shared_ptr<const PTree> idata_, const std::shared_ptr<const Geometry_London>,
             const std::shared_ptr<const Reference>, const bool need_schwarz = false);
    virtual ~SCF_base_London() { }

    virtual void compute() override = 0;

    const std::shared_ptr<const ZCoeff> coeff() const { return coeff_; };
    void set_coeff(const std::shared_ptr<ZCoeff> o) { coeff_ = o; };

    const std::shared_ptr<const ZHcore> hcore() const { return hcore_; };
    const std::vector<double>& schwarz() const { return schwarz_; };

    int nocc() const { return nocc_; };
    int noccB() const { return noccB_; };
    double energy() const { return energy_; };

    double thresh_overlap() const { return thresh_overlap_; }
    double thresh_scf() const { return thresh_scf_; }

    virtual std::shared_ptr<const Reference> conv_to_ref() const override = 0;

    double* eig() { return eig_.data(); };

    std::shared_ptr<ComplexDFHalfDist> half() const { return half_; }
    void discard_half() { half_.reset(); }
};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::SCF_base_London)

#endif
