//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: scf_base.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef __SCF_SCF_BASE_H
#define __SCF_SCF_BASE_H

#include <src/mat1e/overlap.h>
#include <src/mat1e/giao/zoverlap.h>
#include <src/mat1e/hcore.h>
#include <src/mat1e/giao/zhcore.h>
#include <src/wfn/method.h>
#include <src/scf/fmm/fmm.h>

namespace bagel {

template <typename MatType = Matrix, typename OvlType = Overlap, typename HcType = Hcore,
          class Enable = typename std::enable_if<((std::is_same<MatType, Matrix>::value && std::is_same<OvlType, Overlap>::value && std::is_same<HcType, Hcore>::value)
                      || (std::is_same<MatType, ZMatrix>::value && std::is_same<OvlType, ZOverlap>::value && std::is_same<HcType, ZHcore>::value))>::type>
class SCF_base_ : public Method {
  protected:
    std::shared_ptr<const MatType> tildex_;
    std::shared_ptr<const OvlType> overlap_;
    std::shared_ptr<const HcType> hcore_;
    std::shared_ptr<const Coeff_<MatType>> coeff_;

    int max_iter_;

    int diis_start_;
    int diis_size_;

    double thresh_overlap_;
    double thresh_scf_;
    int multipole_print_;
    int dma_print_;

    std::vector<double> schwarz_;
    void init_schwarz();

    VectorB eig_;
    double energy_;
    std::vector<double> scf_dipole_;

    int nocc_;
    int noccB_;

    const std::string indent = "  ";

    // when gradient is requested, we store half-transformed integrals
    // TODO so far only implemented in closed-shell SCF
    bool do_grad_;
    std::shared_ptr<DFHalfDist> half_;

    // FMM
    bool dofmm_;
    std::shared_ptr<const FMM> fmm_;
    std::shared_ptr<const FMM> fmmK_;

    bool restart_;

    void get_coeff(const std::shared_ptr<const Reference> ref) { coeff_ = ref->coeff(); }

  private:
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<Method>(*this);
      ar & tildex_ & overlap_ & hcore_ & coeff_ & max_iter_ & diis_start_ & diis_size_
         & thresh_overlap_ & thresh_scf_ & multipole_print_ & dma_print_ & schwarz_ & eig_ & energy_
         & nocc_ & noccB_ & do_grad_ & restart_ & dofmm_ & fmm_ & fmmK_;
    }

  public:
    SCF_base_() { }
    SCF_base_(std::shared_ptr<const PTree> idata_, std::shared_ptr<const Geometry>, std::shared_ptr<const Reference>, const bool need_schwarz = false);
    virtual ~SCF_base_() { }

    virtual void compute() override = 0;

    const std::shared_ptr<const Coeff_<MatType>> coeff() const { return coeff_; }
    void set_coeff(const std::shared_ptr<Coeff_<MatType>> o) { coeff_ = o; }

    const std::shared_ptr<const HcType> hcore() const { return hcore_; }
    const std::vector<double>& schwarz() const { return schwarz_; }

    int nocc() const { return nocc_; }
    int noccB() const { return noccB_; }
    double energy() const { return energy_; }
    const std::vector<double>& scf_dipole() const { return scf_dipole_; }

    double thresh_overlap() const { return thresh_overlap_; }
    double thresh_scf() const { return thresh_scf_; }

    virtual std::shared_ptr<const Reference> conv_to_ref() const override = 0;

    VectorB& eig() { return eig_; }

    std::shared_ptr<DFHalfDist> half() const { return half_; }
    void discard_half() { half_.reset(); }
};

// specialized for GIAO cases
template <>
void SCF_base_<ZMatrix, ZOverlap, ZHcore, std::enable_if<true>::type>::get_coeff(const std::shared_ptr<const Reference> ref);

using SCF_base = SCF_base_<Matrix, Overlap, Hcore>;
using SCF_base_London = SCF_base_<ZMatrix, ZOverlap, ZHcore>;

}

extern template class bagel::SCF_base_<bagel::Matrix, bagel::Overlap, bagel::Hcore>;
extern template class bagel::SCF_base_<bagel::ZMatrix, bagel::ZOverlap, bagel::ZHcore>;

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::SCF_base)
BOOST_CLASS_EXPORT_KEY(bagel::SCF_base_London)

#endif
