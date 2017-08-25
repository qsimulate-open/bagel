//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: geometry.h
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


#ifndef __SRC_WFN_GEOMETRY_H
#define __SRC_WFN_GEOMETRY_H

#include <src/df/df.h>
#include <src/util/input/input.h>
#include <src/molecule/molecule.h>
#include <src/wfn/hcoreinfo.h>
#include <src/wfn/fmminfo.h>

namespace bagel {

class Geometry : public Molecule {
  protected:
    // integral screening
    double schwarz_thresh_;
    double overlap_thresh_;

    // for DF calculations
    mutable std::shared_ptr<DFDist> df_;
    // small component
    mutable std::shared_ptr<DFDist> dfs_;
    // small-large component
    mutable std::shared_ptr<DFDist> dfsl_;

    // Constructor helpers
    void common_init2(const bool print, const double thresh, const bool nodf = false);
    void compute_integrals(const double thresh) const;
    void get_electric_field(std::shared_ptr<const PTree>& geominfo);
    void set_london(std::shared_ptr<const PTree>& geominfo);
    void init_magnetism();

    // Magnetism-specific parameters
    bool magnetism_;
    bool london_;
    bool use_finite_;

    // Lattice parameters
    bool do_periodic_df_;
    std::vector<std::array<double, 3>> primitive_vectors_;

    // Hcore Information
    std::shared_ptr<const HcoreInfo> hcoreinfo_;

    // FMM
    std::shared_ptr<const FMMInfo> fmm_;

  private:
    // serialization
    friend class boost::serialization::access;

    template<class Archive>
    void save(Archive& ar, const unsigned int) const {
      ar << boost::serialization::base_object<Molecule>(*this);
      ar << schwarz_thresh_ << overlap_thresh_ << magnetism_ << london_ << use_finite_ << do_periodic_df_ << hcoreinfo_ << fmm_;
      const size_t dfindex = !df_ ? 0 : std::hash<DFDist*>()(df_.get());
      ar << dfindex;
      const bool do_rel   = !!dfs_;
      const bool do_gaunt = !!dfsl_;
      ar << do_rel << do_gaunt;
    }

    template<class Archive>
    void load(Archive& ar, const unsigned int) {
      ar >> boost::serialization::base_object<Molecule>(*this);
      ar >> schwarz_thresh_ >> overlap_thresh_ >> magnetism_ >> london_ >> use_finite_ >> do_periodic_df_ >> hcoreinfo_ >> fmm_;
      size_t dfindex;
      ar >> dfindex;
      static std::map<size_t, std::weak_ptr<DFDist>> dfmap;
      if (!fmm_) {
        if (dfmap[dfindex].expired()) {
          compute_integrals(overlap_thresh_);
          dfmap[dfindex] = df_;
        } else {
          df_ = dfmap[dfindex].lock();
        }
      }

      bool do_rel, do_gaunt;
      ar >> do_rel >> do_gaunt;
      if (do_rel)
        compute_relativistic_integrals(do_gaunt);
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned int file_version) {
      boost::serialization::split_member(ar, *this, file_version);
    }

  public:
    Geometry() { }
    Geometry(std::shared_ptr<const PTree> idata);
    Geometry(const std::vector<std::shared_ptr<const Atom>> atoms, std::shared_ptr<const PTree> o);
    Geometry(const Geometry& o, std::shared_ptr<const PTree> idata, const bool discard_prev_df = true);
    Geometry(const Geometry& o, std::shared_ptr<const Matrix> disp, std::shared_ptr<const PTree> geominfo, const bool rotate = true, const bool nodf = false);
    Geometry(const Geometry& o, const std::array<double,3> disp);
    Geometry(std::vector<std::shared_ptr<const Geometry>>, const bool nodf = false);
    Geometry(const Geometry& o, const std::string extent_type);

    // Gradients of the nuclear-nuclear potential energy
    std::shared_ptr<const Matrix> compute_grad_vnuc() const;

    // Thresholds
    double schwarz_thresh() const { return schwarz_thresh_; }
    double overlap_thresh() const { return overlap_thresh_; }
    bool london() const { return london_; }
    bool magnetism() const { return magnetism_; }

    // returns schwarz screening TODO not working for DF yet
    std::vector<double> schwarz() const;

    // Returns DF data
    std::shared_ptr<const DFDist> df() const { return df_; }
    std::shared_ptr<const DFDist> dfs() const { return dfs_; }
    std::shared_ptr<const DFDist> dfsl() const { return dfsl_; }

    // TODO resolve "mutable" issues
    void discard_df() const { df_.reset(); dfs_.reset(); dfsl_.reset(); }

    // type T should be a derived class of DFDist
    template<typename T>
    std::shared_ptr<T> form_fit(const double thr, const bool inverse, const double gam = 0.0, const bool average = false, const std::shared_ptr<Matrix> d2 = nullptr) const {
      return std::make_shared<T>(nbasis(), naux(), atoms(), aux_atoms(), thr, inverse, gam, average, d2);
    }

    // initialize relativistic components
    std::shared_ptr<const Geometry> relativistic(const bool do_gaunt, const bool do_coulomb = true) const;
    void compute_relativistic_integrals(const bool do_gaunt);
    void discard_relativistic() const;

    // Lattice
    std::vector<std::array<double, 3>> primitive_vectors() const { return primitive_vectors_; }
    std::array<double, 3> primitive_vectors(const int i) const { return primitive_vectors_[i]; };
    bool do_periodic_df() const { return do_periodic_df_; }
    std::shared_ptr<const Geometry> periodic(std::vector<std::shared_ptr<const Atom>> new_atoms) const;

    // Hcore Information
    std::shared_ptr<const HcoreInfo> hcoreinfo() const { return hcoreinfo_; }

    // FMM
    std::shared_ptr<const FMMInfo> fmm() const { return fmm_; }
};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::Geometry)

#endif
