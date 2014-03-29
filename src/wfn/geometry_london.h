//
// BAGEL - Parallel electron correlation program.
// Filename: geometry_london.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <rreynoldschem@u.northwestern.edu>
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


#ifndef __SRC_WFN_GEOMETRY_LONDON_H
#define __SRC_WFN_GEOMETRY_LONDON_H

#include <src/df/df_london.h>
#include <src/molecule/molecule.h>
#include <src/input/input.h>

namespace bagel {

class Geometry_London: public Molecule {
  protected:
    // integral screening
    double schwarz_thresh_;
    double overlap_thresh_;

    // for DF calculations
    mutable std::shared_ptr<DFDist_London> df_;
    // small component
    mutable std::shared_ptr<DFDist_London> dfs_;
    // small-large component
    mutable std::shared_ptr<DFDist_London> dfsl_;

    // for R12 calculations
    double gamma_;

    // Constructor helpers
    void common_init2(const bool print, const double thresh, const bool nodf = false);
    void compute_integrals(const double thresh, const bool nodf);

  private:
    // serialization
    friend class boost::serialization::access;

    template<class Archive>
    void save(Archive& ar, const unsigned int) const {
      ar << boost::serialization::base_object<Molecule>(*this);
      ar << spherical_ << aux_merged_ << nbasis_ << nele_ << nfrc_ << naux_ << lmax_ << aux_lmax_
         << offsets_ << aux_offsets_ << basisfile_ << auxfile_ << schwarz_thresh_ << overlap_thresh_ << gamma_;
      const size_t dfindex = !df_ ? 0 : std::hash<DFDist_London*>()(df_.get());
      ar << dfindex;
      const bool do_rel   = !!dfs_;
      const bool do_gaunt = !!dfsl_;
      ar << do_rel << do_gaunt;
    }

    template<class Archive>
    void load(Archive& ar, const unsigned int) {
      ar >> boost::serialization::base_object<Molecule>(*this);
      ar >> spherical_ >> aux_merged_ >> nbasis_ >> nele_ >> nfrc_ >> naux_ >> lmax_ >> aux_lmax_
         >> offsets_ >> aux_offsets_ >> basisfile_ >> auxfile_ >> schwarz_thresh_ >> overlap_thresh_ >> gamma_;

      size_t dfindex;
      ar >> dfindex;
      static std::map<size_t, std::weak_ptr<DFDist_London>> dfmap;
      if (dfmap[dfindex].expired()) {
        compute_integrals(overlap_thresh_, dfindex == 0);
        dfmap[dfindex] = df_;
      } else {
        df_ = dfmap[dfindex].lock();
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
    Geometry_London() { }
    Geometry_London(const std::shared_ptr<const PTree>);
//    Geometry_London(const std::vector<std::shared_ptr<const Atom>> atoms, const std::shared_ptr<const PTree> o);
    Geometry_London(const Geometry_London& o, const std::shared_ptr<const PTree> idata, const bool discard_prev_df = true);
//    Geometry_London(const Geometry_London& o, const std::shared_ptr<const Matrix> disp, const std::shared_ptr<const PTree> geominfo, const bool rotate = true, const bool nodf = false);
//    Geometry_London(const Geometry_London& o, const std::array<double,3> disp);
//    Geometry_London(std::vector<std::shared_ptr<const Geometry_London>>);

    // Returns a constant
    int nirrep() const { return nirrep_; }
    double gamma() const {return gamma_; }
    const std::shared_ptr<const Matrix> compute_grad_vnuc() const;
    double schwarz_thresh() const { return schwarz_thresh_; }
    double overlap_thresh() const { return overlap_thresh_; }

    // The position of the specific function in the basis set.
    const std::vector<std::vector<int>>& offsets() const { return offsets_; }
    const std::vector<std::vector<int>>& aux_offsets() const { return aux_offsets_; }
    const std::vector<int>& offset(const unsigned int i) const { return offsets_.at(i); }
    const std::vector<int>& aux_offset(const unsigned int i) const { return aux_offsets_.at(i); }

    // returns schwarz screening TODO not working for DF yet
    std::vector<double> schwarz() const;

    // Returns the Petite list.
    std::shared_ptr<Petite> plist() const { return plist_; }

    // Returns DF data
    const std::shared_ptr<const DFDist_London> df() const { return df_; }
    const std::shared_ptr<const DFDist_London> dfs() const { return dfs_; }
    const std::shared_ptr<const DFDist_London> dfsl() const { return dfsl_; }
    // TODO resolve "mutable" issues
    void discard_df() const { df_.reset(); dfs_.reset(); dfsl_.reset(); }

    // In R12 methods, we need to construct a union of OBS and CABS.
    // Currently, this is done by creating another object and merge OBS and CABS into atoms_.
    // After this, compute_nuclear_repulsion() should not be called.
    // Not undo-able.
    void merge_obs_aux();

    // type T should be a derived class of DFDist
    template<typename T>
    std::shared_ptr<T> form_fit(const double thr, const bool inverse, const double gam = 0.0, const bool average = false) const {
      return std::make_shared<T>(nbasis(), naux(), atoms(), aux_atoms(), thr, inverse, gam, average);
    }

    // initialize relativistic components
    std::shared_ptr<const Geometry_London> relativistic(const bool do_gaunt) const;
    void compute_relativistic_integrals(const bool do_gaunt);
    void discard_relativistic() const;

};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::Geometry_London)

#endif
