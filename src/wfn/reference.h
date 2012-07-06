//
// Newint - Parallel electron correlation program.
// Filename: reference.h
// Copyright (C) 2012 Toru Shiozaki
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


#ifndef _NEWINT_WFN_REFERENCE_H
#define _NEWINT_WFN_REFERENCE_H

#include <memory>
#include <vector>
#include <src/scf/coeff.h>
#include <src/scf/hcore.h>
#include <src/scf/geometry.h>
#include <src/fci/dvec.h>
#include <src/wfn/rdm.h>

// all the info to construct wave functions

class Reference : public std::enable_shared_from_this<Reference> {

  protected:
    std::shared_ptr<const Geometry> geom_; 
    std::shared_ptr<const Coeff> coeff_;

    const double energy_;

    std::shared_ptr<const Hcore> hcore_;
    std::vector<double> schwarz_;
    std::vector<double> eig_;

    int ncore_;
    const int nclosed_;
    const int nact_;
    const int nvirt_;

    int nstate_;

    std::vector<std::shared_ptr<RDM<1> > >  rdm1_;
    std::vector<std::shared_ptr<RDM<2> > >  rdm2_;
    std::shared_ptr<const RDM<1> > rdm1_av_;
    std::shared_ptr<const RDM<2> > rdm2_av_;

    // this is only for UHF gradient. Somehow I cannot come up with a beautiful design for this.
    std::shared_ptr<const Matrix1e> erdm1_;

  public:
    Reference(std::shared_ptr<const Geometry> g, std::shared_ptr<const Coeff> c,
              const double en, std::shared_ptr<const Hcore> h, const std::vector<double>& s,
              const int& nclo, const int& nact, const int& nvirt,
              const std::vector<std::shared_ptr<RDM<1> > >& rdm1 = std::vector<std::shared_ptr<RDM<1> > >(),
              const std::vector<std::shared_ptr<RDM<2> > >& rdm2 = std::vector<std::shared_ptr<RDM<2> > >(),
              std::shared_ptr<const RDM<1> > rdm1_av = std::shared_ptr<RDM<1> >(),
              std::shared_ptr<const RDM<2> > rdm2_av = std::shared_ptr<RDM<2> >());

    ~Reference() {};

    std::shared_ptr<const Geometry> geom() const { return geom_; };
    const std::vector<double>& schwarz() const { return schwarz_; };
    std::shared_ptr<const Hcore> hcore() const { return hcore_; };
    const std::shared_ptr<const Coeff> coeff() const { return coeff_; };
    void set_coeff(const std::shared_ptr<const Coeff> c) { coeff_ = c; };

    void set_eig(const std::vector<double>& eig) { eig_ = eig; };
    std::vector<double> eig() const { return eig_; };
    void set_erdm1(const std::shared_ptr<const Matrix1e> o) { erdm1_ = o; };
    std::shared_ptr<const Matrix1e> erdm1() const { return erdm1_; };

    int nclosed() const { return nclosed_; };
    int nact() const { return nact_; };
    int nvirt() const { return nvirt_; };
    int nocc() const { return nclosed_ + nact_; };
    int ncore() const { return ncore_; };
    void set_ncore(const int i) { ncore_ = i; };

    // used in SA-CASSCF
    void set_nstate(const int i) { nstate_ = i; };
    int nstate() const { return nstate_; };

    double energy() const { return energy_; };

    std::shared_ptr<const RDM<1> > rdm1(const int irdm) const { return rdm1_.at(irdm); }; 
    std::shared_ptr<const RDM<1> > rdm1_av() const { return rdm1_av_; }; 

    // returns an occ-occ sized 1RDM
    std::shared_ptr<Matrix1e> rdm1_mat(std::shared_ptr<const RDM<1> > o) const; 
    std::shared_ptr<Matrix1e> rdm1_mat(const int irdm) const { return rdm1_mat(rdm1_[irdm]); };
    std::shared_ptr<Matrix1e> rdm1_mat() const { return rdm1_mat(rdm1_av_); };

    std::shared_ptr<const RDM<2> > rdm2(const int irdm) const { return rdm2_.at(irdm); }; 
    std::shared_ptr<const RDM<2> > rdm2_av() const { return rdm2_av_; }; 

    // function to return a CI vectors from orbital info
    std::shared_ptr<Dvec> civectors() const;

};

#endif
