//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: nevpt2.h
// Copyright (C) 2014 Toru Shiozaki
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


#ifndef __SRC_NEVPT2_NEVPT2_H
#define __SRC_NEVPT2_NEVPT2_H

#include <src/pt2/mp2/mp2cache.h>
#include <src/wfn/method.h>
#include <src/df/dfdistt.h>
#include <src/df/reldffullt.h>
#include <src/util/math/quatmatrix.h>
#include <src/wfn/relreference.h>

namespace bagel {

template<typename DataType>
class NEVPT2 : public Method {
  protected:
    using MatType  = typename std::conditional<std::is_same<DataType,double>::value, Matrix, ZMatrix>::type;
    using DiagType = typename std::conditional<std::is_same<DataType,double>::value, Matrix, QuatMatrix>::type;
    using VecType  = typename std::conditional<std::is_same<DataType,double>::value, VectorB, ZVectorB>::type;
    using ViewType = typename std::conditional<std::is_same<DataType,double>::value, MatView, ZMatView>::type;
    using DFType   = typename std::conditional<std::is_same<DataType,double>::value, DFDistT, ListRelDFFullT>::type;

  protected:
    int ncore_;
    int nfrozenvirt_;
    int nclosed_;
    int nact_;
    int nvirt_;
    int istate_;
    double norm_thresh_;

    bool gaunt_;
    bool breit_;

    std::string abasis_;

    double energy_;

    // density matrices to be used
    // particle RDMs
    std::shared_ptr<const MatType> rdm1_;
    std::shared_ptr<const MatType> rdm2_;
    std::shared_ptr<const MatType> rdm3_;
    std::shared_ptr<const MatType> rdm4_;
    // hole RDMs
    std::shared_ptr<const MatType> hrdm1_;
    std::shared_ptr<const MatType> hrdm2_;
    std::shared_ptr<const MatType> hrdm3_;
    // <a+a b+b c+c..>
    std::shared_ptr<const MatType> ardm2_;
    std::shared_ptr<const MatType> ardm3_;
    std::shared_ptr<const MatType> ardm4_;
    // <a+a bb+>
    std::shared_ptr<const MatType> srdm2_;
    // <a+a bb+ c+c>
    std::shared_ptr<const MatType> srdm3_;

    // integrals in physicists notation
    std::shared_ptr<const MatType> ints2_;
    std::shared_ptr<      MatType> fockact_;
    std::shared_ptr<      MatType> fockact_c_;
    std::shared_ptr<      MatType> fockact_h_;
    std::shared_ptr<      MatType> fockact_p_;

    // K and K'mat
    std::shared_ptr<const MatType> qvec_;
    std::shared_ptr<const MatType> kmat_;
    std::shared_ptr<const MatType> kmatp_;
    std::shared_ptr<const MatType> kmat2_;
    std::shared_ptr<const MatType> kmatp2_;

    // A
    std::shared_ptr<const MatType> amat2_;
    std::shared_ptr<const MatType> amat3_;
    std::shared_ptr<const MatType> amat3t_;
    // B
    std::shared_ptr<const MatType> bmat2_;
    std::shared_ptr<const MatType> bmat2t_;
    // C
    std::shared_ptr<const MatType> cmat2_;
    std::shared_ptr<const MatType> cmat2t_;
    // D
    std::shared_ptr<const MatType> dmat2_;
    std::shared_ptr<const MatType> dmat1_;
    std::shared_ptr<const MatType> dmat1t_;

    void init_reference();
    std::shared_ptr<MatType> compute_fock(std::shared_ptr<const Geometry> cgeom, std::shared_ptr<const MatType> hcore,
                                          const ViewType coeff, const double scale_exch = 1.0, const double scale_coulomb = 1.0);
    std::tuple<std::shared_ptr<DFType>,std::shared_ptr<DFType>,std::shared_ptr<DFType>,std::shared_ptr<DFType>>
      compute_full_nevpt2(std::shared_ptr<const Geometry>, std::shared_ptr<const MatType>, std::shared_ptr<const MatType>,
                          std::shared_ptr<const MatType>, std::shared_ptr<const MatType>, const bool gaunt, const bool breit) const;
    std::tuple<std::map<int,std::shared_ptr<const MatType>>,std::map<int,std::shared_ptr<const MatType>>,std::map<int,std::shared_ptr<const MatType>>>
      slice_ax(std::shared_ptr<const DFType> o) const;

    std::tuple<std::map<int,std::shared_ptr<MP2Cache_<DataType>>>,std::map<int,std::shared_ptr<MP2Cache_<DataType>>>,std::map<int,std::shared_ptr<MP2Cache_<DataType>>>>
      set_up_cache(std::shared_ptr<const Geometry>, std::shared_ptr<DFType> fullvi, std::shared_ptr<DFType> fullvi2, std::shared_ptr<DFType> fullvi3,
                   std::vector<std::vector<std::tuple<int,int,MP2Tag<DataType>,MP2Tag<DataType>>>>) const;

    void compute_rdm();
    void compute_hrdm();
    void compute_asrdm();
    void compute_ints();
    void compute_kmat();
    void compute_abcd();

    std::shared_ptr<const MatType> coeff() const;
    std::tuple<std::shared_ptr<MatType>,VectorB> remove_core(std::shared_ptr<const MatType> in, const VectorB& eig) const;
    std::tuple<std::shared_ptr<MatType>,VectorB> remove_frozenvirt(std::shared_ptr<const MatType> in, const VectorB& eig) const;

  public:
    NEVPT2(std::shared_ptr<const PTree>, std::shared_ptr<const Geometry>, std::shared_ptr<const Reference> = nullptr);

    virtual void compute() override;
    virtual std::shared_ptr<const Reference> conv_to_ref() const override { return ref_; }

    double energy() const { return energy_; }
    int ncore() const { return ncore_; }
    std::string abasis() const { return abasis_; }
};

template<> void NEVPT2<double>::init_reference();
template<> void NEVPT2<std::complex<double>>::init_reference();

template<> void NEVPT2<double>::compute_rdm();
template<> void NEVPT2<std::complex<double>>::compute_rdm();

template<> std::shared_ptr<const Matrix> NEVPT2<double>::coeff() const;
template<> std::shared_ptr<const ZMatrix> NEVPT2<std::complex<double>>::coeff() const;
template<> std::tuple<std::shared_ptr<Matrix>,VectorB> NEVPT2<double>::remove_core(std::shared_ptr<const Matrix>, const VectorB&) const;
template<> std::tuple<std::shared_ptr<Matrix>,VectorB> NEVPT2<double>::remove_frozenvirt(std::shared_ptr<const Matrix>, const VectorB&) const;
template<> std::tuple<std::shared_ptr<ZMatrix>,VectorB> NEVPT2<std::complex<double>>::remove_core(std::shared_ptr<const ZMatrix>, const VectorB&) const;
template<> std::tuple<std::shared_ptr<ZMatrix>,VectorB> NEVPT2<std::complex<double>>::remove_frozenvirt(std::shared_ptr<const ZMatrix>, const VectorB&) const;

template<> std::shared_ptr<Matrix> NEVPT2<double>::compute_fock(std::shared_ptr<const Geometry> cgeom, std::shared_ptr<const Matrix> hcore,
                                                                const MatView coeff, const double exch, const double coulomb);
template<> std::shared_ptr<ZMatrix> NEVPT2<std::complex<double>>::compute_fock(std::shared_ptr<const Geometry> cgeom, std::shared_ptr<const ZMatrix> hcore,
                                                                               const ZMatView coeff, const double exch, const double coulomb);

template<> std::tuple<std::shared_ptr<DFDistT>,std::shared_ptr<DFDistT>,std::shared_ptr<DFDistT>,std::shared_ptr<DFDistT>>
  NEVPT2<double>::compute_full_nevpt2(std::shared_ptr<const Geometry>, std::shared_ptr<const Matrix>, std::shared_ptr<const Matrix>,
                                      std::shared_ptr<const Matrix>, std::shared_ptr<const Matrix>, const bool, const bool) const;
template<> std::tuple<std::shared_ptr<ListRelDFFullT>,std::shared_ptr<ListRelDFFullT>,std::shared_ptr<ListRelDFFullT>,std::shared_ptr<ListRelDFFullT>>
  NEVPT2<std::complex<double>>::compute_full_nevpt2(std::shared_ptr<const Geometry>, std::shared_ptr<const ZMatrix>, std::shared_ptr<const ZMatrix>,
                                                    std::shared_ptr<const ZMatrix>, std::shared_ptr<const ZMatrix>, const bool, const bool) const;

template<> std::tuple<std::map<int,std::shared_ptr<const Matrix>>,std::map<int,std::shared_ptr<const Matrix>>,std::map<int,std::shared_ptr<const Matrix>>>
              NEVPT2<double>::slice_ax(std::shared_ptr<const DFDistT> o) const;
template<> std::tuple<std::map<int,std::shared_ptr<const ZMatrix>>,std::map<int,std::shared_ptr<const ZMatrix>>,std::map<int,std::shared_ptr<const ZMatrix>>>
              NEVPT2<std::complex<double>>::slice_ax(std::shared_ptr<const ListRelDFFullT> o) const;

template<>
std::tuple<std::map<int,std::shared_ptr<MP2Cache_<double>>>,
           std::map<int,std::shared_ptr<MP2Cache_<double>>>,
           std::map<int,std::shared_ptr<MP2Cache_<double>>>>
              NEVPT2<double>::set_up_cache(std::shared_ptr<const Geometry>,
                                           std::shared_ptr<DFType> fullvi, std::shared_ptr<DFType> fullvi2, std::shared_ptr<DFType> fullvi3,
                                           std::vector<std::vector<std::tuple<int,int,MP2Tag<double>,MP2Tag<double>>>>) const;
template<>
std::tuple<std::map<int,std::shared_ptr<MP2Cache_<std::complex<double>>>>,
           std::map<int,std::shared_ptr<MP2Cache_<std::complex<double>>>>,
           std::map<int,std::shared_ptr<MP2Cache_<std::complex<double>>>>>
              NEVPT2<std::complex<double>>::set_up_cache(std::shared_ptr<const Geometry>,
                                                         std::shared_ptr<DFType> fullvi, std::shared_ptr<DFType> fullvi2, std::shared_ptr<DFType> fullvi3,
                                                         std::vector<std::vector<std::tuple<int,int,MP2Tag<std::complex<double>>,MP2Tag<std::complex<double>>>>>) const;

extern template class NEVPT2<double>;
extern template class NEVPT2<std::complex<double>>;

}

#endif
