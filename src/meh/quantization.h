//
// BAGEL - Parallel electron correlation program.
// Filename: quantization.h
// Copyright (C) 2013 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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

#ifndef __meh_quantization_h
#define __meh_quantization_h

#include <utility>

#include <src/fci/dvec.h>

namespace bagel {

/************************************************************************************
*  This class computes the action of second quantization operators on a CAS state   *
************************************************************************************/
enum class SQ {
  CreateAlpha = 1,
  AnnihilateAlpha = -1,
  CreateBeta = 2,
  AnnihilateBeta = -2
};

class Quantization { 
  protected:
    static const int Alpha = 0;
    static const int Beta = 1;

    static const int Annihilate = 0;
    static const int Create = 1;

  public:
    Quantization() {}

    virtual const int ij(const int norb) const = 0;
    virtual std::shared_ptr<Dvec> compute(std::shared_ptr<const Civec> ccvec) const = 0;
};

template <SQ oper>
class OneBody : public Quantization {
  public:
    OneBody() {}

    const int ij(const int norb) const override { return norb; }
    std::shared_ptr<Dvec> compute(std::shared_ptr<const Civec> ccvec) const override {
      const int action = (static_cast<int>(oper) > 0 ? Create : Annihilate);
      const int spin = ( (oper==SQ::CreateAlpha || oper==SQ::AnnihilateAlpha) ? Alpha : Beta ); 

      std::shared_ptr<const Determinants> source_det = ccvec->det();
      std::shared_ptr<const Determinants> target_det = ( (spin == Alpha) ? 
        (action == Annihilate ? ccvec->det()->remalpha() : ccvec->det()->addalpha() ) : 
        (action == Annihilate ? ccvec->det()->rembeta() : ccvec->det()->addbeta() )); 

      const int norb = target_det->norb();

      const int source_lena = source_det->lena();
      const int source_lenb = source_det->lenb();

      const int target_lena = target_det->lena();
      const int target_lenb = target_det->lenb();

      const int length = (spin == Alpha ? target_lenb : target_lena);
      const int source_start = ( spin == Alpha ? source_lenb : 1 );
      const int target_start = ( spin == Alpha ? target_lenb : 1 );
      const int source_stride = ( spin == Alpha ? 1 : source_lenb );
      const int target_stride = ( spin == Alpha ? 1 : target_lenb );

      std::shared_ptr<Dvec> out(new Dvec(target_det, norb));

      const double* source_base = ccvec->data();

      for (int i = 0; i < norb; ++i) {
        double* target_base = out->data(i)->data();

        for(auto& iter : (spin == Alpha ? (action == Annihilate ? source_det->phidowna(i) : source_det->phiupa(i)) : (action == Annihilate ? source_det->phidownb(i) : source_det->phiupb(i)) ) ) {
          const double sign = static_cast<double>( iter.sign );
          double* target = target_base + target_start * iter.target;
          const double* source = source_base + source_start * iter.source;
          daxpy_(length, sign, source, source_stride, target, target_stride);
        }   
      }

      return out;
    }
};


template <SQ oper1, SQ oper2>
class TwoBody : public Quantization {
  public:
    TwoBody() {}

    const int ij(const int norb) const override { return norb*norb; }
    std::shared_ptr<Dvec> compute(std::shared_ptr<const Civec> ccvec) const override {
      if ( ((static_cast<int>(oper1) + static_cast<int>(oper2)) == 0) && (static_cast<int>(oper2) < 0) ) {
        const int spin = (oper1==SQ::CreateAlpha ? Alpha : Beta);
        std::shared_ptr<const Determinants> det = ccvec->det();
        const int norb = det->norb();

        const int lena = det->lena();
        const int lenb = det->lenb();

        const int length = (spin == Alpha ? lenb : lena);
        const int start = ( spin == Alpha ? lenb : 1 );
        const int stride = ( spin == Alpha ? 1 : lenb );

        const int sizeij = norb * norb;
        std::shared_ptr<Dvec> out(new Dvec(det, sizeij));

        const double* source_base = ccvec->data();

        for (int ij = 0; ij < sizeij; ++ij) {
          double* target_base = out->data(ij)->data();
          for(auto& iter : (spin==Alpha ? det->phia(ij) : det->phib(ij))) {
            const double sign = static_cast<double>(iter.sign);
            double* target = target_base + start * iter.target;
            const double* source = source_base + start * iter.source;
            daxpy_(length, sign, source, stride, target, stride);
          }
        }

        return out;
      }
      else {
        const int norb = ccvec->det()->norb();

        OneBody<oper2> first;
        OneBody<oper1> second;

        std::shared_ptr<Dvec> intermediate = first.compute(ccvec);

        std::vector<std::shared_ptr<Dvec>> tmp_vec;
        for (int i = 0; i < norb; ++i) {
          tmp_vec.push_back(second.compute(intermediate->data(i)));
        }

        std::shared_ptr<Dvec> out(new Dvec(tmp_vec.front()->det(), norb*norb));
        int ij = 0;
        for (int i = 0; i < norb; ++i) {
          for (int j = 0; j < norb; ++j, ++ij) {
            out->data(ij) = tmp_vec.at(j)->data(i);
          }
        }

        return out;
      }
    }
};

template<SQ oper1, SQ oper2, SQ oper3> class ThreeBody : public Quantization {
  public:
    ThreeBody() {}
  
    const int ij(const int norb) const override { return norb*norb*norb; }
    std::shared_ptr<Dvec> compute(std::shared_ptr<const Civec> ccvec) const override {
      if ( ((static_cast<int>(oper2) + static_cast<int>(oper3)) == 0) && (static_cast<int>(oper3) < 0) ) { // Can start with an optimized TwoBody
        TwoBody<oper2,oper3> two;
        OneBody<oper1> one;

        std::shared_ptr<Dvec> tmp = two.compute(ccvec);
        const int sizejk = tmp->ij();

        std::vector<std::shared_ptr<Dvec>> tmp_vec;
        for (int i = 0; i < sizejk; ++i) {
          tmp_vec.push_back(one.compute(tmp->data(i)));
        }
        const int sizei = tmp_vec.front()->ij();

        std::shared_ptr<Dvec> out(new Dvec(tmp_vec.front()->det(), sizei*sizejk));
        for (int i = 0; i < sizei; ++i) {
          for (int jk = 0; jk < sizejk; ++jk) {
            out->data(i * sizejk + jk) = tmp_vec.at(jk)->data(i);
          }
        }

        return out;
      }
      else {
        OneBody<oper3> one;
        TwoBody<oper1,oper2> two;

        std::shared_ptr<Dvec> tmp = one.compute(ccvec);
        const int sizek = tmp->ij();

        std::vector<std::shared_ptr<Dvec>> tmp_vec;
        for (int k = 0; k < sizek; ++k) {
          tmp_vec.push_back(two.compute(tmp->data(k)));
        }
        const int sizeij = tmp_vec.front()->ij();

        std::shared_ptr<Dvec> out(new Dvec(tmp_vec.front()->det(), sizeij*sizek));
        for (int ij = 0; ij < sizeij; ++ij) {
          for (int k = 0; k < sizek; ++k) {
            out->data(ij * sizek + k) = tmp_vec.at(k)->data(ij);
          }
        }

        return out;
      }
    }
};

}

#endif
