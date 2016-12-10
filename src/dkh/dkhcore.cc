//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dkhcore.cc
// Copyright (C) 2016 Toru Shiozaki
//
// Author: Raymond Wang <raymondwang@u.northwestern.edu> 
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


#include <src/dkh/dkhcore.h>
#include <src/mat1e/kinetic.h>
#include <src/mat1e/nai.h>
#include <src/mat1e/rel/small1e.h>
#include <src/mat1e/overlap.h>
#include <src/integral/os/overlapbatch.h>
#include <src/mat1e/mixedbasis.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(DKHcore_<double>)
BOOST_CLASS_EXPORT_IMPLEMENT(DKHcore_<std::complex<double>>)


template<typename DataType>
DKHcore_<DataType>::DKHcore_(shared_ptr<const Molecule> mol) : MatType(mol->nbasis(), mol->nbasis()) {
  
  cout << "       - Using DKHcore" << endl;
  init(mol);
}

inline Matrix pre_scale(const VectorB& vec, const Matrix& mat) {
//  assert(mat.ndim() == mat.mdim() && mat.ndim() == vec.size());
  Matrix out(mat);
  for (size_t i = 0; i != mat.ndim(); ++i)
    dscal_(mat.ndim(),vec(i),out.element_ptr(i,0),mat.ndim());
  return out;
}

inline Matrix post_scale(const Matrix& mat, const VectorB& vec) {
//  assert(mat.ndim() == mat.mdim() && mat.ndim() == vec.size());
  Matrix out(mat);
  for (size_t i = 0; i != mat.ndim(); ++i) {
    dscal_(mat.ndim(),vec(i),out.element_ptr(0,i),1);
  } 
  return out;
}

template<>
void DKHcore_<double>::init(shared_ptr<const Molecule> mol0) {

  auto mol = make_shared<Molecule>(*mol0);
  mol = mol->uncontract();

  // build DKH Hamiltonian in uncontracted basis
  const Kinetic kinetic(mol); // TODO assert declaration
  const Overlap overlap(mol);
  shared_ptr<const Matrix> tildex = overlap.tildex(1.0e-9);
  assert(kinetic.ndim() == tildex->ndim()); // if it fails, check linear dependency
  Matrix transfer(*tildex % kinetic * *tildex);
  const size_t ndim = transfer.ndim(); 
  VectorB eig(ndim);
  VectorB Ep(ndim);
  transfer.diagonalize(eig);
  transfer = *tildex * transfer;
  const NAI nai(mol);
  const Small1e<NAIBatch> small1e(mol);
  const double c2 = c__ * c__;

  for (size_t i = 0; i != eig.size(); ++i) { 
    Ep(i) = c__ * std::sqrt(2.0 * eig(i) + c2);
  }

  VectorB A(ndim);
  for (size_t i = 0; i != ndim; ++i) {
    A(i) = std::sqrt(0.5 * (Ep(i) + c2) / Ep(i));
  }
  
  Matrix V(transfer % nai * transfer);

  Matrix AVA(V);
  for (size_t i = 0; i != ndim; ++i) {
    for (size_t j = 0; j != ndim; ++j){
      AVA(j,i) = AVA(j,i) *A(j) * A(i) / (Ep(j) + Ep(i));
    }
  }

  VectorB B(ndim);
  for (size_t i = 0; i != ndim; ++i) {
    B(i) = A(i) * c__ / (Ep(i) + c2);
  }

  const Matrix smallnai(transfer % small1e[0] * transfer);

  Matrix BVB(smallnai);
  for(size_t i = 0; i != ndim; ++i) {
    for(size_t j = 0; j != ndim; ++j) {
      BVB(j,i) = BVB(j,i) * B(j) * B(i) / (Ep(j) + Ep(i));
    }
  }

  VectorB RI(ndim);
  VectorB RI_inv(ndim);
  for (size_t i = 0; i != ndim; ++i) {
    RI(i) = 2.0 * eig(i) * pow((c__ / (Ep(i) + c2)), 2);
    RI_inv(i) = 1.0 / RI(i);
  }
 
  Matrix dkh(ndim, ndim);
  for (size_t i = 0; i != ndim; ++i) {
    dkh(i,i) = c2 * 2.0 * eig(i) / (Ep(i) + c2);
  }

  Matrix EAVA(pre_scale(Ep,AVA));
  Matrix AVARI(post_scale(AVA,RI));
  Matrix AVAE(post_scale(AVA,Ep));
  Matrix BVBE(post_scale(BVB,Ep));
  Matrix EBVB(pre_scale(Ep,BVB));
  
  dkh += post_scale(pre_scale(A,V),A)                + post_scale(pre_scale(B,smallnai),B) // Free Particle Projection
       - BVB * (EAVA - 0.5 * pre_scale(RI_inv,BVBE)) - (AVAE + 0.5 * EAVA) * BVB  
       + 0.5 * (post_scale(EAVA,RI) - EBVB) * AVA    + 0.5 * EBVB * pre_scale(RI_inv,BVB)
       + AVARI * (EAVA + 0.5 * AVAE)                 + BVBE * pre_scale(RI_inv,BVB)      
       - 0.5 * BVB * AVAE                            - 0.5 * AVA * BVBE;                   // DKH2
  
/*  dkh += post_scale(pre_scale(A,V),A) + post_scale(pre_scale(B,smallnai),B) 
       - BVB * EAVA                     - AVAE * BVB                         + AVARI * EAVA
       + BVBE * pre_scale(RI_inv,BVB)   - 0.5 * BVB * AVAE                   - 0.5 * AVA * BVBE
       + 0.5 * AVARI * AVAE             + 0.5 * BVB * pre_scale(RI_inv,BVBE) - 0.5 * EBVB * AVA
       - 0.5 * EAVA * BVB               + 0.5 * EAVA * pre_scale(RI,AVA)     + 0.5 * EBVB * pre_scale(RI_inv,BVB); */ 

  const MixedBasis<OverlapBatch> mix(mol0, mol);

  transfer = transfer % mix;
  
  dkh = transfer % dkh * transfer;

  copy_block(0,0,mol0->nbasis(),mol0->nbasis(),dkh);
}


template class DKHcore_<double>;
template class DKHcore_<std::complex<double>>;


// Leave it here if we want do DKH2FULL in the future
/*  
  Matrix smallx(transfer % small1e[2] * transfer);
  Matrix smally(transfer % small1e[3] * transfer);
  Matrix smallz(transfer % small1e[1] * transfer);
  for (size_t i = 0; i != ndim; ++i) {
    for (size_t j = 0; j != ndim; ++j) {
      const size_t denom = Ep(i) + Ep(j);
      smallx(i,j) /= denom; 
      smally(i,j) /= denom;
      smallz(i,j) /= denom;
    }
  }

  const Matrix Xc(B * smallx * B);
  const Matrix Yc(B * smally * B);
  const Matrix Zc(B * smallz * B);
  
     - Xc * Ep * RI_inv * Xc        - Yc * Ep * RI_inv * Yc          - Zc * Ep * RI_inv * Zc
     - 0.5 * Xc * RI_inv * Xc * Ep  - 0.5 * Yc * RI_inv * Yc * Ep    - 0.5 * Zc * RI_inv * Zc * Ep
     - 0.5 * Ep * Xc * RI_inv * Xc  - 0.5 * Ep * Yc * RI_inv * Yc    - 0.5 * Ep * Zc * RI_inv * Zc; // DKH2FULL       
*/
