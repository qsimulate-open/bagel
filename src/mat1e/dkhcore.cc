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


#include <src/mat1e/dkhcore.h>
#include <src/mat1e/kinetic.h>
#include <src/mat1e/nai.h>
#include <src/mat1e/rel/small1e.h>
#include <src/mat1e/overlap.h>
#include <src/integral/os/overlapbatch.h>
#include <src/mat1e/mixedbasis.h>
#include <src/util/math/algo.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(DKHcore)


DKHcore::DKHcore(shared_ptr<const Molecule> mol) : Matrix(mol->nbasis(), mol->nbasis()) {
  init(mol);
}


void DKHcore::init(shared_ptr<const Molecule> mol0) {

  auto pre_scale = [](const VectorB& vec, const Matrix& mat) {
    assert(mat.ndim() == vec.size());
    shared_ptr<Matrix> out = mat.clone();
    for (int i = 0; i != mat.mdim(); ++i)
      for (int j = 0; j != mat.ndim(); ++j)
        (*out)(j, i) = mat(j, i) * vec(j);
    return out;
  };

  auto post_scale = [](const Matrix& mat, const VectorB& vec) {
    assert(mat.mdim() == vec.size());
    shared_ptr<Matrix> out = mat.copy();
    for (int i = 0; i != mat.mdim(); ++i)
      blas::scale_n(vec(i), out->element_ptr(0, i), mat.ndim());
    return out;
  };

  auto mol = make_shared<Molecule>(*mol0);
  mol = mol->uncontract();

  // build DKH Hamiltonian in uncontracted basis
  shared_ptr<const Matrix> transfer;
  VectorB eig;
  {
    const Overlap overlap(mol);
    shared_ptr<const Matrix> tildex = overlap.tildex();
    const Kinetic kinetic(mol);
    auto tmp = make_shared<Matrix>(*tildex % kinetic * *tildex);
    eig = VectorB(tmp->ndim());
    tmp->diagonalize(eig);
    transfer = make_shared<Matrix>(*tildex * *tmp);
  }

  const double c2 = c__ * c__;
  const int nunc = eig.size();
  VectorB Ep(nunc);
  VectorB A(nunc);
  VectorB B(nunc);
  VectorB RI(nunc);
  VectorB RI_inv(nunc);
  for (int i = 0; i != nunc; ++i) {
    Ep(i) = c__ * std::sqrt(2.0 * eig(i) + c2);
    A(i) = std::sqrt(0.5 * (Ep(i) + c2) / Ep(i));
    B(i) = A(i) * c__ / (Ep(i) + c2);
    RI(i) = 2.0 * eig(i) * pow((c__ / (Ep(i) + c2)), 2);
    RI_inv(i) = 1.0 / RI(i);
  }

  const NAI nai(mol);
  const Matrix V(*transfer % nai * *transfer);
  const Small1e<NAIBatch> small1e(mol);
  const Matrix smallnai(*transfer % small1e[0] * *transfer);

  Matrix AVA(nunc, nunc);
  Matrix BVB(nunc, nunc);
  for (int i = 0; i != nunc; ++i)
    for (int j = 0; j != nunc; ++j) {
      const double denom = 1.0 / (Ep(j) + Ep(i));
      AVA(j, i) = V(j, i) * A(j) * A(i) * denom ;
      BVB(j, i) = smallnai(j, i) * B(j) * B(i) * denom;
    }

  Matrix dkh(nunc, nunc);
  for (int i = 0; i != nunc; ++i)
    dkh(i, i) = c2 * 2.0 * eig(i) / (Ep(i) + c2);

  const Matrix EAVA(*pre_scale(Ep, AVA));
  const Matrix AVARI(*post_scale(AVA, RI));
  const Matrix AVAE(*post_scale(AVA, Ep));
  const Matrix BVBE(*post_scale(BVB, Ep));
  const Matrix EBVB(*pre_scale(Ep, BVB));

  dkh += *post_scale(*pre_scale(A, V), A)           + *post_scale(*pre_scale(B, smallnai), B) // Free Particle Projection
       - BVB * (EAVA - 0.5 * (*pre_scale(RI_inv, BVBE)- AVAE))
       - (AVAE + 0.5 * EAVA) * BVB                  + 0.5 * (*post_scale(EAVA, RI) - EBVB) * AVA
       + 0.5 * EBVB * *pre_scale(RI_inv, BVB)       + AVARI * (EAVA + 0.5 * AVAE)
       + BVBE * *pre_scale(RI_inv, BVB)             - 0.5 * AVA * BVBE;                   // DKH2

  const MixedBasis<OverlapBatch> mix(mol0, mol);
  const Matrix transfer2 = *transfer % mix;
  Matrix_base<double>::operator=(transfer2 % dkh * transfer2);
}


/*
   dkh += post_scale(pre_scale(A,V),A) + post_scale(pre_scale(B,smallnai),B)
       - BVB * EAVA                     - AVAE * BVB                         + AVARI * EAVA
       + BVBE * pre_scale(RI_inv,BVB)   - 0.5 * BVB * AVAE                   - 0.5 * AVA * BVBE
       + 0.5 * AVARI * AVAE             + 0.5 * BVB * pre_scale(RI_inv,BVBE) - 0.5 * EBVB * AVA
       - 0.5 * EAVA * BVB               + 0.5 * EAVA * pre_scale(RI,AVA)     + 0.5 * EBVB * pre_scale(RI_inv,BVB);
*/

// Leave it here if we want do DKH2FULL in the future
/*
  Matrix smallx(transfer % small1e[2] * transfer);
  Matrix smally(transfer % small1e[3] * transfer);
  Matrix smallz(transfer % small1e[1] * transfer);
  for (int i = 0; i != ndim; ++i) {
    for (int j = 0; j != ndim; ++j) {
      const int denom = Ep(i) + Ep(j);
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
