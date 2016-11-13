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

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(DKHcore)

DKHcore::DKHcore(shared_ptr<const Molecule> mol) : Matrix(mol->nbasis(), mol->nbasis()) {
  
  init(mol);
  cout<<endl;
  cout<<"  === Using DKHcore === "<<endl;
  cout<<endl;
}


void DKHcore::init(shared_ptr<const Molecule> mol0) {
  
  auto mol = make_shared<Molecule>(*mol0);
  mol = mol->uncontract();
cout<<"mol0->nbasis: "<<mol0->nbasis()<<endl;
cout<<"mol->nbasis: "<<mol->nbasis()<<endl;

  Kinetic kinetic(mol);
  Overlap overlap(mol);
  shared_ptr<const Matrix> tildex = overlap.tildex();
  Matrix transfer(ndim(), mdim(), 0.0);
  transfer = *tildex % kinetic * *tildex;
  VectorB eig(mol->nbasis());
  VectorB ep_vec(mol->nbasis());
  transfer.diagonalize(eig);
  transfer = *tildex * transfer;

  NAI nai(mol);
  Small1e<NAIBatch> small1e(mol);
  
  for (int i = 0; i < eig.size(); ++i) { // Ep vector
    ep_vec(i) = c__ * std::sqrt(2 * eig(i) + c__ * c__);
  }

  Matrix Ep(ndim(), mdim(), 0.0); // Ep matrix
  for (int i = 0; i != ndim(); ++i) {
    Ep(i,i) = ep_vec(i);
  }
  
  Matrix A(ndim(), mdim(), 0.0);
  for (int i = 0; i != ndim(); ++i) {
    A(i,i) = std::sqrt(0.5 * (ep_vec(i) + c__ * c__) / ep_vec(i));
  }
  
  Matrix V(ndim(), mdim(), 0.0); // V in p^2 basis
  V = transfer % nai * transfer;
  
  Matrix tildeV(ndim(), mdim(), 0.0);
  for (int i = 0; i != ndim(); ++i) {
    for (int j = 0; j != mdim(); ++j){
      tildeV(i,j) = V(i,j) / (ep_vec(i) + ep_vec(j));
    }
  }
  
  Matrix AVA(ndim(), mdim(), 0.0); // A * tildeV * A 
  AVA = A * tildeV * A;

  Matrix B(ndim(), mdim(), 0.0);
  for (int i = 0; i != ndim(); ++i) {
    B(i,i) = A(i,i) * c__ / (ep_vec(i) + c__ * c__);
  }

  Matrix smallnai(ndim(), mdim()); // pVp in p^2 basis
  smallnai.copy_block(0, 0, ndim(), mdim(), small1e[0]);
  smallnai = transfer % smallnai * transfer;
  
  Matrix BVB(ndim(), mdim(), 0.0); // A * R * tildeV * R * A
  for(int i = 0; i != ndim(); ++i) {
    for(int j = 0; j != mdim(); ++j) {
      BVB(i,j) = smallnai(i,j) / (ep_vec(i) + ep_vec(j));
    }
  }
  BVB = B * BVB * B;

  Matrix RI(ndim(), mdim(), 0.0); // used for RI terms
  for (int i = 0; i != ndim(); ++i) {
    RI(i,i) = 2 * eig(i) * pow((c__ / (ep_vec(i) + c__ * c__)), 2);
  }
 
  Matrix RI_inv(RI);
  RI_inv.inverse();

  Matrix smallx(ndim(), mdim()); // x component of (p x Vp)
  smallx.copy_block(0, 0, ndim(), mdim(), small1e[2]);
  smallx = transfer % smallx * transfer;
  for (int i = 0; i != ndim(); ++i) {
    for (int j = 0; j != mdim(); ++j) {
      smallx(i,j) /= ep_vec(i) + ep_vec(j);
    }
  }

  Matrix smally(ndim(), mdim()); // y component ...
  smally.copy_block(0, 0, ndim(), mdim(), small1e[3]);
  smally = transfer % smally * transfer;
  for (int i = 0; i != ndim(); ++i) {
    for (int j = 0; j != mdim(); ++j) {
      smally(i,j) /= ep_vec(i) + ep_vec(j);
    }
  }
  
  Matrix smallz(ndim(), mdim()); // z component ...
  smallz.copy_block(0, 0, ndim(), mdim(), small1e[1]);
  smallz = transfer % smallz * transfer;
  for (int i = 0; i != ndim(); ++i) {
    for (int j = 0; j != mdim(); ++j) {
      smallz(i,j) /= ep_vec(i) + ep_vec(j);
    }
  }

  Matrix Xc(ndim(), mdim(), 0.0);
  Xc = B * smallx * B;
  
  Matrix Yc(ndim(), mdim(), 0.0);
  Yc = B * smally * B;
  
  Matrix Zc(ndim(), mdim(), 0.0);
  Zc = B * smallz * B;
   
  Matrix dkh(ndim(), mdim(), 0.0);
  for (int i = 0; i != ndim(); ++i) {
    dkh(i,i) = ep_vec(i) - c__ * c__;  // DKH0
  }
  
  dkh += A * V * A                    + B * smallnai * B  // DKH1
       - BVB * Ep * AVA               - AVA * Ep * BVB                 + AVA * RI * Ep * AVA
       + BVB * Ep * RI_inv * BVB      - 0.5 * BVB * AVA * Ep           - 0.5 * AVA * BVB * Ep
       + 0.5 * AVA * RI * AVA * Ep    + 0.5 * BVB * RI_inv * BVB * Ep  - 0.5 * Ep * BVB * AVA
       - 0.5 * Ep * AVA * BVB         + 0.5 * Ep * AVA * RI * AVA      + 0.5 * Ep * BVB * RI_inv * BVB 
       - Xc * Ep * RI_inv * Xc        - Yc * Ep * RI_inv * Yc          - Zc * Ep * RI_inv * Zc
       - 0.5 * Xc * RI_inv * Xc * Ep  - 0.5 * Yc * RI_inv * Yc * Ep    - 0.5 * Zc * RI_inv * Zc * Ep
       - 0.5 * Ep * Xc * RI_inv * Xc  - 0.5 * Ep * Yc * RI_inv * Yc    - 0.5 * Ep * Zc * RI_inv * Zc; // DKH2

  
  transfer = overlap * transfer;

  dkh = transfer * dkh ^ transfer;
  
  copy_block(0,0,ndim(),mdim(),dkh);
}
