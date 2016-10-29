//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dkhcore.cc
// Copyright (C) 2016 Toru Shiozaki
//
// Author: Raymond Wang <yiqunwang2021@u.northwestern.edu> 
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

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(DKHcore)

DKHcore::DKHcore(shared_ptr<const Molecule> mol) : Matrix(mol->nbasis(), mol->nbasis()) {

  init(mol);
  cout<<endl;
  cout<<"  === Using DKHcore === "<<endl;
  cout<<endl;
}


void DKHcore::init(shared_ptr<const Molecule> mol) {
  
  Kinetic transfer(mol); // we do not use kinetic matrix directly, transfer here is kinetic.
  VectorB eig(mol->nbasis());
  VectorB en(mol->nbasis());
  transfer.diagonalize(eig);
  NAI nai(mol);

  for (int i = 0; i < en.size(); ++i){
    en(i) = c__ * std::sqrt(2*eig(i)+c__*c__);
  }

  Matrix A(ndim(),mdim(),0.0);
  for (int i=0; i<ndim(); i++){
    A(i,i) = std::sqrt(0.5*(en(i)+c__*c__)/en(i));
  }

  Matrix P(ndim(),mdim(),0.0);
  for (int i=0; i<ndim(); i++){
    P(i,i) = std::sqrt(2*eig(i));
  }
  
  Matrix pK2(ndim(),mdim(),0.0);
  for (int i=0; i<ndim(); i++){
    pK2(i,i) = pow((P(i,i)*c__/(en(i)+c__*c__)),2);
  }

  Matrix pK2_inv(ndim(),mdim(),0.0);
  for (int i=0; i<ndim(); i++){
    pK2_inv(i,i) = 1 / pK2(i,i);
  }

  Matrix E(ndim(),mdim(),0.0);
  for (int i=0; i<ndim(); i++){
    E(i,i) = en(i);
  }

  Matrix B(ndim(),mdim(),0.0);
  for (int i=0; i<ndim(); i++){
    B(i,i) = A(i,i) * c__ / (en(i)+c__*c__);
  }

  Matrix V(ndim(),mdim(),0.0);
  V = transfer % nai * transfer;

  Matrix Vp(ndim(),mdim(),0.0);
  for(int i=0; i<ndim(); i++){
    for(int j=0; j<mdim(); j++){
      Vp(i,j) = V(i,j) / (en(i)+en(j));
    }
  }
  
  Matrix dkh(ndim(),mdim(),0.0);
  for (int i=0; i<ndim(); i++){
    dkh(i,i) = en(i) - c__*c__;
  }

  dkh += A*V*A + B*P*V*P*B;

  dkh += (-1)*B*P*Vp*P*B*E*A*Vp*A - A*Vp*A*E*B*P*Vp*P*B + A*Vp*A*pK2*E*A*Vp*A + B*P*Vp*P*B*E*pK2_inv*B*P*Vp*P*B                //NWChem
          -0.5*B*P*Vp*P*B*A*Vp*A*E -0.5*A*Vp*A*B*P*Vp*P*B*E + 0.5*A*Vp*A*pK2*A*Vp*A*E + 0.5*B*P*Vp*P*B*pK2_inv*B*P*Vp*P*B*E
          -0.5*E*B*P*Vp*P*B*A*Vp*A -0.5*E*A*Vp*A*B*P*Vp*P*B + 0.5*E*A*Vp*A*pK2*A*Vp*A + 0.5*E*B*P*Vp*P*B*pK2_inv*B*P*Vp*P*B;

//  dkh += 0.5 * ((-1)*A*P*tildeV*P*A*A*V*A + A*P*tildeV*A*A*V*P*A + A*tildeV*A*P2*A*V*A - A*tildeV*A*A*P*V*P*A               //Reiher 
//                 - A*P*V*P*A*A*tildeV*A + A*P*V*A*A*tildeV*P*A + A*V*A*P2*A*tildeV*A - A*V*A*A*P*tildeV*P*A);

  dkh = transfer * dkh ^ transfer;
  copy_block(0,0,ndim(),mdim(),dkh);
}
