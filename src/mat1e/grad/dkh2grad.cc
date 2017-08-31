//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dkh2grad.cc
// Copyright (C) 2017 Toru Shiozaki
//
// Author: Nils Strand <nilsstrand2022@u.northwestern.edu>
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


#include <src/mat1e/grad/dkh2grad.h>
#include <src/mat1e/kinetic.h>
#include <src/mat1e/nai.h>
#include <src/mat1e/overlap.h>
#include <src/mat1e/rel/small1e.h>
#include <src/util/math/algo.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(DKH2grad)


DKH2grad::DKH2grad(shared_ptr<const Geometry> mol) : Relgrad_base(mol->nbasis(), mol->nbasis()) {
  cout << "       - Using DKH2grad" << endl;
  init();
}


void DKH2grad::init() {

  dkh2grad = vector<vector<Matrix>>(3, vector<Matrix>(natom));
  shared_ptr<Overlap> s = make_shared<Overlap>(molu);
  contracts(s);
  overlapgrad(s);
  const Matrix s_inv12 = *s->tildex();
  Matrix s_inv = *s;
  s_inv.inverse_symmetric();
  const Kinetic T_p(molu);
  kineticgrad(T_p);
  // const Matrix T = U_T % T_p * U_T;
  const Matrix T_pp = s_inv12 % T_p * s_inv12;
  shared_ptr<VectorB> t = make_shared<VectorB>(nunc);
  store_mat(t);
  Matrix W = T_pp;
  W.diagonalize(*t);

  const double c2 = c__ * c__;
  shared_ptr<VectorB> Ep = make_shared<VectorB>(nunc);
  store_mat(Ep);
  shared_ptr<VectorB> A = make_shared<VectorB>(nunc);
  store_mat(A);
  shared_ptr<VectorB> K = make_shared<VectorB>(nunc);
  store_mat(K);
  shared_ptr<VectorB> B = make_shared<VectorB>(nunc);
  store_mat(B);
  shared_ptr<VectorB> R = make_shared<VectorB>(nunc);
  store_mat(R);
  shared_ptr<VectorB> R_inv = make_shared<VectorB>(nunc);
  store_mat(R_inv);
  for (int k = 0; k < nunc; ++k) {
    (*Ep)(k) = c__ * std::sqrt(2 * (*t)(k) + c2);
    (*A)(k) = std::sqrt((c2 + (*Ep)(k)) / (2 * (*Ep)(k)));
    (*K)(k) = c__ / ((*Ep)(k) + c2);
    (*B)(k) = (*A)(k) * (*K)(k);
    (*R)(k) = 2 * (*t)(k) * pow((*K)(k), 2);
    (*R_inv)(k) = 1 / (*R)(k);
  }

  for (int i = 0; i < 3; i++) {

    for (int j = 0; j < natom) {

      const Matrix T_pX = U_T * T_X[i][j] ^ U_T - PU[i][j] * T_p + T_p * PU[i][j];
      const Matrix s_X = U_T * S_X[i][j] ^ U_T - PU[i][j] * s + s * PU[i][j];
      const Matrix T_ppX = s_inv12 % T_pX * s_inv12 - 0.5 * (s_inv * s_X * T_pp + T_pp * s_inv * s_X);
      const Matrix T_ppX_W = W % T_ppX * W;
      Matrix PW(nunc, nunc);
      for (int k = 0; k < nunc; ++k) {
        for (int l = 0; l < nunc; ++l) {
          PW(k, l) = k == l ? 0 : T_ppX_W(k, l) / ((*t)(k) - (*t)(l));
        }
      }
      const Matrix t_X = T_ppX_W - PW * *vec2mat[t] + *vec2mat[t] * PW;

      shared_ptr<VectorB> Ep_X = make_shared<VectorB>(nunc);
      store_mat(Ep_X);
      shared_ptr<VectorB> A_X = make_shared<VectorB>(nunc);
      store_mat(A_X);
      shared_ptr<VectorB> B_X = make_shared<VectorB>(nunc);
      store_mat(B_X);
      shared_ptr<VectorB> R_X = make_shared<VectorB>(nunc);
      store_mat(R_X);
      shared_ptr<VectorB> R_invX = make_shared<VectorB>(nunc);
      store_mat(R_invX);
      for (int k = 0; k < nunc; ++k) {
        (*Ep)(k) = c__ * std::sqrt(2 * (*t)(k) + c2);
        (*A)(k) = std::sqrt((c2 + (*Ep)(k)) / (2 * (*Ep)(k)));
        (*K)(k) = c__ / ((*Ep)(k) + c2);
        (*B)(k) = (*A)(k) * (*K)(k);
        (*Ep_X)(k) = c__ * t_X[i][j](k, k) / std::sqrt(2 * (*t)(k) + c2);
        (*A_X)(k) = c2 * (*Ep_X)(k) / (-4 * pow((*Ep)(k), 2) * (*A)(k));
        (*B_X)(k) = c__ * (*Ep_X)(k) * (c__ * (*K)(k) / (4 * pow((*Ep)(k), 2) * (*A)(k)) - (*A)(k) / pow((*Ep)(k) + c2, 2));
        (*R)(k) = 2 * (*t)(k) * pow((*K)(k), 2);
        (*R_inv)(k) = 1 / (*R)(k);
        (*R_X)(k) = 2 * t_X[i][j](k, k) * (pow((*K)(k), 2) - 2 * c2 / (pow((*Ep)(k) + c2, 2) * std::sqrt(2 * (*t)(k) + c2)));
        (*R_invX)(k) = -1 * (*R_X)(k) / pow((*R)(k), 2);
      }
      
      shared_ptr<VectorB> t_rel = make_shared<VectorB>(nunc);
      store_mat(t_rel);
      for (int k = 0; k < nunc; ++k) {
        (*t_rel)(k) = c__ * std::sqrt(2 * (*t)(k) + c2) - c2;
      }
      const Matrix T_pprel = W * *vec2mat[t_rel] ^ W;
      const Matrix T_prel = s_inv12 * T_pprel ^ s_inv12;
      const Matrix T_rel = *U_T % T_prel * *U_T;
      const Matrix DW = W * PW ^ W;
      const Matrix T_pprelX = W * *vec2mat[Ep] ^ W + DW * T_pprel - T_pprel * DW;
      const Matrix T_prelX = s_inv12 * T_pprelX ^ s_inv12 + 0.5 * (s_inv * s_X * T_prel + T_prel * s_inv * s_X);
      const Matrix DU = *U_T % PU[i][j] * *U_T;
      const Matrix T_relX = *U_T % T_prelX * *U_T + DU * T_rel - T_rel * DU;

      const NAI V_p(molu);
      // const Matrix V = U_T % V_p * U_T;
      const Matrix V_pp = s_inv12 % V_p * s_inv12;
      const Matrix V_ppp = W % V_pp * W;
      const Matrix V_pX = U_T * V_X[i][j] ^ U_T - PU * V_p + V_p * PU;
      const Matrix V_ppX = s_inv12 % V_pX * s_inv12 - 0.5 * (s_inv * s_X * V_pp + V_pp * s_inv * s_X);
      const Matrix V_pppX = W % V_ppX * W - PW * V_ppp + V_ppp * PW;

      const Small1e<NAIBatch> small1e(molu);
      const Matrix O_p = small1e[0];
      // const Matrix O = U_T % O_p * U_T;
      const Matrix O_pp = s_inv12 % O_p * s_inv12;
      const Matrix O_ppp = W % O_pp * W;
      const Matrix O_pX = U_T * O_X[i][j] ^ U_T - PU * O_p + O_p * PU;
      const Matrix O_ppX = s_inv12 % O_pX * s_inv12 - 0.5 * (s_inv * s_X * O_pp + O_pp * s_inv * s_X);
      const Matrix O_pppX = W % O_ppX * W - PW * O_ppp + O_ppp * PW;

      Matrix V_ppprelX = *vec2mat[A] * V_pppX * *vec2mat[A] + *vec2mat[A_X] * V_ppp * *vec2mat[A] + *vec2mat[A] * V_ppp * *vec2mat[A_X]
                        + *vec2mat[B] * O_pppX * *vec2mat[B] + *vec2mat[B_X] * O_ppp * *vec2mat[B] + *vec2mat[B] * O_ppp * *vec2mat[B_X];

      Matrix V_long(nunc, nunc), A_long(nunc, nunc), O_long(nunc, nunc), B_long(nunc, nunc);
      Matrix V_longX(nunc, nunc), A_longX(nunc, nunc), O_longX(nunc, nunc), B_longX(nunc, nunc);
      for (int k = 0; k < nunc; k++) {
        for (int l = 0; l < nunc; l++) {
          V_long(k, l) = V_ppp(k, l) / ((*Ep)(k) + (*Ep)(l));
          A_long(k, l) = (*A)(k) * V_long(k, l) * (*A)(l);
          O_long(k, l) = O_ppp(k, l) / ((*Ep)(k) + (*Ep)(l));
          B_long(k, l) = (*B)(k) * O_long(k, l) * (*B)(l);
          V_longX(k, l) = -1 * pow(*Ep)(k) + (*Ep)(l), -2) * ((*Ep_X)(k) + (*Ep_X)(l)) * V_ppp(k, l) + V_pppX(k, l) / ((*Ep)(k) + (*Ep)(l));
          A_longX(k, l) = (*A)(k) * V_longX(k, l) * (*A)(l) + (*A_X)(k) * V_long(k, l) * (*A)(l) + (*A)(k) * V_long(k, l) * (*A_X)(l);
          O_longX(k, l) = -1 * pow((*Ep)(k) + (*Ep)(l), -2) * ((*Ep_X)(k) + (*Ep_X)(l)) * O_ppp(k, l) + O_pppX(k, l) / ((*Ep)(k) + (*Ep)(l));
          B_longX(k, l) = (*B)(k) * O_longX(k, l) * (*B)(l) + (*B_X)(k) * O_long(k, l) * (*B)(l) + (*B)(k) * O_long(k, l) * (*B_X)(l);
        }
      }

      V_ppprelX += -1 * B_longX * *vec2mat[Ep] * A_long - B_long * *vec2mat[Ep_X] * A_long - B_long * *vec2mat[Ep] * A_longX
                - A_longX *  *vec2mat[Ep] * B_long - A_long * *vec2mat[Ep_X] * B_long - A_long * *vec2mat[Ep] * B_longX
                + A_longX * *vec2mat[Ep % R] * A_long + A_long * *vec2mat[Ep_X % R + Ep % R_X] * A_long + A_long * *vec2mat[Ep % R] * A_longX
                + B_longX * *vec2mat[Ep % R_inv] * B_long + B_long * *vec2mat[Ep_X % R_inv + Ep % R_invX] * B_long + B_long * *vec2mat[Ep % R_inv] * B_longX
                - 0.5 * (B_longX * A_long * *vec2mat[Ep] + B_long * A_longX * *vec2mat[Ep] + B_long * A_long * *vec2mat[Ep_X])
                - 0.5 * (A_longX * B_long * *vec2mat[Ep] + A_long * B_longX * *vec2mat[Ep] + A_long * B_long * *vec2mat[Ep_X])
                + 0.5 * (A_longX * *vec2mat[R] * A_long * *vec2mat[Ep] + A_long * *vec2mat[R_X] * A_long * *vec2mat[Ep])
                + 0.5 * (A_long * *vec2mat[R] * A_longX * *vec2mat[Ep] + A_long * *vec2mat[R] * A_long * *vec2mat[Ep_X])
                + 0.5 * (B_longX * *vec2mat[R_inv] * B_long * *vec2mat[Ep] + B_long * *vec2mat[R_invX] * B_long * *vec2mat[Ep])
                + 0.5 * (B_long * *vec2mat[R_inv] * B_longX * *vec2mat[Ep] + B_long * *vec2mat[R_inv] * B_long * *vec2mat[Ep_X])
                - 0.5 * (*vec2mat[Ep_X] * B_long * A_long + *vec2mat[Ep] * B_longX * A_long + *vec2mat[Ep] * B_long * A_longX)
                - 0.5 * (*vec2mat[Ep_X] * A_long * B_long + *vec2mat[Ep] * A_longX * B_long + *vec2mat[Ep] * A_long * B_longX)
                + 0.5 * (*vec2mat[Ep_X] * A_long * *vec2mat[R] * A_long + *vec2mat[Ep] * A_longX * *vec2mat[R] * A_long)
                + 0.5 * (*vec2mat[Ep] * A_long * *vec2mat[R_X] * A_long + *vec2mat[Ep] * A_long * *vec2mat[R] * A_longX)
                + 0.5 * (*vec2mat[Ep_X] * B_long * *vec2mat[R_inv] * B_long + *vec2mat[Ep] * B_longX * *vec2mat[R_inv] * B_long)
                + 0.5 * (*vec2mat[Ep] * B_long * *vec2mat[R_invX] * B_long + *vec2mat[Ep] * B_long * *vec2mat[R_inv] * B_longX);

      Matrix V_ppprel = *vec2mat[A] * V_ppp * *vec2mat[A] + *vec2mat[B] * O_ppp * *vec2mat[B]
                      - B_long * *vec2mat[Ep] * A_long - A_long * *vec2mat[Ep] * B_long
                      + A_long * *vec2mat[R % Ep] * A_long + B_long * *vec2mat[R_inv % Ep] * B_long
                      - 0.5 * (B_long * A_long * *vec2mat[Ep] + A_long * B_long * *vec2mat[Ep])
                      + 0.5 * (A_long * *vec2mat[R] * A_long * *vec2mat[Ep] + B_long * *vec2mat[R_inv] * B_long * *vec2mat[Ep])
                      - 0.5 * (*vec2mat[Ep] * B_long * A_long + *vec2mat[Ep] * A_long * B_long)
                      + 0.5 * (*vec2mat[Ep] * A_long * *vec2mat[R] * A_long + *vec2mat[Ep] * B_long * *vec2mat[R_inv] * B_long);
      
      const Matrix V_pprel = W * V_ppprel ^ W;
      const Matrix V_prel = s_inv12 * V_pprel ^ s_inv12;
      const Matrix V_rel = *U_T % V_prel * *U_T;
      const Matrix V_pprelX = W * V_ppprelX ^ W + DW * V_pprel - V_pprel * DW;
      const Matrix V_prelX = s_inv12 * V_pprelX ^ s_inv12 + 0.5 * (s_inv * s_X * V_prel + V_prel * s_inv * s_X);
      const Matrix V_relX = *U_T % V_prelX * *U_T + DU * V_rel - V_rel * DU;

      dkh2grad[i][j] = T_relX + V_relX;

    }

  }

  // Matrix_base<double>::operator=(T_relX + V_relX);
}

shared_ptr<vector<vector<Matrix>>> DKH2grad::result() const {
  return make_shared<vector<vector<Matrix>>>(dkh2grad);
}

