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
#include <src/mat1e/rel/small1e.h>
#include <src/mat1e/overlap.h>
#include <src/integral/os/overlapbatch.h>
#include <src/mat1e/mixedbasis.h>
#include <src/util/math/algo.h>
#include <src/mat1e/grad/goverlap.h>
#include <src/mat1e/grad/gkinetic.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(DKH2grad)


DKH2grad::DKH2grad(shared_ptr<const Molecule> mol) : Matrix(mol->nbasis(), mol->nbasis()) {
  cout << "       - Using DKH2grad" << endl;
  init(mol);
}


void DKH2grad::init(shared_ptr<const Molecule> mol) {

  auto vec2mat = [](const VectorB &vec) {
    Matrix out(vec.size(), vec.size());
    for (int i = 0; i < vec.size(); ++i) {
      for (int j = i; j < vec.size(); ++j) {
        out(i, j) = i == j ? vec(i) : 0;
      }
    }
    return out;
  };

  auto naigrad = [](shared_ptr<const Molecule> mol) {
    shared_ptr<Matrix> V_X = make_shared<Matrix>(Matrix(1024, 1024));
    return V_X;
  };

  auto smallnaigrad = [](shared_ptr<const Molecule> mol) {
    shared_ptr<Matrix> O_X = make_shared<Matrix>(Matrix(1024, 1024));
    return O_X;
  };

  shared_ptr<const Molecule> molu = make_shared<Molecule>(*mol);
  molu = molu->uncontract();
  const MixedBasis<OverlapBatch> U_T(mol, molu);
  const Overlap s(molu);
  const Matrix s_inv12 = *s.tildex();
  Matrix s_inv = s;
  s_inv.inverse_symmetric();
  const Kinetic T_p(molu);
  // const Matrix T = U_T % T_p * U_T;
  const Matrix T_pp = s_inv12 % T_p * s_inv12;
  const int nunc = T_pp.ndim();
  VectorB t(nunc);
  Matrix W = T_pp;
  W.diagonalize(t);

  const Matrix T_pX = GKinetic(molu);
  const Matrix S_X = GOverlap(molu);
  const Matrix V_pX = *naigrad(molu);
  const Matrix O_pX = *smallnaigrad(molu);

  const Matrix S_X_U = U_T * S_X ^ U_T;
  Matrix PU(nunc, nunc);
  for (int i = 0; i < nunc; ++i) {
    for (int j = i; j < nunc; ++j) {
      PU(i, j) = i == j ? 0 : S_X_U(i, j) / (s(j, j) - s(i, i));
    }
  }
  // const Matrix T_pX = U_T * T_X ^ U_T - PU * T_p + T_p * PU;

  const Matrix s_X = S_X_U - PU * s + s * PU;
  const Matrix T_ppX = s_inv12 % T_pX * s_inv12 - 0.5 * (s_inv * s_X * T_pp + T_pp * s_inv * s_X);

  const Matrix T_ppX_W = W % T_ppX * W;
  Matrix PW(nunc, nunc);
  for (int i = 0; i < nunc; ++i) {
    for (int j = i; j < nunc; ++j) {
      PW(i, j) = i == j ? 0 : T_ppX_W(i, j) / (t(j) - t(i));
    }
  }
  const Matrix t_X = T_ppX_W - PW * vec2mat(t) + vec2mat(t) * PW;

  const double c2 = c__ * c__;
  VectorB Ep(nunc), A(nunc), K(nunc), B(nunc), Ep_X(nunc), A_X(nunc), K_X(nunc), B_X(nunc), R(nunc), R_inv(nunc), R_X(nunc), R_invX(nunc);
  for (int i = 0; i != nunc; ++i) {
    Ep(i) = c__ * std::sqrt(2 * t(i) + c2);
    A(i) = std::sqrt((c2 + Ep(i)) / (2 * Ep(i)));
    K(i) = c__ / (Ep(i) + c2);
    B(i) = A(i) * K(i);
    Ep_X(i) = c__ * t_X(i, i) / std::sqrt(2 * t(i) + c2);
    A_X(i) = c2 * Ep_X(i) / (-4 * pow(Ep(i), 2) * A(i));
    K_X(i) = -1 * c__ * Ep_X(i) / pow(Ep(i) + c2, 2);
    B_X(i) = c__ * Ep_X(i) * (c__ * K(i) / (4 * pow(Ep(i), 2) * A(i)) - A(i) / pow(Ep(i) + c2, 2));
    R(i) = 2 * t(i) * pow(K(i), 2);
    R_inv(i) = 1 / R(i);
    R_X(i) = 2 * t_X(i, i) * (pow(K(i), 2) - 2 * c2 / (pow(Ep(i) + c2, 2) * std::sqrt(2 * t(i) + c2)));
    R_invX(i) = -1 * R_X(i) / pow(R(i), 2);
  }
  
  VectorB t_rel(nunc);
  for (int i = 0; i < nunc; ++i) {
    t_rel(i) = c__ * std::sqrt(2 * t(i) + c2) - c2;
  }
  const Matrix T_pprel = W * vec2mat(t_rel) ^ W;
  const Matrix T_prel = s_inv12 * T_pprel ^ s_inv12;
  const Matrix T_rel = U_T % T_prel * U_T;
  const Matrix DW = W * PW ^ W;
  const Matrix T_pprelX = W * vec2mat(Ep) ^ W + DW * T_pprel - T_pprel * DW;
  const Matrix T_prelX = s_inv12 * T_pprelX ^ s_inv12 + 0.5 * (s_inv * s_X * T_prel + T_prel * s_inv * s_X);
  const Matrix DU = U_T % PU * U_T;
  const Matrix T_relX = U_T % T_prelX * U_T + DU * T_rel - T_rel * DU;

  const NAI V_p(molu);
  // const Matrix V = U_T % V_p * U_T;
  const Matrix V_pp = s_inv12 % V_p * s_inv12;
  const Matrix V_ppp = W % V_pp * W;
  // const Matrix V_pX = U_T * V_X ^ U_T - PU * V_p + V_p * PU;
  const Matrix V_ppX = s_inv12 % V_pX * s_inv12 - 0.5 * (s_inv * s_X * V_pp + V_pp * s_inv * s_X);
  const Matrix V_pppX = W % V_ppX * W - PW * V_ppp + V_ppp * PW;

  const Small1e<NAIBatch> small1e(molu);
  const Matrix O_p = small1e[0];
  // const Matrix O = U_T % O_p * U_T;
  const Matrix O_pp = s_inv12 % O_p * s_inv12;
  const Matrix O_ppp = W % O_pp * W;
  // const Matrix O_pX = U_T * O_X ^ U_T - PU * O_p + O_p * PU;
  const Matrix O_ppX = s_inv12 % O_pX * s_inv12 - 0.5 * (s_inv * s_X * O_pp + O_pp * s_inv * s_X);
  const Matrix O_pppX = W % O_ppX * W - PW * O_ppp + O_ppp * PW;

  Matrix V_ppprelX = vec2mat(A) * V_pppX * vec2mat(A) + vec2mat(A_X) * V_ppp * vec2mat(A) + vec2mat(A) * V_ppp * vec2mat(A_X)
                    + vec2mat(B) * O_pppX * vec2mat(B) + vec2mat(B_X) * O_ppp * vec2mat(B) + vec2mat(B) * O_ppp * vec2mat(B_X);

  Matrix V_long(nunc, nunc), A_long(nunc, nunc), O_long(nunc, nunc), B_long(nunc, nunc);
  Matrix V_longX(nunc, nunc), A_longX(nunc, nunc), O_longX(nunc, nunc), B_longX(nunc, nunc);
  for (int i = 0; i < nunc; i++) {
    for (int j = i; j < nunc; j++) {
      V_long(i, j) = V_ppp(i, j) / (Ep(i) + Ep(j));
      A_long(i, j) = A(i) * V_long(i, j) * A(j);
      O_long(i, j) = O_ppp(i, j) / (Ep(i) + Ep(j));
      B_long(i, j) = B(i) * O_long(i, j) * B(j);
      V_longX(i, j) = -1 * pow(Ep(i) + Ep(j), -2) * (Ep_X(i) + Ep_X(j)) * V_ppp(i, j) + V_pppX(i, j) / (Ep(i) + Ep(j));
      A_longX(i, j) = A(i) * V_longX(i, j) * A(j) + A_X(i) * V_long(i, j) * A(j) + A(i) * V_long(i, j) * A_X(j);
      O_longX(i, j) = -1 * pow(Ep(i) + Ep(j), -2) * (Ep_X(i) + Ep_X(j)) * O_ppp(i, j) + O_pppX(i, j) / (Ep(i) + Ep(j));
      B_longX(i, j) = B(i) * O_longX(i, j) * B(j) + B_X(i) * O_long(i, j) * B(j) + B(i) * O_long(i, j) * B_X(j);
    }
  }

  V_ppprelX += -1 * B_longX * vec2mat(Ep) * A_long - B_long * vec2mat(Ep_X) * A_long - B_long * vec2mat(Ep) * A_longX
            - A_longX *  vec2mat(Ep) * B_long - A_long * vec2mat(Ep_X) * B_long - A_long * vec2mat(Ep) * B_longX
            + A_longX * vec2mat(Ep % R) * A_long + A_long * vec2mat(Ep_X % R + Ep % R_X) * A_long + A_long * vec2mat(Ep % R) * A_longX
            + B_longX * vec2mat(Ep % R_inv) * B_long + B_long * vec2mat(Ep_X % R_inv + Ep % R_invX) * B_long + B_long * vec2mat(Ep % R_inv) * B_longX
            - 0.5 * (B_longX * A_long * vec2mat(Ep) + B_long * A_longX * vec2mat(Ep) + B_long * A_long * vec2mat(Ep_X))
            - 0.5 * (A_longX * B_long * vec2mat(Ep) + A_long * B_longX * vec2mat(Ep) + A_long * B_long * vec2mat(Ep_X))
            + 0.5 * (A_longX * vec2mat(R) * A_long * vec2mat(Ep) + A_long * vec2mat(R_X) * A_long * vec2mat(Ep))
            + 0.5 * (A_long * vec2mat(R) * A_longX * vec2mat(Ep) + A_long * vec2mat(R) * A_long * vec2mat(Ep_X))
            + 0.5 * (B_longX * vec2mat(R_inv) * B_long * vec2mat(Ep) + B_long * vec2mat(R_invX) * B_long * vec2mat(Ep))
            + 0.5 * (B_long * vec2mat(R_inv) * B_longX * vec2mat(Ep) + B_long * vec2mat(R_inv) * B_long * vec2mat(Ep_X))
            - 0.5 * (vec2mat(Ep_X) * B_long * A_long + vec2mat(Ep) * B_longX * A_long + vec2mat(Ep) * B_long * A_longX)
            - 0.5 * (vec2mat(Ep_X) * A_long * B_long + vec2mat(Ep) * A_longX * B_long + vec2mat(Ep) * A_long * B_longX)
            + 0.5 * (vec2mat(Ep_X) * A_long * vec2mat(R) * A_long + vec2mat(Ep) * A_longX * vec2mat(R) * A_long)
            + 0.5 * (vec2mat(Ep) * A_long * vec2mat(R_X) * A_long + vec2mat(Ep) * A_long * vec2mat(R) * A_longX)
            + 0.5 * (vec2mat(Ep_X) * B_long * vec2mat(R_inv) * B_long + vec2mat(Ep) * B_longX * vec2mat(R_inv) * B_long)
            + 0.5 * (vec2mat(Ep) * B_long * vec2mat(R_invX) * B_long + vec2mat(Ep) * B_long * vec2mat(R_inv) * B_longX);

  Matrix V_ppprel = vec2mat(A) * V_ppp * vec2mat(A) + vec2mat(B) * O_ppp * vec2mat(B)
                  - B_long * vec2mat(Ep) * A_long - A_long * vec2mat(Ep) * B_long
                  + A_long * vec2mat(R % Ep) * A_long + B_long * vec2mat(R_inv % Ep) * B_long
                  - 0.5 * (B_long * A_long * vec2mat(Ep) + A_long * B_long * vec2mat(Ep))
                  + 0.5 * (A_long * vec2mat(R) * A_long * vec2mat(Ep) + B_long * vec2mat(R_inv) * B_long * vec2mat(Ep))
                  - 0.5 * (vec2mat(Ep) * B_long * A_long + vec2mat(Ep) * A_long * B_long)
                  + 0.5 * (vec2mat(Ep) * A_long * vec2mat(R) * A_long + vec2mat(Ep) * B_long * vec2mat(R_inv) * B_long);
  
  const Matrix V_pprel = W * V_ppprel ^ W;
  const Matrix V_prel = s_inv12 * V_pprel ^ s_inv12;
  const Matrix V_rel = U_T % V_prel * U_T;
  const Matrix V_pprelX = W * V_ppprelX ^ W + DW * V_pprel - V_pprel * DW;
  const Matrix V_prelX = s_inv12 * V_pprelX ^ s_inv12 + 0.5 * (s_inv * s_X * V_prel + V_prel * s_inv * s_X);
  const Matrix V_relX = U_T % V_prelX * U_T + DU * V_rel - V_rel * DU;

  Matrix_base<double>::operator=(T_relX + V_relX);
}

