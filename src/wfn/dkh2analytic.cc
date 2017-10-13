//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dkh2analytic.cc
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

#include <src/wfn/dkh2analytic.h>
#include <src/grad/gradeval_base.h>
#include <src/integral/os/overlapbatch.h>
#include <src/mat1e/kinetic.h>
#include <src/mat1e/mixedbasis.h>
#include <src/mat1e/nai.h>
#include <src/mat1e/overlap.h>
#include <src/mat1e/rel/small1e.h>

using namespace std;
using namespace bagel;

vector<shared_ptr<Matrix>> DKH2Analytic::dkh_grad(shared_ptr<const Molecule> current) {
  shared_ptr<const Geometry> mol = dynamic_pointer_cast<const Geometry>(current);
  assert(mol);
  mol = mol->unc_geom();
  gradinit(mol);
  
  vector<shared_ptr<Matrix>> dkh2grad(3 * natom);
  const MixedBasis<OverlapBatch> mix(current, mol);
  contracts(mol);
  overlapgrad(mol);

  DiagVec s_inv12(nunc);
  DiagVec s_inv(nunc);
  for (int k = 0; k < nunc; ++k) {
    s_inv12(k) = 1 / std::sqrt(s(k));
    s_inv(k) = 1 / s(k);
  }

  const Kinetic T(mol);
  auto T_p = make_shared<const Matrix>(U % T * U);
  kineticgrad(mol, T_p);
  Matrix T_pp = s_inv12 * *T_p * s_inv12;
  Matrix W = T_pp;
  VectorB t0(nunc);
  W.diagonalize(t0);
  DiagVec t(t0);

  const double c2 = c__ * c__;
  DiagVec Ep(nunc), A(nunc), K(nunc), B(nunc), R(nunc), R_inv(nunc), EpR(nunc), EpRinv(nunc);
  for (int k = 0; k < nunc; ++k) {
    Ep(k) = c__ * std::sqrt(2 * t(k) + c2);
    A(k) = std::sqrt((c2 + Ep(k)) / (2 * Ep(k)));
    K(k) = c__ / (Ep(k) + c2);
    B(k) = A(k) * K(k);
    R(k) = 2 * t(k) * pow(K(k), 2);
    R_inv(k) = 1 / R(k);
    EpR(k) = Ep(k) * R(k);
    EpRinv(k) = Ep(k) * R_inv(k);
  }

  DiagVec t_rel(nunc);
  for (int k = 0; k < nunc; ++k) {
    t_rel(k) = c__ * std::sqrt(2 * t(k) + c2) - c2;
  }
  const Matrix T_pprel = W * t_rel ^ W;
  const Matrix T_prel = s_inv12 * T_pprel * s_inv12;
  const Matrix T_rel = U * T_prel ^ U;

  const NAI V(mol);
  auto V_p = make_shared<const Matrix>(U % V * U);
  naigrad(mol, V_p);
  const Matrix V_pp = s_inv12 * *V_p * s_inv12;
  const Matrix V_ppp = W % V_pp * W;
  const Small1e<NAIBatch> small1e(mol);
  auto O_p = make_shared<const Matrix>(U % small1e[0] * U);
  smallnaigrad(mol, O_p);
  const Matrix O_pp = s_inv12 * *O_p * s_inv12;
  const Matrix O_ppp = W % O_pp * W;

  Matrix V_long(nunc, nunc), A_long(nunc, nunc), O_long(nunc, nunc), B_long(nunc, nunc);
  for (int k = 0; k < nunc; k++) {
    for (int l = 0; l < nunc; l++) {
      V_long(k, l) = V_ppp(k, l) / (Ep(k) + Ep(l));
      A_long(k, l) = A(k) * V_long(k, l) * A(l);
      O_long(k, l) = O_ppp(k, l) / (Ep(k) + Ep(l));
      B_long(k, l) = B(k) * O_long(k, l) * B(l);
    }
  }
  
  Matrix V_ppprel = A * V_ppp * A + B * O_ppp * B
                  - B_long * Ep * A_long - A_long * Ep * B_long
                  + A_long * EpR * A_long + B_long * EpRinv * B_long
                  - 0.5 * (B_long * A_long * Ep + A_long * B_long * Ep)
                  + 0.5 * (A_long * R * A_long * Ep + B_long * R_inv * B_long * Ep)
                  - 0.5 * (Ep * B_long * A_long + Ep * A_long * B_long)
                  + 0.5 * (Ep * A_long * R * A_long + Ep * B_long * R_inv * B_long);
  const Matrix V_pprel = W * V_ppprel ^ W;
  const Matrix V_prel = s_inv12 * V_pprel * s_inv12;
  const Matrix V_rel = U * V_prel ^ U;

  PU[0].print("PU");
  s_X[0].print("s_X");
  T_pX[0].print("T_pX");
  V_pX[0].print("V_pX");
  O_pX[0].print("O_pX");
  
  for (int i = 0; i < 3 * natom; i++) {
    const Matrix T_ppX = s_inv12 * T_pX[i] * s_inv12 - 0.5 * (s_inv * s_X[i] * T_pp + T_pp *s_inv * s_X[i]);
    T_ppX.print("T_ppX");
    const Matrix T_ppX_W = W % T_ppX * W;
    Matrix PW(nunc, nunc);
    for (int k = 0; k < nunc; ++k) {
      for (int l = 0; l < nunc; ++l) {
        PW(k, l) = k == l || fabs(t(k) - t(l)) < 1.0e-8 ? 0 : T_ppX_W(k, l) / (t(k) - t(l));
      }
    }
    PW.print("PW");
    const Matrix t_X0 = T_ppX_W - PW * t + t * PW;
    t_X0.print("t_X0");
    DiagVec t_X(t_X0);

    DiagVec Ep_X(nunc), A_X(nunc), B_X(nunc), R_X(nunc), R_invX(nunc), EpRX(nunc), EpRinvX(nunc);
    for (int k = 0; k < nunc; ++k) {
      Ep_X(k) = c__ * t_X(k) / std::sqrt(2 * t(k) + c2);
      A_X(k) = c2 * Ep_X(k) / (-4 * pow(Ep(k), 2) * A(k));
      B_X(k) = c__ * Ep_X(k) * (c__ * K(k) / (4 * pow(Ep(k), 2) * A(k)) - A(k) / pow(Ep(k) + c2, 2));
      R_X(k) = 2 * t_X(k) * (pow(K(k), 2) - 2 * c2 / (pow(Ep(k) + c2, 2) * std::sqrt(2 * t(k) + c2)));
      R_invX(k) = -1 * R_X(k) / pow(R(k), 2);
      EpRX(k) = Ep_X(k) * R(k) + Ep(k) * R_X(k);
      EpRinvX(k) = Ep_X(k) * R_inv(k) + Ep(k) * R_invX(k);
    }
    Ep_X.print("Ep_X");
    A_X.print("A_X");
    B_X.print("B_X");
    R_X.print("R_X");
    R_invX.print("R_invX");
    EpRX.print("EpRX");
    EpRinvX.print("EpRinvX");
    
    const Matrix DW = W * PW ^ W;
    const Matrix T_pprelX = W * Ep ^ W + DW * T_pprel - T_pprel * DW;
    const Matrix T_prelX = s_inv12 * T_pprelX * s_inv12 + 0.5 * (s_inv * s_X[i] * T_prel + T_prel * s_inv * s_X[i]);
    const Matrix DU = U * PU[i] ^ U;
    const Matrix T_relX = U * T_prelX ^ U + DU * T_rel - T_rel * DU;

    const Matrix V_ppX = s_inv12 * V_pX[i] * s_inv12 - 0.5 * (s_inv * s_X[i] * V_pp + V_pp * s_inv * s_X[i]);
    const Matrix V_pppX = W % V_ppX * W - PW * V_ppp + V_ppp * PW;
    const Matrix O_ppX = s_inv12 * O_pX[i] * s_inv12 - 0.5 * (s_inv * s_X[i] * O_pp + O_pp * s_inv * s_X[i]);
    const Matrix O_pppX = W % O_ppX * W - PW * O_ppp + O_ppp * PW;

    Matrix V_ppprelX = A * V_pppX * A + A_X * V_ppp * A + A * V_ppp * A_X + B * O_pppX * B + B_X * O_ppp * B + B * O_ppp * B_X;

    Matrix V_longX(nunc, nunc), A_longX(nunc, nunc), O_longX(nunc, nunc), B_longX(nunc, nunc);
    for (int k = 0; k < nunc; k++) {
      for (int l = 0; l < nunc; l++) {
        V_longX(k, l) = -1 * pow(Ep(k) + Ep(l), -2) * (Ep_X(k) + Ep_X(l)) * V_ppp(k, l) + V_pppX(k, l) / (Ep(k) + Ep(l));
        A_longX(k, l) = A(k) * V_longX(k, l) * A(l) + A_X(k) * V_long(k, l) * A(l) + A(k) * V_long(k, l) * A_X(l);
        O_longX(k, l) = -1 * pow(Ep(k) + Ep(l), -2) * (Ep_X(k) + Ep_X(l)) * O_ppp(k, l) + O_pppX(k, l) / (Ep(k) + Ep(l));
        B_longX(k, l) = B(k) * O_longX(k, l) * B(l) + B_X(k) * O_long(k, l) * B(l) + B(k) * O_long(k, l) * B_X(l);
      }
    }
    V_longX.print("V_longX");
    A_longX.print("A_longX");
    O_longX.print("O_longX");
    B_longX.print("B_longX");

    V_ppprelX += -1 * B_longX * Ep * A_long - B_long * Ep_X * A_long - B_long * Ep * A_longX
              - A_longX *  Ep * B_long - A_long * Ep_X * B_long - A_long * Ep * B_longX
              + A_longX * EpR * A_long + A_long * EpRX * A_long + A_long * EpR * A_longX
              + B_longX * EpRinv * B_long + B_long * EpRinvX * B_long + B_long * EpRinv * B_longX
              - 0.5 * (B_longX * A_long * Ep + B_long * A_longX * Ep + B_long * A_long * Ep_X)
              - 0.5 * (A_longX * B_long * Ep + A_long * B_longX * Ep + A_long * B_long * Ep_X)
              + 0.5 * (A_longX * R * A_long * Ep + A_long * R_X * A_long * Ep)
              + 0.5 * (A_long * R * A_longX * Ep + A_long * R * A_long * Ep_X)
              + 0.5 * (B_longX * R_inv * B_long * Ep + B_long * R_invX * B_long * Ep)
              + 0.5 * (B_long * R_inv * B_longX * Ep + B_long * R_inv * B_long * Ep_X)
              - 0.5 * (Ep_X * B_long * A_long + Ep * B_longX * A_long + Ep * B_long * A_longX)
              - 0.5 * (Ep_X * A_long * B_long + Ep * A_longX * B_long + Ep * A_long * B_longX)
              + 0.5 * (Ep_X * A_long * R * A_long + Ep * A_longX * R * A_long)
              + 0.5 * (Ep * A_long * R_X * A_long + Ep * A_long * R * A_longX)
              + 0.5 * (Ep_X * B_long * R_inv * B_long + Ep * B_longX * R_inv * B_long)
              + 0.5 * (Ep * B_long * R_invX * B_long + Ep * B_long * R_inv * B_longX);

    const Matrix V_pprelX = W * V_ppprelX ^ W + DW * V_pprel - V_pprel * DW;
    const Matrix V_prelX = s_inv12 * V_pprelX * s_inv12 + 0.5 * (s_inv * s_X[i] * V_prel + V_prel * s_inv * s_X[i]);
    const Matrix V_relX = U * V_prelX ^ U + DU * V_rel - V_rel * DU;

    T_relX.print("T_relX");
    V_relX.print("V_relX");
    dkh2grad[i] = make_shared<Matrix>(mix % (T_relX + V_relX) * mix);

  }
  return dkh2grad;
}

void DKH2Analytic::gradinit(shared_ptr<const Geometry> mol) {
  natom = mol->natom();
  nunc = mol->nbasis();
  
  const Overlap S(mol);
  U = S;
  VectorB s0(nunc);
  U.diagonalize(s0);
  s = DiagVec(s0);

  PU = vector<Matrix>(3 * natom, Matrix(nunc, nunc));
  s_X = vector<Matrix>(3 * natom, Matrix(nunc, nunc));
  T_pX = vector<Matrix>(3 * natom, Matrix(nunc, nunc));
  V_pX = vector<Matrix>(3 * natom, Matrix(nunc, nunc));
  O_pX = vector<Matrix>(3 * natom, Matrix(nunc, nunc));

  id = Matrix(nunc, nunc);
  for (int i = 0; i < nunc; i++) {
    id(i, i) = 1;
  }
}

void DKH2Analytic::contracts(shared_ptr<const Geometry> mol) {
  const Matrix U_T = U % id;
  for (int k = 0; k < nunc; k++) {
    for (int l = 0; l < nunc; l++) {
      if (k == l) {
        for (int i = 0; i < 3 * natom; i++) {
          PU[i](k, l) = 0;
        }
      }
      else {
        auto grad = make_shared<GradFile>(natom);
        auto den = make_shared<Matrix>(nunc, nunc);

        for (int m = 0; m < nunc; m++) {
          for (int n = 0; n < nunc; n++) {
            (*den)(m, n) = fabs(s(k) - s(l)) < 1.0e-8 ? 0 : U_T(k, m) * U(n, l) / (s(k) - s(l));
          }
        }
        GradEval_base ge(mol);
        grad = ge.contract_overlapgrad(den);
        for (int i = 0; i < natom; i++) {
          for (int j = 0; j < 3; j++) {
            PU[3 * i + j](k, l) += (*grad)(i, j);
          }
        }
      }
    }
  }
}

void DKH2Analytic::overlapgrad(shared_ptr<const Geometry> mol) {
  const Matrix U_T = U % id;
  for (int k = 0; k < nunc; k++) {
    for (int l = 0; l < nunc; l++) {
      auto grad = make_shared<GradFile>(natom);
      auto den = make_shared<Matrix>(nunc, nunc);
      for (int m = 0; m < nunc; m++) {
        for (int n = 0; n < nunc; n++) {
          (*den)(m, n) = U_T(k, m) * U(n, l);
        }
      }
      GradEval_base ge(mol);
      grad = ge.contract_overlapgrad(den);
      for (int i = 0; i < natom; i++) {
        for (int j = 0; j < 3; j++) {
          s_X[3 * i + j](k, l) += (*grad)(i, j);
        }
      }
    }
  }
  for (int i = 0; i < 3 * natom; i++) {
    s_X[i] -= PU[i] * s - s * PU[i];
  }
}

void DKH2Analytic::kineticgrad(shared_ptr<const Geometry> mol, shared_ptr<const Matrix> T_p) {
  const Matrix U_T = U % id;
  for (int k = 0; k < nunc; k++) {
    for (int l = 0; l < nunc; l++) {
      auto grad = make_shared<GradFile>(natom);
      auto den = make_shared<Matrix>(nunc, nunc);
      for (int m = 0; m < nunc; m++) {
        for (int n = 0; n < nunc; n++) {
          (*den)(m, n) = U_T(k, m) * U(n, l);
        }
      }
      GradEval_base ge(mol);
      grad = ge.contract_kineticgrad(den);
      for (int i = 0; i < natom; i++) {
        for (int j = 0; j < 3; j++) {
          T_pX[3 * i + j](k, l) += (*grad)(i, j);
        }
      }
    }
  }
  for (int i = 0; i < 3 * natom; i++) {
    T_pX[i] -= PU[i] * *T_p - *T_p * PU[i];
  }
}

void DKH2Analytic::naigrad(shared_ptr<const Geometry> mol, shared_ptr<const Matrix> V_p) {
  const Matrix U_T = U % id;
  for (int k = 0; k < nunc; k++) {
    for (int l = 0; l < nunc; l++) {
      auto grad = make_shared<GradFile>(natom);
      auto den = make_shared<Matrix>(nunc, nunc);
      for (int m = 0; m < nunc; m++) {
        for (int n = 0; n < nunc; n++) {
          (*den)(m, n) = U_T(k, m) * U(n, l);
        }
      }
      GradEval_base ge(mol);
      grad = ge.contract_naigrad(den);
      for (int i = 0; i < natom; i++) {
        for (int j = 0; j < 3; j++) {
          V_pX[3 * i + j](k, l) += (*grad)(i, j);
        }
      }
    }
  }
  for (int i = 0; i < 3 * natom; i++) {
    V_pX[i] -= PU[i] * *V_p - *V_p * PU[i];
  }
}

void DKH2Analytic::smallnaigrad(shared_ptr<const Geometry> mol, shared_ptr<const Matrix> O_p) {
  mol = make_shared<Geometry>(*mol->relativistic(false));
  const Matrix U_T = U % id;
  for (int k = 0; k < nunc; k++) {
    for (int l = 0; l < nunc; l++) {
      auto grad = make_shared<GradFile>(natom);
      auto den = make_shared<Matrix>(nunc, nunc);
      auto zero = make_shared<Matrix>(*den);
      for (int m = 0; m < nunc; m++) {
        for (int n = 0; n < nunc; n++) {
          (*den)(m, n) = U_T(k, m) * U(n, l);
        }
      }
      GradEval_base ge(mol);
      grad = ge.contract_smallnaigrad({den, zero, zero, den, zero, den});
      for (int i = 0; i < natom; i++) {
        for (int j = 0; j < 3; j++) {
          O_pX[3 * i + j](k, l) += (*grad)(i, j);
        }
      }
    }
  }
  for (int i = 0; i < 3 * natom; i++) {
    O_pX[i] -= PU[i] * *O_p - *O_p * PU[i];
  }
}

