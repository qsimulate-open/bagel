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

  auto s_inv12 = make_shared<VectorB>(nunc);
  auto s_inv = make_shared<VectorB>(nunc);
  for (int k = 0; k < nunc; ++k) {
    (*s_inv12)(k) = 1 / std::sqrt((*s)(k));
    (*s_inv)(k) = 1 / (*s)(k);
  }
  store_mat(s_inv12);
  store_mat(s_inv);

  const Kinetic T(mol);
  auto T_p = make_shared<Matrix>(U % T * U);
  kineticgrad(mol, T_p);
  Matrix T_pp = *vec2mat[s_inv12] % *T_p * *vec2mat[s_inv12];
  Matrix W = T_pp;
  auto t = make_shared<VectorB>(nunc);
  W.diagonalize(*t);
  T_pp.print("T_pp");
  W.print("W");
  store_mat(t);
  cout << "t" << endl;
  for (int k = 0; k < nunc; k++) {
    cout << (*t)(k) << " ";
  }
  cout << endl;

  const double c2 = c__ * c__;
  auto Ep = make_shared<VectorB>(nunc);
  auto A = make_shared<VectorB>(nunc);
  auto K = make_shared<VectorB>(nunc);
  auto B = make_shared<VectorB>(nunc);
  auto R = make_shared<VectorB>(nunc);
  auto R_inv = make_shared<VectorB>(nunc);
  auto EpR = make_shared<VectorB>(nunc);
  auto EpRinv = make_shared<VectorB>(nunc);
  for (int k = 0; k < nunc; ++k) {
    (*Ep)(k) = c__ * std::sqrt(2 * (*t)(k) + c2);
    (*A)(k) = std::sqrt((c2 + (*Ep)(k)) / (2 * (*Ep)(k)));
    (*K)(k) = c__ / ((*Ep)(k) + c2);
    (*B)(k) = (*A)(k) * (*K)(k);
    (*R)(k) = 2 * (*t)(k) * pow((*K)(k), 2);
    (*R_inv)(k) = 1 / (*R)(k);
    (*EpR)(k) = (*Ep)(k) * (*R)(k);
    (*EpRinv)(k) = (*Ep)(k) * (*R_inv)(k);
  }
  store_mat(Ep);
  store_mat(A);
  store_mat(K);
  store_mat(B);
  store_mat(R);
  store_mat(R_inv);
  store_mat(EpR);
  store_mat(EpRinv);

  cout << "Ep" << endl;
  for (int k = 0; k < nunc; k++) {
    cout << (*Ep)(k) << " ";
  }
  cout << endl;
  cout << "A" << endl;
  for (int k = 0; k < nunc; k++) {
    cout << (*A)(k) << " ";
  }
  cout << endl;
  cout << "K" << endl;
  for (int k = 0; k < nunc; k++) {
    cout << (*K)(k) << " ";
  }
  cout << endl;
  cout << "B" << endl;
  for (int k = 0; k < nunc; k++) {
    cout << (*B)(k) << " ";
  }
  cout << endl;
  cout << "R" << endl;
  for (int k = 0; k < nunc; k++) {
    cout << (*R)(k) << " ";
  }
  cout << endl;
  cout << "R_inv" << endl;
  for (int k = 0; k < nunc; k++) {
    cout << (*R_inv)(k) << " ";
  }
  cout << endl;
  cout << "EpR" << endl;
  for (int k = 0; k < nunc; k++) {
    cout << (*EpR)(k) << " ";
  }
  cout << endl;
  cout << "EpRinv" << endl;
  for (int k = 0; k < nunc; k++) {
    cout << (*EpRinv)(k) << " ";
  }
  cout << endl;

  auto t_rel = make_shared<VectorB>(nunc);
  for (int k = 0; k < nunc; ++k) {
    (*t_rel)(k) = c__ * std::sqrt(2 * (*t)(k) + c2) - c2;
  }
  store_mat(t_rel);
  const Matrix T_pprel = W * *vec2mat[t_rel] ^ W;
  const Matrix T_prel = *vec2mat[s_inv12] * T_pprel ^ *vec2mat[s_inv12];
  const Matrix T_rel = U * T_prel ^ U;
  T_rel.print("T_rel");

  const NAI V(mol);
  auto V_p = make_shared<Matrix>(U % V * U);
  naigrad(mol, V_p);
  const Matrix V_pp = *vec2mat[s_inv12] % *V_p * *vec2mat[s_inv12];
  V_pp.print("V_pp");
  const Matrix V_ppp = W % V_pp * W;
  const Small1e<NAIBatch> small1e(mol);
  auto O_p = make_shared<Matrix>(U % small1e[0] * U);
  smallnaigrad(mol, O_p);
  const Matrix O_pp = *vec2mat[s_inv12] % *O_p * *vec2mat[s_inv12];
  O_pp.print("O_pp");
  const Matrix O_ppp = W % O_pp * W;

  Matrix V_long(nunc, nunc), A_long(nunc, nunc), O_long(nunc, nunc), B_long(nunc, nunc);
  for (int k = 0; k < nunc; k++) {
    for (int l = 0; l < nunc; l++) {
      V_long(k, l) = V_ppp(k, l) / ((*Ep)(k) + (*Ep)(l));
      A_long(k, l) = (*A)(k) * V_long(k, l) * (*A)(l);
      O_long(k, l) = O_ppp(k, l) / ((*Ep)(k) + (*Ep)(l));
      B_long(k, l) = (*B)(k) * O_long(k, l) * (*B)(l);
    }
  }
  
  Matrix V_ppprel = *vec2mat[A] * V_ppp * *vec2mat[A] + *vec2mat[B] * O_ppp * *vec2mat[B]
                  - B_long * *vec2mat[Ep] * A_long - A_long * *vec2mat[Ep] * B_long
                  + A_long * *vec2mat[EpR] * A_long + B_long * *vec2mat[EpRinv] * B_long
                  - 0.5 * (B_long * A_long * *vec2mat[Ep] + A_long * B_long * *vec2mat[Ep])
                  + 0.5 * (A_long * *vec2mat[R] * A_long * *vec2mat[Ep] + B_long * *vec2mat[R_inv] * B_long * *vec2mat[Ep])
                  - 0.5 * (*vec2mat[Ep] * B_long * A_long + *vec2mat[Ep] * A_long * B_long)
                  + 0.5 * (*vec2mat[Ep] * A_long * *vec2mat[R] * A_long + *vec2mat[Ep] * B_long * *vec2mat[R_inv] * B_long);
  const Matrix V_pprel = W * V_ppprel ^ W;
  const Matrix V_prel = *vec2mat[s_inv12] * V_pprel ^ *vec2mat[s_inv12];
  const Matrix V_rel = U * V_prel ^ U;
  V_rel.print("V_rel");
  
  for (int i = 0; i < 3 * natom; i++) {
    const Matrix T_ppX = *vec2mat[s_inv12] % T_pX[i] * *vec2mat[s_inv12] - 0.5 * (*vec2mat[s_inv] * s_X[i] * T_pp + T_pp * *vec2mat[s_inv] * s_X[i]);
    T_pX[i].print("MATRIX T_pX");
    s_X[i].print("MATRIX s_X");
    T_ppX.print("MATRIX T_ppX");
    const Matrix T_ppX_W = W % T_ppX * W;
    T_ppX_W.print("MATRIX T_ppX_W");
    Matrix PW(nunc, nunc);
    for (int k = 0; k < nunc; ++k) {
      for (int l = 0; l < nunc; ++l) {
        PW(k, l) = k == l ? 0 : T_ppX_W(k, l) / ((*t)(k) - (*t)(l));
      }
    }
    PW.print("MATRIX PW");
    const Matrix t_X = T_ppX_W - PW * *vec2mat[t] + *vec2mat[t] * PW;

    auto Ep_X = make_shared<VectorB>(nunc);
    auto A_X = make_shared<VectorB>(nunc);
    auto B_X = make_shared<VectorB>(nunc);
    auto R_X = make_shared<VectorB>(nunc);
    auto R_invX = make_shared<VectorB>(nunc);
    auto EpRX = make_shared<VectorB>(nunc);
    auto EpRinvX = make_shared<VectorB>(nunc);
    for (int k = 0; k < nunc; ++k) {
      (*Ep_X)(k) = c__ * t_X(k, k) / std::sqrt(2 * (*t)(k) + c2);
      (*A_X)(k) = c2 * (*Ep_X)(k) / (-4 * pow((*Ep)(k), 2) * (*A)(k));
      (*B_X)(k) = c__ * (*Ep_X)(k) * (c__ * (*K)(k) / (4 * pow((*Ep)(k), 2) * (*A)(k)) - (*A)(k) / pow((*Ep)(k) + c2, 2));
      (*R_X)(k) = 2 * t_X(k, k) * (pow((*K)(k), 2) - 2 * c2 / (pow((*Ep)(k) + c2, 2) * std::sqrt(2 * (*t)(k) + c2)));
      (*R_invX)(k) = -1 * (*R_X)(k) / pow((*R)(k), 2);
      (*EpRX)(k) = (*Ep_X)(k) * (*R)(k) + (*Ep)(k) * (*R_X)(k);
      (*EpRinvX)(k) = (*Ep_X)(k) * (*R_inv)(k) + (*Ep)(k) * (*R_invX)(k);
    }
    store_mat(Ep_X);
    store_mat(A_X);
    store_mat(B_X);
    store_mat(R_X);
    store_mat(R_invX);
    store_mat(EpRX);
    store_mat(EpRinvX);

    cout << "Ep_X" << endl;
    for (int k = 0; k < nunc; k++) {
      cout << (*Ep_X)(k) << " ";
    }
    cout << endl;
    cout << "A_X" << endl;
    for (int k = 0; k < nunc; k++) {
      cout << (*A_X)(k) << " ";
    }
    cout << endl;
    cout << "B_X" << endl;
    for (int k = 0; k < nunc; k++) {
      cout << (*B_X)(k) << " ";
    }
    cout << endl;
    cout << "R_X" << endl;
    for (int k = 0; k < nunc; k++) {
      cout << (*R_X)(k) << " ";
    }
    cout << endl;
    cout << "R_invX" << endl;
    for (int k = 0; k < nunc; k++) {
      cout << (*R_invX)(k) << " ";
    }
    cout << endl;
    cout << "EpRX" << endl;
    for (int k = 0; k < nunc; k++) {
      cout << (*EpRX)(k) << " ";
    }
    cout << endl;
    cout << "EpRinvX" << endl;
    for (int k = 0; k < nunc; k++) {
      cout << (*EpRinvX)(k) << " ";
    }
    cout << endl;
    
    const Matrix DW = W * PW ^ W;
    const Matrix T_pprelX = W * *vec2mat[Ep] ^ W + DW * T_pprel - T_pprel * DW;
    const Matrix T_prelX = *vec2mat[s_inv12] * T_pprelX ^ *vec2mat[s_inv12] + 0.5 * (*vec2mat[s_inv] * s_X[i] * T_prel + T_prel * *vec2mat[s_inv] * s_X[i]);
    const Matrix DU = U * PU[i] ^ U;
    const Matrix T_relX = U * T_prelX ^ U + DU * T_rel - T_rel * DU;

    const Matrix V_ppX = *vec2mat[s_inv12] % V_pX[i] * *vec2mat[s_inv12] - 0.5 * (*vec2mat[s_inv] * s_X[i] * V_pp + V_pp * *vec2mat[s_inv] * s_X[i]);
    const Matrix V_pppX = W % V_ppX * W - PW * V_ppp + V_ppp * PW;
    const Matrix O_ppX = *vec2mat[s_inv12] % O_pX[i] * *vec2mat[s_inv12] - 0.5 * (*vec2mat[s_inv] * s_X[i] * O_pp + O_pp * *vec2mat[s_inv] * s_X[i]);
    const Matrix O_pppX = W % O_ppX * W - PW * O_ppp + O_ppp * PW;

    Matrix V_ppprelX = *vec2mat[A] * V_pppX * *vec2mat[A] + *vec2mat[A_X] * V_ppp * *vec2mat[A] + *vec2mat[A] * V_ppp * *vec2mat[A_X]
                      + *vec2mat[B] * O_pppX * *vec2mat[B] + *vec2mat[B_X] * O_ppp * *vec2mat[B] + *vec2mat[B] * O_ppp * *vec2mat[B_X];

    Matrix V_longX(nunc, nunc), A_longX(nunc, nunc), O_longX(nunc, nunc), B_longX(nunc, nunc);
    for (int k = 0; k < nunc; k++) {
      for (int l = 0; l < nunc; l++) {
        V_longX(k, l) = -1 * pow((*Ep)(k) + (*Ep)(l), -2) * ((*Ep_X)(k) + (*Ep_X)(l)) * V_ppp(k, l) + V_pppX(k, l) / ((*Ep)(k) + (*Ep)(l));
        A_longX(k, l) = (*A)(k) * V_longX(k, l) * (*A)(l) + (*A_X)(k) * V_long(k, l) * (*A)(l) + (*A)(k) * V_long(k, l) * (*A_X)(l);
        O_longX(k, l) = -1 * pow((*Ep)(k) + (*Ep)(l), -2) * ((*Ep_X)(k) + (*Ep_X)(l)) * O_ppp(k, l) + O_pppX(k, l) / ((*Ep)(k) + (*Ep)(l));
        B_longX(k, l) = (*B)(k) * O_longX(k, l) * (*B)(l) + (*B_X)(k) * O_long(k, l) * (*B)(l) + (*B)(k) * O_long(k, l) * (*B_X)(l);
      }
    }

    V_ppprelX += -1 * B_longX * *vec2mat[Ep] * A_long - B_long * *vec2mat[Ep_X] * A_long - B_long * *vec2mat[Ep] * A_longX
              - A_longX *  *vec2mat[Ep] * B_long - A_long * *vec2mat[Ep_X] * B_long - A_long * *vec2mat[Ep] * B_longX
              + A_longX * *vec2mat[EpR] * A_long + A_long * *vec2mat[EpRX] * A_long + A_long * *vec2mat[EpR] * A_longX
              + B_longX * *vec2mat[EpRinv] * B_long + B_long * *vec2mat[EpRinvX] * B_long + B_long * *vec2mat[EpRinv] * B_longX
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

    const Matrix V_pprelX = W * V_ppprelX ^ W + DW * V_pprel - V_pprel * DW;
    const Matrix V_prelX = *vec2mat[s_inv12] * V_pprelX ^ *vec2mat[s_inv12] + 0.5 * (*vec2mat[s_inv] * s_X[i] * V_prel + V_prel * *vec2mat[s_inv] * s_X[i]);
    const Matrix V_relX = U * V_prelX ^ U + DU * V_rel - V_rel * DU;

    dkh2grad[i] = make_shared<Matrix>(mix % (T_relX + V_relX) * mix);

  }
  const Matrix dkh = mix % (T_rel + V_rel) * mix;
  dkh.print();
  return dkh2grad;
}

void DKH2Analytic::gradinit(shared_ptr<const Geometry> mol) {
  natom = mol->natom();
  nunc = mol->nbasis();
  
  const Overlap S(mol);
  U = S;
  s = make_shared<VectorB>(nunc);
  U.diagonalize(*s);
  store_mat(s);

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
            (*den)(m, n) = U_T(k, m) * U(n, l) / ((*s)(k) - (*s)(l));
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
    s_X[i] -= PU[i] * *vec2mat[s] - *vec2mat[s] * PU[i];
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

void DKH2Analytic::store_mat(shared_ptr<const VectorB> vec) {
  auto mat = make_shared<Matrix>(nunc, nunc);
  for (int i = 0; i < nunc; i++) {
    (*mat)(i, i) = (*vec)(i);
  }
  vec2mat[vec] = mat;
}

