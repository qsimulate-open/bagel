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

vector<shared_ptr<Matrix>> DKH2Analytic::dkh_grad(shared_ptr<const Molecule> current) const {
  shared_ptr<const Geometry> mol = dynamic_pointer_cast<const Geometry>(current);
  assert(mol);
  shared_ptr<DKH2AnalyticData> data = gradinit(mol);

  vector<shared_ptr<Matrix>> dkh2grad(3 * data->natom);
  auto s = make_shared<Overlap>(data->molu);
  contracts(mol, s, data);
  overlapgrad(mol, s, data);
  const Matrix s_inv12 = *s->tildex();
  Matrix s_inv = *s;
  s_inv.inverse_symmetric();
  auto T_p = make_shared<Kinetic>(data->molu);
  kineticgrad(mol, T_p, data);
  const Matrix T_pp = s_inv12 % *T_p * s_inv12;
  auto t = make_shared<VectorB>(data->nunc);
  store_mat(t, data);
  Matrix W = T_pp;
  W.diagonalize(*t);

  const double c2 = c__ * c__;
  auto Ep = make_shared<VectorB>(data->nunc);
  store_mat(Ep, data);
  auto A = make_shared<VectorB>(data->nunc);
  store_mat(A, data);
  auto K = make_shared<VectorB>(data->nunc);
  store_mat(K, data);
  auto B = make_shared<VectorB>(data->nunc);
  store_mat(B, data);
  auto R = make_shared<VectorB>(data->nunc);
  store_mat(R, data);
  auto R_inv = make_shared<VectorB>(data->nunc);
  store_mat(R_inv, data);
  for (int k = 0; k < data->nunc; ++k) {
    (*Ep)(k) = c__ * std::sqrt(2 * (*t)(k) + c2);
    (*A)(k) = std::sqrt((c2 + (*Ep)(k)) / (2 * (*Ep)(k)));
    (*K)(k) = c__ / ((*Ep)(k) + c2);
    (*B)(k) = (*A)(k) * (*K)(k);
    (*R)(k) = 2 * (*t)(k) * pow((*K)(k), 2);
    (*R_inv)(k) = 1 / (*R)(k);
  }

  for (int i = 0; i < 3; i++) {
    const Matrix T_ppX = s_inv12 % data->T_pX[i] * s_inv12 - 0.5 * (s_inv * data->s_X[i] * T_pp + T_pp * s_inv * data->s_X[i]);
    const Matrix T_ppX_W = W % T_ppX * W;
    Matrix PW(data->nunc, data->nunc);
    for (int k = 0; k < data->nunc; ++k) {
      for (int l = 0; l < data->nunc; ++l) {
        PW(k, l) = k == l ? 0 : T_ppX_W(k, l) / ((*t)(k) - (*t)(l));
      }
    }
    const Matrix t_X = T_ppX_W - PW * *data->vec2mat[t] + *data->vec2mat[t] * PW;

    auto Ep_X = make_shared<VectorB>(data->nunc);
    store_mat(Ep_X, data);
    auto A_X = make_shared<VectorB>(data->nunc);
    store_mat(A_X, data);
    auto B_X = make_shared<VectorB>(data->nunc);
    store_mat(B_X, data);
    auto R_X = make_shared<VectorB>(data->nunc);
    store_mat(R_X, data);
    auto R_invX = make_shared<VectorB>(data->nunc);
    store_mat(R_invX, data);
    for (int k = 0; k < data->nunc; ++k) {
      (*Ep)(k) = c__ * std::sqrt(2 * (*t)(k) + c2);
      (*A)(k) = std::sqrt((c2 + (*Ep)(k)) / (2 * (*Ep)(k)));
      (*K)(k) = c__ / ((*Ep)(k) + c2);
      (*B)(k) = (*A)(k) * (*K)(k);
      (*Ep_X)(k) = c__ * t_X(k, k) / std::sqrt(2 * (*t)(k) + c2);
      (*A_X)(k) = c2 * (*Ep_X)(k) / (-4 * pow((*Ep)(k), 2) * (*A)(k));
      (*B_X)(k) = c__ * (*Ep_X)(k) * (c__ * (*K)(k) / (4 * pow((*Ep)(k), 2) * (*A)(k)) - (*A)(k) / pow((*Ep)(k) + c2, 2));
      (*R)(k) = 2 * (*t)(k) * pow((*K)(k), 2);
      (*R_inv)(k) = 1 / (*R)(k);
      (*R_X)(k) = 2 * t_X(k, k) * (pow((*K)(k), 2) - 2 * c2 / (pow((*Ep)(k) + c2, 2) * std::sqrt(2 * (*t)(k) + c2)));
      (*R_invX)(k) = -1 * (*R_X)(k) / pow((*R)(k), 2);
    }
    
    auto t_rel = make_shared<VectorB>(data->nunc);
    store_mat(t_rel, data);
    for (int k = 0; k < data->nunc; ++k) {
      (*t_rel)(k) = c__ * std::sqrt(2 * (*t)(k) + c2) - c2;
    }
    const Matrix T_pprel = W * *data->vec2mat[t_rel] ^ W;
    const Matrix T_prel = s_inv12 * T_pprel ^ s_inv12;
    const Matrix T_rel = data->U_T % T_prel * data->U_T;
    const Matrix DW = W * PW ^ W;
    const Matrix T_pprelX = W * *data->vec2mat[Ep] ^ W + DW * T_pprel - T_pprel * DW;
    const Matrix T_prelX = s_inv12 * T_pprelX ^ s_inv12 + 0.5 * (s_inv * data->s_X[i] * T_prel + T_prel * s_inv * data->s_X[i]);
    const Matrix DU = data->U_T % data->PU[i] * data->U_T;
    const Matrix T_relX = data->U_T % T_prelX * data->U_T + DU * T_rel - T_rel * DU;

    auto V_p = make_shared<NAI>(data->molu);
    naigrad(mol, V_p, data);
    const Matrix V_pp = s_inv12 % *V_p * s_inv12;
    const Matrix V_ppp = W % V_pp * W;
    const Matrix V_ppX = s_inv12 % data->V_pX[i] * s_inv12 - 0.5 * (s_inv * data->s_X[i] * V_pp + V_pp * s_inv * data->s_X[i]);
    const Matrix V_pppX = W % V_ppX * W - PW * V_ppp + V_ppp * PW;

    const Small1e<NAIBatch> small1e(data->molu);
    auto O_p = make_shared<Matrix>(small1e[0]);
    smallnaigrad(mol, O_p, data);
    const Matrix O_pp = s_inv12 % *O_p * s_inv12;
    const Matrix O_ppp = W % O_pp * W;
    const Matrix O_ppX = s_inv12 % data->O_pX[i] * s_inv12 - 0.5 * (s_inv * data->s_X[i] * O_pp + O_pp * s_inv * data->s_X[i]);
    const Matrix O_pppX = W % O_ppX * W - PW * O_ppp + O_ppp * PW;

    Matrix V_ppprelX = *data->vec2mat[A] * V_pppX * *data->vec2mat[A] + *data->vec2mat[A_X] * V_ppp * *data->vec2mat[A] + *data->vec2mat[A] * V_ppp * *data->vec2mat[A_X]
                      + *data->vec2mat[B] * O_pppX * *data->vec2mat[B] + *data->vec2mat[B_X] * O_ppp * *data->vec2mat[B] + *data->vec2mat[B] * O_ppp * *data->vec2mat[B_X];

    Matrix V_long(data->nunc, data->nunc), A_long(data->nunc, data->nunc), O_long(data->nunc, data->nunc), B_long(data->nunc, data->nunc);
    Matrix V_longX(data->nunc, data->nunc), A_longX(data->nunc, data->nunc), O_longX(data->nunc, data->nunc), B_longX(data->nunc, data->nunc);
    for (int k = 0; k < data->nunc; k++) {
      for (int l = 0; l < data->nunc; l++) {
        V_long(k, l) = V_ppp(k, l) / ((*Ep)(k) + (*Ep)(l));
        A_long(k, l) = (*A)(k) * V_long(k, l) * (*A)(l);
        O_long(k, l) = O_ppp(k, l) / ((*Ep)(k) + (*Ep)(l));
        B_long(k, l) = (*B)(k) * O_long(k, l) * (*B)(l);
        V_longX(k, l) = -1 * pow((*Ep)(k) + (*Ep)(l), -2) * ((*Ep_X)(k) + (*Ep_X)(l)) * V_ppp(k, l) + V_pppX(k, l) / ((*Ep)(k) + (*Ep)(l));
        A_longX(k, l) = (*A)(k) * V_longX(k, l) * (*A)(l) + (*A_X)(k) * V_long(k, l) * (*A)(l) + (*A)(k) * V_long(k, l) * (*A_X)(l);
        O_longX(k, l) = -1 * pow((*Ep)(k) + (*Ep)(l), -2) * ((*Ep_X)(k) + (*Ep_X)(l)) * O_ppp(k, l) + O_pppX(k, l) / ((*Ep)(k) + (*Ep)(l));
        B_longX(k, l) = (*B)(k) * O_longX(k, l) * (*B)(l) + (*B_X)(k) * O_long(k, l) * (*B)(l) + (*B)(k) * O_long(k, l) * (*B_X)(l);
      }
    }

    auto EpR = make_shared<VectorB>(*Ep % *R);
    store_mat(EpR, data);
    auto EpRX = make_shared<VectorB>(*Ep_X % *R + *Ep % *R_X);
    store_mat(EpRX, data);
    auto EpRinv = make_shared<VectorB>(*Ep % *R_inv);
    store_mat(EpRinv, data);
    auto EpRinvX = make_shared<VectorB>(*Ep_X % *R_inv + *Ep % *R_invX);
    store_mat(EpRinvX, data);
    V_ppprelX += -1 * B_longX * *data->vec2mat[Ep] * A_long - B_long * *data->vec2mat[Ep_X] * A_long - B_long * *data->vec2mat[Ep] * A_longX
              - A_longX *  *data->vec2mat[Ep] * B_long - A_long * *data->vec2mat[Ep_X] * B_long - A_long * *data->vec2mat[Ep] * B_longX
              + A_longX * *data->vec2mat[EpR] * A_long + A_long * *data->vec2mat[EpRX] * A_long + A_long * *data->vec2mat[EpR] * A_longX
              + B_longX * *data->vec2mat[EpRinv] * B_long + B_long * *data->vec2mat[EpRinvX] * B_long + B_long * *data->vec2mat[EpRinv] * B_longX
              - 0.5 * (B_longX * A_long * *data->vec2mat[Ep] + B_long * A_longX * *data->vec2mat[Ep] + B_long * A_long * *data->vec2mat[Ep_X])
              - 0.5 * (A_longX * B_long * *data->vec2mat[Ep] + A_long * B_longX * *data->vec2mat[Ep] + A_long * B_long * *data->vec2mat[Ep_X])
              + 0.5 * (A_longX * *data->vec2mat[R] * A_long * *data->vec2mat[Ep] + A_long * *data->vec2mat[R_X] * A_long * *data->vec2mat[Ep])
              + 0.5 * (A_long * *data->vec2mat[R] * A_longX * *data->vec2mat[Ep] + A_long * *data->vec2mat[R] * A_long * *data->vec2mat[Ep_X])
              + 0.5 * (B_longX * *data->vec2mat[R_inv] * B_long * *data->vec2mat[Ep] + B_long * *data->vec2mat[R_invX] * B_long * *data->vec2mat[Ep])
              + 0.5 * (B_long * *data->vec2mat[R_inv] * B_longX * *data->vec2mat[Ep] + B_long * *data->vec2mat[R_inv] * B_long * *data->vec2mat[Ep_X])
              - 0.5 * (*data->vec2mat[Ep_X] * B_long * A_long + *data->vec2mat[Ep] * B_longX * A_long + *data->vec2mat[Ep] * B_long * A_longX)
              - 0.5 * (*data->vec2mat[Ep_X] * A_long * B_long + *data->vec2mat[Ep] * A_longX * B_long + *data->vec2mat[Ep] * A_long * B_longX)
              + 0.5 * (*data->vec2mat[Ep_X] * A_long * *data->vec2mat[R] * A_long + *data->vec2mat[Ep] * A_longX * *data->vec2mat[R] * A_long)
              + 0.5 * (*data->vec2mat[Ep] * A_long * *data->vec2mat[R_X] * A_long + *data->vec2mat[Ep] * A_long * *data->vec2mat[R] * A_longX)
              + 0.5 * (*data->vec2mat[Ep_X] * B_long * *data->vec2mat[R_inv] * B_long + *data->vec2mat[Ep] * B_longX * *data->vec2mat[R_inv] * B_long)
              + 0.5 * (*data->vec2mat[Ep] * B_long * *data->vec2mat[R_invX] * B_long + *data->vec2mat[Ep] * B_long * *data->vec2mat[R_inv] * B_longX);

    Matrix V_ppprel = *data->vec2mat[A] * V_ppp * *data->vec2mat[A] + *data->vec2mat[B] * O_ppp * *data->vec2mat[B]
                    - B_long * *data->vec2mat[Ep] * A_long - A_long * *data->vec2mat[Ep] * B_long
                    + A_long * *data->vec2mat[EpR] * A_long + B_long * *data->vec2mat[EpRinv] * B_long
                    - 0.5 * (B_long * A_long * *data->vec2mat[Ep] + A_long * B_long * *data->vec2mat[Ep])
                    + 0.5 * (A_long * *data->vec2mat[R] * A_long * *data->vec2mat[Ep] + B_long * *data->vec2mat[R_inv] * B_long * *data->vec2mat[Ep])
                    - 0.5 * (*data->vec2mat[Ep] * B_long * A_long + *data->vec2mat[Ep] * A_long * B_long)
                    + 0.5 * (*data->vec2mat[Ep] * A_long * *data->vec2mat[R] * A_long + *data->vec2mat[Ep] * B_long * *data->vec2mat[R_inv] * B_long);
    
    const Matrix V_pprel = W * V_ppprel ^ W;
    const Matrix V_prel = s_inv12 * V_pprel ^ s_inv12;
    const Matrix V_rel = data->U_T % V_prel * data->U_T;
    const Matrix V_pprelX = W * V_ppprelX ^ W + DW * V_pprel - V_pprel * DW;
    const Matrix V_prelX = s_inv12 * V_pprelX ^ s_inv12 + 0.5 * (s_inv * data->s_X[i] * V_prel + V_prel * s_inv * data->s_X[i]);
    const Matrix V_relX = data->U_T % V_prelX * data->U_T + DU * V_rel - V_rel * DU;

    *dkh2grad[i] = T_relX + V_relX;

  }
  return dkh2grad;
}

shared_ptr<DKH2AnalyticData> DKH2Analytic::gradinit(shared_ptr<const Geometry> mol) const {
  auto data = make_shared<DKH2AnalyticData>();
  data->natom = mol->natom();
  data->nbasis = mol->nbasis();
  data->molu = make_shared<Molecule>(*mol->uncontract());
  data->nunc = data->molu->nbasis();
  data->U_T = MixedBasis<OverlapBatch>(mol, data->molu);

  data->PU = vector<Matrix>(3 * data->natom, Matrix(data->nunc, data->nunc));
  data->s_X = vector<Matrix>(3 * data->natom, Matrix(data->nunc, data->nunc));
  data->T_pX = vector<Matrix>(3 * data->natom, Matrix(data->nunc, data->nunc));
  data->V_pX = vector<Matrix>(3 * data->natom, Matrix(data->nunc, data->nunc));
  data->O_pX = vector<Matrix>(3 * data->natom, Matrix(data->nunc, data->nunc));

  data->id = Matrix(data->nunc, data->nunc);
  for (int i = 0; i < data->nunc; i++) {
    data->id(i, i) = 1;
  }
  return data;
}

void DKH2Analytic::contracts(shared_ptr<const Geometry> mol, shared_ptr<const Matrix> s, shared_ptr<DKH2AnalyticData> data) const {
  const Matrix U = data->U_T % data->id;
  for (int k = 0; k < data->nunc; k++) {
    for (int l = 0; l < data->nunc; l++) {
      if (k == l) {
        for (int i = 0; i < 3 * data->natom; i++) {
          data->PU[i](k, l) = 0;
        }
      }
      else {
        auto grad = make_shared<GradFile>(data->natom);
        auto den = make_shared<Matrix>(data->nbasis, data->nbasis);

        for (int m = 0; m < data->nbasis; m++) {
          for (int n = 0; n < data->nbasis; n++) {
            (*den)(m, n) = data->U_T(k, m) * U(n, l) / ((*s)(k, k) - (*s)(l, l));
          }
        }
        GradEval_base ge(mol);
        grad = ge.contract_overlapgrad(den);
        for (int i = 0; i < data->natom; i++) {
          for (int j = 0; j < 3; j++) {
            data->PU[3 * i + j](k, l) = (*grad)(i, j);
          }
        }
      }
    }
  }
}

void DKH2Analytic::overlapgrad(shared_ptr<const Geometry> mol, shared_ptr<const Matrix> s, shared_ptr<DKH2AnalyticData> data) const {
  const Matrix U = data->U_T % data->id;
  for (int k = 0; k < data->nunc; k++) {
    for (int l = 0; l < data->nunc; l++) {
      auto grad = make_shared<GradFile>(data->natom);
      auto den = make_shared<Matrix>(data->nbasis, data->nbasis);
      for (int m = 0; m < data->nbasis; m++) {
        for (int n = 0; n < data->nbasis; n++) {
          (*den)(m, n) = data->U_T(k, m) * U(n, l);
        }
      }
      GradEval_base ge(mol);
      grad = ge.contract_naigrad(den);
      for (int i = 0; i < data->natom; i++) {
        for (int j = 0; j < 3; j++) {
          data->s_X[3 * i + j](k, l) = (*grad)(i, j);
        }
      }
    }
  }
  for (int i = 0; i < 3 * data->natom; i++) {
    data->s_X[i] -= data->PU[i] * *s - *s * data->PU[i];
  }
}

void DKH2Analytic::kineticgrad(shared_ptr<const Geometry> mol, shared_ptr<const Matrix> T_p, shared_ptr<DKH2AnalyticData> data) const {
  const Matrix U = data->U_T % data->id;
  for (int k = 0; k < data->nunc; k++) {
    for (int l = 0; l < data->nunc; l++) {
      auto grad = make_shared<GradFile>(data->natom);
      auto den = make_shared<Matrix>(data->nbasis, data->nbasis);
      for (int m = 0; m < data->nbasis; m++) {
        for (int n = 0; n < data->nbasis; n++) {
          (*den)(m, n) = data->U_T(k, m) * U(n, l);
        }
      }
      GradEval_base ge(mol);
      grad = ge.contract_kineticgrad(den);
      for (int i = 0; i < data->natom; i++) {
        for (int j = 0; j < 3; j++) {
          data->T_pX[3 * i + j](k, l) = (*grad)(i, j);
        }
      }
    }
  }
  for (int i = 0; i < 3 * data->natom; i++) {
    data->T_pX[i] -= data->PU[i] * *T_p - *T_p * data->PU[i];
  }
}

void DKH2Analytic::naigrad(shared_ptr<const Geometry> mol, shared_ptr<const Matrix> V_p, shared_ptr<DKH2AnalyticData> data) const {
  const Matrix U = data->U_T % data->id;
  for (int k = 0; k < data->nunc; k++) {
    for (int l = 0; l < data->nunc; l++) {
      auto grad = make_shared<GradFile>(data->natom);
      auto den = make_shared<Matrix>(data->nbasis, data->nbasis);
      for (int m = 0; m < data->nbasis; m++) {
        for (int n = 0; n < data->nbasis; n++) {
          (*den)(m, n) = data->U_T(k, m) * U(n, l);
        }
      }
      GradEval_base ge(mol);
      grad = ge.contract_naigrad(den);
      for (int i = 0; i < data->natom; i++) {
        for (int j = 0; j < 3; j++) {
          data->V_pX[3 * i + j](k, l) = (*grad)(i, j);
        }
      }
    }
  }
  for (int i = 0; i < 3 * data->natom; i++) {
    data->V_pX[i] -= data->PU[i] * *V_p - *V_p * data->PU[i];
  }
}

void DKH2Analytic::smallnaigrad(shared_ptr<const Geometry> mol, shared_ptr<const Matrix> O_p, shared_ptr<DKH2AnalyticData> data) const {
  const Matrix U = data->U_T % data->id;
  for (int k = 0; k < data->nunc; k++) {
    for (int l = 0; l < data->nunc; l++) {
      auto grad = make_shared<GradFile>(data->natom);
      auto den = make_shared<Matrix>(data->nbasis, data->nbasis);
      auto zero = make_shared<Matrix>(*den);
      for (int m = 0; m < data->nbasis; m++) {
        for (int n = 0; n < data->nbasis; n++) {
          (*den)(m, n) = data->U_T(k, m) * U(n, l);
        }
      }
      GradEval_base ge(mol);
      grad = ge.contract_smallnaigrad({den, zero, zero, den, zero, den});
      for (int i = 0; i < data->natom; i++) {
        for (int j = 0; j < 3; j++) {
          data->O_pX[3 * i + j](k, l) = (*grad)(i, j);
        }
      }
    }
  }
  for (int i = 0; i < 3 * data->natom; i++) {
    data->O_pX[i] -= data->PU[i] * *O_p - *O_p * data->PU[i];
  }
}

void DKH2Analytic::store_mat(shared_ptr<const VectorB> vec, shared_ptr<DKH2AnalyticData> data) const {
  auto mat = make_shared<Matrix>(data->nunc, data->nunc);
  for (int i = 0; i < data->nunc; i++) {
    (*mat)(i, i) = (*vec)(i);
  }
  data->vec2mat[vec] = mat;
}

