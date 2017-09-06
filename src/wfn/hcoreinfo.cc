//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: hcoreinfo.cc
// Copyright (C) 2017 Toru Shiozaki
//
// Author: Jae Woo Park <jwpk1201@northwestern.edu>
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


#include <src/mat1e/hcore.h>
#include <src/mat1e/dkhcore.h>
#include <src/wfn/hcoreinfo.h>

using namespace std;
using namespace bagel;


HcoreInfo::HcoreInfo(shared_ptr<const PTree> idata) : type_(HcoreType::standard) {
  // DKH
  const bool dkh = idata->get<bool>("dkh", false);
  if (dkh)
    type_ = HcoreType::dkh;
  mat1e_dx_ = idata->get<double>("mat1e_dx", 0.001);
  gradtype_ = idata->get<bool>("gradtype", false);

  // ECP
  const string basisfile = idata->get<string>("basis", "");
  const bool ecp = basisfile.find("ecp") != string::npos;
  if (ecp) {
    if (dkh)
      throw runtime_error("DKH and ECP cannot be used simultaneously");
    type_ = HcoreType::ecp;
  }
} 


vector<shared_ptr<Matrix>> HcoreInfo::dkh_grad(shared_ptr<const Molecule> mol) const {
  int natom = current->natom();
  vector<shared_ptr<Matrix>> dkhgrad;

  for (int i = 0; i != natom; ++i) {
    for (int j = 0; j != 3; ++j) {
      shared_ptr<Matrix> h_plus;
      {
        auto displ = make_shared<XYZFile>(natom);
        displ->element(j,i) = mat1e_dx();
        auto geom_plus = make_shared<Molecule>(*current, displ, false);
        shared_ptr<Matrix> hd_plus = compute_dkh(geom_plus);
        auto ho_plus = make_shared<Hcore>(geom_plus);

        h_plus = make_shared<Matrix>(*hd_plus - *ho_plus);
      }

      shared_ptr<Matrix> h_minus;
      {
        auto displ = make_shared<XYZFile>(natom);
        displ->element(j,i) = -mat1e_dx();
        auto geom_minus = make_shared<Molecule>(*current, displ, false);
        shared_ptr<Matrix> hd_minus = compute_dkh(geom_minus);
        auto ho_minus = make_shared<Hcore>(geom_minus);

        h_minus = make_shared<Matrix>(*hd_minus - *ho_minus);
      }

      dkhgrad.push_back(make_shared<Matrix>(*h_plus - *h_minus));
      dkhgrad[j+i*3]->scale(1.0 / (2.0 * mat1e_dx()));
    }
  }

  return dkhgrad;
}

vector<shared_ptr<Matrix>> HcoreInfo::dkh_analyticgrad(shared_ptr<const Geometry> mol) {
  gradinit(mol);

  vector<shared_ptr<Matrix>> dkh2grad(3 * natom);
  auto s = make_shared<Overlap>(molu);
  contracts(mol, s);
  overlapgrad(mol, s);
  const Matrix s_inv12 = *s->tildex();
  Matrix s_inv = *s;
  s_inv.inverse_symmetric();
  auto T_p = make_shared<Kinetic>(molu);
  kineticgrad(mol, T_p);
  const Matrix T_pp = s_inv12 % *T_p * s_inv12;
  auto t = make_shared<VectorB>(nunc);
  store_mat(t);
  Matrix W = T_pp;
  W.diagonalize(*t);

  const double c2 = c__ * c__;
  auto Ep = make_shared<VectorB>(nunc);
  store_mat(Ep);
  auto A = make_shared<VectorB>(nunc);
  store_mat(A);
  auto K = make_shared<VectorB>(nunc);
  store_mat(K);
  auto B = make_shared<VectorB>(nunc);
  store_mat(B);
  auto R = make_shared<VectorB>(nunc);
  store_mat(R);
  auto R_inv = make_shared<VectorB>(nunc);
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
    const Matrix T_ppX = s_inv12 % T_pX[i] * s_inv12 - 0.5 * (s_inv * s_X[i] * T_pp + T_pp * s_inv * s_X[i]);
    const Matrix T_ppX_W = W % T_ppX * W;
    Matrix PW(nunc, nunc);
    for (int k = 0; k < nunc; ++k) {
      for (int l = 0; l < nunc; ++l) {
        PW(k, l) = k == l ? 0 : T_ppX_W(k, l) / ((*t)(k) - (*t)(l));
      }
    }
    const Matrix t_X = T_ppX_W - PW * *vec2mat[t] + *vec2mat[t] * PW;

    auto Ep_X = make_shared<VectorB>(nunc);
    store_mat(Ep_X);
    auto A_X = make_shared<VectorB>(nunc);
    store_mat(A_X);
    auto B_X = make_shared<VectorB>(nunc);
    store_mat(B_X);
    auto R_X = make_shared<VectorB>(nunc);
    store_mat(R_X);
    auto R_invX = make_shared<VectorB>(nunc);
    store_mat(R_invX);
    for (int k = 0; k < nunc; ++k) {
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
    
    auto t_rel = make_shared<VectorB>(nunc);
    store_mat(t_rel);
    for (int k = 0; k < nunc; ++k) {
      (*t_rel)(k) = c__ * std::sqrt(2 * (*t)(k) + c2) - c2;
    }
    const Matrix T_pprel = W * *vec2mat[t_rel] ^ W;
    const Matrix T_prel = s_inv12 * T_pprel ^ s_inv12;
    const Matrix T_rel = U_T % T_prel * U_T;
    const Matrix DW = W * PW ^ W;
    const Matrix T_pprelX = W * *vec2mat[Ep] ^ W + DW * T_pprel - T_pprel * DW;
    const Matrix T_prelX = s_inv12 * T_pprelX ^ s_inv12 + 0.5 * (s_inv * s_X[i] * T_prel + T_prel * s_inv * s_X[i]);
    const Matrix DU = U_T % PU[i] * U_T;
    const Matrix T_relX = U_T % T_prelX * U_T + DU * T_rel - T_rel * DU;

    auto V_p = make_shared<NAI>(molu);
    naigrad(mol, V_p);
    const Matrix V_pp = s_inv12 % *V_p * s_inv12;
    const Matrix V_ppp = W % V_pp * W;
    const Matrix V_ppX = s_inv12 % V_pX[i] * s_inv12 - 0.5 * (s_inv * s_X[i] * V_pp + V_pp * s_inv * s_X[i]);
    const Matrix V_pppX = W % V_ppX * W - PW * V_ppp + V_ppp * PW;

    const Small1e<NAIBatch> small1e(molu);
    auto O_p = make_shared<Matrix>(small1e[0]);
    smallnaigrad(mol, O_p);
    const Matrix O_pp = s_inv12 % *O_p * s_inv12;
    const Matrix O_ppp = W % O_pp * W;
    const Matrix O_ppX = s_inv12 % O_pX[i] * s_inv12 - 0.5 * (s_inv * s_X[i] * O_pp + O_pp * s_inv * s_X[i]);
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
        V_longX(k, l) = -1 * pow((*Ep)(k) + (*Ep)(l), -2) * ((*Ep_X)(k) + (*Ep_X)(l)) * V_ppp(k, l) + V_pppX(k, l) / ((*Ep)(k) + (*Ep)(l));
        A_longX(k, l) = (*A)(k) * V_longX(k, l) * (*A)(l) + (*A_X)(k) * V_long(k, l) * (*A)(l) + (*A)(k) * V_long(k, l) * (*A_X)(l);
        O_longX(k, l) = -1 * pow((*Ep)(k) + (*Ep)(l), -2) * ((*Ep_X)(k) + (*Ep_X)(l)) * O_ppp(k, l) + O_pppX(k, l) / ((*Ep)(k) + (*Ep)(l));
        B_longX(k, l) = (*B)(k) * O_longX(k, l) * (*B)(l) + (*B_X)(k) * O_long(k, l) * (*B)(l) + (*B)(k) * O_long(k, l) * (*B_X)(l);
      }
    }

    auto EpR = make_shared<VectorB>(*Ep % *R);
    store_mat(EpR);
    auto EpRX = make_shared<VectorB>(*Ep_X % *R + *Ep % *R_X);
    store_mat(EpRX);
    auto EpRinv = make_shared<VectorB>(*Ep % *R_inv);
    store_mat(EpRinv);
    auto EpRinvX = make_shared<VectorB>(*Ep_X % *R_inv + *Ep % *R_invX);
    store_mat(EpRinvX);
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

    Matrix V_ppprel = *vec2mat[A] * V_ppp * *vec2mat[A] + *vec2mat[B] * O_ppp * *vec2mat[B]
                    - B_long * *vec2mat[Ep] * A_long - A_long * *vec2mat[Ep] * B_long
                    + A_long * *vec2mat[EpR] * A_long + B_long * *vec2mat[EpRinv] * B_long
                    - 0.5 * (B_long * A_long * *vec2mat[Ep] + A_long * B_long * *vec2mat[Ep])
                    + 0.5 * (A_long * *vec2mat[R] * A_long * *vec2mat[Ep] + B_long * *vec2mat[R_inv] * B_long * *vec2mat[Ep])
                    - 0.5 * (*vec2mat[Ep] * B_long * A_long + *vec2mat[Ep] * A_long * B_long)
                    + 0.5 * (*vec2mat[Ep] * A_long * *vec2mat[R] * A_long + *vec2mat[Ep] * B_long * *vec2mat[R_inv] * B_long);
    
    const Matrix V_pprel = W * V_ppprel ^ W;
    const Matrix V_prel = s_inv12 * V_pprel ^ s_inv12;
    const Matrix V_rel = U_T % V_prel * U_T;
    const Matrix V_pprelX = W * V_ppprelX ^ W + DW * V_pprel - V_pprel * DW;
    const Matrix V_prelX = s_inv12 * V_pprelX ^ s_inv12 + 0.5 * (s_inv * s_X[i] * V_prel + V_prel * s_inv * s_X[i]);
    const Matrix V_relX = U_T % V_prelX * U_T + DU * V_rel - V_rel * DU;

    *dkh2grad[i] = T_relX + V_relX;
  }
  return dkh2grad;
}

void HcoreInfo::gradinit(shared_ptr<const Geometry> mol) {
  natom = geom->natom();
  nbasis = geom->nbasis();
  molu = make_shared<Molecule>(*mol->uncontract());
  nunc = molu->nbasis();
  U_T = MixedBasis<OverlapBatch>(mol, molu);

  PU = vector<Matrix>(3 * natom, Matrix(nunc, nunc));
  s_X = vector<Matrix>(3 * natom, Matrix(nunc, nunc));
  T_pX = vector<Matrix>(3 * natom, Matrix(nunc, nunc));
  V_pX = vector<Matrix>(3 * natom, Matrix(nunc, nunc));
  O_pX = vector<Matrix>(3 * natom, Matrix(nunc, nunc));

  id = make_shared<Matrix>(nunc, nunc);
  for (int i = 0; i < nunc; i++) {
    (*id)(i, i) = 1;
  }
}

void HcoreInfo::contracts(shared_ptr<const Geometry> mol, shared_ptr<const Matrix> s) {
  const Matrix U = U_T % id;
  for (int k = 0; k < nunc; k++) {
    for (int l = 0; l < nunc; l++) {
      if (k == l) {
        for (int i = 0; i < 3 * natom; i++) {
          PU[i](k, l) = 0;
        }
      }
      else {
        auto grad = make_shared<GradFile>(natom);
        auto den = make_shared<Matrix>(nbasis, nbasis);
        for (int m = 0; m < nbasis; k++) {
          for (int n = 0; n < nbasis; l++) {
            (*den)(m, n) = U_T(k, m) * U(n, l) / ((*s)(k, k) - (*s)(l, l));
          }
        }
        GradEval_base ge(mol);
        grad = ge.contract_overlapgrad(den);
        for (int i = 0; i < natom; i++) {
          for (int j = 0; j < 3; j++) {
            PU[3 * i + j](k, l) = (*grad)(i, j);
          }
        }
      }
    }
  }
}

void HcoreInfo::overlapgrad(shared_ptr<const Geometry> mol, shared_ptr<const Matrix> s) {
  const Matrix U = U_T % id;
  for (int k = 0; k < nunc; k++) {
    for (int l = 0; l < nunc; l++) {
      auto grad = make_shared<GradFile>(natom);
      auto den = make_shared<Matrix>(nbasis, nbasis);
      for (int m = 0; m < nbasis; k++) {
        for (int n = 0; n < nbasis; l++) {
          (*den)(m, n) = U_T(k, m) * U(n, l);
        }
      }
      GradEval_base ge(mol);
      grad = ge.contract_naigrad(den);
      for (int i = 0; i < natom; i++) {
        for (int j = 0; j < 3; j++) {
          s_X[3 * i + j](k, l) = (*grad)(i, j);
        }
      }
    }
  }
  for (int i = 0; i < 3 * natom; i++) {
    s_X[i] -= PU[i] * s - s * PU[i];
  }
}

void HcoreInfo::kineticgrad(shared_ptr<const Geometry> mol, shared_ptr<const Matrix> T_p) {
  const Matrix U = U_T % id;
  for (int k = 0; k < nunc; k++) {
    for (int l = 0; l < nunc; l++) {
      auto grad = make_shared<GradFile>(natom);
      auto den = make_shared<Matrix>(nbasis, nbasis);
      for (int m = 0; m < nbasis; k++) {
        for (int n = 0; n < nbasis; l++) {
          (*den)(m, n) = U_T(k, m) * U(n, l);
        }
      }
      GradEval_base ge(mol);
      grad = ge.contract_kineticgrad(den);
      for (int i = 0; i < natom; i++) {
        for (int j = 0; j < 3; j++) {
          T_pX[3 * i + j](k, l) = (*grad)(i, j);
        }
      }
    }
  }
  for (int i = 0; i < 3 * natom; i++) {
    T_pX[i] -= PU[i] * T_p - T_p * PU[i];
  }
}

void HcoreInfo::naigrad(shared_ptr<const Geometry> mol, shared_ptr<const Matrix> V_p) {
  const Matrix U = U_T % id;
  for (int k = 0; k < nunc; k++) {
    for (int l = 0; l < nunc; l++) {
      auto grad = make_shared<GradFile>(natom);
      auto den = make_shared<Matrix>(nbasis, nbasis);
      for (int m = 0; m < nbasis; k++) {
        for (int n = 0; n < nbasis; l++) {
          (*den)(m, n) = U_T(k, m) * U(n, l);
        }
      }
      GradEval_base ge(mol);
      grad = ge.contract_naigrad(den);
      for (int i = 0; i < natom; i++) {
        for (int j = 0; j < 3; j++) {
          V_pX[3 * i + j](k, l) = (*grad)(i, j);
        }
      }
    }
  }
  for (int i = 0; i < 3 * natom; i++) {
    V_pX[i] -= PU[i] * V_p - V_p * PU[i];
  }
}

void HcoreInfo::smallnaigrad(shared_ptr<const Geometry> mol, shared_ptr<const Matrix> O_p) {
  const Matrix U = U_T % id;
  for (int k = 0; k < nunc; k++) {
    for (int l = 0; l < nunc; l++) {
      auto grad = make_shared<GradFile>(natom);
      auto den = make_shared<Matrix>(nbasis, nbasis);
      for (int m = 0; m < nbasis; k++) {
        for (int n = 0; n < nbasis; l++) {
          (*den)(m, n) = U_T(k, m) * U(n, l);
        }
      }
      GradEval_base ge(mol);
      grad = ge.contract_gradsmall1e({den, den, den, den, den, den});
      for (int i = 0; i < natom; i++) {
        for (int j = 0; j < 3; j++) {
          O_pX[3 * i + j](k, l) = (*grad)(i, j);
        }
      }
    }
  }
  for (int i = 0; i < 3 * natom; i++) {
    O_pX[i] -= PU[i] * O_p - O_p * PU[i];
  }
}

void HcoreInfo::store_mat(shared_ptr<const VectorB> vec) {
  auto mat = make_shared<Matrix>(nunc, nunc);
  for (int i = 0; i < nunc; i++) {
    (*mat)(i, i) = *vec(i);
  }
  vec2mat[vec] = mat;
}

shared_ptr<Matrix> HcoreInfo::compute_grad_dkh(shared_ptr<const Molecule> current, shared_ptr<const Matrix> den) const {
  int natom = current->natom();
  auto out = make_shared<Matrix>(3,natom);
  vector<shared_ptr<Matrix>> dkhg = gradtype_ ? dkh_grad(current) : dkh_analyticgrad(dynamic_cast<const Geometry>(current));

  for (int i = 0; i != natom; ++i)
    for (int j = 0; j != 3; ++j)
      out->element(j,i) += dkhg[j+i*3]->dot_product(den);

  return out;
}


shared_ptr<Matrix> HcoreInfo::compute_grad(shared_ptr<const Molecule> current, shared_ptr<const Matrix> den) const {
  int natom = current->natom();
  auto out = make_shared<Matrix>(3, natom);

  if (dkh())
    out = compute_grad_dkh(current, den);

  return out;
}


shared_ptr<Matrix> HcoreInfo::compute_dkh(shared_ptr<const Molecule> current) const {
  auto out = make_shared<DKHcore>(current);

  return out;
}


shared_ptr<Matrix> HcoreInfo::compute(shared_ptr<const Molecule> current) const {
  shared_ptr<Matrix> out;

  if (dkh())
    out = compute_dkh(current);

  return out;
}


void HcoreInfo::print() const {
  if (dkh())
    cout << "      - Using DKHcore" << endl;
}
