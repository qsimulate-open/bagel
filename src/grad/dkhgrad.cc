//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dkhgrad.cc
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


#include <src/grad/dkhgrad.h>
#include <src/mat1e/kinetic.h>
#include <src/mat1e/nai.h>
#include <src/mat1e/rel/small1e.h>
#include <src/mat1e/overlap.h>
#include <src/wfn/contractmat.h>

using namespace std;
using namespace bagel;


array<shared_ptr<const Matrix>, 4> DKHgrad::compute(shared_ptr<const Matrix> rdm1, shared_ptr<const Matrix> erdm1) {
  // Gradient integrals and RDMs to be computed in uncontracted basis
  auto unc = make_shared<Molecule>(*mol_);
  unc = unc->uncontract();
  // Number of uncontracted basis functions
  const int nbasis = unc->nbasis();

  shared_ptr<const Matrix> wmat, wmat_rev;
  shared_ptr<const DiagMatrix> t;
  {
    // Compute momentum transformation matrix
    const Overlap overlap(unc);
    shared_ptr<const Matrix> gamma = overlap.tildex();
    // Kinetic operator is diagonal in momentum space
    const Kinetic kinetic(unc);
    auto lambda = make_shared<Matrix>(*gamma % kinetic * *gamma);
    VectorB t0(nbasis);
    lambda->diagonalize(t0);
    t = make_shared<const DiagMatrix>(t0);
    // Transformation matrix to momentum space
    wmat = make_shared<const Matrix>(*gamma * *lambda);
    // Reverse transformation from momentum space
    wmat_rev = make_shared<const Matrix>(overlap * *wmat);
  }

  shared_ptr<const Matrix> v, pvp;
  {
    // Evaluate integrals in momentum space
    const NAI nai(unc);
    v = make_shared<const Matrix>(*wmat % nai * *wmat);
    const Small1e<NAIBatch> small1e(unc);
    pvp = make_shared<const Matrix>(*wmat % small1e[0] * *wmat);
  }

  // Projection operator from uncontracted to contracted AO basis
  auto pmat = make_shared<const ContractMat>(mol_, nbasis);

  // Lagrange multiplier for diagonalization of kinetic operator and energy derivative with respect to momentum orbital rotation
  shared_ptr<const Matrix> tden, zpq, ypq;
  tie(tden, zpq, ypq) = compute_tden_(rdm1, pmat, wmat, wmat_rev, t, v, pvp); 

  shared_ptr<const Matrix> sden = compute_sden_(rdm1, erdm1, pmat, wmat, wmat_rev, t, v, pvp, zpq, ypq);
  shared_ptr<const Matrix> vden0, vden1;
  tie(vden0, vden1) = compute_vden_(rdm1, pmat, wmat, wmat_rev, t, v, pvp);
  return { tden, vden0, vden1, sden };
}


// Effective density matrix for kinetic gradient. This function also returns zpq and ypq
tuple<shared_ptr<const Matrix>,shared_ptr<const Matrix>,shared_ptr<const Matrix>>
  DKHgrad::compute_tden_(shared_ptr<const Matrix> rdm1, shared_ptr<const Matrix> pmat, shared_ptr<const Matrix> wmat, shared_ptr<const Matrix> wmat_rev,
                         shared_ptr<const DiagMatrix> t, shared_ptr<const Matrix> v, shared_ptr<const Matrix> pvp) {
  const int nbasis = t->ndim();

  // output
  auto den = make_shared<Matrix>(nbasis, nbasis);
  auto zpq = make_shared<Matrix>(nbasis, nbasis);
  auto ypq = make_shared<Matrix>(nbasis, nbasis);

  const double c2 = c__ * c__;
  DiagMatrix E(nbasis);
  DiagMatrix A(nbasis);
  DiagMatrix B(nbasis);
  DiagMatrix K(nbasis);
  DiagMatrix dE(nbasis);
  DiagMatrix dA(nbasis);
  DiagMatrix dK(nbasis);
  DiagMatrix dB(nbasis);
  for (int p = 0; p != nbasis; ++p) {
    E(p) = c__ * sqrt(2.0 * (*t)(p) + c2);
    A(p) = sqrt((c2 + E(p)) / (2.0 * E(p)));
    K(p) = c__ / (E(p) + c2);
    B(p) = A(p) * K(p);
    dE(p) = c__ / sqrt(2.0 * (*t)(p) + c2);
    dA(p) = -c2 * dE(p) / (4.0 * pow(E(p), 2) * A(p));
    dK(p) = -c__ * dE(p) / pow(E(p) + c2, 2);
    dB(p) = A(p) * dK(p) + dA(p) * K(p);
  }

  const Matrix CPW = (*pmat % *wmat_rev) % *rdm1 * (*pmat % *wmat_rev);
  {
    const Matrix CAN = CPW * A * *v;
    const Matrix NAC = *v * A * CPW;
    const Matrix CBS = CPW * B * *pvp;
    const Matrix SBC = *pvp * B * CPW;
    // DKH1 correction
    *ypq += 2.0 * (CPW * E - c2 * CPW + (CAN + NAC) * A + (CBS + SBC) * B);

    DiagMatrix EC(nbasis);
    for (int p = 0; p != nbasis; ++p) {
      // dE_DKH1/dU_pq
      (*ypq)(p, p) += 2.0 * (CPW(p, p) * dE(p) + ((CAN(p, p) + NAC(p, p)) * dA(p) + (CBS(p, p) + SBC(p, p)) * dB(p))) * (*t)(p);
      EC(p) = dE(p) * CPW(p, p) + 2.0 * (NAC(p, p) * dA(p) + SBC(p, p) * dB(p));
    }
    *den += *wmat * EC ^ *wmat;
  }

  Matrix N(nbasis, nbasis);
  Matrix O(nbasis, nbasis);
  DiagMatrix F(nbasis);
  DiagMatrix G(nbasis);
  DiagMatrix H(nbasis);
  DiagMatrix dF(nbasis);
  DiagMatrix dG(nbasis);
  DiagMatrix dH(nbasis);

  // DKH2 corrections, one for each term
  for (int i = 0; i != 12; ++i) {
    switch (i) {
    case 0:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = B(p);
        G(p) = A(p);
        H(p) = -B(p) * E(p) * A(p);
        dF(p) = dB(p);
        dG(p) = dA(p);
        dH(p) = -dB(p) * E(p) * A(p) - B(p) * dE(p) * A(p) - B(p) * E(p) * dA(p);
      }
      N = *pvp;
      O = *v;
      break;
    case 1:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = A(p);
        G(p) = B(p);
        H(p) = -A(p) * E(p) * B(p);
        dF(p) = dA(p);
        dG(p) = dB(p);
        dH(p) = -dA(p) * E(p) * B(p) - A(p) * dE(p) * B(p) - A(p) * E(p) * dB(p);
      }
      N = *v;
      O = *pvp;
      break;
    case 2:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = G(p) = A(p);
        H(p) = 2.0 * pow(A(p) * K(p), 2) * (*t)(p) * E(p);
        dF(p) = dG(p) = dA(p);
        dH(p) = 4.0 * A(p) * dA(p) * (*t)(p) * pow(K(p), 2) * E(p) + 2.0 * pow(A(p) * K(p), 2) * E(p)
              + 4.0 * pow(A(p), 2) * (*t)(p) * K(p) * dK(p) * E(p) + 2.0 * pow(A(p) * K(p), 2) * (*t)(p) * dE(p);
      }
      N = O = *v;
      break;
    case 3:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = G(p) = B(p);
        H(p) = pow(B(p) / K(p), 2) * E(p) / (2.0 * (*t)(p));
        dF(p) = dG(p) = dB(p);
        dH(p) = ((2.0 * B(p) * dB(p) * E(p) + pow(B(p), 2) * dE(p)) * (*t)(p) * K(p)
              - pow(B(p), 2) * E(p) * (K(p) + 2.0 * (*t)(p) * dK(p))) / (2.0 * pow((*t)(p) * K(p), 2) * K(p));
      }
      N = O = *pvp;
      break;
    case 4:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = B(p);
        G(p) = A(p) * E(p);
        H(p) = -B(p) * A(p) / 2.0;
        dF(p) = dB(p);
        dG(p) = dA(p) * E(p) + A(p) * dE(p);
        dH(p) = -dB(p) * A(p) / 2.0 - B(p) * dA(p) / 2.0;
      }
      N = *pvp;
      O = *v;
      break;
    case 5:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = A(p);
        G(p) = B(p) * E(p);
        H(p) = -A(p) * B(p) / 2.0;
        dF(p) = dA(p);
        dG(p) = dB(p) * E(p) + B(p) * dE(p);
        dH(p) = -dA(p) * B(p) / 2.0 - A(p) * dB(p) / 2.0;
      }
      N = *v;
      O = *pvp;
      break;
    case 6:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = A(p);
        G(p) = A(p) * E(p);
        H(p) = pow(A(p) * K(p), 2) * (*t)(p);
        dF(p) = dA(p);
        dG(p) = dA(p) * E(p) + A(p) * dE(p);
        dH(p) = 2.0 * A(p) * dA(p) * (*t)(p) * pow(K(p), 2) + pow(A(p) * K(p), 2)
              + 2.0 * pow(A(p), 2) * (*t)(p) * K(p) * dK(p);
      }
      N = O = *v;
      break;
    case 7:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = B(p);
        G(p) = B(p) * E(p);
        H(p) = pow(B(p) / K(p), 2) / (4.0 * (*t)(p));
        dF(p) = dB(p);
        dG(p) = dB(p) * E(p) + B(p) * dE(p);
        dH(p) = (2.0 * B(p) * dB(p) * (*t)(p) * K(p)
              - pow(B(p), 2) * (K(p) + 2.0 * (*t)(p) * dK(p))) / (4.0 * pow((*t)(p) * K(p), 2) * K(p));
      }
      N = O = *pvp;
      break;
    case 8:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = E(p) * B(p);
        G(p) = A(p);
        H(p) = -B(p) * A(p) / 2.0;
        dF(p) = dE(p) * B(p) + E(p) * dB(p);
        dG(p) = dA(p);
        dH(p) = -dB(p) * A(p) / 2.0 - B(p) * dA(p) / 2.0;
      }
      N = *pvp;
      O = *v;
      break;
    case 9:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = E(p) * A(p);
        G(p) = B(p);
        H(p) = -A(p) * B(p) / 2.0;
        dF(p) = dE(p) * A(p) + E(p) * dA(p);
        dG(p) = dB(p);
        dH(p) = -dA(p) * B(p) / 2.0 - A(p) * dB(p) / 2.0;
      }
      N = *v;
      O = *pvp;
      break;
    case 10:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = E(p) * A(p);
        G(p) = A(p);
        H(p) = pow(A(p) * K(p), 2) * (*t)(p);
        dF(p) = dE(p) * A(p) + E(p) * dA(p);
        dG(p) = dA(p);
        dH(p) = 2.0 * A(p) * dA(p) * (*t)(p) * pow(K(p), 2) + pow(A(p) * K(p), 2)
              + 2.0 * pow(A(p), 2) * (*t)(p) * K(p) * dK(p);
      }
      N = O = *v;
      break;
    case 11:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = E(p) * B(p);
        G(p) = B(p);
        H(p) = pow(B(p) / K(p), 2) / (4.0 * (*t)(p));
        dF(p) = dE(p) * B(p) + E(p) * dB(p);
        dG(p) = dB(p);
        dH(p) = (2.0 * B(p) * dB(p) * (*t)(p) * K(p)
              - pow(B(p), 2) * (K(p) + 2.0 * (*t)(p) * dK(p))) / (4.0 * pow((*t)(p) * K(p), 2) * K(p));
      }
      N = O = *pvp;
      break;
    }
    Matrix WN(nbasis, nbasis);
    Matrix WO(nbasis, nbasis);
    Matrix WWN(nbasis, nbasis);
    Matrix WWO(nbasis, nbasis);
    for (int p = 0; p != nbasis; ++p) {
      for (int q = 0; q != nbasis; ++q) {
        WN(q, p) = N(q, p) / (E(p) + E(q));
        WO(q, p) = O(q, p) / (E(p) + E(q));
        WWN(q, p) = WN(q, p) / (E(p) + E(q));
        WWO(q, p) = WO(q, p) / (E(p) + E(q));
      }
    }

    Matrix CH(nbasis, nbasis);
    Matrix CH2(nbasis, nbasis);
    for (int p = 0; p != nbasis; ++p) {
      for (int q = 0; q != nbasis; ++q) {
        if (p != q) {
          CH(q, p) = CPW(q, p) * H(p);
          CH2(q, p) = CPW(q, p) * dH(p);
        }
      }
    }

    const Matrix CFN = CPW * F * WN;
    const Matrix CGO = CPW * G * WO;
    {
      const Matrix HNFC = H * WN * F * CPW * G;
      const Matrix HOGC = H * WO * G * CPW * F;
      Matrix WHNFC(nbasis, nbasis);
      Matrix WHOGC(nbasis, nbasis);
      for (int p = 0; p != nbasis; ++p) {
        for (int q = 0; q != nbasis; ++q) {
          WHNFC(q, p) = (HNFC(q, p) + CH(q, p) * G(q) * F(p) * WN(p, p)) / (E(p) + E(q));
          WHOGC(q, p) = (HOGC(q, p) + CH(q, p) * F(q) * G(p) * WO(p, p)) / (E(p) + E(q));
        }
        WHNFC(p, p) += H(p) * CFN(p, p) * G(p) / (2.0 * E(p));
        WHOGC(p, p) += H(p) * CGO(p, p) * F(p) / (2.0 * E(p));
      }
      *ypq += O * WHNFC + N * WHOGC + CFN * H * WO * G + CGO * H * WN * F;
    }

    {
      Matrix FN(nbasis, nbasis);
      Matrix FO(nbasis, nbasis);
      Matrix FNN(nbasis, nbasis);
      Matrix FOO(nbasis, nbasis);
      for (int p = 0; p != nbasis; ++p)
        for (int q = 0; q != nbasis; ++q)
          if (p != q) {
            FN(q, p) = WN(q, p) * F(q);
            FO(q, p) = WO(q, p) * F(q);
            FNN(q, p) = WWN(q, p) * F(q);
            FOO(q, p) = WWO(q, p) * F(q);
          }
      FN = CPW * FN;
      FO = CPW * FO;
      FNN = CPW * FNN;
      FOO = CPW * FOO;

      Matrix GCFN2H(nbasis, nbasis);
      Matrix GCFO2H(nbasis, nbasis);
      Matrix GCFN2H2(nbasis, nbasis);
      Matrix GCFO2H2(nbasis, nbasis);
      for (int p = 0; p != nbasis; ++p)
        for (int q = 0; q != nbasis; ++q)
          if (p != q) {
            const double fac1 = G(q) * H(p) / (E(q) + E(p));
            const double fac2 = G(q) * dH(p) * (*t)(p);
            const double fac3 = G(q) * H(p) * dE(p) * (*t)(p);
            GCFN2H(q, p) = fac1 * FN(q, p);
            GCFO2H(q, p) = fac1 * FO(q, p);
            GCFN2H2(q, p) = fac2 * FN(q, p) - 2.0 * fac3 * FNN(q, p);
            GCFO2H2(q, p) = fac2 * FO(q, p) - 2.0 * fac3 * FOO(q, p);
          }

      // dE_DKH2/dU_pq
      *ypq += O * GCFN2H + N * GCFO2H;

      const Matrix OGCFN2H2 = WO * GCFN2H2 + WN * GCFO2H2;
      for (int p = 0; p != nbasis; ++p)
        (*ypq)(p, p) += OGCFN2H2(p, p);
    }

    {
      const Matrix CFNHOG = (CFN * H * (WO * dG - WWO * G * dE) + CGO * H * (WN * dF - WWN * F * dE)) * *t;
      const Matrix CFNGH = (CFN * G * dH - (CPW * F * WWN * G * H + WWN * F * CH * G) * dE + WN * F * CH2 * G) * *t;
      const Matrix CFNGHE = (CFN * G * H + WN * F * CH * G) * dE * *t;
      const Matrix CGOFH = (CGO * F * dH - (CPW * G * WWO * F * H + WWO * G * CH * F) * dE + WO * G * CH2 * F) * *t;
      const Matrix CGOFHE = (CGO * F * H + WO * G * CH * F) * dE * *t;
      for (int p = 0; p != nbasis; ++p)
        (*ypq)(p, p) += 2.0 * CFNHOG(p, p) + CFNGH(p, p) * WO(p, p) - CFNGHE(p, p) * WWO(p, p)
                      + CGOFH(p, p) * WN(p, p) - CGOFHE(p, p) * WWN(p, p);
    }

    DiagMatrix dkh2(nbasis);
    {
      const Matrix OGCFN = WO * G * CPW * F * WN;
      const Matrix OGCFNN = WO * G * CPW * F * WWN;
      const Matrix OOGCFN = WWO * G * CPW * F * WN;
      const Matrix OHNF = WO * H * (WN * dF - WWN * F * dE);
      const Matrix NHOG = WN * H * (WO * dG - WWO * G * dE);
      Matrix OHNFC(nbasis, nbasis);
      Matrix NHOGC(nbasis, nbasis);
      for (int p = 0; p != nbasis; ++p) {
        for (int q = 0; q != nbasis; ++q) {
          OHNFC(q, p) = OHNF(q, p) * CPW(q, p);
          NHOGC(q, p) = NHOG(q, p) * CPW(q, p);
        }
        dkh2(p) = OGCFN(p, p) * dH(p) - OGCFNN(p, p) * H(p) * dE(p) - OOGCFN(p, p) * H(p) * dE(p);
      }
      dkh2 += OHNFC * G.diag() + NHOGC * F.diag();
    }

    *den += *wmat * dkh2 ^ *wmat;
  }

  // z_pq from Lagrangian
  for (int p = 0; p != nbasis; ++p)
    for (int q = 0; q != nbasis; ++q)
      (*zpq)(q, p) = fabs((*t)(p) - (*t)(q)) > 1.0e-12 ? -0.5 * ((*ypq)(q, p) - (*ypq)(p, q)) / ((*t)(q) - (*t)(p)) : 0.0;

  *den += *wmat * *zpq ^ *wmat;
  return {den, zpq, ypq};
}


// Effective density matrices for NAI/SmallNAI gradients
tuple<shared_ptr<const Matrix>, shared_ptr<const Matrix>>
  DKHgrad::compute_vden_(shared_ptr<const Matrix> rdm1, shared_ptr<const Matrix> pmat, shared_ptr<const Matrix> wmat, shared_ptr<const Matrix> wmat_rev,
                         shared_ptr<const DiagMatrix> t, shared_ptr<const Matrix> v, shared_ptr<const Matrix> pvp) {
  const int nbasis = t->ndim();
  const double c2 = c__ * c__;
  DiagMatrix E(nbasis);
  DiagMatrix A(nbasis);
  DiagMatrix K(nbasis);
  DiagMatrix B(nbasis);
  for (int p = 0; p != nbasis; ++p) {
    E(p) = c__ * sqrt(2.0 * (*t)(p) + c2);
    A(p) = sqrt((E(p) + c2) / (2.0 * E(p)));
    K(p) = c__ / (E(p) + c2);
    B(p) = A(p) * K(p);
  }

  const Matrix CPW = (*pmat % *wmat_rev) % *rdm1 * (*pmat % *wmat_rev);
  auto den = make_shared<Matrix>(nbasis, nbasis);
  auto pvpden = make_shared<Matrix>(nbasis, nbasis);
  // DKH1 correction
  *den += *wmat * A * CPW * A ^ *wmat;
  *pvpden += *wmat * B * CPW * B ^ *wmat;

  Matrix N(nbasis, nbasis);
  Matrix O(nbasis, nbasis);
  DiagMatrix F(nbasis);
  DiagMatrix G(nbasis);
  DiagMatrix H(nbasis);
  pair<int, int> vint;
  // DKH2 corrections, one for each term
  for (int i = 0; i != 12; ++i) {
    switch (i) {
    case 0:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = B(p);
        G(p) = A(p);
        H(p) = -B(p) * E(p) * A(p);
      }
      N = *pvp;
      O = *v;
      vint = make_pair(1, 0);
      break;
    case 1:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = A(p);
        G(p) = B(p);
        H(p) = -A(p) * E(p) * B(p);
      }
      N = *v;
      O = *pvp;
      vint = make_pair(0, 1);
      break;
    case 2:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = G(p) = A(p);
        H(p) = 2.0 * pow(A(p) * K(p), 2) * (*t)(p) * E(p);
      }
      N = O = *v;
      vint = make_pair(0, 0);
      break;
    case 3:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = G(p) = B(p);
        H(p) = pow(B(p) / K(p), 2) * E(p) / (2.0 * (*t)(p));
      }
      N = O = *pvp;
      vint = make_pair(1, 1);
      break;
    case 4:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = B(p);
        G(p) = A(p) * E(p);
        H(p) = -B(p) * A(p) / 2.0;
      }
      N = *pvp;
      O = *v;
      vint = make_pair(1, 0);
      break;
    case 5:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = A(p);
        G(p) = B(p) * E(p);
        H(p) = -A(p) * B(p) / 2.0;
      }
      N = *v;
      O = *pvp;
      vint = make_pair(0, 1);
      break;
    case 6:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = A(p);
        G(p) = A(p) * E(p);
        H(p) = pow(A(p) * K(p), 2) * (*t)(p);
      }
      N = O = *v;
      vint = make_pair(0, 0);
      break;
    case 7:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = B(p);
        G(p) = B(p) * E(p);
        H(p) = pow(B(p) / K(p), 2) / (4.0 * (*t)(p));
      }
      N = O = *pvp;
      vint = make_pair(1, 1);
      break;
    case 8:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = E(p) * B(p);
        G(p) = A(p);
        H(p) = -B(p) * A(p) / 2.0;
      }
      N = *pvp;
      O = *v;
      vint = make_pair(1, 0);
      break;
    case 9:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = E(p) * A(p);
        G(p) = B(p);
        H(p) = -A(p) * B(p) / 2.0;
      }
      N = *v;
      O = *pvp;
      vint = make_pair(0, 1);
      break;
    case 10:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = E(p) * A(p);
        G(p) = A(p);
        H(p) = pow(A(p) * K(p), 2) * (*t)(p);
      }
      N = O = *v;
      vint = make_pair(0, 0);
      break;
    case 11:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = E(p) * B(p);
        G(p) = B(p);
        H(p) = pow(B(p) / K(p), 2) / (4.0 * (*t)(p));
      }
      N = O = *pvp;
      vint = make_pair(1, 1);
      break;
    }

    Matrix WN(nbasis, nbasis);
    Matrix WO(nbasis, nbasis);
    for (int p = 0; p != nbasis; ++p) {
      for (int q = 0; q != nbasis; ++q) {
        WN(q, p) = N(q, p) / (E(p) + E(q));
        WO(q, p) = O(q, p) / (E(p) + E(q));
      }
    }

    Matrix WHOGCF(nbasis, nbasis);
    Matrix WHNFCG(nbasis, nbasis);
    {
      const Matrix HOGCF = H * WO * G * CPW * F;
      const Matrix HNFCG = H * WN * F * CPW * G;
      for (int p = 0; p != nbasis; ++p)
        for (int q = 0; q != nbasis; ++q) {
          WHOGCF(q, p) = HOGCF(q, p) / (E(p) + E(q));
          WHNFCG(q, p) = HNFCG(q, p) / (E(p) + E(q));
        }
    }

    if (vint.first)
      *pvpden += *wmat * WHOGCF ^ *wmat;
    else
      *den += *wmat * WHOGCF ^ *wmat;

    if (vint.second)
      *pvpden += *wmat * WHNFCG ^ *wmat;
    else
      *den += *wmat * WHNFCG ^ *wmat;
  }
  return { den, pvpden };
}


// Effective density matrix for overlap gradient
shared_ptr<const Matrix> DKHgrad::compute_sden_(shared_ptr<const Matrix> rdm1, shared_ptr<const Matrix> erdm1,
                                                shared_ptr<const Matrix> pmat, shared_ptr<const Matrix> wmat,
                                                shared_ptr<const Matrix> wmat_rev, shared_ptr<const DiagMatrix> t,
                                                shared_ptr<const Matrix> v, shared_ptr<const Matrix> pvp,
                                                shared_ptr<const Matrix> zpq, shared_ptr<const Matrix> ypq) {
  // Extra contributions stem for dependence of energy on overlap matrix (due to reverse transformation)
  const int nbasis = t->ndim();
  const double c2 = c__ * c__;
  DiagMatrix E(nbasis);
  DiagMatrix A(nbasis);
  DiagMatrix K(nbasis);
  DiagMatrix B(nbasis);
  for (int p = 0; p != nbasis; ++p) {
    E(p) = c__ * sqrt(2.0 * (*t)(p) + c2);
    A(p) = sqrt((E(p) + c2) / (2.0 * E(p)));
    K(p) = c__ / (E(p) + c2);
    B(p) = A(p) * K(p);
  }

  const Matrix CPW = *pmat * *rdm1 ^ (*wmat_rev % *pmat);
  auto den = make_shared<Matrix>(nbasis, nbasis);
  // DKH1 correction
  *den += 2.0 * ((CPW * E - c2 * CPW) + CPW * (A * *v * A + B * *pvp * B)) ^ *wmat;

  Matrix N(nbasis, nbasis);
  Matrix O(nbasis, nbasis);
  DiagMatrix F(nbasis);
  DiagMatrix G(nbasis);
  DiagMatrix H(nbasis);
  // DKH2 corrections, one for each term
  for (int i = 0; i != 12; ++i) {
    switch (i) {
    case 0:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = B(p);
        G(p) = A(p);
        H(p) = -B(p) * E(p) * A(p);
      }
      N = *pvp;
      O = *v;
      break;
    case 1:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = A(p);
        G(p) = B(p);
        H(p) = -A(p) * E(p) * B(p);
      }
      N = *v;
      O = *pvp;
      break;
    case 2:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = G(p) = A(p);
        H(p) = 2.0 * pow(A(p) * K(p), 2) * (*t)(p) * E(p);
      }
      N = O = *v;
      break;
    case 3:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = G(p) = B(p);
        H(p) = pow(B(p) / K(p), 2) * E(p) / (2.0 * (*t)(p));
      }
      N = O = *pvp;
      break;
    case 4:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = B(p);
        G(p) = A(p) * E(p);
        H(p) = -B(p) * A(p) / 2.0;
      }
      N = *pvp;
      O = *v;
      break;
    case 5:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = A(p);
        G(p) = B(p) * E(p);
        H(p) = -A(p) * B(p) / 2.0;
      }
      N = *v;
      O = *pvp;
      break;
    case 6:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = A(p);
        G(p) = A(p) * E(p);
        H(p) = pow(A(p) * K(p), 2) * (*t)(p);
      }
      N = O = *v;
      break;
    case 7:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = B(p);
        G(p) = B(p) * E(p);
        H(p) = pow(B(p) / K(p), 2) / (4.0 * (*t)(p));
      }
      N = O = *pvp;
      break;
    case 8:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = E(p) * B(p);
        G(p) = A(p);
        H(p) = -B(p) * A(p) / 2.0;
      }
      N = *pvp;
      O = *v;
      break;
    case 9:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = E(p) * A(p);
        G(p) = B(p);
        H(p) = -A(p) * B(p) / 2.0;
      }
      N = *v;
      O = *pvp;
      break;
    case 10:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = E(p) * A(p);
        G(p) = A(p);
        H(p) = pow(A(p) * K(p), 2) * (*t)(p);
      }
      N = O = *v;
      break;
    case 11:
      for (int p = 0; p != nbasis; ++p) {
        F(p) = E(p) * B(p);
        G(p) = B(p);
        H(p) = pow(B(p) / K(p), 2) / (4.0 * (*t)(p));
      }
      N = O = *pvp;
      break;
    }
    Matrix WN(nbasis, nbasis);
    Matrix WO(nbasis, nbasis);
    for (int p = 0; p != nbasis; ++p) {
      for (int q = 0; q != nbasis; ++q) {
        WN(q, p) = N(q, p) / (E(p) + E(q));
        WO(q, p) = O(q, p) / (E(p) + E(q));
      }
    }
    *den += CPW * (F * WN * H * WO * G + G * WO * H * WN * F) ^ *wmat;
  }

  // a_tilde from dL/dU_pq
  const Matrix at = 2.0 * *zpq * *t;
  // X_bar from Lagrangian, satisfies dL/dU_pq = 0
  Matrix xb(nbasis, nbasis);
  for (int p = 0; p != nbasis; ++p) {
    for (int q = 0; q != nbasis; ++q) {
      xb(q, p) = 0.25 * ((*ypq)(q, p) + (*ypq)(p, q) + at(q, p) + at(p, q));
    }
  }

  *den -= (*pmat * *erdm1 ^ *pmat) + (*wmat * xb ^ *wmat);
  return den;
}
