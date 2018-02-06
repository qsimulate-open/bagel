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
#include <src/mat1e/contrcoeff.h>

using namespace std;
using namespace bagel;


DKHgrad::DKHgrad(shared_ptr<const Molecule> current) {
  // Gradient integrals and RDMs to be computed in uncontracted basis
  auto mol = make_shared<Molecule>(*current);
  mol = mol->uncontract();
  nbasis_ = mol->nbasis();

  // Compute momentum transformation matrix
  const Overlap overlap(mol);
  shared_ptr<const Matrix> gamma = overlap.tildex();
  const Kinetic kinetic(mol);
  auto lambda = make_shared<Matrix>(*gamma % kinetic * *gamma);
  VectorB kinetic0(nbasis_);
  // Kinetic operator is diagonal in momentum space
  lambda->diagonalize(kinetic0);
  kinetic_ = DiagVec(kinetic0);
  wtrans_ = *gamma * *lambda;
  wtrans_rev_ = overlap * wtrans_;

  // Evaluate integrals in momentum space
  const NAI nai(mol);
  nai_ = wtrans_ % nai * wtrans_;
  const Small1e<NAIBatch> small1e(mol);
  smallnai_ = wtrans_ % small1e[0] * wtrans_;

  ptrans_ = ContrCoeff(current, nbasis_);

  zmult_ = ederiv_ = Matrix(nbasis_, nbasis_);
}

shared_ptr<const Matrix> DKHgrad::compute_tden(shared_ptr<const Matrix> rdm1) {
  const double c2 = c__ * c__;
  DiagVec E(nbasis_);
  DiagVec A(nbasis_);
  DiagVec B(nbasis_);
  DiagVec K(nbasis_);
  DiagVec dE(nbasis_);
  DiagVec dA(nbasis_);
  DiagVec dK(nbasis_);
  DiagVec dB(nbasis_);
  for (int p = 0; p != nbasis_; ++p) {
    E(p) = c__ * sqrt(2 * kinetic_(p) + c2);
    A(p) = sqrt((c2 + E(p)) / (2 * E(p)));
    K(p) = c__ / (E(p) + c2);
    B(p) = A(p) * K(p);
    dE(p) = c__ / sqrt(2 * kinetic_(p) + c2);
    dA(p) = -c2 * dE(p) / (4 * pow(E(p), 2) * A(p));
    dK(p) = -c__ * dE(p) / pow(E(p) + c2, 2);
    dB(p) = A(p) * dK(p) + dA(p) * K(p);
  }

  const Matrix CPW = (ptrans_ % wtrans_rev_) % *rdm1 * (ptrans_ % wtrans_rev_);
  const Matrix CAN = CPW * A * nai_;
  const Matrix NAC = nai_ * A * CPW;
  const Matrix CBS = CPW * B * smallnai_;
  const Matrix SBC = smallnai_ * B * CPW;
  auto den = make_shared<Matrix>(nbasis_, nbasis_);
  DiagVec EC(nbasis_);
  // DKH1 correction
  ederiv_ += 2 * (CPW * E - c2 * CPW + (CAN + NAC) * A + (CBS + SBC) * B);
  for (int p = 0; p != nbasis_; ++p) {
    // dE_DKH1/dU_pq
    ederiv_(p, p) += 2 * (CPW(p, p) * dE(p) + ((CAN(p, p) + NAC(p, p)) * dA(p) + (CBS(p, p) + SBC(p, p)) * dB(p))) * kinetic_(p);
    EC(p) = dE(p) * CPW(p, p) + 2 * (NAC(p, p) * dA(p) + SBC(p, p) * dB(p));
  }
  *den += wtrans_ * EC ^ wtrans_;
  // den->print("d_tilde dkh1");

  Matrix N(nbasis_, nbasis_);
  Matrix O(nbasis_, nbasis_);
  DiagVec F(nbasis_);
  DiagVec G(nbasis_);
  DiagVec H(nbasis_);
  DiagVec dF(nbasis_);
  DiagVec dG(nbasis_);
  DiagVec dH(nbasis_);
  // DKH2 corrections, one for each term
  for (int i = 0; i != 12; ++i) {
    switch (i) {
    case 0:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = B(p);
        G(p) = A(p);
        H(p) = -B(p) * E(p) * A(p);
        dF(p) = dB(p);
        dG(p) = dA(p);
        dH(p) = -dB(p) * E(p) * A(p) - B(p) * dE(p) * A(p) - B(p) * E(p) * dA(p);
      }
      N = smallnai_;
      O = nai_;
      break;
    case 1:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = A(p);
        G(p) = B(p);
        H(p) = -A(p) * E(p) * B(p);
        dF(p) = dA(p);
        dG(p) = dB(p);
        dH(p) = -dA(p) * E(p) * B(p) - A(p) * dE(p) * B(p) - A(p) * E(p) * dB(p);
      }
      N = nai_;
      O = smallnai_;
      break;
    case 2:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = G(p) = A(p);
        H(p) = 2 * pow(A(p) * K(p), 2) * kinetic_(p) * E(p);
        dF(p) = dG(p) = dA(p);
        dH(p) = 4 * A(p) * dA(p) * kinetic_(p) * pow(K(p), 2) * E(p) + 2 * pow(A(p) * K(p), 2) * E(p)
              + 4 * pow(A(p), 2) * kinetic_(p) * K(p) * dK(p) * E(p) + 2 * pow(A(p) * K(p), 2) * kinetic_(p) * dE(p);
      }
      N = O = nai_;
      break;
    case 3:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = G(p) = B(p);
        H(p) = pow(B(p) / K(p), 2) * E(p) / (2 * kinetic_(p));
        dF(p) = dG(p) = dB(p);
        dH(p) = ((2 * B(p) * dB(p) * E(p) + pow(B(p), 2) * dE(p)) * kinetic_(p) * K(p)
              - pow(B(p), 2) * E(p) * (K(p) + 2 * kinetic_(p) * dK(p))) / (2 * pow(kinetic_(p) * K(p), 2) * K(p));
      }
      N = O = smallnai_;
      break;
    case 4:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = B(p);
        G(p) = A(p) * E(p);
        H(p) = -B(p) * A(p) / 2;
        dF(p) = dB(p);
        dG(p) = dA(p) * E(p) + A(p) * dE(p);
        dH(p) = -dB(p) * A(p) / 2 - B(p) * dA(p) / 2;
      }
      N = smallnai_;
      O = nai_;
      break;
    case 5:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = A(p);
        G(p) = B(p) * E(p);
        H(p) = -A(p) * B(p) / 2;
        dF(p) = dA(p);
        dG(p) = dB(p) * E(p) + B(p) * dE(p);
        dH(p) = -dA(p) * B(p) / 2 - A(p) * dB(p) / 2;
      }
      N = nai_;
      O = smallnai_;
      break;
    case 6:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = A(p);
        G(p) = A(p) * E(p);
        H(p) = pow(A(p) * K(p), 2) * kinetic_(p);
        dF(p) = dA(p);
        dG(p) = dA(p) * E(p) + A(p) * dE(p);
        dH(p) = 2 * A(p) * dA(p) * kinetic_(p) * pow(K(p), 2) + pow(A(p) * K(p), 2)
              + 2 * pow(A(p), 2) * kinetic_(p) * K(p) * dK(p);
      }
      N = O = nai_;
      break;
    case 7:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = B(p);
        G(p) = B(p) * E(p);
        H(p) = pow(B(p) / K(p), 2) / (4 * kinetic_(p));
        dF(p) = dB(p);
        dG(p) = dB(p) * E(p) + B(p) * dE(p);
        dH(p) = (2 * B(p) * dB(p) * kinetic_(p) * K(p)
              - pow(B(p), 2) * (K(p) + 2 * kinetic_(p) * dK(p))) / (4 * pow(kinetic_(p) * K(p), 2) * K(p));
      }
      N = O = smallnai_;
      break;
    case 8:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = E(p) * B(p);
        G(p) = A(p);
        H(p) = -B(p) * A(p) / 2;
        dF(p) = dE(p) * B(p) + E(p) * dB(p);
        dG(p) = dA(p);
        dH(p) = -dB(p) * A(p) / 2 - B(p) * dA(p) / 2;
      }
      N = smallnai_;
      O = nai_;
      break;
    case 9:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = E(p) * A(p);
        G(p) = B(p);
        H(p) = -A(p) * B(p) / 2;
        dF(p) = dE(p) * A(p) + E(p) * dA(p);
        dG(p) = dB(p);
        dH(p) = -dA(p) * B(p) / 2 - A(p) * dB(p) / 2;
      }
      N = nai_;
      O = smallnai_;
      break;
    case 10:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = E(p) * A(p);
        G(p) = A(p);
        H(p) = pow(A(p) * K(p), 2) * kinetic_(p);
        dF(p) = dE(p) * A(p) + E(p) * dA(p);
        dG(p) = dA(p);
        dH(p) = 2 * A(p) * dA(p) * kinetic_(p) * pow(K(p), 2) + pow(A(p) * K(p), 2)
              + 2 * pow(A(p), 2) * kinetic_(p) * K(p) * dK(p);
      }
      N = O = nai_;
      break;
    case 11:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = E(p) * B(p);
        G(p) = B(p);
        H(p) = pow(B(p) / K(p), 2) / (4 * kinetic_(p));
        dF(p) = dE(p) * B(p) + E(p) * dB(p);
        dG(p) = dB(p);
        dH(p) = (2 * B(p) * dB(p) * kinetic_(p) * K(p)
              - pow(B(p), 2) * (K(p) + 2 * kinetic_(p) * dK(p))) / (4 * pow(kinetic_(p) * K(p), 2) * K(p));
      }
      N = O = smallnai_;
      break;
    }
    Matrix WN(nbasis_, nbasis_);
    Matrix WO(nbasis_, nbasis_);
    Matrix WWN(nbasis_, nbasis_);
    Matrix WWO(nbasis_, nbasis_);
    for (int p = 0; p != nbasis_; ++p) {
      for (int q = 0; q != nbasis_; ++q) {
        WN(q, p) = N(q, p) / (E(p) + E(q));
        WO(q, p) = O(q, p) / (E(p) + E(q));
        WWN(q, p) = WN(q, p) / (E(p) + E(q));
        WWO(q, p) = WO(q, p) / (E(p) + E(q));
      }
    }
    const Matrix CFN = CPW * F * WN;
    const Matrix CGO = CPW * G * WO;
    const Matrix HNFC = H * WN * F * CPW * G;
    const Matrix HOGC = H * WO * G * CPW * F;
    const Matrix CFNHOG = CFN * H * WO * dG * kinetic_;
    const Matrix CFNHOOGE = CFN * H * WWO * G * dE * kinetic_;
    const Matrix CGOHNF = CGO * H * WN * dF * kinetic_;
    const Matrix CGOHNNFE = CGO * H * WWN * F * dE * kinetic_;
    const Matrix CFNGH = (CFN * G * dH - CPW * F * WWN * G * H * dE) * kinetic_;
    const Matrix CFNGHE = CFN * G * H * dE * kinetic_;
    const Matrix CGOFH = (CGO * F * dH - CPW * G * WWO * F * H * dE) * kinetic_;
    const Matrix CGOFHE = CGO * F * H * dE * kinetic_;
    Matrix CH(nbasis_, nbasis_);
    Matrix CH2(nbasis_, nbasis_);
    for (int p = 0; p != nbasis_; ++p) {
      for (int q = 0; q != nbasis_; ++q) {
        if (p != q) {
          CH(q, p) = CPW(q, p) * H(p);
          CH2(q, p) = CPW(q, p) * dH(p);
        }
      }
    }
    Matrix WHNFC(nbasis_, nbasis_);
    Matrix WHOGC(nbasis_, nbasis_);
    Matrix CHFGO(nbasis_, nbasis_);
    Matrix CHGFN(nbasis_, nbasis_);
    DiagVec HCFNG(nbasis_);
    DiagVec HCGOF(nbasis_);
    for (int p = 0; p != nbasis_; ++p) {
      for (int q = 0; q != nbasis_; ++q) {
        WHNFC(q, p) = HNFC(q, p) / (E(p) + E(q));
        WHOGC(q, p) = HOGC(q, p) / (E(p) + E(q));
        CHFGO(q, p) = CH(q, p) * F(q) * G(p) * WO(p, p) / (E(q) + E(p));
        CHGFN(q, p) = CH(q, p) * G(q) * F(p) * WN(p, p) / (E(q) + E(p));
      }
      HCFNG(p) = H(p) * CFN(p, p) * G(p) / (2 * E(p));
      HCGOF(p) = H(p) * CGO(p, p) * F(p) / (2 * E(p));
    }
    const Matrix CH2FNG = WN * F * CH2 * G * kinetic_;
    const Matrix CHFNGE = WN * F * CH * G * dE * kinetic_;
    const Matrix CHFNNGE = WWN * F * CH * G * dE * kinetic_;
    const Matrix CH2GOF = WO * G * CH2 * F * kinetic_;
    const Matrix CHGOFE = WO * G * CH * F * dE * kinetic_;
    const Matrix CHGOOFE = WWO * G * CH * F * dE * kinetic_;
    Matrix FN(nbasis_, nbasis_);
    Matrix FO(nbasis_, nbasis_);
    Matrix FNN(nbasis_, nbasis_);
    Matrix FOO(nbasis_, nbasis_);
    for (int p = 0; p != nbasis_; ++p) {
      for (int q = 0; q != nbasis_; ++q) {
        if (p != q) {
          FN(q, p) = WN(q, p) * F(q);
          FO(q, p) = WO(q, p) * F(q);
          FNN(q, p) = WWN(q, p) * F(q);
          FOO(q, p) = WWO(q, p) * F(q);
        }
      }
    }
    const Matrix CFN2 = CPW * FN;
    const Matrix CFO2 = CPW * FO;
    const Matrix CFNN2 = CPW * FNN;
    const Matrix CFOO2 = CPW * FOO;
    Matrix GCFN2H(nbasis_, nbasis_);
    Matrix GCFO2H(nbasis_, nbasis_);
    Matrix GCFN2H2(nbasis_, nbasis_);
    Matrix GCFO2H2(nbasis_, nbasis_);
    Matrix GCFNN2H2(nbasis_, nbasis_);
    Matrix GCFOO2H2(nbasis_, nbasis_);
    for (int p = 0; p != nbasis_; ++p) {
      for (int q = 0; q != nbasis_; ++q) {
        if (p != q) {
          GCFN2H(q, p) = G(q) * CFN2(q, p) * H(p) / (E(q) + E(p));
          GCFO2H(q, p) = G(q) * CFO2(q, p) * H(p) / (E(q) + E(p));
          GCFN2H2(q, p) = G(q) * CFN2(q, p) * dH(p) * kinetic_(p);
          GCFO2H2(q, p) = G(q) * CFO2(q, p) * dH(p) * kinetic_(p);
          GCFNN2H2(q, p) = G(q) * CFNN2(q, p) * H(p) * dE(p) * kinetic_(p);
          GCFOO2H2(q, p) = G(q) * CFOO2(q, p) * H(p) * dE(p) * kinetic_(p);
        }
      }
    }
    const Matrix OGCFN2H2 = WO * (GCFN2H2 - 2 * GCFNN2H2);
    const Matrix NGCFO2H2 = WN * (GCFO2H2 - 2 * GCFOO2H2);
    // dE_DKH2/dU_pq
    ederiv_ += CFN * H * WO * G + CGO * H * WN * F + O * (WHNFC + CHGFN + GCFN2H) + N * (WHOGC + CHFGO + GCFO2H)
            + O * HCFNG + N * HCGOF;

    for (int p = 0; p != nbasis_; ++p) {
      ederiv_(p, p) += 2 * (CFNHOG(p, p) + CGOHNF(p, p) - CFNHOOGE(p, p) - CGOHNNFE(p, p))
                    + (CFNGH(p, p) + CH2FNG(p, p) - CHFNNGE(p, p)) * WO(p, p) - (CFNGHE(p, p) + CHFNGE(p, p)) * WWO(p, p)
                    + (CGOFH(p, p) + CH2GOF(p, p) - CHGOOFE(p, p)) * WN(p, p) - (CGOFHE(p, p) + CHGOFE(p, p)) * WWN(p, p)
                    + OGCFN2H2(p, p) + NGCFO2H2(p, p);
    }
    const Matrix OHNF = WO * H * (WN * dF - WWN * F * dE);
    const Matrix NHOG = WN * H * (WO * dG - WWO * G * dE);
    const Matrix OGCFN = WO * G * CPW * F * WN;
    const Matrix OGCFNN = WO * G * CPW * F * WWN;
    const Matrix OOGCFN = WWO * G * CPW * F * WN;
    Matrix OHNFC(nbasis_, nbasis_);
    Matrix NHOGC(nbasis_, nbasis_);
    DiagVec OGCFNH(nbasis_);
    DiagVec OGCFNNH(nbasis_);
    DiagVec OOGCFNH(nbasis_);
    for (int p = 0; p != nbasis_; ++p) {
      for (int q = 0; q != nbasis_; ++q) {
        OHNFC(q, p) = OHNF(q, p) * CPW(q, p);
        NHOGC(q, p) = NHOG(q, p) * CPW(q, p);
      }
      OGCFNH(p) = OGCFN(p, p) * dH(p);
      OGCFNNH(p) = OGCFNN(p, p) * H(p) * dE(p);
      OOGCFNH(p) = OOGCFN(p, p) * H(p) * dE(p);
    }
    const DiagVec dkh2(OHNFC * *G.data() + NHOGC * *F.data() + *OGCFNH.data() - *OGCFNNH.data() - *OOGCFNH.data());
    *den += wtrans_ * dkh2 ^ wtrans_;
    // den->print("d_tilde " + to_string(i));
  }

  // z_pq from Lagrangian
  for (int p = 0; p != nbasis_; ++p) {
    for (int q = 0; q != nbasis_; ++q) {
      zmult_(q, p) = p == q ? 0 : -0.5 * (ederiv_(q, p) - ederiv_(p, q)) / (kinetic_(q) - kinetic_(p));
    }
  }

  *den += wtrans_ * zmult_ ^ wtrans_;
  // zmult_.print("z_pq");
  // den->print("d_tilde full");
  return den;
}

array<shared_ptr<const Matrix>, 2> DKHgrad::compute_vden(shared_ptr<const Matrix> rdm1) {
  const double c2 = c__ * c__;
  DiagVec E(nbasis_);
  DiagVec A(nbasis_);
  DiagVec K(nbasis_);
  DiagVec B(nbasis_);
  for (int p = 0; p != nbasis_; ++p) {
    E(p) = c__ * sqrt(2 * kinetic_(p) + c2);
    A(p) = sqrt((E(p) + c2) / (2 * E(p)));
    K(p) = c__ / (E(p) + c2);
    B(p) = A(p) * K(p);
  }

  const Matrix CPW = (ptrans_ % wtrans_rev_) % *rdm1 * (ptrans_ % wtrans_rev_);
  auto den = make_shared<Matrix>(nbasis_, nbasis_);
  auto pvpden = make_shared<Matrix>(nbasis_, nbasis_);
  // DKH1 correction
  *den += wtrans_ * A * CPW * A ^ wtrans_;
  *pvpden += wtrans_ * B * CPW * B ^ wtrans_;

  Matrix N(nbasis_, nbasis_);
  Matrix O(nbasis_, nbasis_);
  DiagVec F(nbasis_);
  DiagVec G(nbasis_);
  DiagVec H(nbasis_);
  pair<int, int> vint;
  // DKH2 corrections, one for each term
  for (int i = 0; i != 12; ++i) {
    switch (i) {
    case 0:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = B(p);
        G(p) = A(p);
        H(p) = -B(p) * E(p) * A(p);
      }
      N = smallnai_;
      O = nai_;
      vint = make_pair(1, 0);
      break;
    case 1:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = A(p);
        G(p) = B(p);
        H(p) = -A(p) * E(p) * B(p);
      }
      N = nai_;
      O = smallnai_;
      vint = make_pair(0, 1);
      break;
    case 2:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = G(p) = A(p);
        H(p) = 2 * pow(A(p) * K(p), 2) * kinetic_(p) * E(p);
      }
      N = O = nai_;
      vint = make_pair(0, 0);
      break;
    case 3:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = G(p) = B(p);
        H(p) = pow(B(p) / K(p), 2) * E(p) / (2 * kinetic_(p));
      }
      N = O = smallnai_;
      vint = make_pair(1, 1);
      break;
    case 4:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = B(p);
        G(p) = A(p) * E(p);
        H(p) = -B(p) * A(p) / 2;
      }
      N = smallnai_;
      O = nai_;
      vint = make_pair(1, 0);
      break;
    case 5:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = A(p);
        G(p) = B(p) * E(p);
        H(p) = -A(p) * B(p) / 2;
      }
      N = nai_;
      O = smallnai_;
      vint = make_pair(0, 1);
      break;
    case 6:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = A(p);
        G(p) = A(p) * E(p);
        H(p) = pow(A(p) * K(p), 2) * kinetic_(p);
      }
      N = O = nai_;
      vint = make_pair(0, 0);
      break;
    case 7:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = B(p);
        G(p) = B(p) * E(p);
        H(p) = pow(B(p) / K(p), 2) / (4 * kinetic_(p));
      }
      N = O = smallnai_;
      vint = make_pair(1, 1);
      break;
    case 8:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = E(p) * B(p);
        G(p) = A(p);
        H(p) = -B(p) * A(p) / 2;
      }
      N = smallnai_;
      O = nai_;
      vint = make_pair(1, 0);
      break;
    case 9:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = E(p) * A(p);
        G(p) = B(p);
        H(p) = -A(p) * B(p) / 2;
      }
      N = nai_;
      O = smallnai_;
      vint = make_pair(0, 1);
      break;
    case 10:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = E(p) * A(p);
        G(p) = A(p);
        H(p) = pow(A(p) * K(p), 2) * kinetic_(p);
      }
      N = O = nai_;
      vint = make_pair(0, 0);
      break;
    case 11:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = E(p) * B(p);
        G(p) = B(p);
        H(p) = pow(B(p) / K(p), 2) / (4 * kinetic_(p));
      }
      N = O = smallnai_;
      vint = make_pair(1, 1);
      break;
    }
    Matrix WN(nbasis_, nbasis_);
    Matrix WO(nbasis_, nbasis_);
    for (int p = 0; p != nbasis_; ++p) {
      for (int q = 0; q != nbasis_; ++q) {
        WN(q, p) = N(q, p) / (E(p) + E(q));
        WO(q, p) = O(q, p) / (E(p) + E(q));
      }
    }
    const Matrix HOGCF = H * WO * G * CPW * F;
    const Matrix HNFCG = H * WN * F * CPW * G;
    Matrix WHOGCF(nbasis_, nbasis_);
    Matrix WHNFCG(nbasis_, nbasis_);
    for (int p = 0; p != nbasis_; ++p) {
      for (int q = 0; q != nbasis_; ++q) {
        WHOGCF(q, p) = HOGCF(q, p) / (E(p) + E(q));
        WHNFCG(q, p) = HNFCG(q, p) / (E(p) + E(q));
      }
    }
    if (vint.first) {
      *pvpden += wtrans_ * WHOGCF ^ wtrans_;
    } else {
      *den += wtrans_ * WHOGCF ^ wtrans_;
    }
    if (vint.second) {
      *pvpden += wtrans_ * WHNFCG ^ wtrans_;
    } else {
      *den += wtrans_ * WHNFCG ^ wtrans_;
    }
  }
  return { den, pvpden };
}

shared_ptr<const Matrix> DKHgrad::compute_sden(shared_ptr<const Matrix> rdm1, shared_ptr<const Matrix> erdm1) {
  // Extra contributions stem for dependence of energy on overlap matrix (due to reverse transformation)
  const double c2 = c__ * c__;
  DiagVec E(nbasis_);
  DiagVec A(nbasis_);
  DiagVec K(nbasis_);
  DiagVec B(nbasis_);
  for (int p = 0; p != nbasis_; ++p) {
    E(p) = c__ * sqrt(2 * kinetic_(p) + c2);
    A(p) = sqrt((E(p) + c2) / (2 * E(p)));
    K(p) = c__ / (E(p) + c2);
    B(p) = A(p) * K(p);
  }

  const Matrix CPW = ptrans_ * *rdm1 ^ (wtrans_rev_ % ptrans_);
  auto den = make_shared<Matrix>(nbasis_, nbasis_);
  // DKH1 correction
  *den += 2 * ((CPW * E - c2 * CPW) + CPW * (A * nai_ * A + B * smallnai_ * B)) ^ wtrans_;
  // den->print("X_tilde dkh1");

  Matrix N(nbasis_, nbasis_);
  Matrix O(nbasis_, nbasis_);
  DiagVec F(nbasis_);
  DiagVec G(nbasis_);
  DiagVec H(nbasis_);
  // DKH2 corrections, one for each term
  for (int i = 0; i != 12; ++i) {
    switch (i) {
    case 0:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = B(p);
        G(p) = A(p);
        H(p) = -B(p) * E(p) * A(p);
      }
      N = smallnai_;
      O = nai_;
      break;
    case 1:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = A(p);
        G(p) = B(p);
        H(p) = -A(p) * E(p) * B(p);
      }
      N = nai_;
      O = smallnai_;
      break;
    case 2:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = G(p) = A(p);
        H(p) = 2 * pow(A(p) * K(p), 2) * kinetic_(p) * E(p);
      }
      N = O = nai_;
      break;
    case 3:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = G(p) = B(p);
        H(p) = pow(B(p) / K(p), 2) * E(p) / (2 * kinetic_(p));
      }
      N = O = smallnai_;
      break;
    case 4:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = B(p);
        G(p) = A(p) * E(p);
        H(p) = -B(p) * A(p) / 2;
      }
      N = smallnai_;
      O = nai_;
      break;
    case 5:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = A(p);
        G(p) = B(p) * E(p);
        H(p) = -A(p) * B(p) / 2;
      }
      N = nai_;
      O = smallnai_;
      break;
    case 6:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = A(p);
        G(p) = A(p) * E(p);
        H(p) = pow(A(p) * K(p), 2) * kinetic_(p);
      }
      N = O = nai_;
      break;
    case 7:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = B(p);
        G(p) = B(p) * E(p);
        H(p) = pow(B(p) / K(p), 2) / (4 * kinetic_(p));
      }
      N = O = smallnai_;
      break;
    case 8:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = E(p) * B(p);
        G(p) = A(p);
        H(p) = -B(p) * A(p) / 2;
      }
      N = smallnai_;
      O = nai_;
      break;
    case 9:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = E(p) * A(p);
        G(p) = B(p);
        H(p) = -A(p) * B(p) / 2;
      }
      N = nai_;
      O = smallnai_;
      break;
    case 10:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = E(p) * A(p);
        G(p) = A(p);
        H(p) = pow(A(p) * K(p), 2) * kinetic_(p);
      }
      N = O = nai_;
      break;
    case 11:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = E(p) * B(p);
        G(p) = B(p);
        H(p) = pow(B(p) / K(p), 2) / (4 * kinetic_(p));
      }
      N = O = smallnai_;
      break;
    }
    Matrix WN(nbasis_, nbasis_);
    Matrix WO(nbasis_, nbasis_);
    for (int p = 0; p != nbasis_; ++p) {
      for (int q = 0; q != nbasis_; ++q) {
        WN(q, p) = N(q, p) / (E(p) + E(q));
        WO(q, p) = O(q, p) / (E(p) + E(q));
      }
    }
    *den += CPW * (F * WN * H * WO * G + G * WO * H * WN * F) ^ wtrans_;
    // den->print("X_tilde " + to_string(i));
  }

  // a_tilde from dL/dU_pq
  const Matrix at = 2 * zmult_ * kinetic_;
  // X_bar from Lagrangian, satisfies dL/dU_pq = 0
  Matrix xb(nbasis_, nbasis_);
  for (int p = 0; p != nbasis_; ++p) {
    for (int q = 0; q != nbasis_; ++q) {
      xb(q, p) = 0.25 * (ederiv_(q, p) + ederiv_(p, q) + at(q, p) + at(p, q));
    }
  }

  *den -= (ptrans_ * *erdm1 ^ ptrans_) + (wtrans_ * xb ^ wtrans_);
  // den->print("X_tilde full");
  return den;
}
