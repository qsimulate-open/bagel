//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dkhcoreinfo.cc
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


#include <src/wfn/dkhcoreinfo.h>
#include <src/mat1e/kinetic.h>
#include <src/mat1e/nai.h>
#include <src/mat1e/rel/small1e.h>
#include <src/mat1e/overlap.h>
#include <src/mat1e/contrcoeff.h>

using namespace std;
using namespace bagel;


DKHcoreInfo::DKHcoreInfo(shared_ptr<const Molecule> current) {
  auto mol = make_shared<Molecule>(*current);
  mol = mol->uncontract();
  nbasis_ = mol->nbasis();

  const Overlap overlap(mol);
  shared_ptr<const Matrix> gamma = overlap.tildex();
  const Kinetic kinetic(mol);
  auto lambda = make_shared<Matrix>(*gamma % kinetic * *gamma);
  VectorB kinetic0(nbasis_);
  lambda->diagonalize(kinetic0);
  kinetic_ = DiagVec(kinetic0);
  wtrans_ = *gamma * *lambda;
  wtrans_rev_ = overlap * wtrans_;

  const NAI nai(mol);
  nai_ = wtrans_ % nai * wtrans_;
  const Small1e<NAIBatch> small1e(mol);
  smallnai_ = wtrans_ % small1e[0] * wtrans_;

  ptrans_ = ContrCoeff(current, nbasis_);

  wtrans_.print("W");
  wtrans_rev_.print("S^unc W");
  ptrans_.print("P");

  zmult_ = ederiv_ = Matrix(nbasis_, nbasis_);
}

shared_ptr<const Matrix> DKHcoreInfo::compute_tden(shared_ptr<const Matrix> rdm1) {
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
  ederiv_ += 2 * (CPW * E - c2 * CPW + (CAN + NAC) * A + (CBS + SBC) * B);
  for (int p = 0; p != nbasis_; ++p) {
    ederiv_(p, p) += 2 * (CPW(p, p) * dE(p) + ((CAN(p, p) + NAC(p, p)) * dA(p) + (CBS(p, p) + SBC(p, p)) * dB(p))) * kinetic_(p);
    for (int a = 0; a != nbasis_; ++a) {
      for (int b = 0; b != nbasis_; ++b) {
        (*den)(b, a) += (dE(p) * CPW(p, p) + 2 * (NAC(p, p) * dA(p) + SBC(p, p) * dB(p))) * wtrans_(a, p) * wtrans_(b, p);
      }
    }
  }

  Matrix N(nbasis_, nbasis_);
  Matrix O(nbasis_, nbasis_);
  DiagVec F(nbasis_);
  DiagVec G(nbasis_);
  DiagVec H(nbasis_);
  DiagVec dF(nbasis_);
  DiagVec dG(nbasis_);
  DiagVec dH(nbasis_);
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
    const Matrix HSGC = H * WO * G * CPW * F;
    const Matrix CFNHO = CFN * H * WO;
    const Matrix CFNHOO = CFN * H * WWO;
    const Matrix CGOHN = CGO * H * WN;
    const Matrix CGOHNN = CGO * H * WWN;
    const Matrix CFNGH = (CFN * G * dH - CPW * F * WWN * G * H * dE) * kinetic_;
    const Matrix CFNGHE = CFN * G * H * dE * kinetic_;
    const Matrix CGOFH = (CGO * F * dH - CPW * G * WWO * F * H * dE) * kinetic_;
    const Matrix CGOFHE = CGO * F * H * dE * kinetic_;
    ederiv_ += CFN * H * WO * G + CGO * H * WN * F;
    for (int p = 0; p != nbasis_; ++p) {
      for (int q = 0; q != nbasis_; ++q) {
        for (int r = 0; r != nbasis_; ++r) {
          ederiv_(q, p) += (1 / (E(p) + E(r))) * (O(r, q) * HNFC(r, p) + N(r, q) * HSGC(r, p));
          if (p != r) {
            ederiv_(q, p) += (CPW(r, p) * H(p) / (E(r) + E(p))) * (F(r) * G(p) * N(r, q) * WO(p, p) + F(p) * G(r) * WN(p, p) * O(r, q));
            if (p == q) {
              ederiv_(p, p) += CPW(r, p) * (F(r) * (G(p) * dH(p) * WN(r, p) * WO(p, p) - G(p) * H(p) * dE(p) * (WN(r, p) * WWO(p, p)
                            + WWN(r, p) * WO(p, p))) + G(r) * (F(p) * dH(p) * WN(p, p) * WO(r, p) - F(p) * H(p) * dE(p) * (WWN(p, p) * WO(r, p)
                            + WN(p, p) * WWO(r, p)))) * kinetic_(p);
            }
            for (int s = 0; s != nbasis_; ++s) {
              if (p != s) {
                ederiv_(q, p) += CPW(s, r) * H(p) * (F(s) * G(r) * WN(s, p) * O(r, q) / (E(r) + E(p)) + F(r) * G(s) * N(s, q) * WO(r, p) / (E(s) + E(p)));
                if (p == q) {
                  ederiv_(p, p) += CPW(s, r) * (F(s) * G(r) * (dH(p) * WN(s, p) - 2 * H(p) * dE(p) * WWN(s, p)) * WO(r, p)
                                + F(r) * G(s) * (dH(p) * WO(r, p) - 2 * H(p) * dE(p) * WWO(r, p)) * WN(s, p)) * kinetic_(p);
                }
              }
            }
          }
          // for (int a = 0; a != nbasis_; ++a) {
          //   for (int b = 0; b != nbasis_; ++b) {
          //     (*den)(b, a) += (WO(r, q) * (WN(r, p) * G(q) * dF(p) - WWN(r, p) * F(p) * G(q) * dE(p)) * H(r) * wtrans_(a, p) * wtrans_(b, p)
          //                 + WN(r, p) * (WO(r, q) * F(p) * dG(q) - WWO(r, q) * F(p) * G(q) * dE(q)) * H(r) * wtrans_(a, q) * wtrans_(b, q)
          //                 + (WN(r, p) * WO(r, q) * F(p) * G(q) * dH(r) - WWN(r, p) * WO(r, q) * F(p) * G(q) * H(r) * dE(r)
          //                 - WN(r, p) * WWO(r, q) * F(p) * G(q) * H(r) * dE(r)) * wtrans_(a, r) * wtrans_(b, r)) * CPW(q, p);
          //   }
          // }
        }
        ederiv_(q, p) += (H(p) / (2 * E(p))) * (CFN(p, p) * G(p) * O(q, p) + CGO(p, p) * F(p) * N(q, p));
      }
      ederiv_(p, p) += 2 * (CFNHO(p, p) * dG(p) + CGOHN(p, p) * dF(p) - (CFNHOO(p, p) * G(p) + CGOHNN(p, p) * F(p)) * dE(p)) * kinetic_(p)
                    + CFNGH(p, p) * WO(p, p) - CFNGHE(p, p) * WWO(p, p) + CGOFH(p, p) * WN(p, p) - CGOFHE(p, p) * WWN(p, p);
    }
    const Matrix GOHNF = G * WO * H * WN * dF;
    const Matrix GOHNFE = G * WO * H * WWN * F * dE;
    for (int p = 0; p != nbasis_; ++p) {
      for (int q = 0; q != nbasis_; ++q) {
        for (int r = 0; r != nbasis_; ++r) {
          for (int a = 0; a != nbasis_; ++a) {
            for (int b = 0; b != nbasis_; ++b) {
              // how to modify the line below
              (*den)(b, a) += WO(r, q) * WN(r, p) * dF(p) * G(q) * H(r) * wtrans_(a, p) * wtrans_(b, p) * CPW(q, p)
                          - WO(r, q) * WWN(r, p) * F(p) * dE(p) * G(q) * H(r) * wtrans_(a, p) * wtrans_(b, p) * CPW(q, p)
                          + WN(r, p) * WO(r, q) * F(p) * dG(q) * H(r) * wtrans_(a, q) * wtrans_(b, q) * CPW(q, p)
                          - WN(r, p) * WWO(r, q) * F(p) * G(q) * dE(q) * H(r) * wtrans_(a, q) * wtrans_(b, q) * CPW(q, p)
              // how to modify the line below
                          + WN(r, p) * WO(r, q) * F(p) * G(q) * dH(r) * wtrans_(a, r) * wtrans_(b, r) * CPW(q, p)
                          - WWN(r, p) * WO(r, q) * F(p) * G(q) * H(r) * dE(r) * wtrans_(a, r) * wtrans_(b, r) * CPW(q, p)
                          - WN(r, p) * WWO(r, q) * F(p) * G(q) * H(r) * dE(r) * wtrans_(a, r) * wtrans_(b, r) * CPW(q, p);
            }
          }
        }
      }
    }
  }

  for (int p = 0; p != nbasis_; ++p) {
    for (int q = 0; q != nbasis_; ++q) {
      zmult_(q, p) = p == q ? 0 : -0.5 * (ederiv_(q, p) - ederiv_(p, q)) / (kinetic_(q) - kinetic_(p));
    }
  }

  *den += wtrans_ * zmult_ ^ wtrans_;
  ederiv_.print("Y_pq");
  zmult_.print("z_pq");
  den->print("d_tilde");
  return den;
}

array<shared_ptr<const Matrix>, 2> DKHcoreInfo::compute_vden(shared_ptr<const Matrix> rdm1) {
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
  *den += wtrans_ * A * CPW * A ^ wtrans_;
  *pvpden += wtrans_ * B * CPW * B ^ wtrans_;

  Matrix N(nbasis_, nbasis_);
  Matrix O(nbasis_, nbasis_);
  DiagVec F(nbasis_);
  DiagVec G(nbasis_);
  DiagVec H(nbasis_);
  pair<int, int> vint;
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
    const Matrix WHNF = wtrans_ * H * WN * F;
    const Matrix WHOG = wtrans_ * H * WO * G;
    for (int p = 0; p != nbasis_; ++p) {
      for (int q = 0; q != nbasis_; ++q) {
        for (int r = 0; r != nbasis_; ++r) {
          for (int a = 0; a != nbasis_; ++a) {
            for (int b = 0; b != nbasis_; ++b) {
              double nden = (H(r) / (E(p) + E(r))) * F(p) * G(q) * WO(r, q) * wtrans_(a, p) * wtrans_(b, r) * CPW(q, p);
              if (vint.first) {
                (*pvpden)(b, a) += nden;
              }
              else {
                (*den)(b, a) += nden;
              }
              double oden = (H(r) / (E(q) + E(r))) * F(p) * G(q) * WN(r, p) * wtrans_(a, q) * wtrans_(b, r) * CPW(p, q);
              if (vint.second) {
                (*pvpden)(b, a) += oden;
              }
              else {
                (*den)(b, a) += oden;
              }
            }
          }
        }
      }
    }
  }
  den->print("d_bar");
  pvpden->print("d_check");
  return { den, pvpden };
}

shared_ptr<const Matrix> DKHcoreInfo::compute_sden(shared_ptr<const Matrix> rdm1, shared_ptr<const Matrix> erdm1) {
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
  *den += 2 * ((CPW * E - c2 * CPW) + CPW * (A * nai_ * A + B * smallnai_ * B)) ^ wtrans_;

  Matrix N(nbasis_, nbasis_);
  Matrix O(nbasis_, nbasis_);
  DiagVec F(nbasis_);
  DiagVec G(nbasis_);
  DiagVec H(nbasis_);
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
  }

  Matrix at(nbasis_, nbasis_);
  for (int p = 0; p != nbasis_; ++p) {
    for (int q = 0; q != nbasis_; ++q) {
      at(q, p) = 2 * zmult_(q, p) * kinetic_(q);
    }
  }

  Matrix xb(nbasis_, nbasis_);
  for (int p = 0; p != nbasis_; ++p) {
    for (int q = 0; q != nbasis_; ++q) {
      xb(q, p) = 0.25 * (ederiv_(q, p) + ederiv_(p, q) + at(q, p) + at(p, q));
    }
  }

  *den -= (ptrans_ * *erdm1 ^ ptrans_) + (wtrans_ * xb ^ wtrans_);
  at.print("a_tilde");
  xb.print("X_bar");
  den->print("X_tilde");
  return den;
}
