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
  kinetic_ = VectorB(nbasis_);
  lambda->diagonalize(kinetic_);
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
  VectorB E(nbasis_), A(nbasis_), B(nbasis_), K(nbasis_), dE(nbasis_), dA(nbasis_), dK(nbasis_), dB(nbasis_);
  for (int p = 0; p != nbasis_; ++p) {
    E(p) = c__ * sqrt(2 * kinetic_(p) + c2);
    A(p) = sqrt(0.5 * (E(p) + c2) / E(p));
    K(p) = c__ / (E(p) + c2);
    B(p) = A(p) * K(p);
    dE(p) = c__ / sqrt(2 * kinetic_(p) + c2);
    dA(p) = -c2 * c__ / (4 * E(p) * E(p) * A(p) * sqrt(2 * kinetic_(p) + c2));
    dK(p) = -c2 / (pow(E(p) + c2, 2) * sqrt(2 * kinetic_(p) + c2));
    dB(p) = A(p) * dK(p) + dA(p) * K(p);
  }

  const Matrix CPW = (ptrans_ % wtrans_rev_) % *rdm1 * (ptrans_ % wtrans_rev_);
  for (int q = 0; q != nbasis_; ++q) {
    for (int p = 0; p != nbasis_; ++p) {
      ederiv_(p, q) += 2 * CPW(p, q) * (E(q) - c2);
      if (p == q) {
        ederiv_(p, q) += 2 * CPW(p, p) * dE(p) * kinetic_(p);
      }
      for (int r = 0; r != nbasis_; ++r) {
        ederiv_(p, q) += 2 * CPW(r, p) * (A(q) * A(r) * nai_(r, q) + B(q) * B(r) * smallnai_(r, q));
        ederiv_(p, q) += 2 * CPW(r, q) * (A(q) * A(r) * nai_(r, p) + B(q) * B(r) * smallnai_(r, p));
        if (p == q) {
          ederiv_(p, q) += 2 * CPW(r, p) * (dA(p) * kinetic_(p) * A(r) * nai_(r, q) + dB(p) * kinetic_(p) * B(r) * smallnai_(r, q));
          ederiv_(p, q) += 2 * CPW(r, q) * (dA(p) * kinetic_(p) * A(r) * nai_(r, p) + dB(p) * kinetic_(p) * B(r) * smallnai_(r, p));
        }
      }
    }
  }

  shared_ptr<Matrix> den = make_shared<Matrix>(nbasis_, nbasis_);
  for (int b = 0; b != nbasis_; ++b) {
    for (int a = 0; a != nbasis_; ++a) {
      for (int p = 0; p != nbasis_; ++p) {
        (*den)(a, b) += dE(p) * wtrans_(a, p) * wtrans_(b, p) * CPW(p, p);
      }
      for (int q = 0; q != nbasis_; ++q) {
        for (int p = 0; p != nbasis_; ++p) {
          (*den)(a, b) += (nai_(p, q) * (A(q) * dA(p) * wtrans_(a, p) * wtrans_(b, p) + A(p) * dA(q) * wtrans_(a, q) * wtrans_(b, q))
                      + smallnai_(p, q) * (B(q) * dB(p) * wtrans_(a, p) * wtrans_(b, p) + B(p) * dB(q) * wtrans_(a, q) * wtrans_(b, q))) * CPW(p, q);
        }
      }
    }
  }

  Matrix N(nbasis_, nbasis_), O(nbasis_, nbasis_);
  VectorB F(nbasis_), G(nbasis_), H(nbasis_), dF(nbasis_), dG(nbasis_), dH(nbasis_);
  for (int i = 0; i != 7; ++i) {
    switch (i) {
    case 0:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = B(p);
        G(p) = A(p);
        H(p) = -2 * B(p) * E(p) * A(p);
        dF(p) = dB(p);
        dG(p) = dA(p);
        dH(p) = -2 * dB(p) * E(p) * A(p) - 2 * B(p) * dE(p) * A(p) - 2 * B(p) * E(p) * dA(p);
      }
      N = smallnai_;
      O = nai_;
      break;
    case 1:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = G(p) = A(p);
        H(p) = 2 * pow(A(p) * K(p), 2) * kinetic_(p) * E(p);
        dF(p) = dG(p) = dA(p);
        dH(p) = 4 * A(p) * dA(p) * kinetic_(p) * pow(K(p), 2) * E(p) + 2 * pow(A(p) * K(p), 2) * E(p)
              + 4 * pow(A(p), 2) * kinetic_(p) * K(p) * dK(p) * E(p) + 2 * pow(A(p) * K(p), 2) * kinetic_(p) * dE(p);
      }
      N = O = nai_;
      break;
    case 2:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = G(p) = B(p);
        H(p) = pow(B(p) / K(p), 2) * E(p) / (2 * kinetic_(p));
        dF(p) = dG(p) = dB(p);
        dH(p) = ((2 * B(p) * dB(p) * E(p) + pow(B(p), 2) * dE(p)) * kinetic_(p) * K(p)
              + pow(B(p), 2) * E(p) * (K(p) + 2 * kinetic_(p) * dK(p))) / (2 * pow(kinetic_(p) * K(p), 2) * K(p));
      }
      N = O = smallnai_;
      break;
    case 3:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = B(p);
        G(p) = A(p) * E(p);
        H(p) = -B(p) * A(p);
        dF(p) = dB(p);
        dG(p) = dA(p) * E(p) + A(p) * dE(p);
        dH(p) = -dB(p) * A(p) - B(p) * dA(p);
      }
      N = smallnai_;
      O = nai_;
      break;
    case 4:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = A(p);
        G(p) = B(p) * E(p);
        H(p) = -A(p) * B(p);
        dF(p) = dA(p);
        dG(p) = dB(p) * E(p) + B(p) * dE(p);
        dH(p) = -dA(p) * B(p) - A(p) * dB(p);
      }
      N = nai_;
      O = smallnai_;
      break;
    case 5:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = A(p);
        G(p) = A(p) * E(p);
        H(p) = 2 * pow(A(p) * K(p), 2) * kinetic_(p);
        dF(p) = dA(p);
        dG(p) = dA(p) * E(p) + A(p) * dE(p);
        dH(p) = 4 * A(p) * dA(p) * kinetic_(p) * pow(K(p), 2) + 2 * pow(A(p) * K(p), 2)
              + 4 * pow(A(p), 2) * kinetic_(p) * K(p) * dK(p);
      }
      N = O = nai_;
      break;
    case 6:
      for (int p = 0; p != nbasis_; ++p) {
        F(p) = B(p);
        G(p) = B(p) * E(p);
        H(p) = pow(B(p) / K(p), 2) / (2 * kinetic_(p));
        dF(p) = dB(p);
        dG(p) = dB(p) * E(p) + B(p) * dE(p);
        dH(p) = (2 * B(p) * dB(p) * kinetic_(p) * K(p)
              + pow(B(p), 2) * (K(p) + 2 * kinetic_(p) * dK(p))) / (2 * pow(kinetic_(p) * K(p), 2) * K(p));
      }
      N = O = smallnai_;
      break;
    }
    for (int q = 0; q != nbasis_; ++q) {
      for (int p = 0; p != nbasis_; ++p) {
        for (int s = 0; s != nbasis_; ++s) {
          for (int r = 0; r != nbasis_; ++r) {
            ederiv_(p, q) += (CPW(r, p) * H(s) / ((E(q) + E(s)) * (E(r) + E(s)))) * (F(r) * G(q) * N(r, s) * O(s, q) + F(q) * G(r) * N(s, q) * O(r, s));
            ederiv_(p, q) += (CPW(r, q) * H(s) / ((E(q) + E(s)) * (E(r) + E(s)))) * (F(r) * G(q) * N(r, s) * O(s, p) + F(q) * G(r) * N(s, p) * O(r, s));
            if (p == q) {
              ederiv_(p, q) += (CPW(r, p) * H(s) / ((E(p) + E(s)) * (E(r) + E(s)))) * (F(r) * (2 * dG(p) * kinetic_(p)
                            - 2 * G(p) * dE(p) * kinetic_(p) / (E(p) + E(s))) * N(r, s) * O(s, p) + G(r) * (2 * dF(p) * kinetic_(p)
                            - 2 * F(p) * dE(p) * kinetic_(p) / (E(p) + E(s))) * N(s, p) * O(r, s));              
            }
            if (q == s) {
              ederiv_(p, q) += (CPW(r, q) * H(q) / (2 * E(q) * (E(r) + E(q)))) * (F(r) * G(q) * N(r, q) * O(p, q) + F(q) * G(r) * N(p, q) * O(r, q));
              if (p == q) {
                ederiv_(p, q) += (CPW(r, p) / (2 * E(p) * (E(r) + E(p)))) * (F(r) * (G(p) * dH(p) * kinetic_(p)
                              - G(p) * H(p) * dE(p) * kinetic_(p) / (2 * E(p)) - G(p) * H(p) * dE(p) * kinetic_(p) / (E(p)
                              + E(r))) * N(r, p) * O(p, p) + G(r) * (F(p) * dH(p) * kinetic_(p) - F(p) * H(p) * dE(p) * kinetic_(p) / (2 * E(p))
                              - F(p) * H(p) * dE(p) * kinetic_(p) / (E(p) + E(r))) * N(p, p) * O(r, p));
              }
              if (q != r) {
                ederiv_(p, q) += (CPW(r, q) * H(q) / (2 * E(q) * (E(r) + E(q)))) * (F(r) * G(q) * N(r, p) * O(q, q) + F(q) * G(r) * N(q, q) * O(r, p));
                if (p == q) {
                  ederiv_(p, q) += (CPW(r, p) / (2 * E(p) * (E(r) + E(p)))) * (F(r) * (G(p) * dH(p) * kinetic_(p)
                                - G(p) * H(p) * dE(p) * kinetic_(p) / (2 * E(p)) - G(p) * H(p) * dE(p) * kinetic_(p) / (E(p)
                                + E(r))) * N(r, p) * O(p, p) + G(r) * (F(p) * dH(p) * kinetic_(p) - F(p) * H(p) * dE(p) * kinetic_(p) / (2 * E(p))
                                - F(p) * H(p) * dE(p) * kinetic_(p) / (E(p) + E(r))) * N(p, p) * O(r, p));
                }
              }
            }
            if (q != r && q != s) {
              ederiv_(p, q) += (CPW(r, s) * H(q) / ((E(r) + E(q)) * (E(s) + E(q)))) * (F(r) * G(s) * N(r, q) * O(s, p) + F(s) * G(r) * N(r, p) * O(s, q));
              if (p == q) {
                ederiv_(p, q) += (CPW(r, s) / ((E(r) + E(p)) * (E(s) + E(p)))) * (F(r) * G(s) * (dH(p) * kinetic_(p)
                              - 2 * H(p) * dE(p) * kinetic_(p) / (E(p) + E(r))) * N(r, p) * O(s, p) + F(s) * G(r) * (dH(p) * kinetic_(p)
                              - 2 * H(p) * dE(p) * kinetic_(p) / (E(p) + E(s))) * N(r, p) * O(s, p));
              }
            }
          }
        }
      }
    }

    for (int b = 0; b != nbasis_; ++b) {
      for (int a = 0; a != nbasis_; ++a) {
        for (int r = 0; r != nbasis_; ++r) {
          for (int q = 0; q != nbasis_; ++q) {
            for (int p = 0; p != nbasis_; ++p) {
              (*den)(a, b) += (N(p, r) * O(q, r) / ((2 * E(p) + E(r)) * (E(q) + E(r)))) * ((G(q) * dF(p) - F(p) * G(q) * dE(p) / (E(p)
                          + E(r))) * H(r) * wtrans_(a, p) * wtrans_(b, p) + (F(p) * dG(q) - F(p) * G(q) * dE(q) / (E(q)
                          + E(r))) * H(r) * wtrans_(a, q) * wtrans_(b, q) + (F(p) * G(q) * dH(r) - F(p) * G(q) * H(r) * dE(r) / (E(p) + E(r))
                          - F(p) * G(q) * H(r) * dE(r) / (E(q) + E(r))) * wtrans_(a, r) * wtrans_(b, r)) * CPW(p, q);
              (*den)(a, b) += (N(q, r) * O(p, r) / (2 * (E(p) + E(r)) * (E(q) + E(r)))) * ((F(q) * dG(p) - F(q) * G(p) * dE(p) / (E(p)
                          + E(r))) * H(r) * wtrans_(a, p) * wtrans_(b, p) + (G(p) * dF(q) - F(q) * G(p) * dE(q) / (E(q)
                          + E(r))) * H(r) * wtrans_(a, q) * wtrans_(b, q) + (F(q) * G(p) * dH(r) - F(q) * G(p) * H(r) * dE(r) / (E(p) + E(r))
                          - F(q) * G(p) * H(r) * dE(r) / (E(q) + E(r))) * wtrans_(a, r) * wtrans_(b, r)) * CPW(p, q);
            }
          }
        }
      }
    }
  }

  for (int q = 0; q != nbasis_; ++q) {
    for (int p = 0; p != nbasis_; ++p) {
      zmult_(p, q) = p == q ? 0 : -0.5 * (ederiv_(p, q) - ederiv_(q, p)) / (kinetic_(p) - kinetic_(q));
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
  VectorB E(nbasis_), A(nbasis_), K(nbasis_), B(nbasis_);
  for (int p = 0; p != nbasis_; ++p) {
    E(p) = c__ * sqrt(2 * kinetic_(p) + c2);
    A(p) = sqrt(0.5 * (E(p) + c2) / E(p));
    K(p) = c__ / (E(p) + c2);
    B(p) = A(p) * K(p);
  }

  const Matrix CPW = (ptrans_ % wtrans_rev_) % *rdm1 * (ptrans_ % wtrans_rev_);
  shared_ptr<Matrix> den = make_shared<Matrix>(nbasis_, nbasis_), pvpden = make_shared<Matrix>(nbasis_, nbasis_);
  for (int b = 0; b != nbasis_; ++b) {
    for (int a = 0; a != nbasis_; ++a) {
      for (int q = 0; q != nbasis_; ++q) {
        for (int p = 0; p != nbasis_; ++p) {
          (*den)(a, b) += A(p) * A(q) * wtrans_(a, p) * wtrans_(b, q) * CPW(p, q);
          (*pvpden)(a, b) += B(p) * B(q) * wtrans_(a, p) * wtrans_(b, q) * CPW(p, q);
          for (int r = 0; r != nbasis_; ++r) {
            (*den)(a, b) -= (B(r) * E(r) * A(r) / ((E(p) + E(r)) * (E(q) + E(r)))) * (B(p) * A(q) * smallnai_(r, p) * wtrans_(a, q) * wtrans_(b, r)
                        + B(q) * A(p) * smallnai_(r, q) * wtrans_(a, p) * wtrans_(b, r)) * CPW(p, q);
            (*pvpden)(a, b) -= (B(r) * E(r) * A(r) / ((E(p) + E(r)) * (E(q) + E(r)))) * (B(p) * A(q) * nai_(r, q) * wtrans_(a, p) * wtrans_(b, r)
                            + B(q) * A(p) * nai_(r, p) * wtrans_(a, q) * wtrans_(b, r)) * CPW(p, q);
            (*den)(a, b) += (pow(A(r) * K(r), 2) * kinetic_(r) * E(r) / ((E(p) + E(r)) * (E(q)
                        + E(r)))) * (A(p) * A(q) * nai_(r, p) * wtrans_(a, q) * wtrans_(b, r)
                        + A(q) * A(p) * nai_(r, q) * wtrans_(a, p) * wtrans_(b, r)) * CPW(p, q);
            (*den)(a, b) += (pow(A(r) * K(r), 2) * kinetic_(r) * E(r) / ((E(p) + E(r)) * (E(q)
                        + E(r)))) * (A(p) * A(q) * nai_(r, q) * wtrans_(a, p) * wtrans_(b, r)
                        + A(q) * A(p) * nai_(r, p) * wtrans_(a, q) * wtrans_(b, r)) * CPW(p, q);
            (*pvpden)(a, b) += ((pow(B(r) / K(r), 2) * E(r) / kinetic_(r)) / (4 * (E(p) + E(r)) * (E(q)
                            + E(r)))) * (B(p) * B(q) * smallnai_(r, p) * wtrans_(a, q) * wtrans_(b, r)
                            + B(q) * B(p) * smallnai_(r, q) * wtrans_(a, p) * wtrans_(b, r)) * CPW(p, q);
            (*pvpden)(a, b) += ((pow(B(r) / K(r), 2) * E(r) / kinetic_(r)) / (4 * (E(p) + E(r)) * (E(q)
                            + E(r)))) * (B(p) * B(q) * smallnai_(r, q) * wtrans_(a, p) * wtrans_(b, r)
                            + B(q) * B(p) * smallnai_(r, p) * wtrans_(a, q) * wtrans_(b, r)) * CPW(p, q);
            (*den)(a, b) -= (B(r) * A(r) / (2 * (E(p) + E(r)) * (E(q) + E(r)))) * (B(p) * A(q) * E(q) * smallnai_(r, p) * wtrans_(a, q) * wtrans_(b, r)
                        + B(q) * A(p) * E(p) * smallnai_(r, q) * wtrans_(a, p) * wtrans_(b, r)) * CPW(p, q);
            (*pvpden)(a, b) -= (B(r) * A(r) / (2 * (E(p) + E(r)) * (E(q) + E(r)))) * (B(p) * A(q) * E(q) * nai_(r, q) * wtrans_(a, p) * wtrans_(b, r)
                            + B(q) * A(p) * E(p) * nai_(r, p) * wtrans_(a, q) * wtrans_(b, r)) * CPW(p, q);
            (*pvpden)(a, b) -= (A(r) * B(r) / (2 * (E(p) + E(r)) * (E(q) + E(r)))) * (A(p) * B(q) * E(q) * nai_(r, p) * wtrans_(a, q) * wtrans_(b, r)
                            + A(q) * B(p) * E(p) * nai_(r, q) * wtrans_(a, p) * wtrans_(b, r)) * CPW(p, q);
            (*den)(a, b) -= (A(r) * B(r) / (2 * (E(p) + E(r)) * (E(q) + E(r)))) * (A(p) * B(q) * E(q) * smallnai_(r, q) * wtrans_(a, p) * wtrans_(b, r)
                        + A(q) * B(p) * E(p) * nai_(r, p) * wtrans_(a, q) * wtrans_(b, r)) * CPW(p, q);
            (*den)(a, b) += (pow(A(r) * K(r), 2) * kinetic_(r) / ((E(p) + E(r)) * (E(q)
                        + E(r)))) * (A(p) * A(q) * E(q) * nai_(r, p) * wtrans_(a, q) * wtrans_(b, r)
                        + A(q) * A(p) * E(p) * nai_(r, q) * wtrans_(a, p) * wtrans_(b, r)) * CPW(p, q);
            (*den)(a, b) += (pow(A(r) * K(r), 2) * kinetic_(r) / ((E(p) + E(r)) * (E(q)
                        + E(r)))) * (A(p) * A(q) * E(q) * nai_(r, q) * wtrans_(a, p) * wtrans_(b, r)
                        + A(q) * A(p) * E(p) * nai_(r, p) * wtrans_(a, q) * wtrans_(b, r)) * CPW(p, q);
            (*pvpden)(a, b) += ((pow(B(r) / K(r), 2) / kinetic_(r)) / (4 * (E(p) + E(r)) * (E(q)
                            + E(r)))) * (B(p) * B(q) * E(q) * smallnai_(r, p) * wtrans_(a, q) * wtrans_(b, r)
                            + B(q) * B(p) * E(p) * smallnai_(r, q) * wtrans_(a, p) * wtrans_(b, r)) * CPW(p, q);
            (*pvpden)(a, b) += ((pow(B(r) / K(r), 2) / kinetic_(r)) / (4 * (E(p) + E(r)) * (E(q)
                            + E(r)))) * (B(p) * B(q) * E(q) * smallnai_(r, q) * wtrans_(a, p) * wtrans_(b, r)
                            + B(q) * B(p) * E(p) * smallnai_(r, p) * wtrans_(a, q) * wtrans_(b, r)) * CPW(p, q);
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
  VectorB E(nbasis_), A(nbasis_), K(nbasis_), B(nbasis_);
  for (int p = 0; p != nbasis_; ++p) {
    E(p) = c__ * sqrt(2 * kinetic_(p) + c2);
    A(p) = sqrt(0.5 * (E(p) + c2) / E(p));
    K(p) = c__ / (E(p) + c2);
    B(p) = A(p) * K(p);
  }

  Matrix at(nbasis_, nbasis_);
  for (int q = 0; q != nbasis_; ++q) {
    for (int p = 0; p != nbasis_; ++p) {
      at(p, q) = (zmult_(p, q) + zmult_(q, p)) * kinetic_(p);
    }
  }

  Matrix xb(nbasis_, nbasis_);
  for (int q = 0; q != nbasis_; ++q) {
    for (int p = 0; p != nbasis_; ++p) {
      xb(p, q) = 0.25 * (ederiv_(p, q) + ederiv_(q, p) + at(p, q) + at(q, p));
    }
  }

  shared_ptr<Matrix> den = make_shared<Matrix>(nbasis_, nbasis_);
  const Matrix CPW = ptrans_ * *rdm1 ^ (wtrans_rev_ % ptrans_);
  for (int b = 0; b != nbasis_; ++b) {
    for (int a = 0; a != nbasis_; ++a) {
      for (int p = 0; p != nbasis_; ++p) {
        (*den)(a, b) += 2 * (E(p) - c2) * wtrans_(b, p) * CPW(a, p);
      }
      for (int q = 0; q != nbasis_; ++q) {
        for (int p = 0; p != nbasis_; ++p) {
          (*den)(a, b) += (A(p) * nai_(p, q) * A(q) + B(p) * smallnai_(p, q) * B(q)) * (wtrans_(b, q) * CPW(a, p) + wtrans_(b, p) * CPW(a, q));
          for (int r = 0; r != nbasis_; ++r) {
            (*den)(a, b) += (1 / ((E(p) + E(r)) * (E(q) + E(r)))) * (-B(p) * smallnai_(r, p) * B(r) * E(r) * A(r) * nai_(r, q) * A(q)
                        - A(p) * nai_(r, p) * A(r) * E(r) * B(r) * smallnai_(r, q) * B(q)
                        + 2 * A(p) * nai_(r, p) * pow(A(r) * K(r), 2) * kinetic_(r) * E(r) * nai_(r, q) * A(q)
                        + 0.5 * B(p) * smallnai_(r, p) * (pow(B(r) / K(r), 2) * E(p) / kinetic_(r)) * smallnai_(r, q) * B(q)
                        - 0.5 * B(p) * smallnai_(r, p) * B(r) * A(r) * nai_(r, q) * A(q) * E(q)
                        - 0.5 * A(p) * nai_(r, p) * A(r) * B(r) * smallnai_(r, q) * B(q) * E(q)
                        + A(p) * nai_(r, p) * pow(A(r) * K(r), 2) * kinetic_(r) * nai_(r, q) * A(q) * E(q)
                        + 0.25 * B(p) * smallnai_(r, p) * (pow(B(r) / K(r), 2) / kinetic_(r)) * smallnai_(r, q) * B(q) * E(q)
                        - 0.5 * E(p) * B(p) * smallnai_(r, p) * B(r) * A(r) * nai_(r, q) * A(q)
                        - 0.5 * E(p) * A(p) * nai_(r, p) * A(r) * B(r) * smallnai_(r, q) * B(q)
                        + E(p) * A(p) * nai_(r, p) * pow(A(r) * K(r), 2) * kinetic_(r) * nai_(r, q) * A(q)
                        + 0.25 * E(p) * B(p) * smallnai_(r, p) * (pow(B(r) / K(r), 2) / kinetic_(r)) * smallnai_(r, q) * B(q))
                        * (wtrans_(b, q) * CPW(a, p) + wtrans_(b, p) * CPW(a, q));
          }
        }
      }
    }
  }

  *den -= (ptrans_ * *erdm1 ^ ptrans_) + (wtrans_ * xb ^ wtrans_);
  at.print("a_tilde");
  xb.print("X_bar");
  den->print("X_tilde");
  return den;
}
