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
  const double c2 = c__ * c__, c3 = c2 * c__;
  VectorB E(nbasis_), A(nbasis_), B(nbasis_), K(nbasis_), dE(nbasis_), dA(nbasis_), dB(nbasis_);
  for (int p = 0; p != nbasis_; ++p) {
    E(p) = c__ * sqrt(2.0 * kinetic_(p) + c2);
    A(p) = sqrt(0.5 * (E(p) + c2) / E(p));
    K(p) = c__ / (E(p) + c2);
    B(p) = A(p) * K(p);
    dE(p) = c__ / sqrt(2.0 * kinetic_(p) + c2);
    dA(p) = -c3 / (4.0 * E(p) * E(p) * A(p) * sqrt(2.0 * kinetic_(p) + c2));
    dB(p) = (c__ * K(p) / (4.0 * E(p) * E(p) * A(p)) - A(p) / pow(E(p) + c2, 2)) * c2 / sqrt(2.0 * kinetic_(p) + c2);
  }

  const Matrix CPW = (ptrans_ % wtrans_rev_) % *rdm1 * (ptrans_ % wtrans_rev_);
  for (int q = 0; q != nbasis_; ++q) {
    for (int p = 0; p != nbasis_; ++p) {
      ederiv_(p, q) += 2.0 * CPW(p, q) * (E(q) - c2);
      if (p == q) {
        ederiv_(p, q) += 2.0 * CPW(p, p) * dE(p) * kinetic_(p);
      }
      for (int r = 0; r != nbasis_; ++r) {
        ederiv_(p, q) += 2.0 * CPW(r, p) * (A(q) * A(r) * nai_(r, q) + B(q) * B(r) * smallnai_(r, q));
        ederiv_(p, q) += 2.0 * CPW(r, q) * (A(q) * A(r) * nai_(r, p) + B(q) * B(r) * smallnai_(r, p));
        if (p == q) {
          ederiv_(p, q) += 2.0 * CPW(r, p) * (dA(p) * kinetic_(p) * A(r) * nai_(r, q) + dB(p) * kinetic_(p) * B(r) * smallnai_(r, q));
          ederiv_(p, q) += 2.0 * CPW(r, q) * (dA(p) * kinetic_(p) * A(r) * nai_(r, p) + dB(p) * kinetic_(p) * B(r) * smallnai_(r, p));
        }
      }
    }
  }

  for (int q = 0; q != nbasis_; ++q) {
    for (int p = 0; p != nbasis_; ++p) {
      zmult_(p, q) = p == q ? 0.0 : -0.5 * (ederiv_(p, q) - ederiv_(q, p)) / (kinetic_(p) - kinetic_(q));
    }
  }

  shared_ptr<Matrix> den = make_shared<Matrix>(wtrans_ * zmult_ ^ wtrans_);
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
  ederiv_.print("Y_pq");
  zmult_.print("z_pq");
  den->print("d_tilde");
  return den;

  // const double c2 = c__ * c__, c3 = c2 * c__;
  // VectorB E(nbasis_), dE(nbasis_);
  // for (int p = 0; p != nbasis_; ++p) {
  //   E(p) = c__ * sqrt(2.0 * kinetic_(p) + c2);
  //   dE(p) = c__ / sqrt(2.0 * kinetic_(p) + c2);
  // }

  // const Matrix CPW = (ptrans_ % wtrans_rev_) % *rdm1 * (ptrans_ % wtrans_rev_);
  // for (int q = 0; q != nbasis_; ++q) {
  //   for (int p = 0; p != nbasis_; ++p) {
  //     ederiv_(p, q) = 2.0 * CPW(p, q) * (E(q) - c2);
  //     if (p == q) {
  //       ederiv_(p, q) += 2.0 * CPW(p, p) * dE(p) * kinetic_(p);
  //     }
  //   }
  // }

  // for (int q = 0; q != nbasis_; ++q) {
  //   for (int p = 0; p != nbasis_; ++p) {
  //     zmult_(p, q) = p == q ? 0.0 : -0.5 * (ederiv_(p, q) - ederiv_(q, p)) / (kinetic_(p) - kinetic_(q));
  //   }
  // }

  // shared_ptr<Matrix> den = make_shared<Matrix>(wtrans_ * zmult_ ^ wtrans_);
  // for (int b = 0; b != nbasis_; ++b) {
  //   for (int a = 0; a != nbasis_; ++a) {
  //     for (int p = 0; p != nbasis_; ++p) {
  //       (*den)(a, b) += dE(p) * wtrans_(a, p) * wtrans_(b, p) * CPW(p, p);
  //     }
  //   }
  // }
  // ederiv_.print("Y_pq");
  // zmult_.print("z_pq");
  // den->print("d_tilde");
  // return den;
}

array<shared_ptr<const Matrix>, 2> DKHcoreInfo::compute_vden(shared_ptr<const Matrix> rdm1) {
  const double c2 = c__ * c__;
  VectorB E(nbasis_), A(nbasis_), K(nbasis_), B(nbasis_);
  for (int p = 0; p != nbasis_; ++p) {
    E(p) = c__ * sqrt(2.0 * kinetic_(p) + c2);
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
        }
      }
    }
  }
  den->print("d_bar");
  pvpden->print("d_check");
  return { den, pvpden };

  // return { make_shared<Matrix>(ptrans_ * *rdm1 ^ ptrans_), make_shared<Matrix>(nbasis_, nbasis_) };
}

shared_ptr<const Matrix> DKHcoreInfo::compute_sden(shared_ptr<const Matrix> rdm1, shared_ptr<const Matrix> erdm1) {
  const double c2 = c__ * c__;
  VectorB E(nbasis_), A(nbasis_), K(nbasis_), B(nbasis_);
  for (int p = 0; p != nbasis_; ++p) {
    E(p) = c__ * sqrt(2.0 * kinetic_(p) + c2);
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

  shared_ptr<Matrix> den = make_shared<Matrix>((ptrans_ * *erdm1 ^ ptrans_) + (wtrans_ * xb ^ wtrans_));
  const Matrix CPW = ptrans_ * *rdm1 ^ (wtrans_rev_ % ptrans_);
  for (int b = 0; b != nbasis_; ++b) {
    for (int a = 0; a != nbasis_; ++a) {
      for (int p = 0; p != nbasis_; ++p) {
        (*den)(a, b) -= 2 * (E(p) - c2) * wtrans_(b, p) * CPW(a, p);
      }
      for (int q = 0; q != nbasis_; ++q) {
        for (int p = 0; p != nbasis_; ++p) {
          (*den)(a, b) -= (A(p) * nai_(p, q) * A(q) + B(p) * smallnai_(p, q) * B(q)) * (wtrans_(b, q) * CPW(a, p) + wtrans_(b, p) * CPW(a, q));
        }
      }
    }
  }
  at.print("a_tilde");
  xb.print("X_bar");
  den->print("X_tilde");
  return den;

  // const double c2 = c__ * c__;
  // VectorB E(nbasis_);
  // for (int p = 0; p != nbasis_; ++p) {
  //   E(p) = c__ * sqrt(2.0 * kinetic_(p) + c2);
  // }

  // Matrix at(nbasis_, nbasis_);
  // for (int q = 0; q != nbasis_; ++q) {
  //   for (int p = 0; p != nbasis_; ++p) {
  //     at(p, q) = (zmult_(p, q) + zmult_(q, p)) * kinetic_(p);
  //   }
  // }

  // Matrix xb(nbasis_, nbasis_);
  // for (int q = 0; q != nbasis_; ++q) {
  //   for (int p = 0; p != nbasis_; ++p) {
  //     xb(p, q) = 0.25 * (ederiv_(p, q) + ederiv_(q, p) + at(p, q) + at(q, p));
  //   }
  // }

  // shared_ptr<Matrix> den = make_shared<Matrix>((ptrans_ * *erdm1 ^ ptrans_) + (wtrans_ * xb ^ wtrans_));
  // const Matrix CPW = ptrans_ * *rdm1 ^ (wtrans_rev_ % ptrans_);
  // for (int b = 0; b != nbasis_; ++b) {
  //   for (int a = 0; a != nbasis_; ++a) {
  //     for (int p = 0; p != nbasis_; ++p) {
  //       (*den)(a, b) -= 2 * (E(p) - c2) * wtrans_(b, p) * CPW(a, p);
  //     }
  //   }
  // }
  // at.print("a_tilde");
  // xb.print("X_bar");
  // den->print("X_tilde");
  // return den;
}
