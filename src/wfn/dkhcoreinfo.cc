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
#include <src/integral/os/overlapbatch.h>
#include <src/mat1e/mixedbasis.h>

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
  kinetic_ = VectorB(mol->nbasis());
  lambda->diagonalize(kinetic_);
  wtrans_ = *gamma * *lambda;

  const NAI nai(mol);
  nai_ = wtrans_ % nai * wtrans_;
  const Small1e<NAIBatch> small1e(mol);
  smallnai_ = wtrans_ % small1e[0] * wtrans_;

  ptrans_ = MixedBasis<OverlapBatch>(current, mol);

  zmult_ = ederiv_ = Matrix(nbasis_, nbasis_);
}

shared_ptr<const Matrix> DKHcoreInfo::compute_tden(shared_ptr<const Matrix> rdm1) {
  // const double c2 = c__ * c__, c3 = c2 * c__;
  // VectorB E(nbasis_), A(nbasis_), B(nbasis_), K(nbasis_), dE(nbasis_), dA(nbasis_), dB(nbasis_);
  // for (int p = 0; p != nbasis_; ++p) {
  //   E(p) = c__ * sqrt(2.0 * kinetic_(p) + c2);
  //   A(p) = sqrt(0.5 * (E(p) + c2) / E(p));
  //   K(p) = c__ / (E(p) + c2);
  //   B(p) = A(p) * K(p);
  //   dE(p) = c__ / sqrt(2.0 * kinetic_(p) + c2);
  //   dA(p) = -c3 / (4.0 * E(p) * E(p) * A(p) * sqrt(2.0 * kinetic_(p) + c2));
  //   dB(p) = (c__ * K(p) / (4.0 * E(p) * E(p) * A(p)) - A(p) / pow(E(p) + c2, 2)) * c2 / sqrt(2.0 * kinetic_(p) + c2);
  // }

  // const Matrix CPW = (ptrans_ % wtrans_) % *rdm1 * (ptrans_ % wtrans_);
  // for (int p = 0; p != nbasis_; ++p) {
  //   for (int q = 0; q != nbasis_; ++q) {
  //     ederiv_(p, q) += 2.0 * CPW(p, q) * (E(q) - c2);
  //     if (p == q) {
  //       ederiv_(p, q) += 2.0 * CPW(p, p) * dE(p) * kinetic_(p);
  //     }
  //     for (int r = 0; r != nbasis_; ++r) {
  //       ederiv_(p, q) += 2.0 * CPW(p, r) * (A(q) * A(r) * (nai_(p, r) + nai_(r, q)) + B(q) * B(r) * (smallnai_(p, r) + smallnai_(r, q)));
  //       if (p == q) {
  //         ederiv_(p, q) += 2.0 * CPW(p, r) * (dA(p) * kinetic_(p) * A(r) * (nai_(p, r) + nai_(r, q)) + dB(p) * kinetic_(p) * B(r) * (smallnai_(p, r) + smallnai_(r, q)));
  //       }
  //     }
  //   }
  // }

  // for (int p = 0; p != nbasis_; ++p) {
  //   for (int q = 0; q != nbasis_; ++q) {
  //     zmult_(p, q) = p == q ? 0.0 : -0.5 * (ederiv_(p, q) - ederiv_(q, p)) / (kinetic_(p) - kinetic_(q));
  //   }
  // }

  // shared_ptr<Matrix> den = make_shared<Matrix>(wtrans_ * zmult_ ^ wtrans_);
  // for (int a = 0; a != nbasis_; ++a) {
  //   for (int b = 0; b != nbasis_; ++b) {
  //     for (int p = 0; p != nbasis_; ++p) {
  //       (*den)(a, b) += dE(p) * wtrans_(a, p) * wtrans_(b, p) * CPW(p, p);
  //     }
  //     for (int p = 0; p != nbasis_; ++p) {
  //       for (int q = 0; q != nbasis_; ++q) {
  //         (*den)(a, b) += (nai_(p, q) * (A(q) * dA(p) * wtrans_(a, p) * wtrans_(b, p) + A(p) * dA(q) * wtrans_(a, q) * wtrans_(b, q))
  //                     + smallnai_(p, q) * (B(q) * dB(p) * wtrans_(a, p) * wtrans_(b, p) + B(p) * dB(q) * wtrans_(a, q) * wtrans_(b, q))) * CPW(p, q);
  //       }
  //     }
  //   }
  // }
  // ederiv_.print("Y_pq");
  // zmult_.print("z_pq");
  // den->print("d_tilde");
  // return den;

  // return make_shared<Matrix>(ptrans_ * *rdm1 ^ ptrans_);

  const double c2 = c__ * c__, c3 = c2 * c__;
  VectorB E(nbasis_), A(nbasis_), B(nbasis_), K(nbasis_), dE(nbasis_), dA(nbasis_), dB(nbasis_);
  for (int p = 0; p != nbasis_; ++p) {
    E(p) = c__ * sqrt(2.0 * kinetic_(p) + c2);
    dE(p) = c__ / sqrt(2.0 * kinetic_(p) + c2);
  }

  const Matrix CPW = (ptrans_ % wtrans_) % *rdm1 * (ptrans_ % wtrans_);
  for (int q = 0; q != nbasis_; ++q) {
    for (int p = 0; p != nbasis_; ++p) {
      ederiv_(p, q) = 2.0 * CPW(p, q) * (E(q) - c2);
      if (p == q) {
        ederiv_(p, q) += 2.0 * CPW(p, p) * dE(p) * kinetic_(p);
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
    }
  }
  ederiv_.print("Y_pq");
  zmult_.print("z_pq");
  den->print("d_tilde");
  return den;
}

shared_ptr<const Matrix> DKHcoreInfo::compute_vden(shared_ptr<const Matrix> rdm1) {
  // const double c2 = c__ * c__;
  // VectorB E(nbasis_), A(nbasis_);
  // for (int p = 0; p != nbasis_; ++p) {
  //   E(p) = c__ * sqrt(2.0 * kinetic_(p) + c2);
  //   A(p) = sqrt(0.5 * (E(p) + c2) / E(p));
  // }

  // const Matrix CPW = (ptrans_ % wtrans_) % *rdm1 * (ptrans_ % wtrans_);
  // shared_ptr<Matrix> den = make_shared<Matrix>(nbasis_, nbasis_);
  // for (int a = 0; a != nbasis_; ++a) {
  //   for (int b = 0; b != nbasis_; ++b) {
  //     for (int p = 0; p != nbasis_; ++p) {
  //       for (int q = 0; q != nbasis_; ++q) {
  //         (*den)(a, b) += A(p) * A(q) * wtrans_(a, p) * wtrans_(b, q) * CPW(p, q);
  //       }
  //     }
  //   }
  // }
  // den->print("d_bar");
  // return den;

  return make_shared<Matrix>(ptrans_ * *rdm1 ^ ptrans_);
}

shared_ptr<const Matrix> DKHcoreInfo::compute_pvpden(shared_ptr<const Matrix> rdm1) {
  // const double c2 = c__ * c__;
  // VectorB E(nbasis_), A(nbasis_), B(nbasis_), K(nbasis_);
  // for (int p = 0; p != nbasis_; ++p) {
  //   E(p) = c__ * sqrt(2.0 * kinetic_(p) + c2);
  //   A(p) = sqrt(0.5 * (E(p) + c2) / E(p));
  //   K(p) = c__ / (E(p) + c2);
  //   B(p) = A(p) * K(p);
  // }

  // const Matrix CPW = (ptrans_ % wtrans_) % *rdm1 * (ptrans_ % wtrans_);
  // shared_ptr<Matrix> den = make_shared<Matrix>(nbasis_, nbasis_);
  // for (int a = 0; a != nbasis_; ++a) {
  //   for (int b = 0; b != nbasis_; ++b) {
  //     for (int p = 0; p != nbasis_; ++p) {
  //       for (int q = 0; q != nbasis_; ++q) {
  //         (*den)(a, b) += B(p) * B(q) * wtrans_(a, p) * wtrans_(b, q) * CPW(p, q);
  //       }
  //     }
  //   }
  // }
  // den->print("d_check");
  // return den;

  return make_shared<Matrix>(nbasis_, nbasis_);
}

shared_ptr<const Matrix> DKHcoreInfo::compute_sden(shared_ptr<const Matrix> erdm1) {
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

  at.print("a_tilde");
  xb.print("X_bar");
  Matrix dldu = ederiv_ + at - 2.0 * xb;
  dldu.print("dL/dU_pq");
  return make_shared<const Matrix>((ptrans_ * *erdm1 ^ ptrans_) + (wtrans_ * xb ^ wtrans_));

  // return make_shared<Matrix>(ptrans_ * *erdm1 ^ ptrans_);
}
