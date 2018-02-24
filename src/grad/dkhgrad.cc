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

  shared_ptr<const Matrix> wmat, wmat_rev;
  shared_ptr<const DiagMatrix> t;
  {
    // Compute momentum transformation matrix
    const Overlap overlap(unc);
    shared_ptr<const Matrix> xmat = overlap.tildex();
    // Kinetic operator is diagonal in momentum space
    const Kinetic kinetic(unc);
    auto gamma = make_shared<Matrix>(*xmat % kinetic * *xmat);
    VectorB t0(gamma->ndim());
    gamma->diagonalize(t0);
    t = make_shared<const DiagMatrix>(t0);
    // Transformation matrix to momentum space
    wmat = make_shared<const Matrix>(*xmat * *gamma);
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
  auto pmat = make_shared<const ContractMat>(mol_, unc->nbasis());

  // Lagrange multiplier for diagonalization of kinetic operator and energy derivative with respect to momentum orbital rotation
  shared_ptr<const Matrix> tden, zpq, ypq;
  tie(tden, zpq, ypq) = compute_tden_(rdm1, pmat, wmat, wmat_rev, t, v, pvp); 

  shared_ptr<const Matrix> sden = compute_sden_(rdm1, erdm1, pmat, wmat, wmat_rev, t, v, pvp, zpq, ypq);
  shared_ptr<const Matrix> vden0, vden1;
  tie(vden0, vden1) = compute_vden_(rdm1, pmat, wmat, wmat_rev, t, v, pvp);
  return { tden, vden0, vden1, sden };
}


namespace {
struct TempArrays { 
  array<shared_ptr<const Matrix>,12> N;
  array<shared_ptr<const Matrix>,12> O;
  array<DiagMatrix,12> F;
  array<DiagMatrix,12> G;
  array<DiagMatrix,12> H;
  array<DiagMatrix,12> dF;
  array<DiagMatrix,12> dG;
  array<DiagMatrix,12> dH;
  array<pair<int,int>,12> vint;
  DiagMatrix E;
  DiagMatrix A;
  DiagMatrix B;
  DiagMatrix K;
  DiagMatrix dE;
  DiagMatrix dA;
  DiagMatrix dB;
  DiagMatrix dK;

  TempArrays(shared_ptr<const DiagMatrix> t, shared_ptr<const Matrix> v, shared_ptr<const Matrix> pvp)
    : E(v->ndim()), A(v->ndim()), B(v->ndim()), K(v->ndim()), dE(v->ndim()), dA(v->ndim()), dB(v->ndim()), dK(v->ndim()) {

    // Number of uncontracted basis functions
    const int nbasis = t->ndim();
    const double c2 = c__ * c__;
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

    for (int i = 0; i != 12; ++i) { // loop over terms
      F[i] = DiagMatrix(nbasis);
      G[i] = DiagMatrix(nbasis);
      H[i] = DiagMatrix(nbasis);
      dF[i] = DiagMatrix(nbasis);
      dG[i] = DiagMatrix(nbasis);
      dH[i] = DiagMatrix(nbasis);
      switch (i) {
      case 0:
        for (int p = 0; p != nbasis; ++p) {
          F[i](p) = B(p);
          G[i](p) = A(p);
          H[i](p) = -B(p) * E(p) * A(p);
          dF[i](p) = dB(p);
          dG[i](p) = dA(p);
          dH[i](p) = -dB(p) * E(p) * A(p) - B(p) * dE(p) * A(p) - B(p) * E(p) * dA(p);
        }
        N[i] = pvp;
        O[i] = v;
        vint[i] = make_pair(1, 0);
        break;
      case 1:
        for (int p = 0; p != nbasis; ++p) {
          F[i](p) = A(p);
          G[i](p) = B(p);
          H[i](p) = -A(p) * E(p) * B(p);
          dF[i](p) = dA(p);
          dG[i](p) = dB(p);
          dH[i](p) = -dA(p) * E(p) * B(p) - A(p) * dE(p) * B(p) - A(p) * E(p) * dB(p);
        }
        N[i] = v;
        O[i] = pvp;
        vint[i] = make_pair(0, 1);
        break;
      case 2:
        for (int p = 0; p != nbasis; ++p) {
          F[i](p) = G[i](p) = A(p);
          H[i](p) = 2.0 * pow(A(p) * K(p), 2) * (*t)(p) * E(p);
          dF[i](p) = dG[i](p) = dA(p);
          dH[i](p) = 4.0 * A(p) * dA(p) * (*t)(p) * pow(K(p), 2) * E(p) + 2.0 * pow(A(p) * K(p), 2) * E(p)
                + 4.0 * pow(A(p), 2) * (*t)(p) * K(p) * dK(p) * E(p) + 2.0 * pow(A(p) * K(p), 2) * (*t)(p) * dE(p);
        }
        N[i] = O[i] = v;
        vint[i] = make_pair(0, 0);
        break;
      case 3:
        for (int p = 0; p != nbasis; ++p) {
          F[i](p) = G[i](p) = B(p);
          H[i](p) = pow(B(p) / K(p), 2) * E(p) / (2.0 * (*t)(p));
          dF[i](p) = dG[i](p) = dB(p);
          dH[i](p) = ((2.0 * B(p) * dB(p) * E(p) + pow(B(p), 2) * dE(p)) * (*t)(p) * K(p)
                - pow(B(p), 2) * E(p) * (K(p) + 2.0 * (*t)(p) * dK(p))) / (2.0 * pow((*t)(p) * K(p), 2) * K(p));
        }
        N[i] = O[i] = pvp;
        vint[i] = make_pair(1, 1);
        break;
      case 4:
        for (int p = 0; p != nbasis; ++p) {
          F[i](p) = B(p);
          G[i](p) = A(p) * E(p);
          H[i](p) = -B(p) * A(p) / 2.0;
          dF[i](p) = dB(p);
          dG[i](p) = dA(p) * E(p) + A(p) * dE(p);
          dH[i](p) = -dB(p) * A(p) / 2.0 - B(p) * dA(p) / 2.0;
        }
        N[i] = pvp;
        O[i] = v;
        vint[i] = make_pair(1, 0);
        break;
      case 5:
        for (int p = 0; p != nbasis; ++p) {
          F[i](p) = A(p);
          G[i](p) = B(p) * E(p);
          H[i](p) = -A(p) * B(p) / 2.0;
          dF[i](p) = dA(p);
          dG[i](p) = dB(p) * E(p) + B(p) * dE(p);
          dH[i](p) = -dA(p) * B(p) / 2.0 - A(p) * dB(p) / 2.0;
        }
        N[i] = v;
        O[i] = pvp;
        vint[i] = make_pair(0, 1);
        break;
      case 6:
        for (int p = 0; p != nbasis; ++p) {
          F[i](p) = A(p);
          G[i](p) = A(p) * E(p);
          H[i](p) = pow(A(p) * K(p), 2) * (*t)(p);
          dF[i](p) = dA(p);
          dG[i](p) = dA(p) * E(p) + A(p) * dE(p);
          dH[i](p) = 2.0 * A(p) * dA(p) * (*t)(p) * pow(K(p), 2) + pow(A(p) * K(p), 2)
                + 2.0 * pow(A(p), 2) * (*t)(p) * K(p) * dK(p);
        }
        N[i] = O[i] = v;
        vint[i] = make_pair(0, 0);
        break;
      case 7:
        for (int p = 0; p != nbasis; ++p) {
          F[i](p) = B(p);
          G[i](p) = B(p) * E(p);
          H[i](p) = pow(B(p) / K(p), 2) / (4.0 * (*t)(p));
          dF[i](p) = dB(p);
          dG[i](p) = dB(p) * E(p) + B(p) * dE(p);
          dH[i](p) = (2.0 * B(p) * dB(p) * (*t)(p) * K(p)
                - pow(B(p), 2) * (K(p) + 2.0 * (*t)(p) * dK(p))) / (4.0 * pow((*t)(p) * K(p), 2) * K(p));
        }
        N[i] = O[i] = pvp;
        vint[i] = make_pair(1, 1);
        break;
      case 8:
        for (int p = 0; p != nbasis; ++p) {
          F[i](p) = E(p) * B(p);
          G[i](p) = A(p);
          H[i](p) = -B(p) * A(p) / 2.0;
          dF[i](p) = dE(p) * B(p) + E(p) * dB(p);
          dG[i](p) = dA(p);
          dH[i](p) = -dB(p) * A(p) / 2.0 - B(p) * dA(p) / 2.0;
        }
        N[i] = pvp;
        O[i] = v;
        vint[i] = make_pair(1, 0);
        break;
      case 9:
        for (int p = 0; p != nbasis; ++p) {
          F[i](p) = E(p) * A(p);
          G[i](p) = B(p);
          H[i](p) = -A(p) * B(p) / 2.0;
          dF[i](p) = dE(p) * A(p) + E(p) * dA(p);
          dG[i](p) = dB(p);
          dH[i](p) = -dA(p) * B(p) / 2.0 - A(p) * dB(p) / 2.0;
        }
        N[i] = v;
        O[i] = pvp;
        vint[i] = make_pair(0, 1);
        break;
      case 10:
        for (int p = 0; p != nbasis; ++p) {
          F[i](p) = E(p) * A(p);
          G[i](p) = A(p);
          H[i](p) = pow(A(p) * K(p), 2) * (*t)(p);
          dF[i](p) = dE(p) * A(p) + E(p) * dA(p);
          dG[i](p) = dA(p);
          dH[i](p) = 2.0 * A(p) * dA(p) * (*t)(p) * pow(K(p), 2) + pow(A(p) * K(p), 2)
                + 2.0 * pow(A(p), 2) * (*t)(p) * K(p) * dK(p);
        }
        N[i] = O[i] = v;
        vint[i] = make_pair(0, 0);
        break;
      case 11:
        for (int p = 0; p != nbasis; ++p) {
          F[i](p) = E(p) * B(p);
          G[i](p) = B(p);
          H[i](p) = pow(B(p) / K(p), 2) / (4.0 * (*t)(p));
          dF[i](p) = dE(p) * B(p) + E(p) * dB(p);
          dG[i](p) = dB(p);
          dH[i](p) = (2.0 * B(p) * dB(p) * (*t)(p) * K(p)
                - pow(B(p), 2) * (K(p) + 2.0 * (*t)(p) * dK(p))) / (4.0 * pow((*t)(p) * K(p), 2) * K(p));
        }
        N[i] = O[i] = pvp;
        vint[i] = make_pair(1, 1);
        break;
      }
    }
  }
};
}

// Effective density matrix for kinetic gradient. This function also returns zpq and ypq
tuple<shared_ptr<const Matrix>,shared_ptr<const Matrix>,shared_ptr<const Matrix>>
  DKHgrad::compute_tden_(shared_ptr<const Matrix> rdm1, shared_ptr<const Matrix> pmat, shared_ptr<const Matrix> wmat, shared_ptr<const Matrix> wmat_rev,
                         shared_ptr<const DiagMatrix> t, shared_ptr<const Matrix> v, shared_ptr<const Matrix> pvp) {
  const double c2 = c__ * c__;
  // Number of uncontracted basis functions
  const int nbasis = t->ndim();

  const TempArrays ta(t, v, pvp);

  const DiagMatrix& A = ta.A;
  const DiagMatrix& B = ta.B;
  const DiagMatrix& E = ta.E;
  const DiagMatrix& dA = ta.dA;
  const DiagMatrix& dB = ta.dB;
  const DiagMatrix& dE = ta.dE;

  const Matrix mmat = (*pmat % *wmat_rev) % *rdm1 * (*pmat % *wmat_rev);
  shared_ptr<Matrix> den, ypq;
  {
    const Matrix man = mmat * A * *v;
    const Matrix nam = *v * A * mmat;
    const Matrix mbs = mmat * B * *pvp;
    const Matrix sbm = *pvp * B * mmat;
    // DKH1 correction
    ypq = make_shared<Matrix>(2.0 * (mmat * E - c2 * mmat + (man + nam) * A + (mbs + sbm) * B));

    DiagMatrix em(nbasis);
    for (int p = 0; p != nbasis; ++p) {
      // dE_DKH1/dU_pq
      (*ypq)(p, p) += 2.0 * (mmat(p, p) * dE(p) + ((man(p, p) + nam(p, p)) * dA(p) + (mbs(p, p) + sbm(p, p)) * dB(p))) * (*t)(p);
      em(p) = dE(p) * mmat(p, p) + 2.0 * (nam(p, p) * dA(p) + sbm(p, p) * dB(p));
    }
    den = make_shared<Matrix>(*wmat * em ^ *wmat);
  }

  auto divide_Epq = [&E](const Matrix& input) {
    shared_ptr<Matrix> out = input.clone();
    for (int p = 0; p != input.mdim(); ++p)
      for (int q = 0; q != input.ndim(); ++q)
        (*out)(q, p) = input(q, p) / (E(p) + E(q));
    return *out;
  };

  auto mult_offdiag = [](const Matrix& m, const DiagMatrix& d, const bool rev) {
    shared_ptr<Matrix> out = m.clone(); 
    for (int p = 0; p != m.mdim(); ++p)
      for (int q = 0; q != m.ndim(); ++q)
        if (p != q)
          (*out)(q, p) = m(q, p) * d(!rev ? p : q);
    return *out;
  };

  for (int i = 0; i != 12; ++i) {
    const DiagMatrix& F = ta.F.at(i); 
    const DiagMatrix& G = ta.G.at(i); 
    const DiagMatrix& H = ta.H.at(i); 
    const DiagMatrix& dF = ta.dF.at(i); 
    const DiagMatrix& dG = ta.dG.at(i); 
    const DiagMatrix& dH = ta.dH.at(i); 

    const Matrix& N = *ta.N.at(i);
    const Matrix& O = *ta.O.at(i);
    const Matrix wn = divide_Epq(N);
    const Matrix wo = divide_Epq(O);
    const Matrix wwn = divide_Epq(wn);
    const Matrix wwo = divide_Epq(wo);

    const Matrix mh = mult_offdiag(mmat, H, false);
    const Matrix mh2 = mult_offdiag(mmat, dH, false);

    const Matrix mfn = mmat * F * wn;
    const Matrix mgo = mmat * G * wo;
    {
      Matrix whnfm = divide_Epq(H * wn * F * mmat * G);
      Matrix whogm = divide_Epq(H * wo * G * mmat * F);
      for (int p = 0; p != nbasis; ++p) {
        for (int q = 0; q != nbasis; ++q) {
          whnfm(q, p) += mh(q, p) * G(q) * F(p) * wn(p, p) / (E(p) + E(q));
          whogm(q, p) += mh(q, p) * F(q) * G(p) * wo(p, p) / (E(p) + E(q));
        }
        whnfm(p, p) += H(p) * mfn(p, p) * G(p) / (2.0 * E(p));
        whogm(p, p) += H(p) * mgo(p, p) * F(p) / (2.0 * E(p));
      }
      *ypq += O * whnfm + N * whogm + mfn * H * wo * G + mgo * H * wn * F;
    }

    {
      const Matrix fn = mmat * mult_offdiag(wn, F, true);
      const Matrix fo = mmat * mult_offdiag(wo, F, true);
      const Matrix fnn = mmat * mult_offdiag(wwn, F, true);
      const Matrix foo = mmat * mult_offdiag(wwo, F, true);

      Matrix gmfn2h(nbasis, nbasis);
      Matrix gmfo2h(nbasis, nbasis);
      Matrix gmfn2h2(nbasis, nbasis);
      Matrix gmfo2h2(nbasis, nbasis);
      for (int p = 0; p != nbasis; ++p)
        for (int q = 0; q != nbasis; ++q)
          if (p != q) {
            const double fac1 = G(q) * H(p) / (E(q) + E(p));
            const double fac2 = G(q) * dH(p) * (*t)(p);
            const double fac3 = G(q) * H(p) * dE(p) * (*t)(p);
            gmfn2h(q, p) = fac1 * fn(q, p);
            gmfo2h(q, p) = fac1 * fo(q, p);
            gmfn2h2(q, p) = fac2 * fn(q, p) - 2.0 * fac3 * fnn(q, p);
            gmfo2h2(q, p) = fac2 * fo(q, p) - 2.0 * fac3 * foo(q, p);
          }

      // dE_DKH2/dU_pq
      *ypq += O * gmfn2h + N * gmfo2h;

      const Matrix ogmfn2h2 = wo * gmfn2h2 + wn * gmfo2h2;
      for (int p = 0; p != nbasis; ++p)
        (*ypq)(p, p) += ogmfn2h2(p, p);
    }

    {
      const Matrix mfnhog = (mfn * H * (wo * dG - wwo * G * dE) + mgo * H * (wn * dF - wwn * F * dE)) * *t;
      const Matrix mfngh = (mfn * G * dH - (mmat * F * wwn * G * H + wwn * F * mh * G) * dE + wn * F * mh2 * G) * *t;
      const Matrix mfnghe = (mfn * G * H + wn * F * mh * G) * dE * *t;
      const Matrix cgofh = (mgo * F * dH - (mmat * G * wwo * F * H + wwo * G * mh * F) * dE + wo * G * mh2 * F) * *t;
      const Matrix cgofhe = (mgo * F * H + wo * G * mh * F) * dE * *t;
      for (int p = 0; p != nbasis; ++p)
        (*ypq)(p, p) += 2.0 * mfnhog(p, p) + mfngh(p, p) * wo(p, p) - mfnghe(p, p) * wwo(p, p)
                      + cgofh(p, p) * wn(p, p) - cgofhe(p, p) * wwn(p, p);
    }

    DiagMatrix dkh2(nbasis);
    {
      const Matrix ogmfn = wo * G * mmat * F * wn;
      const Matrix ogmfnn = wo * G * mmat * F * wwn;
      const Matrix oogmfn = wwo * G * mmat * F * wn;
      for (int p = 0; p != nbasis; ++p)
        dkh2(p) = ogmfn(p, p) * dH(p) - ogmfnn(p, p) * H(p) * dE(p) - oogmfn(p, p) * H(p) * dE(p);

      const Matrix ohnf = *mmat.hadamard_product(wo * H * (wn * dF - wwn * F * dE));
      const Matrix nhog = *mmat.hadamard_product(wn * H * (wo * dG - wwo * G * dE));
      dkh2 += DiagMatrix(ohnf * G.diag() + nhog * F.diag());
    }

    *den += *wmat * dkh2 ^ *wmat;
  }

  // z_pq from Lagrangian
  auto zpq = make_shared<Matrix>(nbasis, nbasis);
  for (int p = 0; p != nbasis; ++p)
    for (int q = 0; q != nbasis; ++q)
      (*zpq)(q, p) = fabs((*t)(p) - (*t)(q)) > 1.0e-12 ? -0.5 * ((*ypq)(q, p) - (*ypq)(p, q)) / ((*t)(q) - (*t)(p)) : 0.0;

  *den += *wmat * *zpq ^ *wmat;
  return make_tuple(den, zpq, ypq);
}


// Effective density matrices for NAI/SmallNAI gradients
tuple<shared_ptr<const Matrix>, shared_ptr<const Matrix>>
  DKHgrad::compute_vden_(shared_ptr<const Matrix> rdm1, shared_ptr<const Matrix> pmat, shared_ptr<const Matrix> wmat, shared_ptr<const Matrix> wmat_rev,
                         shared_ptr<const DiagMatrix> t, shared_ptr<const Matrix> v, shared_ptr<const Matrix> pvp) {
  const TempArrays ta(t, v, pvp);
  const DiagMatrix& A = ta.A;
  const DiagMatrix& B = ta.B;
  const DiagMatrix& E = ta.E;

  const Matrix mmat = (*pmat % *wmat_rev) % *rdm1 * (*pmat % *wmat_rev);
  // DKH1 correction
  auto den = make_shared<Matrix>(*wmat * A * mmat * A ^ *wmat);
  auto pvpden = make_shared<Matrix>(*wmat * B * mmat * B ^ *wmat);

  auto divide_Epq = [&E](const Matrix& input) {
    shared_ptr<Matrix> out = input.clone();
    for (int p = 0; p != input.mdim(); ++p)
      for (int q = 0; q != input.ndim(); ++q)
        (*out)(q, p) = input(q, p) / (E(p) + E(q));
    return *out;
  };

  for (int i = 0; i != 12; ++i) {
    const DiagMatrix& F = ta.F.at(i); 
    const DiagMatrix& G = ta.G.at(i); 
    const DiagMatrix& H = ta.H.at(i); 

    const Matrix wn = divide_Epq(*ta.N.at(i));
    const Matrix wo = divide_Epq(*ta.O.at(i));

    const Matrix whogmf = divide_Epq(H * wo * G * mmat * F);
    const Matrix whnfmg = divide_Epq(H * wn * F * mmat * G);

    const pair<int,int>& vint = ta.vint.at(i);
    (vint.first ? *pvpden : *den) += *wmat * whogmf ^ *wmat;
    (vint.second ? *pvpden : *den) += *wmat * whnfmg ^ *wmat;
  }
  return make_tuple(den, pvpden);
}


// Effective density matrix for overlap gradient
shared_ptr<const Matrix> DKHgrad::compute_sden_(shared_ptr<const Matrix> rdm1, shared_ptr<const Matrix> erdm1,
                                                shared_ptr<const Matrix> pmat, shared_ptr<const Matrix> wmat,
                                                shared_ptr<const Matrix> wmat_rev, shared_ptr<const DiagMatrix> t,
                                                shared_ptr<const Matrix> v, shared_ptr<const Matrix> pvp,
                                                shared_ptr<const Matrix> zpq, shared_ptr<const Matrix> ypq) {
  // Extra contributions stem for dependence of energy on overlap matrix (due to reverse transformation)
  const TempArrays ta(t, v, pvp);
  const DiagMatrix& E = ta.E; 
  const DiagMatrix& A = ta.A; 
  const DiagMatrix& B = ta.B; 

  const Matrix mup = *pmat * *rdm1 ^ (*wmat_rev % *pmat);
  // DKH1 correction
  const double c2 = c__ * c__;
  auto den = make_shared<Matrix>(2.0 * ((mup * E - c2 * mup) + mup * (A * *v * A + B * *pvp * B)) ^ *wmat);

  auto divide_Epq = [&E](const Matrix& input) {
    shared_ptr<Matrix> out = input.clone();
    for (int p = 0; p != input.mdim(); ++p)
      for (int q = 0; q != input.ndim(); ++q)
        (*out)(q, p) = input(q, p) / (E(p) + E(q));
    return *out;
  };

  // DKH2 corrections, one for each term
  for (int i = 0; i != 12; ++i) {
    const DiagMatrix& F = ta.F.at(i); 
    const DiagMatrix& G = ta.G.at(i); 
    const DiagMatrix& H = ta.H.at(i); 
    const Matrix wn = divide_Epq(*ta.N.at(i));
    const Matrix wo = divide_Epq(*ta.O.at(i));
    *den += mup * (F * wn * H * wo * G + G * wo * H * wn * F) ^ *wmat;
  }

  // a_tilde from dL/dU_pq
  // X_bar from Lagrangian, satisfies dL/dU_pq = 0
  Matrix xb(0.5 * *ypq + *zpq * *t);
  xb.symmetrize();

  *den -= (*pmat * *erdm1 ^ *pmat) + (*wmat * xb ^ *wmat);
  return den;
}
