//
// BAGEL - Parallel electron correlation program.
// Filename: f12int4.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


// carbon copy of what I wrote in the orz package
// meant to be standalone

#include <src/scf/scf.h>
#include <src/integral/rys/eribatch.h>
#include <src/integral/rys/slaterbatch.h>
#include <src/wfn/geometry.h>
#include <src/mp2/f12int4.h>

using namespace std;
using namespace bagel;

extern "C" { void start_up_slater_(); };

#ifdef HAVE_LIBSLATER
void F12Ref::compute() {
    // TODO this 2 is bad
    const size_t ncore = ncore_;
    const size_t nval = geom_->nele()/2 - ncore;
    const size_t nocc = nval + ncore;
    const size_t nmobasis = ref_->coeff()->mdim();
    const size_t naobasis = geom_->nbasis();
    const size_t nvirt = nmobasis - nocc;
    const double* const coeff = ref_->coeff()->data();
    const double* const ocoeff = coeff + ncore*naobasis;
    const double* const vcoeff = ocoeff + nval*naobasis;

    start_up_slater_();

    shared_ptr<Matrix> caom, cxxm, callm;
    int ncabs;
    tie(caom, cxxm, callm, ncabs) = generate_cabs();
    double *cao, *cxx, *call;

    // first compute half-transformed integrals
    shared_ptr<File2> eri, slater, yukawa, slater2;
    shared_ptr<File2> eric, slaterc, slater2c, slatercc;

    // integrals with 2gamma.
    {
    cout << "  * AO integrals..." << endl;
    AOInt<ERIBatch>    erib    (geom_);
    AOInt<SlaterBatch> slaterb (geom_, gamma_, true);
    AOInt<SlaterBatch> slater2b(geom_, gamma_*2.0);
    eri =       erib.data()->half_transform(ocoeff, nval);
    slater = slaterb.data()->half_transform(ocoeff, nval);
    yukawa = slaterb.data2()->half_transform(ocoeff, nval);
    slater2 = slater2b.data()->half_transform(ocoeff, nval);
    } {
    cout << "  * AO integrals with 1 aux" << endl;
    AOInt<ERIBatch>    erib    (geom_, 0.0,        false, false, true);
    AOInt<SlaterBatch> slaterb (geom_, gamma_,     false, false, true);
    AOInt<SlaterBatch> slater2b(geom_, gamma_*2.0, false, false, true);
    eric =       erib.data()->half_transform(ocoeff, nval);
    slaterc = slaterb.data()->half_transform(ocoeff, nval);
    slater2c = slater2b.data()->half_transform(ocoeff, nval);
    } {
    cout << "  * AO integrals with 2 aux" << endl;
    AOInt<SlaterBatch> slaterb (geom_, gamma_,     false, true, true);
    slatercc = slaterb.data()->half_transform(ocoeff, nval);
    }

    // X intermediate
    shared_ptr<F12Mat> x2 = slater2->f12mat(ocoeff);
    shared_ptr<F12Mat> x(new F12Mat(*x2));
#if 0
    shared_ptr<F12Ten> s = slater->f12ten(coeff, coeff, nbasis, nbasis);
    shared_ptr<F12Ten> sox =   slater->f12ten(coeff, cao, nocc, ncabs);
                      *sox += *slaterc->f12ten(coeff, cxx, nocc, ncabs);
    {
      *x -= *s->contract(s) + *sox->contract(sox);
    }
x->print();

    // V intermediate
    shared_ptr<F12Mat> v = yukawa->f12mat(ocoeff);
    {
      shared_ptr<F12Ten> r = eri->f12ten(coeff, coeff, nbasis, nbasis);
      *v -= *r->contract(s);
    }
v->print();
#endif

    // B intermediate
    //  - Y term
    shared_ptr<F12Mat> y(new F12Mat(*x2 * gamma_*gamma_));
#if 0
    DTensor xterm = gtrans_int2(ctinp, orbSymIRs, aote2_sl1, C_aoOBS, "_a__", C_aoOBS, "_a__", C_aoOBS, "_a__", C_aoOBS, "_a__");
    {
      DTensor frs = gtrans_int2(ctinp, orbSymIRs, aotei_sl1, C_aoOBS, "fav_", C_aoOBS, "_a__", C_aoOBS, "fav_", C_aoOBS, "_a__");
      DTensor xrs(ngen,ngen,ngen,ngen);
      contraction(xrs, frs, frs, nmoc, nmoc, ngen);
      xterm -= xrs * 0.5; // 0.5 because of the definition of contraction
    }
    {
      DTensor fxm =  gtrans_int2(ctinp, orbSymIRs, aotei_sl1, C_aoOBS, "___r", C_aoOBS, "_a__", C_aoOBS, "fa__", C_aoOBS, "_a__")
                        + gtrans_int2(ctinp, orbSymIRs, aotei_sl2, C_aoABS, "___r", C_aoOBS, "_a__", C_aoOBS, "fa__", C_aoOBS, "_a__");
      DTensor xxm(ngen,ngen,ngen,ngen);
      contraction(xxm, fxm, fxm, ncabs, ngen+nfrozen, ngen);
      xterm -= xxm;
    }
    // symmetrize
    {
      xterm += xterm.swapdim(2,3,0,1);
      xterm += xterm.swapdim(1,0,3,2);
      xterm *= 0.25;
    }
    print_local(xterm);
    ddot_to_amp(xterm, ngen, slater_param, "X term total (just for debugging)");
#endif

#if 0
    const int ngen = ctinp.nocc()+ctinp.nclosed();
    const int nmoc = nmo + nfrozen;

    // B intermediate
    DTensor bterm;
    {
      // T = gamma^2 X
      const double sl2 = slater_param*slater_param;
      DTensor tterm = gtrans_int2(ctinp, orbSymIRs, aote2_sl1, C_aoOBS, "_a__", C_aoOBS, "_a__", C_aoOBS, "_a__", C_aoOBS, "_a__") * sl2;
      ddot_to_amp(tterm, ngen, slater_param, "T term");
      bterm = tterm;
    }

    {
      // Q = X*h
      DTensor tmp1 = fock_aoOBS_aoOBS + exch_aoOBS_aoOBS; // hartree operator
      DTensor tmp2 = fock_aoABS_aoOBS + exch_aoABS_aoOBS;
      DTensor tmp3 = gtrans_int1(ctinp, orbSymIRs, tmp1, C_aoOBS, "favr", C_aoOBS, "fa__")
                        + gtrans_int1(ctinp, orbSymIRs, tmp2, C_aoABS, "favr", C_aoOBS, "fa__");
      DTensor hartree_weighted = tmp3 % C_ao;
      DTensor tmp4 = hartree_weighted(Slice(), Slice(0,nao)).copy();
      DTensor tmp5 = hartree_weighted(Slice(), Slice(nao,nao+nao_RI_F12)).copy();
      DTensor qterm = (gtrans_int2(ctinp, orbSymIRs, aote2_sl1, tmp4, "_a__", C_aoOBS, "_a__", C_aoOBS, "_a__", C_aoOBS, "_a__")
                          + gtrans_int2(ctinp, orbSymIRs, aote2_sl2, tmp5, "_a__", C_aoOBS, "_a__", C_aoOBS, "_a__", C_aoOBS, "_a__")) * 2.0;
      ddot_to_amp(qterm, ngen, slater_param, "Q term");
      bterm += qterm;
    }
    {
      // P1 = R^PQ_ij K^R_P R^kl_RQ (approximated by HY2)
      DTensor exch_obs_obs = gtrans_int1(ctinp, orbSymIRs, exch_aoOBS_aoOBS, C_aoOBS, "fav_", C_aoOBS, "fav_");
      DTensor exch_abs_obs = gtrans_int1(ctinp, orbSymIRs, exch_aoOBS_aoOBS, C_aoOBS, "___r", C_aoOBS, "fav_")
                           + gtrans_int1(ctinp, orbSymIRs, exch_aoABS_aoOBS, C_aoABS, "___r", C_aoOBS, "fav_");
      DTensor exch_abs_abs = gtrans_int1(ctinp, orbSymIRs, exch_aoOBS_aoOBS, C_aoOBS, "___r", C_aoOBS, "___r")
                           + gtrans_int1(ctinp, orbSymIRs, exch_aoABS_aoABS, C_aoABS, "___r", C_aoABS, "___r")
                           + gtrans_int1(ctinp, orbSymIRs, exch_aoABS_aoOBS, C_aoABS, "___r", C_aoOBS, "___r")
                           + gtrans_int1(ctinp, orbSymIRs, exch_aoOBS_aoABS, C_aoOBS, "___r", C_aoABS, "___r"); // cheap anyway
      DTensor exch_all = combine_matrices(exch_obs_obs, exch_abs_obs, exch_abs_abs, nmoc, ncabs);
      DTensor exch_weighted = exch_all % C_ao;
      DTensor tmp1 = exch_weighted(Slice(),Slice(0,nao)).copy();
      DTensor tmp2 = exch_weighted(Slice(),Slice(nao,nao+nao_RI_F12)).copy();

      // R^PQ_ij K^R_P R^kl_RQ without any approximation
      DTensor left = gtrans_int2(ctinp, orbSymIRs, aotei_sl1, tmp1   , "favr", C_aoOBS, "_a__", C_aoOBS, "favr", C_aoOBS, "_a__")
                   + gtrans_int2(ctinp, orbSymIRs, aotei_sl2, tmp2   , "favr", C_aoOBS, "_a__", C_aoOBS, "favr", C_aoOBS, "_a__")
                   + gtrans_int2(ctinp, orbSymIRs, aotei_sl2, C_aoABS, "favr", C_aoOBS, "_a__", tmp1   , "favr", C_aoOBS, "_a__").swapdim(2,3,0,1)
                   + gtrans_int2(ctinp, orbSymIRs, aotei_sl3, tmp2   , "favr", C_aoOBS, "_a__", C_aoABS, "favr", C_aoOBS, "_a__");
      DTensor tmp3 = gtrans_int2(ctinp, orbSymIRs, aotei_sl2, C_aoABS, "favr", C_aoOBS, "_a__", C_aoOBS, "favr", C_aoOBS, "_a__");
      DTensor right= gtrans_int2(ctinp, orbSymIRs, aotei_sl1, C_aoOBS, "favr", C_aoOBS, "_a__", C_aoOBS, "favr", C_aoOBS, "_a__")
                   + gtrans_int2(ctinp, orbSymIRs, aotei_sl3, C_aoABS, "favr", C_aoOBS, "_a__", C_aoABS, "favr", C_aoOBS, "_a__")
                   + tmp3 + tmp3.swapdim(2,3,0,1);
      DTensor p1(ngen,ngen,ngen,ngen);
      contraction(p1, left, right, nmoc+ncabs, nmoc+ncabs, ngen);
      ddot_to_amp(p1, ngen, slater_param, "P1 term");
      bterm -= p1;
    }
    DTensor fock_weighted;
    {
      // P2 = R^Pm_ij f^Q_P R^kl_Qm
      DTensor fock_obs_obs = gtrans_int1(ctinp, orbSymIRs, fock_aoOBS_aoOBS, C_aoOBS, "fav_", C_aoOBS, "fav_");
      DTensor fock_abs_obs = gtrans_int1(ctinp, orbSymIRs, fock_aoOBS_aoOBS, C_aoOBS, "___r", C_aoOBS, "fav_")
                           + gtrans_int1(ctinp, orbSymIRs, fock_aoABS_aoOBS, C_aoABS, "___r", C_aoOBS, "fav_");
      DTensor fock_abs_abs = gtrans_int1(ctinp, orbSymIRs, fock_aoOBS_aoOBS, C_aoOBS, "___r", C_aoOBS, "___r")
                           + gtrans_int1(ctinp, orbSymIRs, fock_aoABS_aoABS, C_aoABS, "___r", C_aoABS, "___r")
                           + gtrans_int1(ctinp, orbSymIRs, fock_aoABS_aoOBS, C_aoABS, "___r", C_aoOBS, "___r")
                           + gtrans_int1(ctinp, orbSymIRs, fock_aoOBS_aoABS, C_aoOBS, "___r", C_aoABS, "___r"); // cheap anyway
      DTensor fock_all = combine_matrices(fock_obs_obs, fock_abs_obs, fock_abs_abs, nmoc, ncabs);
      fock_weighted = fock_all % C_ao;
      DTensor tmp1 = fock_weighted(Slice(),Slice(0,nao)).copy();
      DTensor tmp2 = fock_weighted(Slice(),Slice(nao,nao+nao_RI_F12)).copy();
      DTensor left = gtrans_int2(ctinp, orbSymIRs, aotei_sl1, tmp1   , "favr", C_aoOBS, "_a__", C_aoOBS, "fa__", C_aoOBS, "_a__")
                   + gtrans_int2(ctinp, orbSymIRs, aotei_sl2, tmp2   , "favr", C_aoOBS, "_a__", C_aoOBS, "fa__", C_aoOBS, "_a__");
      DTensor right= gtrans_int2(ctinp, orbSymIRs, aotei_sl1, C_aoOBS, "favr", C_aoOBS, "_a__", C_aoOBS, "fa__", C_aoOBS, "_a__")
                   + gtrans_int2(ctinp, orbSymIRs, aotei_sl2, C_aoABS, "favr", C_aoOBS, "_a__", C_aoOBS, "fa__", C_aoOBS, "_a__");
      DTensor p2(ngen,ngen,ngen,ngen);
      contraction(p2, left, right, nmoc+ncabs, ngen+ctinp.nfrozen(), ngen);
      ddot_to_amp(p2, ngen, slater_param, "P2 term");
      bterm -= p2;
    }
    {
      // P3  = R^mA_ij f^n_m R^kl_nA
      DTensor fock_occ_occ = gtrans_int1(ctinp, orbSymIRs, fock_aoOBS_aoOBS, C_aoOBS, "fa__", C_aoOBS, "fa__");
      DTensor tmp0 = C_aoOBS(Slice(0,ngen+ctinp.nfrozen()),Slice(0,nao)).copy();
      DTensor tmp1 = fock_occ_occ % tmp0;
      DTensor left = gtrans_int2(ctinp, orbSymIRs, aotei_sl1, C_aoOBS, "___r", C_aoOBS, "_a__", tmp1   , "fa__", C_aoOBS, "_a__")
                   + gtrans_int2(ctinp, orbSymIRs, aotei_sl2, C_aoABS, "___r", C_aoOBS, "_a__", tmp1   , "fa__", C_aoOBS, "_a__");
      DTensor right= gtrans_int2(ctinp, orbSymIRs, aotei_sl1, C_aoOBS, "___r", C_aoOBS, "_a__", C_aoOBS, "fa__", C_aoOBS, "_a__")
                   + gtrans_int2(ctinp, orbSymIRs, aotei_sl2, C_aoABS, "___r", C_aoOBS, "_a__", C_aoOBS, "fa__", C_aoOBS, "_a__");
      // P5A = R^mA_ij f^P_m R^kl_PA
      DTensor tmp2 = fock_weighted(Slice(0,ngen+nfrozen),Slice(0,nao)).copy();
      DTensor tmp3 = fock_weighted(Slice(0,ngen+nfrozen),Slice(nao,nao+nao_RI_F12)).copy();
      DTensor tmp4 = gtrans_int2(ctinp, orbSymIRs, aotei_sl1, C_aoOBS, "___r", C_aoOBS, "_a__", tmp2   , "fa__", C_aoOBS, "_a__")
                   + gtrans_int2(ctinp, orbSymIRs, aotei_sl2, C_aoABS, "___r", C_aoOBS, "_a__", tmp2   , "fa__", C_aoOBS, "_a__")
                   + gtrans_int2(ctinp, orbSymIRs, aotei_sl2, tmp3   , "fa__", C_aoOBS, "_a__", C_aoOBS, "___r", C_aoOBS, "_a__").swapdim(2,3,0,1)
                   + gtrans_int2(ctinp, orbSymIRs, aotei_sl3, C_aoABS, "___r", C_aoOBS, "_a__", tmp3   , "fa__", C_aoOBS, "_a__");
      // -P3+2P5A
      left = tmp4*2.0 - left;
      DTensor p3(ngen,ngen,ngen,ngen);
      contraction(p3, left, right, ncabs, ngen+ctinp.nfrozen(), ngen);
      ddot_to_amp(p3, ngen, slater_param, "P3+P5A term");
      bterm -= p3;
    }
    {
      // P4  = R^rb_ij f^p_r R^kl_pb
      DTensor fock_obs_obs = gtrans_int1(ctinp, orbSymIRs, fock_aoOBS_aoOBS, C_aoOBS, "fav_", C_aoOBS, "fav_");
      DTensor tmp0 = C_aoOBS(Slice(0,nmoc),Slice(0,nao)).copy();
      DTensor tmp1 = fock_obs_obs % tmp0;
      DTensor left = gtrans_int2(ctinp, orbSymIRs, aotei_sl1, tmp1   , "fav_", C_aoOBS, "_a__", C_aoOBS, "__v_", C_aoOBS, "_a__");
      DTensor right= gtrans_int2(ctinp, orbSymIRs, aotei_sl1, C_aoOBS, "fav_", C_aoOBS, "_a__", C_aoOBS, "__v_", C_aoOBS, "_a__");
      // P5B = R^Ab_ij f^p_A R^kl_pb
      DTensor fock_abs_obs = gtrans_int1(ctinp, orbSymIRs, fock_aoOBS_aoOBS, C_aoOBS, "___r", C_aoOBS, "fav_")
                           + gtrans_int1(ctinp, orbSymIRs, fock_aoABS_aoOBS, C_aoABS, "___r", C_aoOBS, "fav_");
      DTensor tmp2 = fock_abs_obs % C_aoOBS(Slice(nmoc,nmoc+ncabs),Slice()).copy();
      DTensor tmp3 = fock_abs_obs % C_aoABS(Slice(nmoc,nmoc+ncabs),Slice()).copy();
      left += (gtrans_int2(ctinp, orbSymIRs, aotei_sl1, tmp2, "fav_", C_aoOBS, "_a__", C_aoOBS, "__v_", C_aoOBS, "_a__")
             + gtrans_int2(ctinp, orbSymIRs, aotei_sl2, tmp3, "fav_", C_aoOBS, "_a__", C_aoOBS, "__v_", C_aoOBS, "_a__")) * 2.0; // 2 from H.C.
      DTensor p5(ngen,ngen,ngen,ngen);
      contraction(p5, left, right, nmoc, ctinp.nvir(), ngen);
      ddot_to_amp(p5, ngen, slater_param, "P4+P5b term");
      bterm -= p5;
    }
    // symmetrize
    {
      bterm += bterm.swapdim(2,3,0,1);
      bterm += bterm.swapdim(1,0,3,2);
      bterm *= 0.25;
    }
    print_local(bterm);
    ddot_to_amp(bterm, ngen, slater_param, "B term total");

#endif

}

#define thresh 1.0e-8

tuple<shared_ptr<Matrix>, shared_ptr<Matrix>, shared_ptr<Matrix>, int> F12Ref::generate_cabs() const {

  // Form RI space which is a union of OBS and CABS.
  shared_ptr<Geometry> newgeom(new Geometry(*geom_));
  newgeom->merge_obs_aux();

  // TODO this does not work with linearly dependent orbital basis
  shared_ptr<Overlap> union_overlap(new Overlap(newgeom));
  shared_ptr<Matrix> ri_coeff = union_overlap->tildex(thresh);
  shared_ptr<Matrix> ri_reshaped = ref_->coeff()->cut(0,ri_coeff->ndim());

  // SVD to project out OBS component. Note singular values are all 1 as OBS is a subset of RI space.
  shared_ptr<Matrix> tmp(new Matrix(*ri_coeff % *union_overlap * *ri_reshaped));

  const int tmdim = tmp->mdim();
  const int tndim = tmp->ndim();

  shared_ptr<Matrix> U, V;
  tie(U, V) = tmp->svd();

  shared_ptr<Matrix> Ured = U->slice(tmdim, tndim); //(new Matrix(U, make_pair(tmdim, tndim)));
  shared_ptr<Coeff> coeff_cabs = shared_ptr<Coeff>(new Coeff(*ri_coeff * *Ured));

  shared_ptr<Matrix> coeff_entire = ri_reshaped->merge(coeff_cabs);

  // TODO note geom_->nbasis() here
  pair<shared_ptr<Matrix>, shared_ptr<Matrix>> t = coeff_cabs->split(geom_->nbasis(), geom_->naux());

  // TODO check
  int ncabs = ri_coeff->mdim();
  return make_tuple(t.first, t.second, coeff_entire, ncabs);
}

#endif
