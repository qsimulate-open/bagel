//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: finite.cc
// Copyright (C) 2012 Toru Shiozaki
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

#include <src/scf/hf/rohf.h>
#include <src/scf/ks/ks.h>
#include <src/scf/sohf/soscf.h>
#include <src/scf/dhf/dirac.h>
#include <src/scf/giaohf/rhf_london.h>
#include <src/ci/fci/distfci.h>
#include <src/ci/fci/harrison.h>
#include <src/ci/fci/knowles.h>
#include <src/ci/ras/rasci.h>
#include <src/ci/zfci/zharrison.h>
#include <src/pt2/nevpt2/nevpt2.h>
#include <src/pt2/mp2/mp2.h>
#include <src/pt2/dmp2/dmp2.h>
#include <src/multi/casscf/cassecond.h>
#include <src/multi/casscf/casnoopt.h>
#include <src/multi/zcasscf/zcassecond.h>
#include <src/multi/zcasscf/zcasnoopt.h>
#include <src/smith/smith.h>
#include <src/periodic/pscf.h>
#include <src/grad/gradeval.h>
#include <src/grad/finite.h>
#include <src/util/timer.h>
#include <src/wfn/construct_method.h>

using namespace std;
using namespace bagel;

tuple<double, shared_ptr<const Reference>>
FiniteGrad::get_energy(string title, shared_ptr<const PTree> itree, shared_ptr<const Geometry> geom, shared_ptr<const Reference> ref) const {
  double out = 0.0;

  if (!geom || !geom->magnetism()) {
    if (title == "hf")           { auto m = make_shared<RHF>(itree, geom, ref);       m->compute();   out = m->energy();                        ref = m->conv_to_ref(); }
    else if (title == "ks")      { auto m = make_shared<KS>(itree, geom, ref);        m->compute();   out = m->energy();                        ref = m->conv_to_ref(); }
    else if (title == "uhf")     { auto m = make_shared<UHF>(itree, geom, ref);       m->compute();   out = m->energy();                        ref = m->conv_to_ref(); }
    else if (title == "rohf")    { auto m = make_shared<ROHF>(itree, geom, ref);      m->compute();   out = m->energy();                        ref = m->conv_to_ref(); }
    else if (title == "soscf")   { auto m = make_shared<SOSCF>(itree, geom, ref);     m->compute();   out = m->energy();                        ref = m->conv_to_ref(); }
    else if (title == "mp2")     { auto m = make_shared<MP2>(itree, geom, ref);       m->compute();   out = m->energy();                        ref = m->conv_to_ref(); }
    else if (title == "dhf")     { auto m = make_shared<Dirac>(itree, geom, ref);     m->compute();   out = m->energy();                        ref = m->conv_to_ref(); }
    else if (title == "dmp2")    { auto m = make_shared<DMP2>(itree, geom, ref);      m->compute();   out = m->energy();                        ref = m->conv_to_ref(); }
    else if (title == "smith")   { auto m = make_shared<Smith>(itree, geom, ref);     m->compute();   out = m->algo()->energy(target_state_);   ref = m->conv_to_ref(); }
    else if (title == "relsmith"){ auto m = make_shared<RelSmith>(itree, geom, ref);  m->compute();   out = m->algo()->energy(target_state_);   ref = m->conv_to_ref(); }
    else if (title == "zfci")    { auto m = make_shared<ZHarrison>(itree, geom, ref); m->compute();   out = m->energy(target_state_);           ref = m->conv_to_ref(); }
    else if (title == "ras") {
      const string algorithm = itree->get<string>("algorithm", "");
      if ( algorithm == "local" || algorithm == "" ) { auto m = make_shared<RASCI>(itree, geom, ref); m->compute(); out = m->energy(target_state_); ref = m->conv_to_ref(); }
#ifdef HAVE_MPI_H
      else if ( algorithm == "dist" || algorithm == "parallel" )
        throw logic_error("Parallel RASCI not implemented");
#endif
      else
        throw runtime_error("unknown RASCI algorithm specified. " + algorithm);
    }
    else if (title == "fci") {
      const string algorithm = itree->get<string>("algorithm", "");
      const bool dokh = (algorithm == "" || algorithm == "auto") && geom->nele() > geom->nbasis();
      if (dokh || algorithm == "kh" || algorithm == "knowles" || algorithm == "handy") {
        auto m = make_shared<KnowlesHandy>(itree, geom, ref);        m->compute();   out = m->energy(target_state_);   ref = m->conv_to_ref();
      } else if (algorithm == "hz" || algorithm == "harrison" || algorithm == "zarrabian" || algorithm == "") {
        auto m = make_shared<HarrisonZarrabian>(itree, geom, ref);   m->compute();   out = m->energy(target_state_);   ref = m->conv_to_ref();
#ifdef HAVE_MPI_H
      } else if (algorithm == "parallel" || algorithm == "dist") {
        auto m = make_shared<DistFCI>(itree, geom, ref);             m->compute();   out = m->energy(target_state_);   ref = m->conv_to_ref();
#endif
      } else
        throw runtime_error("unknown FCI algorithm specified. " + algorithm);
    }
    else if (title == "casscf") {
      string algorithm = itree->get<string>("algorithm", "");
      if (algorithm == "second" || algorithm == "") {
        auto m = make_shared<CASSecond>(itree, geom, ref);           m->compute();   out = m->energy(target_state_);   ref = m->conv_to_ref();
      } else if (algorithm == "noopt") {
        auto m = make_shared<CASNoopt>(itree, geom, ref);            m->compute();   out = m->energy(target_state_);   ref = m->conv_to_ref();
      } else
        throw runtime_error("unknown CASSCF algorithm specified: " + algorithm);
    }
    else if (title == "nevpt2")  { auto m = make_shared<NEVPT2<double>>(itree, geom, ref);          m->compute();   out = m->energy();    ref = m->conv_to_ref(); }
    else if (title == "dnevpt2") { auto m = make_shared<NEVPT2<complex<double>>>(itree, geom, ref); m->compute();   out = m->energy();    ref = m->conv_to_ref(); }
    else if (title == "zcasscf") {
      string algorithm = itree->get<string>("algorithm", "");
      if (algorithm == "second" || algorithm == "") {
        auto m = make_shared<ZCASSecond>(itree, geom, ref);          m->compute();   out = m->energy(target_state_);   ref = m->conv_to_ref();
      } else if (algorithm == "noopt") {
        auto m = make_shared<ZCASNoopt>(itree, geom, ref);           m->compute();   out = m->energy(target_state_);   ref = m->conv_to_ref();
      } else
        cout << " Optimization algorithm " << algorithm << " is not compatible with ZCASSCF " << endl;
    } else if (title == "pscf")    { auto m = make_shared<PSCF>(itree, geom, ref);             m->compute();   out = m->energy();    ref = m->conv_to_ref(); }

  // now the versions to use with magnetic fields
  } else {
    if (title == "hf")           { auto m = make_shared<RHF_London>(itree, geom, ref);       m->compute();   out = m->energy();                        ref = m->conv_to_ref(); }
    else if (title == "dhf")     { auto m = make_shared<Dirac>(itree, geom, ref);            m->compute();   out = m->energy();                        ref = m->conv_to_ref(); }
    else if (title == "zfci")    { auto m = make_shared<ZHarrison>(itree, geom, ref);        m->compute();   out = m->energy(target_state_);           ref = m->conv_to_ref(); }
    else if (title == "zcasscf") {
      string algorithm = itree->get<string>("algorithm", "");
      if (algorithm == "second" || algorithm == "") {
        auto m = make_shared<ZCASSecond>(itree, geom, ref);          m->compute();   out = m->energy(target_state_);   ref = m->conv_to_ref();
      } else if (algorithm == "noopt") {
        auto m = make_shared<ZCASNoopt>(itree, geom, ref);           m->compute();   out = m->energy(target_state_);   ref = m->conv_to_ref();
      } else
        cout << " Optimization algorithm " << algorithm << " is not compatible with ZCASSCF " << endl;
    } else
      throw runtime_error(to_upper(title) + " method has not been implemented with an applied magnetic field.");
  }

  return tie(out, ref);
}

shared_ptr<GradFile> FiniteGrad::compute() {
  for (auto m = idata_->begin(); m != idata_->end(); ++m) {
    const string title = to_lower((*m)->get<string>("title", ""));
    tie(energy_, ref_) = get_energy(title, *m, geom_, ref_);
  }

  cout << "  Gradient evaluation with respect to " << geom_->natom() * 3 << " DOFs" << endl;
  cout << "  Finite difference size (dx) is " << setprecision(8) << dx_ << " Bohr" << endl;

  Timer timer;
  muffle_ = make_shared<Muffle>("finite.log");

  const int natom = geom_->natom();
  auto grad = make_shared<GradFile>(natom);

  for (int i = 0; i != natom; ++i) {    // atom i
    for (int j = 0; j != 3; ++j) {      // xyz
      muffle_->mute();

      double energy_plus = 0.0;
      {
        auto displ = make_shared<XYZFile>(natom);
        displ->element(j,i) = dx_;
        auto geom_plus = make_shared<Geometry>(*geom_, displ, make_shared<PTree>(), false, false);
        geom_plus->print_atoms();

        shared_ptr<const Reference> ref_plus;
        if (ref_)
          ref_plus = ref_->project_coeff(geom_plus);

        for (auto m = idata_->begin(); m != idata_->end(); ++m) {
          const string title = to_lower((*m)->get<string>("title", ""));
          tie(energy_plus, ref_plus) = get_energy(title, *m, geom_plus, ref_plus);
        }
      }

      double energy_minus = 0.0;
      {
        auto displ = make_shared<XYZFile>(natom);
        displ->element(j,i) = -dx_;
        auto geom_minus = make_shared<Geometry>(*geom_, displ, make_shared<PTree>(), false, false);
        geom_minus->print_atoms();

        shared_ptr<const Reference> ref_minus;
        if (ref_)
          ref_minus = ref_->project_coeff(geom_minus);

        for (auto m = idata_->begin(); m != idata_->end(); ++m) {
          const string title = to_lower((*m)->get<string>("title", ""));
          tie(energy_minus, ref_minus) = get_energy(title, *m, geom_minus, ref_minus);
        }
      }
      grad->element(j,i) = (energy_plus - energy_minus) / (2.0 * dx_);
      muffle_->unmute();
      stringstream ss; ss << "Finite difference evaluation (" << setw(2) << i*3+j+1 << " / " << geom_->natom() * 3 << ")";
      timer.tick_print(ss.str());
    }
  }

  grad->print(": Calculated with finite difference", 0);
  return grad;
}

template<>
shared_ptr<GradFile> FiniteNacm<CASSCF>::compute() {
  cout << "  NACME evaluation with respect to " << geom_->natom() * 3 << " DOFs" << endl;
  cout << "  Finite difference size (dx) is " << setprecision(8) << dx_ << " Bohr(s)" << endl;

  const int natom = geom_->natom();

  shared_ptr<Dvec> civ_ref = ref_->civectors()->copy();
  civ_ref->print (/*sort=*/false);

  const int nclosed = ref_->nclosed();
  const int nocc = ref_->nocc();

  Timer timer;
  muffle_ = make_shared<Muffle>("finite.log");

  auto acoeff_ref = make_shared<Matrix>(ref_->coeff()->slice(nclosed, nocc));
  auto grad = make_shared<GradFile>(natom);
  const int norb = civ_ref->det()->norb();
  auto gmo = make_shared<Matrix>(norb, norb);
  gmo->zero();

  for (int i = 0; i != natom; ++i) {      // atom i
    for (int j = 0; j != 3; ++j) {        // xyz
      muffle_->mute();

      shared_ptr<Matrix> acoeff_plus;
      shared_ptr<Dvec> civ_plus;
      {
        auto displ = make_shared<XYZFile>(natom);
        displ->element(j,i) = dx_;
        auto geom_plus = make_shared<Geometry>(*geom_, displ, make_shared<PTree>(), false, false);
        geom_plus->print_atoms();

        shared_ptr<const Reference> ref_plus;
        if (ref_)
          ref_plus = ref_->project_coeff(geom_plus);

        auto energy_method = construct_method(method_, idata_, geom_plus, ref_plus);
        energy_method->compute();
        auto refgrad_plus = energy_method->conv_to_ref();
        acoeff_plus = make_shared<Matrix>(refgrad_plus->coeff()->slice(nclosed, nocc));
        for (int im = 0; im != acoeff_ref->mdim(); ++im) {
          double dmatch = blas::dot_product(acoeff_ref->element_ptr(0, im), acoeff_ref->ndim(), acoeff_plus->element_ptr(0,im));
          if (dmatch < 0.0)
            blas::scale_n(-1.0, acoeff_plus->element_ptr(0, im), acoeff_ref->ndim());
        }
        civ_plus = refgrad_plus->civectors()->copy();
        civ_plus->match(civ_ref);
      }

      shared_ptr<Matrix> acoeff_minus;
      shared_ptr<Dvec> civ_minus;
      {
        auto displ = make_shared<XYZFile>(natom);
        displ->element(j,i) = -dx_;
        auto geom_minus = make_shared<Geometry>(*geom_, displ, make_shared<PTree>(), false, false);
        geom_minus->print_atoms();

        shared_ptr<const Reference> ref_minus;
        if (ref_)
          ref_minus = ref_->project_coeff(geom_minus);

        auto energy_method = construct_method(method_, idata_, geom_minus, ref_minus);
        energy_method->compute();
        auto refgrad_minus = energy_method->conv_to_ref();
        acoeff_minus = make_shared<Matrix>(refgrad_minus->coeff()->slice(nclosed, nocc));
        for (int im = 0; im != acoeff_ref->mdim(); ++im) {
          double dmatch = blas::dot_product(acoeff_ref->element_ptr(0, im), acoeff_ref->ndim(), acoeff_minus->element_ptr(0,im));
          if (dmatch < 0.0)
            blas::scale_n(-1.0, acoeff_minus->element_ptr(0, im), acoeff_ref->ndim());
        }
        civ_minus = refgrad_minus->civectors()->copy();
        civ_minus->match(civ_ref);
      }

      auto civ_diff = civ_plus->copy();
      *civ_diff -= *civ_minus;
      civ_diff->scale(1.0 / (2.0 * dx_));
      grad->element(j,i) = civ_ref->data(target_state1_)->dot_product(civ_diff->data(target_state2_));

      auto acoeff_diff = make_shared<Matrix>(*acoeff_plus - *acoeff_minus);
      acoeff_diff->scale(1.0 / (2.0 * dx_));

      auto Smn = make_shared<Overlap>(geom_);
      auto Uij = make_shared<Matrix>(*acoeff_ref % *Smn * *acoeff_diff);
      for (int ii = 0; ii != norb; ++ii) {
        for (int ij = 0; ij != norb; ++ij) {
          const int lena = civ_ref->det()->lena();
          const int lenb = civ_ref->det()->lenb();
          if (ii != ij) {
            for (auto& iter : civ_ref->det()->phia(ii, ij)) {
              size_t iaA = iter.source;
              size_t iaB = iter.target;
              double sign = static_cast<double>(iter.sign);

              for (size_t ib = 0; ib != lenb; ++ib) {
                double factor = civ_ref->data(target_state1_)->data(ib+iaB*lenb) * civ_ref->data(target_state2_)->data(ib+iaA*lenb) * sign;
                grad->element(j,i) += factor * (Uij->element(ij, ii) - Uij->element(ii, ij)) * .5;
                if ((i + j * 3) == 0) {
                  gmo->element(ij, ii) += factor * .5;
                  gmo->element(ii, ij) -= factor * .5;
                }
              }
            }
            for (size_t ia = 0; ia != lena; ++ia) {
              for (auto& iter : civ_ref->det()->phib(ii, ij)) {
                size_t ibA = iter.source;
                size_t ibB = iter.target;
                double sign = static_cast<double>(iter.sign);
                double factor = civ_ref->data(target_state1_)->data(ibB+ia*lenb) * civ_ref->data(target_state2_)->data(ibA+ia*lenb) * sign;
                grad->element(j,i) += factor * (Uij->element(ij, ii) - Uij->element(ii, ij)) * .5;
                if ((i + j * 3) == 0) {
                  gmo->element(ij, ii) += factor * .5;
                  gmo->element(ii, ij) -= factor * .5;
                }
              }
            }
          }
        }
      }
      muffle_->unmute();
      stringstream ss; ss << "Finite difference evaluation (" << setw(2) << i*3+j+1 << " / " << geom_->natom() * 3 << ")";
      timer.tick_print(ss.str());
    }
  }
  auto gfin = make_shared<Matrix>(*acoeff_ref * *gmo ^ *acoeff_ref);
  auto grad_basis = make_shared<GradFile>(natom);
  grad_basis = contract_nacme(nullptr, nullptr, nullptr, nullptr, gfin, /*numerical=*/true);
  *grad += *grad_basis;

  grad->print(": NACME calculated with finite difference", 0);

  return grad;
}
