//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: hess.cc
// Copyright (C) 2016 Toru Shiozaki
//
// Author: Bess Vlaisavljevich <bess.vlaisavljevich@northwestern.edu>
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

#include <string>
#include <src/grad/hess.h>
#include <src/grad/force.h>
#include <src/grad/finite.h>
#include <src/wfn/get_energy.h>
#include <src/grad/gradeval.h>
#include <src/util/atommap.h>
#include <src/util/constants.h>
#include <src/util/timer.h>
#include <src/prop/multipole.h>

using namespace std;
using namespace bagel;

static const AtomMap atommap;

Hess::Hess(shared_ptr<const PTree> idata, shared_ptr<const Geometry> g, shared_ptr<const Reference> r) : idata_(idata), geom_(g), ref_(r) {
  numhess_ = idata_->get<bool>("numhess", true);
  numforce_ = idata_->get<bool>("numforce", false);
  if (numhess_) {
    if (!numforce_)
      cout << "  The Hessian will be computed with central gradient differences (analytical gradients)" << endl;
    else
      throw logic_error("The code to compute the Hessian with central finite differences is not implemented"); // TODO
  } else {
    throw logic_error("Analytical Hessian has not been implemented");
  }

  auto input = idata_->get_child("method");
  auto m = input->begin();
  for ( ; m != --input->end(); ++m) {
    const string title = to_lower((*m)->get<string>("title", ""));
    if (title != "molecule") {
      tie(energy_, r) = get_energy(title, *m, geom_, r);
    } else {
      geom_ = make_shared<Geometry>(*geom_, *m);
      if (r) r = r->project_coeff(geom_);
    }
  }
  ref_ = r;

  dx_ = idata_->get<double>("dx", 1.0e-3);
  cout << "  Finite difference displacement (dx) is " << setprecision(8) << dx_ << " bohr" << endl;

  nproc_ = idata_->get<int>("nproc", 1);

  const int natom = geom_->natom();
  const int ndispl = natom * 3;
  hess_      = make_shared<Matrix>(ndispl, ndispl);
  mw_hess_   = make_shared<Matrix>(ndispl, ndispl);
  cartesian_ = make_shared<Matrix>(3, ndispl); //matrix of dmu/dR
}


void Hess::compute() {

  const int natom = geom_->natom();
  const int ndispl = natom * 3;

  muffle_ = make_shared<Muffle>("freq.log");

  // compute Hessian and dipole derivatives using finite difference
  compute_finite_diff_();

  // symmetrize mass weighted hessian
  hess_->print("Hessian");
  mw_hess_->print("Mass Weighted Hessian", ndispl);
  mw_hess_->symmetrize();

  // check if all of the mass are equal to the averaged mass
  bool averaged = true;
  for (auto& i : geom_->atoms())
    averaged &= fabs(i->mass() - atommap.averaged_mass(i->name())) < 1.0e-8;
  if (averaged)
    cout << "    (masses averaged over the natural occurance of isotopes)" << endl << endl;
  else
    cout << "    (custom masses were specified in the input)" << endl << endl;

  mw_hess_->print("Symmetrized Mass Weighted Hessian", ndispl);

  // compute projected Hessian
  project_zero_freq_();

  // diagonalize hessian; eig(i) in Hartree/bohr^2*amu
  VectorB eig(ndispl);
  proj_hess_->diagonalize(eig);

  cout << endl << " Mass Weighted Hessian Eigenvalues" << endl; // units Hartree/bohr^2*amu
  for (int i = 0; i != ndispl; ++i )
    cout << setw(10) << setprecision(5) << eig(i);
  cout << endl;

  proj_hess_->print("Mass Weighted Hessian Eigenvectors", ndispl);

  // convert mw eigenvectors to normalized cartesian modes
  eigvec_cart_ = make_shared<Matrix>(ndispl,ndispl);

  for (int i = 0, counter = 0; i != natom; ++i)
    for (int j = 0; j != 3; ++j, ++counter)
      for (int k = 0, step = 0; k != natom; ++k)
        for (int l = 0; l != 3; ++l, ++step)
          eigvec_cart_->element(step, counter) =  proj_hess_->element(step,counter) / sqrt(geom_->atoms(k)->mass());

  // calculate IR intensities:
  auto normal = make_shared<Matrix>(*cartesian_ * *eigvec_cart_); // dipole derivatives for the normal modes dmu/dQ; units (e bohr / bohr) * 1/sqrt(amu)

  VectorB dmudq2(ndispl); //square of the dipole derivative in hartree*bohr/amu
  for (int i = 0; i != ndispl; ++i)
    dmudq2(i) = blas::dot_product(normal->element_ptr(0, i), 3, normal->element_ptr(0, i));

  ir_ = vector<double>(ndispl, 0.0);
  freq_ = vector<double>(ndispl, 0.0);

  //frequences and IR intensity ( N*pi/3*c^2 * (dmu/dQ)^2 )
  for (int i = 0; i != ndispl; ++i) {
    freq_[i] = (fabs(eig(i)) > 1.0e-6 ? (eig(i) > 0.0 ? sqrt((eig(i) * au2joule__) / amu2kilogram__) / (100.0 * au2meter__ * 2.0 * pi__ * csi__) :
      -sqrt((-eig(i)     * au2joule__) / amu2kilogram__) / (100.0 * au2meter__ * 2.0 * pi__ * csi__)) : 0) ;
    ir_[i] = (fabs(eig(i)) > 1.0e-6 ? (eig(i) > 0.0 ? ((avogadro__ * pi__ *au2meter__ * au2joule__)/ (3.0 * 1000.0 * csi__ * csi__ * amu2kilogram__)) * dmudq2(i) :
      ((avogadro__ * pi__ *au2meter__ * au2joule__)/ (3.0 * 1000.0 * csi__ * csi__ * amu2kilogram__)) * dmudq2(i)) : 0.0);
  }

  print_ir_();
}


void Hess::compute_finite_diff_() {
  Timer timer;
  const int natom = geom_->natom();
  const int ncomm = mpi__->world_size() / nproc_;
  const int npass = (natom * 3 - 1) / ncomm + 1;

  for (int ipass = 0; ipass != npass; ++ipass) {
    const int ncolor = (ipass == (npass-1)) ? (natom * 3) % ncomm : ncomm;
    const int icomm = mpi__->world_rank() % ncolor;
    if (ncolor != 0) {
      mpi__->split(ncolor);
    } else {
      continue;
    }

    const int counter = icomm + ncomm * ipass;
    const int i = counter / 3;
    const int j = counter % 3;

    muffle_->mute();

    vector<double> dipole_plus;
    shared_ptr<const GradFile> outplus;
    //displace +dx
    {
      auto displ = make_shared<XYZFile>(natom);
      displ->element(j,i) = dx_;
      auto geom_plus = make_shared<Geometry>(*geom_, displ, make_shared<PTree>(), false, false);
      geom_plus->print_atoms();

      shared_ptr<const Reference> ref_plus;
      if (ref_)
        ref_plus = ref_->project_coeff(geom_plus);

      auto plus = make_shared<Force>(idata_, geom_plus, ref_plus);
      outplus = plus->compute();
      dipole_plus = plus->force_dipole();
    }

    // displace -dx
    vector<double> dipole_minus;
    shared_ptr<const GradFile> outminus;
    {
      auto displ = make_shared<XYZFile>(natom);
      displ->element(j,i) = -dx_;
      auto geom_minus = make_shared<Geometry>(*geom_, displ, make_shared<PTree>(), false, false);
      geom_minus->print_atoms();

      shared_ptr<const Reference> ref_minus;
      if (ref_)
        ref_minus = ref_->project_coeff(geom_minus);

      auto minus = make_shared<Force>(idata_, geom_minus, ref_minus);
      outminus = minus->compute();
      dipole_minus = minus->force_dipole();
    }

    if (mpi__->rank() == 0) {
      for (int k = 0, step = 0; k != natom; ++k) { // atom j
        for (int l = 0; l != 3; ++l, ++step) { //xyz
          (*hess_)(counter,step) = (outplus->element(l,k) - outminus->element(l,k)) / (2*dx_);
          (*mw_hess_)(counter,step) =  (*hess_)(counter,step) / sqrt(geom_->atoms(i)->mass() * geom_->atoms(k)->mass());
          (*cartesian_)(l,counter) = (dipole_plus[l] - dipole_minus[l]) / (2*dx_);
        }
      }
    }
    muffle_->unmute();
    stringstream ss; ss << "Hessian evaluation (" << setw(2) << i*3+j+1 << " / " << natom * 3 << ")";
    timer.tick_print(ss.str());

    mpi__->merge();
  }

  hess_->allreduce();
  mw_hess_->allreduce();
  cartesian_->allreduce();
}


void Hess::project_zero_freq_() {
  const int natom = geom_->natom();
  const int ndispl = natom * 3;

  // calculate center of mass
  VectorB cmass(3); // values needed to calc center of mass. mi*xi, mi*yi, mi*zi, and total mass
  double total_mass = 0.0;
  // compute center of mass
  for (auto& atom : geom_->atoms()) {
    for (int i = 0; i != 3; ++i)
      cmass(i) += atom->mass() * atom->position(i);
    total_mass += atom->mass();
  }
  blas::scale_n(1.0/total_mass, cmass.data(), 3);
  cout << "    * Projecting out translational and rotational degrees of freedom " << endl;

  Matrix proj(6, ndispl);
  for (int i = 0; i != natom; ++i) {
    const double imass = sqrt(geom_->atoms(i)->mass());
    const array<double,3> pos {{geom_->atoms(i)->position(0) - cmass(0),
                                geom_->atoms(i)->position(1) - cmass(1),
                                geom_->atoms(i)->position(2) - cmass(2)}};
    for (int j = 0; j != 3; ++ j)
      proj(j, 3*i+j) = imass;

    proj(3,3*i)   =  0.0;
    proj(3,3*i+1) = -imass * pos[2];
    proj(3,3*i+2) =  imass * pos[1];

    proj(4,3*i)   =  imass * pos[2];
    proj(4,3*i+1) =  0.0;
    proj(4,3*i+2) = -imass * pos[0];

    proj(5,3*i)   = -imass * pos[1];
    proj(5,3*i+1) =  imass * pos[0];
    proj(5,3*i+2) =  0.0;
  }

  //normalize the set of six orthogonal vectors
  vector<double> norm(6, 0.0);
  for (int i = 0; i != ndispl; ++i)
    for (int j = 0; j != 6; ++j)
      norm[j] += proj(j, i) * proj(j, i);

  for (int i = 0; i != ndispl; ++i)
    for (int j = 0; j != 6; ++j)
      if (fabs(norm[j]) > 1.0e-15)
        proj(j, i) /= sqrt(norm[j]);

  Matrix p(ndispl, ndispl);
  p.unit();
  p -= proj % proj;

  proj_hess_ = make_shared<Matrix>(p % *mw_hess_ * p);
}


void Hess::print_ir_() const {
  cout << "    * Vibrational frequencies, IR intensities, and corresponding cartesian eigenvectors" << endl << endl;
  const int len_n = eigvec_cart_->ndim();
  const int len_m = eigvec_cart_->mdim();

  for (int i = 0; i < len_m; i += 6) {
    const int kmax = min(6, len_m-i);
    cout << setw(17) << " ";
    for (int k = 0; k != kmax; ++k)
      cout << setw(20) << i + k;
    cout << endl << setw(17) << "Freq (cm-1)";
    for (int k = 0; k != kmax; ++k)
      cout << setw(20) << setprecision(2) << freq_[i+k] ;

    cout << endl << endl << setw(17) << "IR Int. (km/mol)";
    for (int k = 0; k != kmax; ++k)
      cout << setw(20) << setprecision(2) << ir_[i+k] ;

    cout << endl << setw(17) << "Rel. IR Int.";
    for (int k = 0; k != kmax; ++k)
      cout << setw(20) << setprecision(2) << (ir_[i+k]/(*max_element(ir_.begin(), ir_.end())))*100.0;
    cout << endl << endl;

    for (int j = 0; j != len_n; ++j) {
      cout << setw(17) << j;
      for (int k = 0; k != kmax; ++k)
        cout << setw(20) << setprecision(5) << eigvec_cart_->element(j, i+k);
      cout << endl;
    }
    cout << endl;
  }
}


shared_ptr<const Reference> Hess::conv_to_ref() const {
  shared_ptr<Reference> out;
  if (ref_) {
    out = make_shared<Reference>(*ref_);
  } else {
    out = make_shared<Reference>(geom_, nullptr, 0, 0, 0);
    cout << "  ** CAUTION ** Reference object being created by Hessian is only valid for printing!" << endl;
  }
  out->set_prop_freq(freq_);
  out->set_prop_ir(ir_);
  out->set_prop_eig(eigvec_cart_);
  return out;
}
