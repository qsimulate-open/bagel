//
// BAGEL - Parallel electron correlation program.
// Filename: pseudospin.cc
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Ryan D. Reynolds <rreynolds2018@u.northwestern.edu>
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

// A function to compute the coefficients of Extended Stevens Operators, for the pseudospin Hamiltonian
// Notation follows I. D. Ryabov, Appl. Magn. Reson. (2009) 35, 481-494.
// Some equations also come from I. D. Ryabov, J. Magn. Reson. (1999) 140, 141-145.

#include <src/prop/pseudospin/pseudospin.h>
#include <src/mat1e/rel/spinint.h>
#include <src/mat1e/angmom.h>
#include <src/integral/os/angmombatch.h>
#include <src/mat1e/rel/small1e.h>

using namespace std;
using namespace bagel;


namespace {
// simple local function to convert N electrons into S string - for printout
string spin_val(const int nspin) {
  string out = to_string(nspin / 2);
  if (nspin % 2 == 1)
    out += " 1/2";
  if (nspin == 1)
    out.erase(0, 2);
  return out;
}
}


Stevens_Operator::Stevens_Operator(shared_ptr<const ZMatrix> _mat, const int _ord, const int _ind)
 : nspin_(_mat->ndim() - 1), matrix_(_mat), order_(_ord), index_(_ind) {
  coeff_ = nan("");
  assert(matrix_->ndim() == matrix_->mdim());
  assert(abs(index_) <= order_);
  assert(order_ >= 0);
  assert((order_ <= nspin_ || matrix_->rms() < 1.0e-8)); // High-order contributions should be zero for low spin
}

string Stevens_Operator::operator_name() const {
  string out = "O_" + to_string(order_) + "^" + to_string(index_);
  if (index_ >= 0) out += " ";
  return out;
}

string Stevens_Operator::coeff_name() const {
  string out = operator_name();
  out[0] = 'B';
  return out;
}


shared_ptr<const Matrix> Pseudospin::read_axes(shared_ptr<const Matrix> default_axes) const {
  array<double,3> default_x, default_z;
  const array<int,3> fwd = {{ 1, 2, 0 }};
  const array<int,3> bck = {{ 2, 0, 1 }};

  for (int i = 0; i != 3; ++i) {
    default_x[i] = default_axes->element(i, 0);
    default_z[i] = default_axes->element(i, 2);
  }

  array<array<double, 3>, 3> new_axes;
  new_axes[2] = idata_->get_array<double,3>("zaxis", default_z);
  new_axes[0] = idata_->get_array<double,3>("xaxis", default_x);
  if (default_z != new_axes[2] && default_x == new_axes[0]) {
    // No x-axis was given to us, so generate one by taking cross product of default y and input z axes
    for (int i = 0; i != 3; ++i)
      new_axes[0][i] = ((default_axes->element(fwd[i], 1) * new_axes[2][bck[i]]) - (new_axes[2][fwd[i]] * default_axes->element(bck[i], 1)));
  } else {
    const double dotprod = new_axes[0][0] * new_axes[2][0] + new_axes[0][1] * new_axes[2][1] + new_axes[0][2] * new_axes[2][2];
    if (std::abs(dotprod) > 1.0e-6)
      throw runtime_error("Axes defining the quantization of spin must be orthogonal.");
  }
  for (int i = 0; i != 3; ++i)
    new_axes[1][i] = ((new_axes[2][fwd[i]] * new_axes[0][bck[i]]) - (new_axes[0][fwd[i]] * new_axes[2][bck[i]]));

  array<double, 3> factors;
  for (int i = 0; i != 3; ++i) {
    double tmp = 0.0;
    for (int j = 0; j != 3; ++j) {
      tmp += new_axes[i][j] * new_axes[i][j];
    }
    factors[i] = 1.0 / std::sqrt(tmp);
  }

  shared_ptr<Matrix> out = default_axes->clone();
  for (int i = 0; i != 3; ++i)
    for (int j = 0; j != 3; ++j)
      out->element(i, j) = factors[j] * new_axes[j][i];

  out->print("New spin quantization axes");
#ifndef NDEBUG
  // Ensure the rotation matrix defining spin quantization axes is unitary
  auto iden = out->clone();
  iden->unit();
  assert(((*out ^ *out) - *iden).rms() < 1.0e-8);
#endif

  return out;
}


Pseudospin::Pseudospin(const int _nspin, shared_ptr<const PTree> _idata) : nspin_(_nspin), nspin1_(_nspin + 1), idata_(_idata) {

  VectorB spinvals(nspin1_);
  for (int i = 0; i != nspin1_; ++i)
    spinvals[i] = (nspin_ / 2.0) - i;
  update_spin_matrices(spinvals);
}


void Pseudospin::compute(const ZHarrison& zfci) {

  // Which ranks of extended Stevens operators to use
  // Default should grab the nonzero time-reversal symmetric orders, but can be specified in input
  vector<int> ranks = {};
  for (int i = 2; i <= nspin_; i += 2)
    ranks.push_back(i);
  const shared_ptr<const PTree> eso_ranks = idata_->get_child_optional("ranks");
  if (eso_ranks) {
    ranks = {};
    for (auto& i : *eso_ranks)
      ranks.push_back(lexical_cast<int>(i->data()));
  }

  cout << setprecision(8);
  cout << endl << "    ********      " << endl;
  cout << endl << "    Modeling Pseudospin Hamiltonian for S = " << spin_val(nspin_) << endl;

  vector<Stevens_Operator> ESO = build_extended_stevens_operators(ranks);

  if (idata_->get<bool>("print_operators", false)) {
    cout << "Number of Stevens operators = " << ESO.size() << endl;
    for (int i = 0; i != ESO.size(); ++i)
      ESO[i].print();
  }

  if (nspin_ > 0) {

    compute_numerical_hamiltonian(zfci, zfci.jop()->coeff_input()->active_part());

    shared_ptr<const Matrix> mag_axes = identify_magnetic_axes();
    spin_axes_ = read_axes(mag_axes);

    for (int i = 0; i != 3; ++i) {
      zfci2_mu_[i] = make_shared<ZMatrix>(nspin1_, nspin1_);
      zfci2_spin_[i] = make_shared<ZMatrix>(nspin1_, nspin1_);
      zfci2_orbang_[i] = make_shared<ZMatrix>(nspin1_, nspin1_);
      for (int j = 0; j != 3; ++j) {
        *zfci2_mu_[i] += spin_axes_->element(j, i) * *zfci_mu_[j];
        *zfci2_spin_[i] += spin_axes_->element(j, i) * *zfci_spin_[j];
        *zfci2_orbang_[i] += spin_axes_->element(j, i) * *zfci_orbang_[j];
      }
    }

    shared_ptr<const ZMatrix> spinham_s = compute_spin_eigenvalues();

    if (nspin_ > 1) {
      ESO = extract_hamiltonian_parameters(ESO, spinham_s);
      shared_ptr<Matrix> dtens = compute_Dtensor(ESO);
    }
  } else {
    cout << "    There is no zero-field splitting or g-tensor to compute for an S = 0 system." << endl;
  }
}


// Compute S_x, S_y, and S_z plus the raising and lowering operators operators in pseudospin basis
void Pseudospin::update_spin_matrices(VectorB spinvals) {
  assert(spinvals.size() == nspin1_);
  for (int i = 0; i != nspin1_ / 2; ++i) {
    assert(std::abs(spinvals[i] + spinvals[nspin_ - i]) < 1.0e-6);
  }

  for (int i = 0; i != 3; ++i)
    spin_xyz_[i] = make_shared<ZMatrix>(nspin1_, nspin1_);
  spin_plus_ = make_shared<ZMatrix>(nspin1_, nspin1_);
  spin_minus_ = make_shared<ZMatrix>(nspin1_, nspin1_);

  const double sval = spinvals[0];
  const double ssp1 = sval * (sval + 1.0);

  for (int i = 0; i != nspin1_; ++i) {
    spin_xyz_[2]->element(i,i) = spinvals[i];
    if (i < nspin_)
      spin_plus_->element(i,i+1) = std::sqrt(ssp1 - spinvals[i]*spinvals[i+1]);
    if (i > 0)
      spin_minus_->element(i,i-1) = std::sqrt(ssp1 - spinvals[i]*spinvals[i-1]);
  }

  spin_xyz_[0]->add_block( 0.5, 0, 0, nspin1_, nspin1_, spin_plus_);
  spin_xyz_[0]->add_block( 0.5, 0, 0, nspin1_, nspin1_, spin_minus_);
  spin_xyz_[1]->add_block( complex<double>( 0.0, -0.5), 0, 0, nspin1_, nspin1_, spin_plus_);
  spin_xyz_[1]->add_block( complex<double>( 0.0,  0.5), 0, 0, nspin1_, nspin1_, spin_minus_);
}


// Compute numerical pseudospin Hamiltonian by diagonalizing S_z matrix
void Pseudospin::compute_numerical_hamiltonian(const ZHarrison& zfci, shared_ptr<const RelCoeff_Block> active_coeff) {
  const complex<double> imag(0.0, 1.0);
  const int norb = zfci.norb();

  // First, we create matrices of the magnetic moment in atomic orbital basis
  array<shared_ptr<ZMatrix>,3> magnetic_moment;
  RelSpinInt ao_spin(zfci.geom());

#ifndef NDEBUG
  { // Calculation of time-reversal matrix assumes it has this form in MO basis, so let's verify...
    auto ao_trev = make_shared<const RelTRevInt>(zfci.geom());
    auto mo_trev = make_shared<ZMatrix>(*active_coeff % *ao_trev * *active_coeff->get_conjg());
    auto mo_trev_exp = make_shared<ZMatrix>(2*norb, 2*norb);
    auto identity = make_shared<ZMatrix>(norb, norb);
    identity->unit();
    mo_trev_exp->add_block( 1.0, norb, 0, norb, norb, identity);
    mo_trev_exp->add_block(-1.0, 0, norb, norb, norb, identity);
    assert((*mo_trev - *mo_trev_exp).rms() < 1.0e-10);
  }
#endif

  { // Compute the matrix representation of the time-reversal operator   (This matrix + complex conjugation)
    trev_h_ = make_shared<ZMatrix>(nspin1_, nspin1_);
    trev_h_->zero();
    const int maxa = std::min(zfci.nele(), zfci.norb());
    const int mina = std::max(zfci.nele() - zfci.norb(), 0);

    /********************/
    // Redundant, but need this info to be available
    vector<int> aniso_state;
    aniso_state.resize(nspin1_);
    for (int i = 0; i != nspin1_; ++i)
      aniso_state[i] = i;

    // aniso_state can be used to request mapping excited states instead
    const shared_ptr<const PTree> exstates = idata_->get_child_optional("states");
    if (exstates) {
      aniso_state = {};
      for (auto& i : *exstates)
        aniso_state.push_back(lexical_cast<int>(i->data()) - 1);
      if (aniso_state.size() != nspin1_)
        throw runtime_error("Aniso:  Wrong number of states requested for this S value (should be " + to_string(nspin1_) + ")");
      for (int i = 0; i != nspin1_; ++i)
        if (aniso_state[i] < 0 || aniso_state[i] >= zfci.nstate())
         throw runtime_error("Aniso:  Invalid state requested (should be between 1 and " + to_string(zfci.nstate()) + ")");
    }
    /********************/

    vector<array<int,2>> ab = {};
    for (int j = maxa; j >= mina; --j)
      ab.push_back({{j, zfci.nele()-j}});

    // Loop over spin sectors
    for (int k = 0; k != ab.size(); ++k) {

      // bra
      auto dvec_i = zfci.cc()->find(ab[k][0], ab[k][1]);

      // ket (just flip alpha and beta - all other contributions will be zero)
      auto dvec_j = zfci.cc()->find(ab[k][1], ab[k][0]);

      auto det_i = dvec_i->data(0)->det();
      auto det_j = dvec_j->data(0)->det();

      // Loop over determinants
      for (auto& ia : det_i->string_bits_a()) {
        for (auto& ib : det_i->string_bits_b()) {
          const int pos1_i = det_i->lexical<0>(ia);
          const int pos2_i = det_i->lexical<1>(ib);
          const int fullpos_i = pos1_i*det_i->lenb() + pos2_i;

          const int pos1_j = det_j->lexical<0>(ib);
          const int pos2_j = det_j->lexical<1>(ia);
          const int fullpos_j = pos1_j*det_j->lenb() + pos2_j;

          // We get a -1 from each beta electron in the ket
          const double sign1 = (ia.count() % 2 == 0) ? 1.0 : -1.0;

          // We can also have a sign change due to the reordering of orbitals within the Slater Determinant
          const double sign2 = ((ia.count() * ib.count()) % 2 == 0) ? 1.0 : -1.0;

          // Loop over pairs of ZFCI eigenstates
          for (int i = 0; i != nspin1_; ++i) {
            for (int j = 0; j != nspin1_; ++j) {

              // Since time-reversal includes the complex conjugation operator, both the bra and ket use conjugate here
              const complex<double> cival_ist = std::conj(dvec_i->data(aniso_state[i])->data(fullpos_i));
              const complex<double> cival_jst = std::conj(dvec_j->data(aniso_state[j])->data(fullpos_j));

              trev_h_->element(i, j) += sign1 * sign2 * cival_ist * cival_jst;
            }
          }
        }
      }
    }
    trev_h_->print(" Time-reversal matrix in ZFCI states", 24);
  }

  { // spin angular momentum
    for (int i = 0; i != 3; ++i) {
      magnetic_moment[i] = ao_spin(i)->copy();
      magnetic_moment[i]->scale(g_elec__);
    }
  }

  array<shared_ptr<ZMatrix>,3> ao_orbang;
  { // orbital angular momentum
    // TODO For geometries with only one metal atom, use that atom's position as default mcoord
    const array<double, 3> mcoord = idata_->get_array<double,3>("center", array<double, 3>({{0.0, 0.0, 0.0}}));
    const int n = zfci.geom()->nbasis();

    array<shared_ptr<ZMatrix>,3> angmom_large;
    array<array<shared_ptr<ZMatrix>,4>,3> angmom_small;

    {
      AngMom angmom(zfci.geom(), mcoord);
      array<shared_ptr<Matrix>,3> mom = angmom.compute();
      for (int i = 0; i != 3; ++i)
        angmom_large[i] = make_shared<ZMatrix>(*mom[i], imag);
    }
    {
      auto smallmom = make_shared<Small1e<AngMomBatch, array<double,3>>>(zfci.geom(), mcoord);
      for (int i = 0; i != 3; ++i)
        for (int j = 0; j != 4; ++j)
          angmom_small[i][j] = make_shared<ZMatrix>((*smallmom)[4*i+j], imag);
    }

    const complex<double> w(0.25/(c__*c__));
    const complex<double> wi(0.0, w.real());
    for (int i = 0; i != 3; ++i) {
      ao_orbang[i] = magnetic_moment[i]->clone();

      ao_orbang[i]->add_block(1.0, 0, 0, n, n, angmom_large[i]);
      ao_orbang[i]->add_block(1.0, n, n, n, n, angmom_large[i]);

      ao_orbang[i]->add_block(  w, 2*n, 2*n, n, n, angmom_small[i][0]);
      ao_orbang[i]->add_block(  w, 3*n, 3*n, n, n, angmom_small[i][0]);
      ao_orbang[i]->add_block( wi, 2*n, 2*n, n, n, angmom_small[i][1]);
      ao_orbang[i]->add_block(-wi, 3*n, 3*n, n, n, angmom_small[i][1]);
      ao_orbang[i]->add_block( wi, 2*n, 3*n, n, n, angmom_small[i][2]);
      ao_orbang[i]->add_block( wi, 3*n, 2*n, n, n, angmom_small[i][2]);
      ao_orbang[i]->add_block(  w, 2*n, 3*n, n, n, angmom_small[i][3]);
      ao_orbang[i]->add_block( -w, 3*n, 2*n, n, n, angmom_small[i][3]);

      (*magnetic_moment[i]).add_block(1.0, 0, 0, 4*n, 4*n, ao_orbang[i]);
    }
  }

  // By default, just use the ground states
  vector<int> aniso_state;
  aniso_state.resize(nspin1_);
  ref_energy_.resize(nspin1_);
  for (int i = 0; i != nspin1_; ++i)
    aniso_state[i] = i;

  // aniso_state can be used to request mapping excited states instead
  const shared_ptr<const PTree> exstates = idata_->get_child_optional("states");
  if (exstates) {
    aniso_state = {};
    for (auto& i : *exstates)
      aniso_state.push_back(lexical_cast<int>(i->data()) - 1);
    if (aniso_state.size() != nspin1_)
      throw runtime_error("Aniso:  Wrong number of states requested for this S value (should be " + to_string(nspin1_) + ")");
    for (int i = 0; i != nspin1_; ++i)
      if (aniso_state[i] < 0 || aniso_state[i] >= zfci.nstate())
        throw runtime_error("Aniso:  Invalid state requested (should be between 1 and " + to_string(zfci.nstate()) + ")");
    cout << "    For the following states:  ";
    for (int i = 0; i != nspin1_; ++i)
      cout << aniso_state[i] << "  ";
    cout << endl;
  } else {
    cout << "    For the ground spin-manifold" << endl;
  }
  cout << endl;

  for (int i = 0; i != nspin1_; ++i)
    ref_energy_[i] = zfci.energy()[aniso_state[i]];

  // Compute spin matrices in the basis of ZFCI Hamiltonian eigenstates
  for (int i = 0; i != 3; ++i) {
    zfci_mu_[i] = make_shared<ZMatrix>(nspin1_, nspin1_);
    zfci_spin_[i] = make_shared<ZMatrix>(nspin1_, nspin1_);
    zfci_orbang_[i] = make_shared<ZMatrix>(nspin1_, nspin1_);
  }
  for (int i = 0; i != nspin1_; ++i) {
    for (int j = 0; j != nspin1_; ++j) {
      shared_ptr<Kramers<2,ZRDM<1>>> temprdm = zfci.rdm1(aniso_state[i], aniso_state[j]);
      if (!temprdm->exist({1,0})) {
        cout << " * Need to generate an off-diagonal rdm of zeroes." << endl;
        temprdm->add({1,0}, temprdm->at({0,0})->clone());
      }
      shared_ptr<const ZRDM<1>> tmp = expand_kramers<1,complex<double>>(temprdm, norb);

      auto rdmmat = make_shared<ZMatrix>(norb * 2, norb * 2);
      copy_n(tmp->data(), tmp->size(), rdmmat->data());

      ZMatrix modensity (2 * norb, 2 * norb);
      modensity.copy_block(0, 0, 2 * norb, 2 * norb, rdmmat);
      ZMatrix aodenconj = (*active_coeff * *modensity.get_conjg() ^ *active_coeff);

      for (int k = 0; k != 3; ++k) {
        zfci_mu_[k]->element(i,j) = aodenconj.dot_product(*magnetic_moment[k]);
        zfci_spin_[k]->element(i,j) = aodenconj.dot_product(*ao_spin(k));
        zfci_orbang_[k]->element(i,j) = aodenconj.dot_product(*ao_orbang[k]);
      }
    }
  }

#if 0
  for (int k = 0; k != 3; ++k) {
    zfci_mu_[k]->print("Magnetic moment in ZFCI basis - " + to_string(k), 24);
    zfci_spin_[k]->print("Spin angular momentum in ZFCI basis - " + to_string(k), 24);
    zfci_orbang_[k]->print("Orbital angular momentum in ZFCI basis - " + to_string(k), 24);
  }
#endif

  // We will subtract out average energy so the pseudospin Hamiltonian is traceless
  complex<double> energy_avg = 0.0;
  for (int i = 0; i != nspin1_; ++i)
    energy_avg += ref_energy_[i];
  energy_avg /= nspin1_;

  // Now build up the numerical pseudospin Hamiltonian!
  spinham_h_ = make_shared<ZMatrix>(nspin1_, nspin1_);
  for (int i = 0; i != nspin1_; ++i) {
    spinham_h_->element(i,i) = ref_energy_[i] - energy_avg;
  }

}


shared_ptr<const Matrix> Pseudospin::identify_magnetic_axes() const {
  auto Atensor = make_shared<const Matrix>(3, 3);
  shared_ptr<Matrix> Atransform;
  VectorB Aeig(3);
  {
    auto temp = make_shared<ZMatrix>(3, 3);
    for (int i = 0; i != 3; ++i)
      for (int j = 0; j != 3; ++j)
        for (int k = 0; k != nspin1_; ++k)
          for (int l = 0; l != nspin1_; ++l)
            temp->element(i, j) += 0.5 * zfci_mu_[i]->element(k, l) * zfci_mu_[j]->element(l, k);

    Atensor = temp->get_real_part();
    assert(temp->get_imag_part()->rms() < 1.0e-10);
    Atransform = Atensor->copy();
    Atransform->diagonalize(Aeig);

    // All eigenvalues of A should be positive, since they are proportional to squares of the principle g-values
    assert(Aeig[0] > 0.0 && Aeig[1] > 0.0 && Aeig[2] > 0.0);

    // Reorder eigenvectors so we quantize spin along the most anisotropic g-axis, rather than just the greatest g
    const double Asqrt_avg = (std::sqrt(Aeig[0]) + std::sqrt(Aeig[1]) + std::sqrt(Aeig[2])) / 3.0;
    if (std::sqrt(Aeig[1]) - Asqrt_avg > 0.0) {
      const double temp = Aeig[0];
      Aeig[0] = Aeig[2];
      Aeig[2] = temp;
      auto tmp = Atransform->copy();
      for (int i = 0; i != 3; ++i) {
        tmp->element(i, 0) = Atransform->element(i, 2);
        tmp->element(i, 2) = Atransform->element(i, 0);
      }
      Atransform = tmp;
    }
#ifndef NDEBUG
    auto Adiag = Atransform->clone();
    for (int i = 0; i != 3; ++i)
      Adiag->element(i, i) = Aeig[i];
    assert((*Atensor - (*Atransform * *Adiag ^ *Atransform)).rms() < 1.0e-10);
    assert(std::abs(std::sqrt(Aeig[2]) - Asqrt_avg) > std::abs(std::sqrt(Aeig[1]) - Asqrt_avg));
    assert(std::abs(std::sqrt(Aeig[2]) - Asqrt_avg) > std::abs(std::sqrt(Aeig[0]) - Asqrt_avg));
#endif

    Atensor->print("A tensor");
    cout << endl;
    for (int i = 0; i != 3; ++i)
      cout << " *** A tensor eigenvalue " << i << " = " << Aeig[i] << endl;
    cout << endl;
  }

  {
    auto gtensor = make_shared<Matrix>(3, 3);
    gtensor->zero();
    array<double,3> gval;
    const double factor = 12.0 / (nspin_ * (0.5 * nspin_ + 1.0) * (nspin_ + 1.0)); //  6.0 / ( S * (S+1) * (2S+1) )
    if (nspin_ > 2)
      cout << "  **  Use caution:  This mapping to the pseudospin Hamiltonian is approximate.  (3rd-order and above terms in A are neglected.)" << endl;
    for (int i = 0; i != 3; ++i) {
      gval[i] = 2.0 * std::sqrt(factor * Aeig[i]);
      gtensor->element(i, i) = gval[i];
    }

    *gtensor = (*Atransform * *gtensor ^ *Atransform);
    gtensor->print("g-tensor in the original coordinate system");

    cout << endl;
    auto Gtensor = make_shared<Matrix>(*gtensor ^ *gtensor);
    Gtensor->print("G-tensor in the original coordinate system");
    cout << endl;

    assert((*Gtensor - 4.0 * factor * *Atensor).rms() < 1.0e-8);

    cout << "  Main axes of magnetic anisotropy:" << endl;
    for (int i = 0; i != 3; ++i) {
      cout << "   " << i << " |g_" << i << "| = " << setw(12) << gval[i] << ",  axis = ( ";
      cout << setw(12) << Atransform->element(0, i) << ", ";
      cout << setw(12) << Atransform->element(1, i) << ", ";
      cout << setw(12) << Atransform->element(2, i) << " )" << endl;
    }
    cout << endl;
  }

  return Atransform;
}


shared_ptr<const ZMatrix> Pseudospin::compute_spin_eigenvalues() const {

  // Diagonalize S_z to get pseudospin eigenstates as combinations of ZFCI Hamiltonian eigenstates
  ZMatrix transform(nspin1_, nspin1_);
  const string diagset = idata_->get<string>("diagop", "Mu");
  if (diagset != "Mu" && diagset != "J" && diagset != "S" && diagset != "L")
    throw runtime_error("Sorry, the only options for which angular momentum to diagonalize are S, L, J and Mu for the magnetic moment");

  for (int i = 0; i != 3; ++i) {
    if (diagset == "Mu")
      transform += spin_axes_->element(i, 2) * *zfci_mu_[i];
    if (diagset == "S" || diagset == "J")
      transform += spin_axes_->element(i, 2) * *zfci_spin_[i];
    if (diagset == "L" || diagset == "J")
      transform += spin_axes_->element(i, 2) * *zfci_orbang_[i];
  }
  VectorB zeig(nspin1_);
#ifndef NDEBUG
  auto spinmat_to_diag = transform.copy();
#endif
  transform.diagonalize(zeig);

  { // Reorder eigenvectors so positive M_s come first
    ZMatrix tempm(nspin1_, nspin1_);
    VectorB tempv = *zeig.clone();
    for (int i = 0; i != nspin1_; ++i) {
      tempv[i] = zeig[nspin_ - i];
      tempm.copy_block(0, i, nspin1_, 1, transform.slice(nspin_ - i, nspin_ - i + 1));
    }
    transform = tempm;
    zeig = tempv;
  }

  { // Adjust the phases of eigenvectors to ensure proper time-reversal symmetry (using the matrix form of the time-reversal operator)
    auto trev_s = make_shared<ZMatrix>(transform % *trev_h_ * *transform.get_conjg());
    for (int k = 0; k <= nspin_ / 2; ++k) {
      const double target_phase = (k % 2 == 0) ? pi__ : 0.0;
      const double phase_error = std::arg(trev_s->element(k, nspin_ - k)) - target_phase;
      const complex<double> adjust = std::polar(1.0, 0.5 * phase_error);
      for (int i = 0; i != nspin1_; ++i)
        transform.element(i, k) = adjust * transform.element(i, k);
      if (nspin_ % 2 == 1 || k < nspin_ / 2) {
        for (int i = 0; i != nspin1_; ++i)
          transform.element(i, nspin_ - k) = adjust * transform.element(i, nspin_ - k);
      }
    }
  }

  { // Adjust the phase to make the (M+1, M) elements of the raising operator real (default choice)
    const string diagset = idata_->get<string>("diagop", "Mu");
    ZMatrix raising_op(nspin1_, nspin1_);

    if (diagset == "Mu") {
      raising_op.add_block(complex<double>(1.0, 0.0), 0, 0, nspin1_, nspin1_, transform % *zfci2_mu_[0] * transform);
      raising_op.add_block(complex<double>(0.0, 1.0), 0, 0, nspin1_, nspin1_, transform % *zfci2_mu_[1] * transform);
    }
    if (diagset == "S" || diagset == "J") {
      raising_op.add_block(complex<double>(1.0, 0.0), 0, 0, nspin1_, nspin1_, transform % *zfci2_spin_[0] * transform);
      raising_op.add_block(complex<double>(0.0, 1.0), 0, 0, nspin1_, nspin1_, transform % *zfci2_spin_[1] * transform);
    }
    if (diagset == "L" || diagset == "J") {
      raising_op.add_block(complex<double>(1.0, 0.0), 0, 0, nspin1_, nspin1_, transform % *zfci2_orbang_[0] * transform);
      raising_op.add_block(complex<double>(0.0, 1.0), 0, 0, nspin1_, nspin1_, transform % *zfci2_orbang_[1] * transform);
    }

    complex<double> adjust = 1.0;
    for (int i = nspin1_ / 2 - 1; i >= 0; --i) {

      double phase_check = std::arg(raising_op.element(i, i+1));
      if (i + 1 == nspin1_ / 2 && nspin1_ % 2 == 0)
        phase_check *= 0.5;
      adjust *= std::polar(1.0, phase_check);

      for (int j = 0; j != nspin1_; ++j) {
        transform.element(j, i) = adjust * transform.element(j, i);
        transform.element(j, nspin_ - i) = std::conj(adjust) * transform.element(j, nspin_ - i);
      }
    }
  }
  transform.print("Pseudospin eigenvector matrix with time-reversal symmetry fixed");

  // For testing arbitrary phase shifts applied to pseudospin eigenfunctions
  const shared_ptr<const PTree> phase_input = idata_->get_child_optional("phases");
  if (phase_input) {
    vector<double> phase_adjust = {};
    for (auto& i : *phase_input)
      phase_adjust.push_back(lexical_cast<double>(i->data()));

    if (phase_adjust.size() != nspin1_ / 2)
      throw runtime_error("Sorry, you seem to be trying to adjust the phase of pseudospin eigenfunctions.  We expect " + to_string(nspin1_/2) + " phases, one for each pair of time-reversal symmetric states.");

    for (int i = 0; i != nspin1_ / 2; ++i) {
      const complex<double> adjustment = std::polar(1.0, phase_adjust[i]);
      cout << "  **  The phase of the m_s = " << " " + spin_val(nspin_ - 2*i) << " pseudospin function will be shifted by " << std::setw(11) <<  phase_adjust[i] << " radians." << endl;
      cout << "  **  The phase of the m_s = " << "-" + spin_val(nspin_ - 2*i) << " pseudospin function will be shifted by " << std::setw(11) << -phase_adjust[i] << " radians." << endl;
      for (int j = 0; j != nspin1_; ++j) {
        transform.element(j, i) = adjustment * transform.element(j, i);
        transform.element(j, nspin_ - i) = std::conj(adjustment) * transform.element(j, nspin_ - i);
      }
    }
    cout << endl;
  }

  const double phase_input_2 = idata_->get<double>("phase_full", 0.0);
  if (phase_input_2 != 0.0) {
    transform.scale(std::polar(1.0, phase_input_2));
    cout << "  **  The phase of all pseudospin functions will be shifted by " << setw(4) <<  phase_input_2 << " radians.  (This should have no effect.)" << endl << endl;
  }

#ifndef NDEBUG
  {
    ZMatrix diagspin(nspin1_, nspin1_);
    diagspin.zero();
    for (int i = 0; i != nspin1_; ++i)
      diagspin.element(i, i) = zeig[i];
    diagspin = transform * diagspin ^ transform;
    assert((diagspin - *spinmat_to_diag).rms() < 1.0e-6);
  }
#endif

  cout << "    The z-axis is set to (" << spin_axes_->element(0, 2) << ", " << spin_axes_->element(1, 2) << ", " << spin_axes_->element(2, 2) << ")." << endl << endl;
  for (int i = 0; i != nspin1_; ++i)
    cout << "    Pseudospin eigenvalue " << i+1 << " = " << setw(12) << zeig[i] << endl;

  // We can no longer use this option, since I made this function const...
  //if (numerical_eig) {
  //  cout << "  **  By request, we compute the Hamiltonian using eigenvalues of S_z, rather than the canonical m_s values" << endl;
  //  update_spin_matrices(zeig);
  //}

  shared_ptr<ZMatrix> spinham_s = make_shared<ZMatrix>(transform % *spinham_h_ * transform);
  array<shared_ptr<ZMatrix>, 3> mu_s;
  array<shared_ptr<ZMatrix>, 3> spin_s;
  array<shared_ptr<ZMatrix>, 3> orbang_s;
  auto trev_s = make_shared<ZMatrix>(transform % *trev_h_ * *transform.get_conjg());
  for (int i = 0; i != 3; ++i) {
    mu_s[i] = make_shared<ZMatrix>(transform % *zfci2_mu_[i] * transform);
    spin_s[i] = make_shared<ZMatrix>(transform % *zfci2_spin_[i] * transform);
    orbang_s[i] = make_shared<ZMatrix>(transform % *zfci2_orbang_[i] * transform);
  }

  spinham_s->print("Pseudospin Hamiltonian", 24);
  cout << endl;
  trev_s->print("Time-reversal operator in pseudospin states", 24);
  cout << endl;
  for (int i = 0; i != 3; ++i) {
    mu_s[i]->print("Magnetic moment matrix, component " + to_string(i), 24);
    spin_s[i]->print("Spin matrix, component " + to_string(i), 24);
    orbang_s[i]->print("Orbital moment matrix, component " + to_string(i), 24);
    cout << endl;
  }
  cout << endl;

  if (!is_t_symmetric(*spinham_s, /*hermitian*/true, /*time reversal*/true))
    throw runtime_error("The spin Hamiltonian seems to not have proper time-reversal symmetry.  Check that your spin value and states mapped are reasonable.");

  cout << endl;
  for (int i = 0; i != 3; ++i) {
    assert(is_t_symmetric(*mu_s[i], /*hermitian*/true, /*time reversal*/false));
  }
  assert(is_t_symmetric(*spinham_s, /*hermitian*/true, /*time reversal*/true));

#ifndef NDEBUG
  { // Verify that the time-reversal operator has taken the correct form
    // The numerical threshold is arbitrary and might need adjustment
    auto trev_s = make_shared<ZMatrix>(transform % *trev_h_ * *transform.get_conjg());
    auto trev_target = make_shared<ZMatrix>(nspin1_, nspin1_);
    for (int k = 0; k <= nspin_; ++k) {
      const double val = (k % 2 == 0) ? -1.0 : 1.0;
      trev_target->element(k, nspin_ - k) = val;
    }
    if (phase_input_2 != 0.0) {
      trev_target->scale(std::polar(1.0, -2.0 * phase_input_2));
      assert((*trev_s - *trev_target).rms() < 1.0e-10);
    }
  }
#endif

  return spinham_s;
}


// Extract Stevens coefficients from the numerical pseudospin Hamiltonian
vector<Stevens_Operator> Pseudospin::extract_hamiltonian_parameters(const vector<Stevens_Operator> param, shared_ptr<const ZMatrix> spinham_s) const {
  const int nop = param.size();
  vector<Stevens_Operator> out = param;

  // s2h is the transformation matrix to build pseudospin Hamiltonian from extended Stevens operators
  // Rows correspond to pairs of pseudospins (SS, S-1S, S-2S...)
  // Columns correspond to Stevens coefficients (B22, B2-2, B21...)
  // Note that we look over the first indices before the second, so we can copy data between vectors and matrices and have it come out in the right order
  auto s2h = make_shared<ZMatrix>(nspin1_ * nspin1_, nop);

  for (int i = 0; i != nop; ++i)
    s2h->copy_block(0, i, nspin1_ * nspin1_, 1, param[i].matrix()->data());

  auto checkham = make_shared<ZMatrix>(nspin1_, nspin1_);

  { // Convert from the pseudospin Hamiltonian to the Stevens coefficients using the left-inverse of s2h
    // We use complex algebra, but the result should be purely real
    auto h2s = make_shared<ZMatrix>(nop, nspin1_ * nspin1_);
    ZMatrix s2h_sqinv = *s2h % *s2h;
    s2h_sqinv.inverse();
    *h2s = s2h_sqinv ^ *s2h;
    assert((*h2s * *s2h).is_identity());

    auto spinham_vec = make_shared<ZMatrix>(nspin1_ * nspin1_,1);
    spinham_vec->copy_block(0, 0, nspin1_ * nspin1_, 1, spinham_s->element_ptr(0,0));
    ZMatrix stevop_vec = *h2s * *spinham_vec;

    ZMatrix check_spinham_vec = *s2h * stevop_vec;
    checkham->copy_block(0, 0, nspin1_, nspin1_, check_spinham_vec.element_ptr(0,0));

    for (int i = 0; i != nop; ++i) {
      out[i].set_coeff(std::real(stevop_vec.element(i, 0)));
      if (std::abs(std::imag(stevop_vec.element(i, 0))) > 1.0e-8)
        throw runtime_error("For some reason, we have obtained a complex coefficient for an extended Stevens operator...  It should be real.");
    }
  }

  cout << endl << "    Stevens coefficients:  " << endl << endl;
  for (int i = 0; i != out.size(); ++i)
    cout << "    " << setw(8) << out[i].coeff_name() << " = " << setw(12) << out[i].coeff() << endl;
  cout << endl;

  const double checkham_error = (*checkham - *spinham_s).rms();
  VectorB shenergies(nspin1_);
  checkham->diagonalize(shenergies);

  if (checkham_error > 1.0e-8) {
    cout << "  **** CAUTION ****  The pseudospin Hamiltonian does not fully reproduce the ab initio Hamiltonian.  RMS error = " << checkham_error << endl;

    cout << "  ** Relative energies from the pseudospin Hamiltonian: " << endl;
    for (int i = nspin_; i >= 0; --i)
      cout << "     " << i << "  " << setw(12) << shenergies[i] - shenergies[0] << " E_h  =  " << setw(14) << (shenergies[i] - shenergies[0])*au2wavenumber__ << " cm-1" << endl;
    cout << endl;

    cout << "  ** Relative energies from the ab initio (relativistic configuration interaction) Hamiltonian: " << endl;
    for (int i = nspin_; i >= 0; --i)
      cout << "     " << i << "  " << setw(12) << ref_energy_[i] - ref_energy_[0] << " E_h  =  " << setw(14) << (ref_energy_[i] - ref_energy_[0])*au2wavenumber__ << " cm-1" << endl;
    cout << endl;
  } else {
    cout << "  ** Relative energies: " << endl;
    for (int i = nspin_; i >= 0; --i) {
      cout << "     " << i << "  " << setw(12) << shenergies[i] - shenergies[0] << " E_h  =  " << setw(14) << (shenergies[i] - shenergies[0])*au2wavenumber__ << " cm-1" << endl;
      assert(std::abs(shenergies[i] - shenergies[0] - ref_energy_[i] + ref_energy_[0]) < 1.0e-7);
    }
    cout << endl;
  }

  return out;
}


shared_ptr<Matrix> Pseudospin::compute_Dtensor(const vector<Stevens_Operator> input) {
  auto out = make_shared<Matrix>(3, 3);
  out->zero();

  // Get D from second-order extended Stevens Operators
  for (int i = 0; i != input.size(); ++i) {
    if (input[i].order() == 2) {
      switch (input[i].index()) {
        case  0:
          out->element(0, 0) -= 1.0 * input[i].coeff();
          out->element(1, 1) -= 1.0 * input[i].coeff();
          out->element(2, 2) += 2.0 * input[i].coeff();
          break;
        case -1:
          out->element(1, 2) += 0.5 * input[i].coeff();
          out->element(2, 1) += 0.5 * input[i].coeff();
          break;
        case  1:
          out->element(2, 0) += 0.5 * input[i].coeff();
          out->element(0, 2) += 0.5 * input[i].coeff();
          break;
        case -2:
          out->element(0, 1) += 1.0 * input[i].coeff();
          out->element(1, 0) += 1.0 * input[i].coeff();
          break;
        case  2:
          out->element(0, 0) += 1.0 * input[i].coeff();
          out->element(1, 1) -= 1.0 * input[i].coeff();
          break;
        default:
          throw logic_error("Some invalid operator was found in Pseudospin::compute_Dtensor(...)");
      }
    }
  }

  /**** PRINTOUT ***/

  shared_ptr<Matrix> Dtensor_diag = out->copy();
  Dtensor_diag->print("D tensor");
  cout << setprecision(8);
  VectorB Ddiag(3);
  Dtensor_diag->diagonalize(Ddiag);

  // Compute Davg so that it works even if D is not traceless (which shouldn't happen on accident)
  const double Davg = 1.0 / 3.0 * (Ddiag[0] + Ddiag[1] + Ddiag[2]);

  int jmax = 0;
  const array<int,3> fwd = {{ 1, 2, 0 }};
  const array<int,3> bck = {{ 2, 0, 1 }};
  if (std::abs(Ddiag[1]-Davg) > std::abs(Ddiag[jmax]-Davg)) jmax = 1;
  if (std::abs(Ddiag[2]-Davg) > std::abs(Ddiag[jmax]-Davg)) jmax = 2;

  cout << endl << "    Upon diagonalization," << endl;
  cout << "      Dxx = " << setw(12) << Ddiag[fwd[jmax]] << endl;
  cout << "      Dyy = " << setw(12) << Ddiag[bck[jmax]] << endl;
  cout << "      Dzz = " << setw(12) << Ddiag[jmax] << endl << endl;
  const double Dval = Ddiag[jmax] - 0.5*(Ddiag[fwd[jmax]] + Ddiag[bck[jmax]]);
  const double Eval = 0.5*(Ddiag[fwd[jmax]] - Ddiag[bck[jmax]]);
  cout << " ** D = " << setw(12) << Dval << " E_h = " << setw(14) << Dval * au2wavenumber__ << " cm-1" << endl;
  cout << " ** E = " << setw(12) << std::abs(Eval) << " E_h = " << setw(14) << std::abs(Eval * au2wavenumber__) << " cm-1" << endl;
  cout << " ** E / D = " << std::abs(Eval / Dval) << endl;

  cout << endl << " ** Axis of principle D-value (relative to spin quant. axes) = (" << Dtensor_diag->element(0, jmax) << ", " << Dtensor_diag->element(1, jmax) << ", " << Dtensor_diag->element(2, jmax) << ")" << endl;

  Matrix full_rotation = *spin_axes_ * *Dtensor_diag;
  cout << endl << " ** Axis of principle D-value (relative to input geometry)  =  (" << full_rotation.element(0, jmax) << ", " << full_rotation.element(1, jmax) << ", " << full_rotation.element(2, jmax) << ")" << endl;
  return out;
}


bool Pseudospin::is_t_symmetric(const ZMatrix& in, const bool hermitian, const bool t_symmetric, const double thresh) const {
  // The matrix must be either Hermitian or skew-Hermitian
  if (hermitian) {
    assert(in.is_hermitian(thresh));
  } else {
    assert((in + *in.transpose_conjg()).rms() < thresh);
  }
  assert(in.ndim() == in.mdim());
  assert(in.ndim() == nspin1_);
  const double t = t_symmetric ? 1.0 : -1.0;

  // The time-reversal operator is kramersop + complex conjugation
  ZMatrix kramersop(nspin1_, nspin1_);
  double val = -1.0;
  for (int i = 0; i != nspin1_; ++i) {
    kramersop.element(nspin_ - i, i) = val;
    val *= -1.0;
  }
  ZMatrix in_reversed = kramersop * *in.get_conjg() ^ kramersop;  // K O K^-1
  const double error = (in - (t * in_reversed)).rms();
  //(in - (t * in_reversed)).print("Time-reversal symmetry error");

  const bool out = error < thresh;
  return out;
}

