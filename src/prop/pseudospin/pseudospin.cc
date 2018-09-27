//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: pseudospin.cc
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Ryan D. Reynolds <rreynolds2018@u.northwestern.edu>
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


Pseudospin::Pseudospin(const int _nspin, shared_ptr<const Geometry> _geom, shared_ptr<const RelCIWfn> _ciwfn, shared_ptr<const PTree> _idata)
 : nspin_(_nspin), nspin1_(_nspin + 1), geom_(_geom), idata_(_idata), ciwfn_(_ciwfn) {
  norb_ = ciwfn_->nact();
  if (nspin_ >= ciwfn_->nstates())
    throw runtime_error("Error in pseudospin module:  Not enough states for the requested spin multiplicity.");

  VectorB spinvals(nspin1_);
  for (int i = 0; i != nspin1_; ++i)
    spinvals[i] = (nspin_ / 2.0) - i;
  update_spin_matrices(spinvals);
}


void Pseudospin::compute(const vector<double> energy_in, shared_ptr<const ZCoeff_Block> active_coeff) {

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

  eso_ = build_extended_stevens_operators(ranks);

  if (idata_->get<bool>("print_operators", false)) {
    cout << endl << "    Number of Stevens operators = " << eso_.size() << endl;
    for (int i = 0; i != eso_.size(); ++i)
      eso_[i].print();
  }
  cout << endl;

  if (nspin_ > 0) {

    // Computes spin, orbital angular momentum, Hamiltonian, and time-reversal operators in the basis of ZFCI eigenstates
    compute_numerical_hamiltonian(energy_in, active_coeff);

    // Compute G and diagonalize to give main magnetic axes, but allow the user to quantify spin along some other axis
    pair<shared_ptr<const Matrix>, array<double,3>> mag_info = identify_magnetic_axes();
    spin_axes_ = read_axes(mag_info.first);
    gval_ = mag_info.second;

    // Rotate spin & related matrices to our new set of axes
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

    // Diagonalize an angular momentum matrix, fix phase, and use the eigenvectors to transform Hamiltonian pseudospin basis
    shared_ptr<const ZMatrix> spinham_s = compute_spin_eigenvalues();

    if (nspin_ > 1) {
      // Decompose the pseudospin Hamiltonian to find the coefficient B_k^q that goes with each Extended Stevens Operator
      eso_ = extract_hamiltonian_parameters(eso_, spinham_s);

      // Then just extract the D-tensor!
      tuple<shared_ptr<const Matrix>, double, double> d_params = compute_dtensor(eso_);
      dtensor_ = get<0>(d_params);
      dval_ = get<1>(d_params);
      eval_ = get<2>(d_params);
    }
  } else {
    cout << endl << "    There is no zero-field splitting or g-tensor to compute for an S = 0 system." << endl;
  }
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
    if (abs(dotprod) > 1.0e-6)
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
    factors[i] = 1.0 / sqrt(tmp);
  }

  shared_ptr<Matrix> out = default_axes->clone();
  for (int i = 0; i != 3; ++i)
    for (int j = 0; j != 3; ++j)
      out->element(i, j) = factors[j] * new_axes[j][i];

  //out->print("New spin quantization axes");
#ifndef NDEBUG
  // Ensure the rotation matrix defining spin quantization axes is unitary
  auto iden = out->clone();
  iden->unit();
  assert(((*out ^ *out) - *iden).rms() < 1.0e-8);
#endif

  return out;
}


// Compute S_x, S_y, and S_z plus the raising and lowering operators operators in pseudospin basis
void Pseudospin::update_spin_matrices(VectorB spinvals) {
  assert(spinvals.size() == nspin1_);
  for (int i = 0; i != nspin1_ / 2; ++i) {
    assert(abs(spinvals[i] + spinvals[nspin_ - i]) < 1.0e-6);
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
      spin_plus_->element(i,i+1) = sqrt(ssp1 - spinvals[i]*spinvals[i+1]);
    if (i > 0)
      spin_minus_->element(i,i-1) = sqrt(ssp1 - spinvals[i]*spinvals[i-1]);
  }

  spin_xyz_[0]->add_block( 0.5, 0, 0, nspin1_, nspin1_, spin_plus_);
  spin_xyz_[0]->add_block( 0.5, 0, 0, nspin1_, nspin1_, spin_minus_);
  spin_xyz_[1]->add_block( complex<double>( 0.0, -0.5), 0, 0, nspin1_, nspin1_, spin_plus_);
  spin_xyz_[1]->add_block( complex<double>( 0.0,  0.5), 0, 0, nspin1_, nspin1_, spin_minus_);
}


// Compute numerical pseudospin Hamiltonian by diagonalizing S_z matrix
void Pseudospin::compute_numerical_hamiltonian(const vector<double> energy_in, shared_ptr<const ZCoeff_Block> active_coeff) {
  const complex<double> imag(0.0, 1.0);

  // First, we create matrices of the magnetic moment in atomic orbital basis
  array<shared_ptr<ZMatrix>,3> magnetic_moment;
  RelSpinInt ao_spin(geom_);

#ifndef NDEBUG
  { // Calculation of time-reversal matrix assumes it has this form in MO basis, so let's verify...
    auto ao_trev = make_shared<const RelTRevInt>(geom_);
    auto mo_trev = make_shared<ZMatrix>(*active_coeff % *ao_trev * *active_coeff->get_conjg());
    auto mo_trev_exp = make_shared<ZMatrix>(2*norb_, 2*norb_);
    auto identity = make_shared<ZMatrix>(norb_, norb_);
    identity->unit();
    mo_trev_exp->add_block( 1.0, norb_, 0, norb_, norb_, identity);
    mo_trev_exp->add_block(-1.0, 0, norb_, norb_, norb_, identity);
    assert((*mo_trev - *mo_trev_exp).rms() < 1.0e-10);
  }
#endif

  vector<int> aniso_state;
  { // By default, just use the ground states
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
        if (aniso_state[i] < 0 || aniso_state[i] >= energy_in.size())
         throw runtime_error("Aniso:  Invalid state requested (should be between 1 and " + to_string(energy_in.size()) + ")");
    }
  }

  { // Compute the matrix representation of the time-reversal operator   (This matrix + complex conjugation)
    trev_h_ = make_shared<ZMatrix>(nspin1_, nspin1_);
    trev_h_->zero();
    const int nele = ciwfn_->det()->first->nele();
    const int maxa = min(nele, norb_);
    const int mina = max(nele - norb_, 0);

    vector<array<int,2>> ab = {};
    for (int j = maxa; j >= mina; --j)
      ab.push_back({{j, nele-j}});

    // Loop over spin sectors
    for (int k = 0; k != ab.size(); ++k) {

      // bra
      auto dvec_i = ciwfn_->civectors()->find(ab[k][0], ab[k][1]);

      // ket (just flip alpha and beta - all other contributions will be zero)
      auto dvec_j = ciwfn_->civectors()->find(ab[k][1], ab[k][0]);

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
              const complex<double> cival_ist = conj(dvec_i->data(aniso_state[i])->data(fullpos_i));
              const complex<double> cival_jst = conj(dvec_j->data(aniso_state[j])->data(fullpos_j));

              trev_h_->element(i, j) += sign1 * sign2 * cival_ist * cival_jst;
            }
          }
        }
      }
    }
    //trev_h_->print(" Time-reversal matrix in ZFCI states", 24);
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
    array<double, 3> mcoord = idata_->get_array<double,3>("center", array<double, 3>({{0.0, 0.0, 0.0}}));

    if (idata_->get<bool>("anstrom", false)) {
      for (auto& i : mcoord) {
        i /= au2angstrom__;
      }
    }

    const int n = geom_->nbasis();

    array<shared_ptr<ZMatrix>,3> angmom_large;
    array<array<shared_ptr<ZMatrix>,4>,3> angmom_small;

    {
      AngMom angmom(geom_, mcoord);
      array<shared_ptr<Matrix>,3> mom = angmom.compute();
      for (int i = 0; i != 3; ++i)
        angmom_large[i] = make_shared<ZMatrix>(*mom[i], imag);
    }
    {
      auto smallmom = make_shared<Small1e<AngMomBatch, array<double,3>>>(geom_, mcoord);
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
      magnetic_moment[i]->scale(-0.5);
    }
  }

  // Load up the energies of the nspin1_ states
  // Default is to take from ZFCI, but we allow custom input so we can use correlated energies, such as from Dirac--NEVPT2
  ref_energy_.resize(nspin1_);
  const shared_ptr<const PTree> input_energy = idata_->get_child_optional("energies");
  if (input_energy) {
    if (input_energy->size() != nspin1_)
      throw runtime_error("Wrong number of energies given; one is needed for each of " + to_string(nspin1_) + " states.");
    auto en = input_energy->begin();
    for (int i = 0; i != nspin1_; ++i) {
      ref_energy_[i] = lexical_cast<double>((*en)->data());
      en++;
    }
    cout << "  *  Energies of Hamiltonian eigenstates are taken from input rather than relativistic FCI." << endl;
  } else {
    for (int i = 0; i != nspin1_; ++i)
      ref_energy_[i] = energy_in[aniso_state[i]];
  }

  // Compute spin matrices in the basis of ZFCI Hamiltonian eigenstates
  for (int i = 0; i != 3; ++i) {
    zfci_mu_[i] = make_shared<ZMatrix>(nspin1_, nspin1_);
    zfci_spin_[i] = make_shared<ZMatrix>(nspin1_, nspin1_);
    zfci_orbang_[i] = make_shared<ZMatrix>(nspin1_, nspin1_);
  }
  ZFCI_bare fci(ciwfn_);
  for (int i = 0; i != nspin1_; ++i) {
    for (int j = 0; j != nspin1_; ++j) {
      shared_ptr<Kramers<2,ZRDM<1>>> temprdm = fci.rdm1(aniso_state[i], aniso_state[j]);
      if (!temprdm->exist({1,0})) {
        cout << " * Need to generate an off-diagonal rdm of zeroes." << endl;
        temprdm->add({1,0}, temprdm->at({0,0})->clone());
      }
      shared_ptr<const ZRDM<1>> tmp = expand_kramers<1,complex<double>>(temprdm, norb_);

      auto rdmmat = make_shared<ZMatrix>(norb_ * 2, norb_ * 2);
      copy_n(tmp->data(), tmp->size(), rdmmat->data());

      ZMatrix modensity (2 * norb_, 2 * norb_);
      modensity.copy_block(0, 0, 2 * norb_, 2 * norb_, rdmmat);
      ZMatrix aodenconj = (*active_coeff * *modensity.get_conjg() ^ *active_coeff);

      for (int k = 0; k != 3; ++k) {
        zfci_mu_[k]->element(i,j) = aodenconj.dot_product(*magnetic_moment[k]);
        zfci_spin_[k]->element(i,j) = aodenconj.dot_product(*ao_spin(k));
        zfci_orbang_[k]->element(i,j) = aodenconj.dot_product(*ao_orbang[k]);
      }
    }
  }

  if (idata_->get<bool>("print_operators", false)) {
    for (int k = 0; k != 3; ++k) {
      zfci_mu_[k]->print("Magnetic moment in ZFCI basis - " + to_string(k), 24);
      zfci_spin_[k]->print("Spin angular momentum in ZFCI basis - " + to_string(k), 24);
      zfci_orbang_[k]->print("Orbital angular momentum in ZFCI basis - " + to_string(k), 24);
    }
  }

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


pair<shared_ptr<const Matrix>, array<double,3>> Pseudospin::identify_magnetic_axes() const {
  auto atensor = make_shared<const Matrix>(3, 3);
  array<double,3> gval;
  shared_ptr<Matrix> atransform;
  VectorB aeig(3);
  {
    auto temp = make_shared<ZMatrix>(3, 3);
    for (int i = 0; i != 3; ++i)
      for (int j = 0; j != 3; ++j)
        for (int k = 0; k != nspin1_; ++k)
          for (int l = 0; l != nspin1_; ++l)
            temp->element(i, j) += 0.5 * zfci_mu_[i]->element(k, l) * zfci_mu_[j]->element(l, k);

    atensor = temp->get_real_part();
    assert(temp->get_imag_part()->rms() < 1.0e-10);
    atransform = atensor->copy();
    atransform->diagonalize(aeig);

    // Zero out any numerically zero values, because -1.0e-16 would cause problems...
    for (int i = 0; i != 3; ++i)
      if (std::abs(aeig[i]) < 1.0e-12)
        aeig[i] = 0.0;

    // All eigenvalues of A should be positive, since they are proportional to squares of the principle g-values
    assert(aeig[0] >= 0.0 && aeig[1] >= 0.0 && aeig[2] >= 0.0);

    // Reorder eigenvectors so we quantize spin along the most anisotropic g-axis, rather than just the greatest g
    const double asqrt_avg = (sqrt(aeig[0]) + sqrt(aeig[1]) + sqrt(aeig[2])) / 3.0;
    if (sqrt(aeig[1]) - asqrt_avg > 0.0) {
      const double temp = aeig[0];
      aeig[0] = aeig[2];
      aeig[2] = temp;
      auto tmp = atransform->copy();
      for (int i = 0; i != 3; ++i) {
        tmp->element(i, 0) = atransform->element(i, 2);
        tmp->element(i, 2) = atransform->element(i, 0);
      }
      atransform = tmp;
    }
#ifndef NDEBUG
    auto adiag = atransform->clone();
    for (int i = 0; i != 3; ++i)
      adiag->element(i, i) = aeig[i];
    assert((*atensor - (*atransform * *adiag ^ *atransform)).rms() < 1.0e-10);
    assert(abs(sqrt(aeig[2]) - asqrt_avg) > abs(sqrt(aeig[1]) - asqrt_avg));
    assert(abs(sqrt(aeig[2]) - asqrt_avg) > abs(sqrt(aeig[0]) - asqrt_avg));
#endif

    //atensor->print("A tensor");
    //cout << endl;
    //for (int i = 0; i != 3; ++i)
    //  cout << " *** A tensor eigenvalue " << i << " = " << aeig[i] << endl;
    //cout << endl;
  }

  {
    auto gtensor = make_shared<Matrix>(3, 3);
    gtensor->zero();
    const double factor = 12.0 / (nspin_ * (0.5 * nspin_ + 1.0) * (nspin_ + 1.0)); //  6.0 / ( S * (S+1) * (2S+1) )
    if (nspin_ > 2) {
      cout << "  **  Note that there is some approximation in the determination of the g-tensor, " << endl;
      cout << "  **  since we are not separating the first-order and higher-order contributions to the magnetic moment." << endl << endl;
    }
    for (int i = 0; i != 3; ++i) {
      gval[i] = 2.0 * sqrt(factor * aeig[i]);
      gtensor->element(i, i) = gval[i];
    }

    *gtensor = (*atransform * *gtensor ^ *atransform);
    gtensor->print("g-tensor");

    // This would have the spatial and spin axes not aligned, which is weird and there's probably no reason to do
    //*gtensor = (*atransform * *gtensor);
    //gtensor->print("g-tensor (rotating only left side by A-diag form)");

    cout << endl;
    auto ggtensor = make_shared<Matrix>(*gtensor ^ *gtensor);
    ggtensor->print("G-tensor");
    cout << endl;

    assert((*ggtensor - 4.0 * factor * *atensor).rms() < 1.0e-8);

    cout << "  Main axes of magnetic anisotropy:" << endl;
    for (int i = 0; i != 3; ++i) {
      cout << "   " << i << " |g_" << i << "| = " << setprecision(5) << setw(8) << gval[i] << ",  axis = (";
      cout << setw(8) << atransform->element(0, i) << ", ";
      cout << setw(8) << atransform->element(1, i) << ", ";
      cout << setw(8) << atransform->element(2, i) << ")" << endl;
    }
    cout << endl;
  }

  pair<shared_ptr<const Matrix>, array<double,3>> out(atransform, gval);
  return out;
}


shared_ptr<const ZMatrix> Pseudospin::compute_spin_eigenvalues() const {

  // Diagonalize S_z to get pseudospin eigenstates as combinations of ZFCI Hamiltonian eigenstates
  ZMatrix transform(nspin1_, nspin1_);
  const string diagset = to_lower(idata_->get<string>("diagop", "Mu"));
  if (diagset != "mu" && diagset != "j" && diagset != "s" && diagset != "l")
    throw runtime_error("Sorry, the only options for which angular momentum to diagonalize are S, L, J and Mu for the magnetic moment");

  for (int i = 0; i != 3; ++i) {
    if (diagset == "mu")
      transform += spin_axes_->element(i, 2) * *zfci_mu_[i];
    if (diagset == "s" || diagset == "j")
      transform += spin_axes_->element(i, 2) * *zfci_spin_[i];
    if (diagset == "l" || diagset == "j")
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
    //trev_s->print("Time-reversal matrix before phase adjustment");
    for (int k = 0; k <= nspin_ / 2; ++k) {
      const double target_phase = (k % 2 == 0) ? pi__ : 0.0;
      const double phase_error = arg(trev_s->element(k, nspin_ - k)) - target_phase;
      const complex<double> adjust = polar(1.0, 0.5 * phase_error);
      for (int i = 0; i != nspin1_; ++i)
        transform.element(i, k) = adjust * transform.element(i, k);
      if (nspin_ % 2 == 1 || k < nspin_ / 2) {
        for (int i = 0; i != nspin1_; ++i)
          transform.element(i, nspin_ - k) = adjust * transform.element(i, nspin_ - k);
      }
    }
  }

  { // Adjust the phase to make the (M+1, M) elements of the raising operator real
    ZMatrix raising_op(nspin1_, nspin1_);

    if (diagset == "mu") {
      raising_op.add_block(complex<double>(1.0, 0.0), 0, 0, nspin1_, nspin1_, transform % *zfci2_mu_[0] * transform);
      raising_op.add_block(complex<double>(0.0, 1.0), 0, 0, nspin1_, nspin1_, transform % *zfci2_mu_[1] * transform);
    }
    if (diagset == "s" || diagset == "j") {
      raising_op.add_block(complex<double>(1.0, 0.0), 0, 0, nspin1_, nspin1_, transform % *zfci2_spin_[0] * transform);
      raising_op.add_block(complex<double>(0.0, 1.0), 0, 0, nspin1_, nspin1_, transform % *zfci2_spin_[1] * transform);
    }
    if (diagset == "l" || diagset == "j") {
      raising_op.add_block(complex<double>(1.0, 0.0), 0, 0, nspin1_, nspin1_, transform % *zfci2_orbang_[0] * transform);
      raising_op.add_block(complex<double>(0.0, 1.0), 0, 0, nspin1_, nspin1_, transform % *zfci2_orbang_[1] * transform);
    }

    complex<double> adjust = 1.0;
    for (int i = nspin1_ / 2 - 1; i >= 0; --i) {

      double phase_check = arg(raising_op.element(i, i+1));
      if (i + 1 == nspin1_ / 2 && nspin1_ % 2 == 0)
        phase_check *= 0.5;
      adjust *= polar(1.0, phase_check);

      //cout << " *** Eigenvector " << i << " time-reversal-allowed phased adjustment = " << arg(adjust) << endl;
      for (int j = 0; j != nspin1_; ++j) {
        transform.element(j, i) = adjust * transform.element(j, i);
        transform.element(j, nspin_ - i) = conj(adjust) * transform.element(j, nspin_ - i);
      }
    }
  }
  cout << endl;

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

  transform.print("Transformation matrix from Relativistic Full CI eigenstates to Pseudospin eigenfunction");

  cout << fixed << setprecision(5) << endl;
  cout << "    The z-axis is set to (";
  cout << setw(8) << spin_axes_->element(0, 2) << ", " << setw(8) << spin_axes_->element(1, 2) << ", " << setw(8) << spin_axes_->element(2, 2) << ")." << endl << endl;
  for (int i = 0; i != nspin1_; ++i)
    cout << "    " << to_upper(diagset) << " diagonal element " << i+1 << " = " << setw(12) << zeig[i] << endl;

  // We can no longer use this option, since I made this function const...
  //if (numerical_eig) {
  //  cout << "  **  By request, we compute the Hamiltonian using eigenvalues of S_z, rather than the canonical m_s values" << endl;
  //  update_spin_matrices(zeig);
  //}

  shared_ptr<ZMatrix> spinham_s = make_shared<ZMatrix>(transform % *spinham_h_ * transform);
  array<shared_ptr<ZMatrix>, 3> mu_s;
#ifdef LOCAL_DEBUG
  array<shared_ptr<ZMatrix>, 3> spin_s;
  array<shared_ptr<ZMatrix>, 3> orbang_s;
  auto trev_s = make_shared<ZMatrix>(transform % *trev_h_ * *transform.get_conjg());
#endif
  for (int i = 0; i != 3; ++i) {
    mu_s[i] = make_shared<ZMatrix>(transform % *zfci2_mu_[i] * transform);
#ifdef LOCAL_DEBUG
    spin_s[i] = make_shared<ZMatrix>(transform % *zfci2_spin_[i] * transform);
    orbang_s[i] = make_shared<ZMatrix>(transform % *zfci2_orbang_[i] * transform);
#endif
  }

  cout << endl;
  spinham_s->print("Pseudospin Hamiltonian", 24);
  cout << endl;

#ifdef LOCAL_DEBUG
  trev_s->print("Time-reversal operator in pseudospin states", 24);
  cout << endl;
  for (int i = 0; i != 3; ++i) {
    mu_s[i]->print("Magnetic moment matrix, component " + to_string(i), 24);
    spin_s[i]->print("Spin matrix, component " + to_string(i), 24);
    orbang_s[i]->print("Orbital moment matrix, component " + to_string(i), 24);
    cout << endl;
  }
  cout << endl;
#endif

  if (!is_t_symmetric(*spinham_s, /*hermitian*/true, /*time reversal*/true))
    throw runtime_error("The spin Hamiltonian seems to not have proper time-reversal symmetry.  Check that your spin value and states mapped are reasonable.");

  // Failures here can sometimes be fixed by using a tighter convergence threshold in the FCI part
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
    assert((*trev_s - *trev_target).rms() < 1.0e-10);
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
      out[i].set_coeff(real(stevop_vec.element(i, 0)));
      if (abs(imag(stevop_vec.element(i, 0))) > 1.0e-8)
        throw runtime_error("For some reason, we have obtained a complex coefficient for an extended Stevens operator...  It should be real.");
    }
  }

  cout << "    Stevens coefficients:  " << endl << endl;
  for (int i = 0; i != out.size(); ++i)
    cout << "    " << setw(8) << out[i].coeff_name() << " = " << setw(12) << out[i].coeff() << endl;
  cout << endl;

  const double checkham_error = (*checkham - *spinham_s).rms();
  VectorB shenergies(nspin1_);
  checkham->diagonalize(shenergies);

  vector<double> ord_energy = ref_energy_;
  sort(ord_energy.begin(), ord_energy.end());

  if (checkham_error > 1.0e-8) {
    cout << "  **** CAUTION ****  The pseudospin Hamiltonian does not fully reproduce the ab initio Hamiltonian.  RMS error = " << checkham_error << endl;

    cout << "  ** Relative energies from the pseudospin Hamiltonian: " << endl;
    for (int i = nspin_; i >= 0; --i)
      cout << "     " << i << "  " << setprecision(8) << setw(12) << shenergies[i] - shenergies[0] << " E_h  =  " << setprecision(4) << setw(8) << (shenergies[i] - shenergies[0])*au2wavenumber__ << " cm-1" << endl;
    cout << endl;

    cout << "  ** Relative energies from the ab initio Hamiltonian: " << endl;
    for (int i = nspin_; i >= 0; --i)
      cout << "     " << i << "  " << setprecision(8) << setw(12) << ord_energy[i] - ord_energy[0] << " E_h  =  " << setprecision(4) << setw(8) << (ord_energy[i] - ord_energy[0])*au2wavenumber__ << " cm-1" << endl;
    cout << endl;
  } else {

    cout << "  ** Relative energies: " << endl;
    for (int i = nspin_; i >= 0; --i) {
      cout << "     " << i << "  " << setprecision(8) << setw(12) << shenergies[i] - shenergies[0] << " E_h  =  " << setprecision(4) << setw(8) << (shenergies[i] - shenergies[0])*au2wavenumber__ << " cm-1" << endl;
      assert(abs(shenergies[i] - shenergies[0] - ord_energy[i] + ord_energy[0]) < 1.0e-7);
    }
    cout << endl;
  }

  return out;
}


tuple<shared_ptr<const Matrix>, double, double> Pseudospin::compute_dtensor(const vector<Stevens_Operator> input) const {
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
          throw logic_error("Some invalid operator was found in Pseudospin::compute_dtensor(...)");
      }
    }
  }

  /**** PRINTOUT ***/
  shared_ptr<Matrix> dtensor_diag = out->copy();
  dtensor_diag->print("D tensor");
  cout << setprecision(8);
  VectorB ddiag(3);
  dtensor_diag->diagonalize(ddiag);

  // Compute Davg so that it works even if D is not traceless (which shouldn't happen on accident)
  const double Davg = 1.0 / 3.0 * (ddiag[0] + ddiag[1] + ddiag[2]);

  int jmax = 0;
  const array<int,3> fwd = {{ 1, 2, 0 }};
  const array<int,3> bck = {{ 2, 0, 1 }};
  if (abs(ddiag[1]-Davg) > abs(ddiag[jmax]-Davg)) jmax = 1;
  if (abs(ddiag[2]-Davg) > abs(ddiag[jmax]-Davg)) jmax = 2;

  cout << endl << "    Upon diagonalization," << endl;
  cout << "      Dxx = " << setw(12) << ddiag[fwd[jmax]] << endl;
  cout << "      Dyy = " << setw(12) << ddiag[bck[jmax]] << endl;
  cout << "      Dzz = " << setw(12) << ddiag[jmax] << endl << endl;
  const double dval = ddiag[jmax] - 0.5*(ddiag[fwd[jmax]] + ddiag[bck[jmax]]);
  const double eval = abs(0.5*(ddiag[fwd[jmax]] - ddiag[bck[jmax]]));
  cout << " ** D = " << setw(12) << setprecision(8) << dval << " E_h = " << setprecision(4) << setw(8) << dval * au2wavenumber__ << " cm-1" << endl;
  cout << " ** E = " << setw(12) << setprecision(8) << eval << " E_h = " << setprecision(4) << setw(8) << eval * au2wavenumber__ << " cm-1" << endl;
  cout << " ** |E / D| = " << abs(eval / dval) << endl;

  //dtensor_diag->print("Transformation matrix of D-tensor");

  Matrix full_rotation = *spin_axes_ * *dtensor_diag;
  cout << fixed << setprecision(5) << endl;
  cout << endl << " ** Axis of principle D-value (relative to spin quant. axes) = (";
  cout << setw(8) << dtensor_diag->element(0, jmax) << ", " << setw(8) << dtensor_diag->element(1, jmax) << ", " << setw(8) << dtensor_diag->element(2, jmax) << ")" << endl;
  cout << endl << " ** Axis of principle D-value (relative to input geometry)  =  (";
  cout << setw(8) << full_rotation.element(0, jmax) << ", " << setw(8) << full_rotation.element(1, jmax) << ", " << setw(8) << full_rotation.element(2, jmax) << ")" << endl;

  tuple<shared_ptr<Matrix>, double, double> results(out, dval, eval);
  return results;
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

