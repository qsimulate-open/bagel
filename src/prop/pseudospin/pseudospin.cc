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

using namespace std;
using namespace bagel;


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


Pseudospin::Pseudospin(const int _nspin) : nspin_(_nspin), nspin1_(_nspin + 1) {

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
  const shared_ptr<const PTree> eso_ranks = zfci.idata()->get_child_optional("aniso_rank");
  if (eso_ranks) {
    ranks = {};
    for (auto& i : *eso_ranks)
      ranks.push_back(lexical_cast<int>(i->data()));
  }

  array<double, 3> rotin = zfci.idata()->get_array<double,3>("aniso_axis", array<double, 3>({{0.0, 0.0, 1.0}}));

  cout << setprecision(8);
  cout << endl << "    ********      " << endl;
  cout << endl << "    Modeling Pseudospin Hamiltonian for S = " << nspin_ / 2 << (nspin_ % 2 == 0 ? "" : " 1/2") << endl;

  vector<Stevens_Operator> ESO = build_extended_stevens_operators(ranks);

  compute_numerical_hamiltonian(zfci, zfci.jop()->coeff_input()->active_part());

  shared_ptr<ZMatrix> spinham_s = compute_spin_eigegenvalues(rotin);
  ESO = extract_hamiltonian_parameters(ESO, spinham_s);
  shared_ptr<Matrix> dtens = compute_Dtensor(ESO);
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

  // First, we create spin matrices in the atomic orbital basis
  RelSpinInt aospin(zfci.geom());

  const int norb = zfci.norb();

  // By default, just use the ground states
  vector<int> aniso_state;
  aniso_state.resize(nspin1_);
  ref_energy_.resize(nspin1_);
  for (int i = 0; i != nspin1_; ++i)
    aniso_state[i] = i;

  // aniso_state can be used to request mapping excited states instead
  const shared_ptr<const PTree> exstates = zfci.idata()->get_child_optional("aniso_state");
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
    spinop_h_[i] = make_shared<ZMatrix>(nspin1_, nspin1_);
  }
  for (int i = 0; i != nspin1_; ++i) {
    for (int j = 0; j != nspin1_; ++j) {
//      shared_ptr<Kramers<2,ZRDM<1>>> temprdm = (*rdm1)(aniso_state[i], aniso_state[j]);
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

      ZMatrix aodensity = (*active_coeff * modensity ^ *active_coeff);

      for (int k = 0; k != 3; ++k)
        spinop_h_[k]->element(i,j) = aodensity.dot_product(*aospin(k));
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


shared_ptr<ZMatrix> Pseudospin::compute_spin_eigegenvalues(const array<double, 3> rotation) const {

  // Diagonalize S_z to get pseudospin eigenstates as combinations of ZFCI Hamiltonian eigenstates
  auto transform = make_shared<ZMatrix>(nspin1_, nspin1_);
  const complex<double> scale = 1.0 / std::sqrt(rotation[0]*rotation[0] + rotation[1]*rotation[1] + rotation[2]*rotation[2]);
  assert(std::abs(std::imag(scale)) < 1.0e-10);
  for (int i = 0; i != 3; ++i)
    *transform += scale * rotation[i] * *spinop_h_[i];
  VectorB zeig(nspin1_);
  transform->diagonalize(zeig);

  { // Reorder eigenvectors so positive M_s come first
    shared_ptr<ZMatrix> tempm = transform->clone();
    VectorB tempv = *zeig.clone();
    for (int i = 0; i != nspin1_; ++i) {
      tempv[i] = zeig[nspin_ - i];
      tempm->copy_block(0, i, nspin1_, 1, transform->slice(nspin_ - i, nspin_ - i + 1));
    }
    transform = tempm;
    zeig = tempv;
  }

  { // Adjust the phases of eigenvectors to ensure proper time-reversal symmetry
    ZMatrix spinham_s = *transform % *spinham_h_ * *transform;
    complex<double> adjust = 1.0;
    for (int k = nspin_ / 2; k > 0; --k) {
      const double phase_error = std::arg(spinham_s.element(nspin_ - k, nspin1_ - k)) - std::arg(spinham_s.element(k - 1, k)) - pi__;
      adjust *= std::polar(1.0, -1.0 * phase_error);
      for (int i = 0; i != nspin1_; ++i)
        transform->element(i, nspin1_ - k) = adjust * transform->element(i, nspin1_ - k);
    }
  }

  { // Check a spin matrix as well with basically the same procedure, since sometimes we miss a phase due to numerically zero entries in the Hamiltonian
    const ZMatrix spinop_x = *transform % *spinop_h_[0] * *transform;
    if (!is_t_symmetric(spinop_x, /*hermitian*/ true, /*t_symmetric*/ false, 1.0e-8)) {
      complex<double> adjust = 1.0;
      for (int k = nspin_ / 2; k > 0; --k) {
        const double phase_error = std::arg(spinop_x.element(nspin_ - k, nspin1_ - k)) - std::arg(spinop_x.element(k - 1, k));
        assert(std::abs(phase_error) < 1.0e-2); // should be a small correction
        adjust *= std::polar(1.0, -1.0 * phase_error);
        for (int i = 0; i != nspin1_; ++i)
          transform->element(i, nspin1_ - k) = adjust * transform->element(i, nspin1_ - k);
      }
    }
  }

  cout << "    The z-axis is set to (" << rotation[0] << ", " << rotation[1] << ", " << rotation[2] << ")." << endl << endl;
  for (int i = 0; i != nspin1_; ++i)
    cout << "    Pseudospin eigenvalue " << i+1 << " = " << setw(12) << zeig[i] << endl;

  // We can no longer use this option, since I made this function const...
  //if (numerical_eig) {
  //  cout << "  **  By request, we compute the Hamiltonian using eigenvalues of S_z, rather than the canonical m_s values" << endl;
  //  update_spin_matrices(zeig);
  //}

  shared_ptr<ZMatrix> spinham_s = make_shared<ZMatrix>(*transform % *spinham_h_ * *transform);
  if (!is_t_symmetric(*spinham_s, /*hermitian*/true, /*time reversal*/true))
    throw runtime_error("The spin Hamiltonian seems to not have proper time-reversal symmetry.  Check that your spin value and states mapped are reasonable.");

  array<shared_ptr<ZMatrix>, 3> spinop_s;
  for (int i = 0; i != 3; ++i) {
    spinop_s[i] = make_shared<ZMatrix>(*transform % *spinop_h_[i] * *transform);
    assert(is_t_symmetric(*spinop_s[i], /*hermitian*/true, /*time reversal*/false));
  }

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

  return out;
}


bool Pseudospin::is_t_symmetric(const ZMatrix& in, const bool hermitian, const bool t_symmetric, const double thresh) const {
  // The matrix must be either Hermitian or skew-Hermitian
  if (hermitian) {
    assert(in.is_hermitian());
  } else {
    assert((in + *in.transpose_conjg()).rms() < 1.0e-8);
  }

  assert(in.ndim() == in.mdim());
  assert(in.ndim() == nspin1_);
  const double h = hermitian ? 1.0 : -1.0;
  const double t = t_symmetric ? 1.0 : -1.0;

  // We explicitly check the upper-triangular part element by element
  bool out = true;
  for (int i = 0; i != nspin1_; ++i) {
    for (int j = i; j != nspin1_; ++j) {
      // fac is +/- 1.0 depending on how close we are to the diagonal
      const int evendiag = (j - i) % 2;
      const double fac = evendiag ? -1.0 : 1.0;

      const complex<double> val = in.element(i, j) - t * h * fac * in.element(nspin_ - j, nspin_ - i);
      if (std::abs(val) > thresh)
        out = false;
    }
  }

#ifndef NDEBUG
  { // Ensure that the element-by-element check for observed pattern agrees with the actual time-reversal operator
    // The time-reversal operator is kramersop + complex conjugation
    ZMatrix kramersop(nspin1_, nspin1_);
    double val = -1.0;
    for (int i = 0; i != nspin1_; ++i) {
      kramersop.element(nspin_ - i, i) = val;
      val *= -1.0;
    }
    ZMatrix in_reversed = kramersop * *in.get_conjg() * *kramersop.transpose_conjg();  // K O K^-1
    const double error = (in - (t * in_reversed)).rms();
    assert((out == true && error < thresh) || (out == false && error >= thresh));
  }
#endif

  return out;
}

