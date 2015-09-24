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


Spin_Operator::Spin_Operator(shared_ptr<const ZMatrix> _mat, const bool _stev, const int _ord, const int _ind)
 : nspin_(_mat->ndim() - 1), matrix_(_mat), stevens_(_stev), order_(_ord), index_(_ind) {
  coeff_ = nan("");
  assert(matrix_->ndim() == matrix_->mdim());                                                         // Matrix representation must be square
  assert((stevens_ && abs(index_) <= order_) || (!stevens_ && index_ >= 0 && index_ <= 8));           // Validity of index_
  assert((stevens_ && order_ >= 0) || (!stevens_ && order_ == 2));                                    // Validity of order_
  assert((order_ <= nspin_ || matrix_->rms() < 1.0e-8));                                              // High-order contributions should be zero for low spin
}

string Spin_Operator::operator_name() const {
  string out = "";
  if (stevens_) {
    out += "O_" + to_string(order_) + "^" + to_string(index_);
    if (index_ >= 0) out += " ";
  } else {
    const array<string, 3> dim = {{ "x", "y", "z" }};
    out += "S" + dim[index_ % 3] + "S" + dim[index_ / 3];
  }
  return out;
}

string Spin_Operator::coeff_name() const {
  string out = operator_name();
  if (stevens_) {
    out[0] = 'B';
  } else {
    out[0] = 'D';
    out.erase(2, 1);
  }
  return out;
}


Pseudospin::Pseudospin(const int _nspin) : nspin_(_nspin), nspin1_(_nspin + 1) {

  VectorB spinvals(nspin1_);
  for (int i = 0; i != nspin1_; ++i)
    spinvals[i] = (nspin_ / 2.0) - i;
  update_spin_matrices(spinvals);
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


vector<Spin_Operator> Pseudospin::build_2ndorder_zfs_operators() const {
  vector<Spin_Operator> out;
  for (int j = 0; j != 3; ++j) {
    for (int i = 0; i != 3; ++i) {
      auto mat = make_shared<ZMatrix>(*spin_xyz(i) * *spin_xyz(j));
      Spin_Operator tmp(mat, false, 2, 3 * j + i);
      out.push_back(tmp);
    }
  }
  return out;
}


// Compute numerical pseudospin Hamiltonian by diagonalizing S_z matrix
void Pseudospin::compute_numerical_hamiltonian(const ZHarrison& zfci, shared_ptr<const RelCoeff_Block> active_coeff) {

  // First, we create spin matrices in the atomic orbital basis
  RelSpinInt aospin(zfci.geom());

  const int norb = zfci.norb();

  // S value of spin manifold to be mapped
  cout << endl << endl;
  cout << "    Modeling Pseudospin Hamiltonian for S = " << nspin_ / 2 << (nspin_ % 2 == 0 ? "" : " 1/2") << endl;

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

      /*******/
      // Test i j symmetry
      if (i > j && zfci.idata()->get<bool>("aniso_test_symm", false)) {
        temprdm = zfci.rdm1(aniso_state[j], aniso_state[i]);
        if (!temprdm->exist({1,0})) {
          cout << " * Need to generate an off-diagonal rdm of zeroes." << endl;
          temprdm->add({1,0}, temprdm->at({0,0})->clone());
        }
        tmp = expand_kramers<1,complex<double>>(temprdm, norb);
        rdmmat->zero();
        copy_n(tmp->data(), tmp->size(), rdmmat->data());
        rdmmat = rdmmat->transpose_conjg();
      }
      /*******/

      ZMatrix modensity (2 * norb, 2 * norb);
      modensity.copy_block(0, 0, 2 * norb, 2 * norb, rdmmat);

      ZMatrix aodensity = (*active_coeff * modensity ^ *active_coeff);

      for (int k = 0; k != 3; ++k)
        spinop_h_[k]->element(i,j) = aodensity.dot_product(*aospin(k));
    }
  }

  for (int i = 0; i != 3; ++i)
    spinop_h_[i]->print("Spin matrix, for component " + to_string(i) + " over ZFCI states");

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


shared_ptr<ZMatrix> Pseudospin::compute_spin_eigegenvalues(const bool symmetrize, const array<complex<double>, 3> rotation) const {

  // Diagonalize S_z to get pseudospin eigenstates as combinations of ZFCI Hamiltonian eigenstates
  auto transform = make_shared<ZMatrix>(nspin1_, nspin1_);
  const complex<double> scale = 1.0 / std::sqrt(std::conj(rotation[0])*rotation[0] + std::conj(rotation[1])*rotation[1] + std::conj(rotation[2])*rotation[2]);
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

  for (int i = 0; i != nspin1_; ++i)
    cout << "    Spin-z eigenvalue " << i+1 << " = " << zeig[i] << endl;

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

  cout << endl;

  cout << endl;
  spinham_s->print("Pseudospin Hamiltonian!");

  spinop_s[0]->print("Spin matrix - x-component");
  spinop_s[1]->print("Spin matrix - y-component");
  spinop_s[2]->print("Spin matrix - z-component");


  (*spinop_s[0] * *spinop_s[0]).print("Spin matrix, x^2");
  (*spinop_s[1] * *spinop_s[1]).print("Spin matrix, y^2");
  (*spinop_s[2] * *spinop_s[2]).print("Spin matrix, z^2");

  /************/
  // To average out broken symmetry and obtain a consistent set of linear equations
  if (symmetrize) {
    if (nspin_ == 3) {
      cout << "  **  By request, we will average out time-reversal asymmetry for an S = 3/2 Hamiltonian to ensure a consistent set of equations" << endl;
      const complex<double> offdiag1 = 0.25 * (spinham_s->element(0,1) + std::conj(spinham_s->element(1,0)) - spinham_s->element(2,3) - std::conj(spinham_s->element(3,2)));
      const complex<double> offdiag2 = 0.25 * (spinham_s->element(0,2) + std::conj(spinham_s->element(2,0)) + spinham_s->element(1,3) + std::conj(spinham_s->element(3,1)));
      spinham_s->element(0,1) = offdiag1;
      spinham_s->element(1,0) = std::conj(offdiag1);
      spinham_s->element(2,3) = -offdiag1;
      spinham_s->element(3,2) = -std::conj(offdiag1);
      spinham_s->element(0,2) = offdiag2;
      spinham_s->element(1,3) = offdiag2;
      spinham_s->element(2,0) = std::conj(offdiag2);
      spinham_s->element(3,1) = std::conj(offdiag2);
      spinham_s->print("Pseudospin Hamiltonian!...  Forcibly symmetrized");
    }
  }

  return spinham_s;
}


// Extract D-tensor or Stevens coefficients from the numerical pseudospin Hamiltonian
vector<Spin_Operator> Pseudospin::extract_hamiltonian_parameters(const bool real, const vector<Spin_Operator> param, shared_ptr<const ZMatrix> spinham_s) const {
  const int nop = param.size();
  vector<Spin_Operator> out = param;

  // d2h is the transformation matrix to build pseudospin Hamiltonian from D-tensor or Extended Stevens Operators
  // Rows correspond to pairs of pseudospins (SS, S-1S, S-2S...)
  // Columns correspond to spin Hamiltonian parameters (e.g., Dxx, Dyx, Dzx, Dxy...)
  // Note that we look over the first indices before the second, so we can copy data between vectors and matrices and have it come out in the right order
  auto d2h = make_shared<ZMatrix>(nspin1_ * nspin1_, nop);

  for (int i = 0; i != nop; ++i)
    d2h->copy_block(0, i, nspin1_ * nspin1_, 1, param[i].matrix()->data());

  auto Dtensor = make_shared<ZMatrix>(3,3);
  auto checkham = make_shared<ZMatrix>(nspin1_, nspin1_);


  // Convert from the pseudospin Hamiltonian to the D-tensor using the left-inverse of d2h
  if (real) {
    // By default, force the D tensor to be real

    // Separate out real and imaginary parts
    auto d2h_real = make_shared<Matrix>(nspin1_ * nspin1_ * 2, nop);
    auto spinham_vec_real = make_shared<Matrix>(nspin1_ * nspin1_ * 2, 1);
    d2h_real->copy_block(            0, 0, nspin1_ * nspin1_, nop, d2h->get_real_part());
    d2h_real->copy_block(nspin1_ * nspin1_, 0, nspin1_ * nspin1_, nop, d2h->get_imag_part());
    spinham_vec_real->copy_block(            0, 0, nspin1_ * nspin1_, 1, spinham_s->get_real_part()->element_ptr(0,0));
    spinham_vec_real->copy_block(nspin1_ * nspin1_, 0, nspin1_ * nspin1_, 1, spinham_s->get_imag_part()->element_ptr(0,0));

    // Compute left-inverse as  (A^T A)^-1 A^T
    Matrix d2h_sqinv_real = *d2h_real % *d2h_real;
    d2h_sqinv_real.inverse();
    Matrix h2d_real = d2h_sqinv_real ^ *d2h_real;
    assert((h2d_real * *d2h_real).is_identity());

    // Extract D-tensor from it
    Matrix Dtensor_vec_real = h2d_real * *spinham_vec_real;
/*
    Matrix Dtensor_real(3, 3);
    Dtensor_real.copy_block(0, 0, 3, 3, Dtensor_vec_real.element_ptr(0,0));

    Dtensor->copy_real_block(1.0, 0, 0, 3, 3, Dtensor_real);
*/

    // Recompute Hamiltonian from D so we can check the fit
    ZMatrix Dtensor_vec(nop, 1);
    //Dtensor_vec.copy_block(0, 0, nop, 1, Dtensor->element_ptr(0,0));
    Dtensor_vec.add_real_block(1.0, 0, 0, nop, 1, Dtensor_vec_real);
    ZMatrix checkham_vec = *d2h * Dtensor_vec;
    checkham->copy_block(0, 0, nspin1_, nspin1_, checkham_vec.element_ptr(0,0));

    for (int i = 0; i != nop; ++i) {
      cout << "  Pseudospin Hamiltonian parameter: " << setw(7) << param[i].coeff_name() << " = " << Dtensor_vec.element(i, 0) << endl;
      out[i].set_coeff(Dtensor_vec.element(i, 0));
    }

  } else {
    // On request, allow complex ZFS parameters
    // Same algorithm, working directly with complex matrices
    auto h2d = make_shared<ZMatrix>(nop, nspin1_ * nspin1_);
    ZMatrix d2h_sqinv = *d2h % *d2h;
    d2h_sqinv.inverse();
    *h2d = d2h_sqinv ^ *d2h;
    assert((*h2d * *d2h).is_identity());

    auto spinham_vec = make_shared<ZMatrix>(nspin1_ * nspin1_,1);
    spinham_vec->copy_block(0, 0, nspin1_ * nspin1_, 1, spinham_s->element_ptr(0,0));
    ZMatrix Dtensor_vec = *h2d * *spinham_vec;
//    Dtensor->copy_block(0, 0, 3, 3, Dtensor_vec.element_ptr(0,0));

    ZMatrix check_spinham_vec = *d2h * Dtensor_vec;
    checkham->copy_block(0, 0, nspin1_, nspin1_, check_spinham_vec.element_ptr(0,0));

    for (int i = 0; i != nop; ++i) {
      cout << "  Pseudospin Hamiltonian parameter: " << setw(7) << param[i].coeff_name() << " = " << Dtensor_vec.element(i, 0) << endl;
      out[i].set_coeff(Dtensor_vec.element(i, 0));
    }
  }

  checkham->print("Pseudospin Hamiltonian, recomputed", 30);
  cout << "  Error in recomputation of spin Hamiltonian from D = " << (*checkham - *spinham_s).rms() << endl << endl;

  VectorB shenergies(nspin1_);
  checkham->diagonalize(shenergies);

  cout << "  ** Relative energies expected from the recomputed Pseudospin Hamiltonian: " << endl;
  for (int i = nspin_; i >= 0; --i)
    cout << "     " << i << "  " << shenergies[i] - shenergies[0] << " E_h  =  " << (shenergies[i] - shenergies[0])*au2wavenumber__ << " cm-1" << endl;
  cout << endl;

  cout << "  ** Relative energies observed by relativistic configuration interaction: " << endl;
  for (int i = nspin_; i >= 0; --i)
    cout << "     " << i << "  " << ref_energy_[i] - ref_energy_[0] << " E_h  =  " << (ref_energy_[i] - ref_energy_[0])*au2wavenumber__ << " cm-1" << endl;
  cout << endl;

  return out;
}


// Working with complex algebra, although it should be fully real...
shared_ptr<ZMatrix> Pseudospin::compute_Dtensor(const vector<Spin_Operator> input) {
  auto out = make_shared<ZMatrix>(3, 3);
  out->zero();

  if (input[0].stevens()) {
    // Get D from second-order extended Stevens Operators
    for (int i = 0; i != input.size(); ++i) {
      assert(input[i].stevens());
      if (input[i].order() == 2) {
        switch (input[i].index()) {
          case  0:
            assert(std::abs(std::imag(input[i].coeff())) < 1.0e-8);
            out->element(0, 0) -= 1.0 * input[i].coeff();
            out->element(1, 1) -= 1.0 * input[i].coeff();
            out->element(2, 2) += 2.0 * input[i].coeff();
            break;
          case -1:
            out->element(1, 2) += 0.5 * input[i].coeff();
            out->element(2, 1) += 0.5 * std::conj(input[i].coeff());
            break;
          case  1:
            out->element(2, 0) += 0.5 * input[i].coeff();
            out->element(0, 2) += 0.5 * std::conj(input[i].coeff());
            break;
          case -2:
            out->element(0, 1) += 1.0 * input[i].coeff();
            out->element(1, 0) += 1.0 * std::conj(input[i].coeff());
            break;
          case  2:
            assert(std::abs(std::imag(input[i].coeff())) < 1.0e-8);
            out->element(0, 0) += 1.0 * input[i].coeff();
            out->element(1, 1) -= 1.0 * input[i].coeff();
            break;
          default:
            throw logic_error("Some invalid operator was found in Pseudospin::compute_Dtensor(...)");
        }
      }
    }
  } else {
    // Just have to copy the data if we have fit the D-tensor directly
    assert(input.size() == 9);
    for (int i = 0; i != 9; ++i) {
      assert(!input[i].stevens() && input[i].order() == 2);
      out->element(i % 3, i / 3) = input[i].coeff();
    }
  }

  /**** PRINTOUT ***/
  shared_ptr<ZMatrix> Dtensor_diag = out->copy();
  //Dtensor_diag->print("D tensor");
  VectorB Ddiag(3);
  Dtensor_diag->diagonalize(Ddiag);
  for (int i = 0; i != 3; ++i)
    cout << "Diagonalized D-tensor value " << i << " = " << Ddiag[i] << endl;

  // Compute Davg so that it works even if D is not traceless (which shouldn't happen on accident)
  const double Davg = 1.0 / 3.0 * (Ddiag[0] + Ddiag[1] + Ddiag[2]);

  int jmax = 0;
  const array<int,3> fwd = {{ 1, 2, 0 }};
  const array<int,3> bck = {{ 2, 0, 1 }};
  if (std::abs(Ddiag[1]-Davg) > std::abs(Ddiag[jmax]-Davg)) jmax = 1;
  if (std::abs(Ddiag[2]-Davg) > std::abs(Ddiag[jmax]-Davg)) jmax = 2;
  const double Dval = Ddiag[jmax] - 0.5*(Ddiag[fwd[jmax]] + Ddiag[bck[jmax]]);
  const double Eval = 0.5*(Ddiag[fwd[jmax]] - Ddiag[bck[jmax]]);
  cout << " ** D = " << Dval << " E_h, or " << Dval * au2wavenumber__ << " cm-1" << endl;
  cout << " ** E = " << std::abs(Eval) << " E_h, or " << std::abs(Eval * au2wavenumber__) << " cm-1" << endl;
  cout << " ** E / D = " << std::abs(Eval / Dval) << endl;

  if (input[0].nspin() == 2) {
    cout << "  ** Relative energies expected from diagonalized D parameters: " << endl;
    if (Dval > 0.0) {
      cout << "     2  " << Dval + std::abs(Eval) << " E_h  =  " << (Dval + std::abs(Eval))*au2wavenumber__ << " cm-1" << endl;
      cout << "     1  " << Dval - std::abs(Eval) << " E_h  =  " << (Dval - std::abs(Eval))*au2wavenumber__ << " cm-1" << endl;
      cout << "     0  " << 0.0 << " E_h  =  " << 0.0 << " cm-1" << endl << endl;
    } else {
      cout << "     2  " << -Dval + 0.5*std::abs(Eval) << " E_h  =  " << (-Dval + 0.5*std::abs(Eval))*au2wavenumber__ << " cm-1" << endl;
      cout << "     1  " << std::abs(Eval) << " E_h  =  " << std::abs(Eval)*au2wavenumber__ << " cm-1" << endl;
      cout << "     0  " << 0.0 << " E_h  =  " << 0.0 << " cm-1" << endl << endl;
    }
  } else if (input[0].nspin() == 3) {
    cout << "  ** Relative energies expected from diagonalized D parameters: " << endl;
    const double energy32 = 2.0*std::sqrt(Dval*Dval + 3.0*Eval*Eval);
    cout << "     3  " << energy32 << " E_h  =  " << energy32*au2wavenumber__ << " cm-1" << endl;
    cout << "     2  " << energy32 << " E_h  =  " << energy32*au2wavenumber__ << " cm-1" << endl;
    cout << "     1  " << 0.0 << " E_h  =  " << 0.0 << " cm-1" << endl;
    cout << "     0  " << 0.0 << " E_h  =  " << 0.0 << " cm-1" << endl << endl;
  }

  return out;
}


bool Pseudospin::is_t_symmetric(const ZMatrix& in, const bool hermitian, const bool t_symmetric, const double thresh) const {
  in.print("Checking t-symmetry of this matrix; hermitian = " + to_string(hermitian) + ", t_symmetric = " + to_string(t_symmetric));
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
  return out;
}

