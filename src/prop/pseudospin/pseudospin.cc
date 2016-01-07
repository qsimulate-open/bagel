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

// for debugging:
#include <src/util/math/comb.h>

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

  cout << setprecision(8);
  cout << endl << "    ********      " << endl;
  cout << endl << "    Modeling Pseudospin Hamiltonian for S = " << nspin_ / 2 << (nspin_ % 2 == 0 ? "" : " 1/2") << endl;

  vector<Stevens_Operator> ESO = build_extended_stevens_operators(ranks);

#if 0
  cout << "Number of Stevens operators = " << ESO.size() << endl;
  for (int i = 0; i != ESO.size(); ++i)
    ESO[i].print();
#endif

  compute_numerical_hamiltonian(zfci, zfci.jop()->coeff_input()->active_part());

  shared_ptr<const Matrix> mag_axes = identify_magnetic_axes();
  array<double,3> rotation;
    for (int i = 0; i != 3; ++i)
      rotation[i] = mag_axes->element(2, i);

  array<double, 3> rotin = zfci.idata()->get_array<double,3>("aniso_axis", rotation);

  shared_ptr<const ZMatrix> spinham_s = compute_spin_eigenvalues(rotin, zfci);

  if (nspin_ > 1) {
    ESO = extract_hamiltonian_parameters(ESO, spinham_s);
    shared_ptr<Matrix> dtens = compute_Dtensor(ESO);
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
  auto ao_trev = make_shared<const RelTRevInt>(zfci.geom());
  auto mo_trev = make_shared<ZMatrix>(*active_coeff % *ao_trev * *active_coeff->get_conjg());

# if 1
  { // Attempt to compute time-reversal matrix in ZFCI states from mo_trev
    trev_h_ = make_shared<ZMatrix>(nspin1_, nspin1_);
    trev_h_->zero();
    const int maxa = std::min(zfci.nele(), zfci.norb());
    const int mina = std::max(zfci.nele() - zfci.norb(), 0);

    /********************/
    // MAKE IT RETURN THE SPIN-Y Matrix rather than K, as a debug check
    if (zfci.idata()->get<bool>("aniso_testconj", false)) {
      mo_trev = make_shared<ZMatrix>(*active_coeff % *ao_spin(1) * *active_coeff);
    }
    const bool conjop = !zfci.idata()->get<bool>("aniso_testconj", false);
    /********************/

    /********************/
    // Redundant, but need this info to be available
    vector<int> aniso_state;
    aniso_state.resize(nspin1_);
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
    }
    /********************/

    vector<array<int,2>> ab = {};
    for (int j = maxa; j >= mina; --j)
      ab.push_back({{j, zfci.nele()-j}});

    // Loop over pairs of spin sectors
    for (int k = 0; k != ab.size(); ++k) {
      for (int l = 0; l != ab.size(); ++l) {
        auto det_i = zfci.cc()->find(ab[k][0], ab[k][1])->data(0)->det();
        auto det_j = zfci.cc()->find(ab[l][0], ab[l][1])->data(0)->det();

        // Loop over pairs of determinants
        for (auto& ia : det_i->string_bits_a()) {
          for (auto& ib : det_i->string_bits_b()) {
            const int pos1_i = det_i->lexical<0>(ia);
            const int pos2_i = det_i->lexical<1>(ib);
            const int fullpos_i = pos1_i*det_i->lenb() + pos2_i;

            for (auto& ja : det_j->string_bits_a()) {
              for (auto& jb : det_j->string_bits_b()) {
                const int pos1_j = det_j->lexical<0>(ja);
                const int pos2_j = det_j->lexical<1>(jb);
                const int fullpos_j = pos1_j*det_j->lenb() + pos2_j;

                // Loop over pairs of ZFCI eigenstates
                for (int i = 0; i != nspin1_; ++i) {
                  for (int j = 0; j != nspin1_; ++j) {

                    auto civec_i = zfci.cc()->find(ab[k][0], ab[k][1])->data(aniso_state[i]);
                    auto civec_j = zfci.cc()->find(ab[l][0], ab[l][1])->data(aniso_state[j]);

                    // 1 electron version - Loop over pairs of MOs
                    if (zfci.nele() == 1) {
                      for (int Ap = 0; Ap != 2*norb; ++Ap) {
                        const int A = Ap % norb;
                        if ((ia[A] && Ap < norb) || (ib[A] && Ap >= norb)) {
                          for (int Bp = 0; Bp != 2*norb; ++Bp) {
                            const int B = Bp % norb;
                            if ((ja[B] && Bp < norb) || (jb[B] && Bp >= norb)) {
                              if (conjop) {
                                trev_h_->element(i, j) += std::conj(civec_i->data(fullpos_i)) * std::conj(civec_j->data(fullpos_j)) * mo_trev->element(Ap, Bp);
                              } else {
                                trev_h_->element(i, j) += std::conj(civec_i->data(fullpos_i)) * civec_j->data(fullpos_j) * mo_trev->element(Ap, Bp);
                              }
                            }
                          }
                        }
                      }
                    }

                    // 2 electron version - Loop over sets of 4 MOs
                    if (zfci.nele() == 2) {
                      for (int Ap = 0; Ap != 2*norb; ++Ap) {
                        const int A = Ap % norb;
                        if ((ia[A] && Ap < norb) || (ib[A] && Ap >= norb)) {
                          for (int Bp = 0; Bp != 2*norb; ++Bp) {
                            if (Ap == Bp) continue;
                            const int B = Bp % norb;
                            if ((ia[B] && Bp < norb) || (ib[B] && Bp >= norb)) {
                              for (int Cp = 0; Cp != 2*norb; ++Cp) {
                                const int C = Cp % norb;
                                if ((ja[C] && Cp < norb) || (jb[C] && Cp >= norb)) {
                                  for (int Dp = 0; Dp != 2*norb; ++Dp) {
                                    if (Cp == Dp) continue;
                                    const int D = Dp % norb;
                                    if ((ja[D] && Dp < norb) || (jb[D] && Dp >= norb)) {
                                      const double sign1 = Ap < Bp ? 1.0 : -1.0;
                                      const double sign2 = Cp < Dp ? 1.0 : -1.0;
                                      const double sign3 = sign1 * sign2 / 2.0;
                                      if (conjop) {
                                        trev_h_->element(i, j) += sign3 * std::conj(civec_i->data(fullpos_i)) * std::conj(civec_j->data(fullpos_j)) * mo_trev->element(Bp, Dp) * mo_trev->element(Ap, Cp);
                                      } else {
                                        if (Ap == Cp)
                                          trev_h_->element(i, j) += sign3 * std::conj(civec_i->data(fullpos_i)) * civec_j->data(fullpos_j) * mo_trev->element(Bp, Dp);
                                        if (Bp == Dp)
                                          trev_h_->element(i, j) += sign3 * std::conj(civec_i->data(fullpos_i)) * civec_j->data(fullpos_j) * mo_trev->element(Ap, Cp);
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }

                    // 3 electron version - Loop over sets of 6 MOs
                    if (zfci.nele() == 3) {
                      for (int Ap = 0; Ap != 2*norb; ++Ap) {
                        const int A = Ap % norb;
                        if ((ia[A] && Ap < norb) || (ib[A] && Ap >= norb)) {

                          for (int Bp = 0; Bp != 2*norb; ++Bp) {
                            if (Ap == Bp) continue;
                            const int B = Bp % norb;
                            if ((ia[B] && Bp < norb) || (ib[B] && Bp >= norb)) {

                              for (int Cp = 0; Cp != 2*norb; ++Cp) {
                                if ((Ap == Cp) || (Bp == Cp)) continue;
                                const int C = Cp % norb;
                                if ((ia[C] && Cp < norb) || (ib[C] && Cp >= norb)) {

                                  for (int Dp = 0; Dp != 2*norb; ++Dp) {
                                    const int D = Dp % norb;
                                    if ((ja[D] && Dp < norb) || (jb[D] && Dp >= norb)) {

                                      for (int Ep = 0; Ep != 2*norb; ++Ep) {
                                        if (Dp == Ep) continue;
                                        const int E = Ep % norb;
                                        if ((ja[E] && Ep < norb) || (jb[E] && Ep >= norb)) {

                                          for (int Fp = 0; Fp != 2*norb; ++Fp) {
                                            if ((Dp == Fp) || (Ep == Fp)) continue;
                                            const int F = Fp % norb;
                                            if ((ja[F] && Fp < norb) || (jb[F] && Fp >= norb)) {

                                              const double sign1 = Ap < Bp ? 1.0 : -1.0;
                                              const double sign2 = Ap < Cp ? 1.0 : -1.0;
                                              const double sign3 = Bp < Cp ? 1.0 : -1.0;
                                              const double sign4 = Dp < Ep ? 1.0 : -1.0;
                                              const double sign5 = Dp < Fp ? 1.0 : -1.0;
                                              const double sign6 = Ep < Fp ? 1.0 : -1.0;
                                              const double sign7 = sign1 * sign2 * sign3 * sign4 * sign5 * sign6 / 6.0;

                                              if (conjop) {
                                                trev_h_->element(i, j) += sign7 * std::conj(civec_i->data(fullpos_i)) * std::conj(civec_j->data(fullpos_j)) * mo_trev->element(Cp, Fp) * mo_trev->element(Ap, Dp) * mo_trev->element(Bp, Ep);
                                              } else {
                                                if ((Ap == Dp) && (Bp == Ep))
                                                  trev_h_->element(i, j) += sign7 * std::conj(civec_i->data(fullpos_i)) * civec_j->data(fullpos_j) * mo_trev->element(Cp, Fp);
                                                if ((Bp == Ep) && (Cp == Fp))
                                                  trev_h_->element(i, j) += sign7 * std::conj(civec_i->data(fullpos_i)) * civec_j->data(fullpos_j) * mo_trev->element(Ap, Dp);
                                                if ((Ap == Dp) && (Cp == Fp))
                                                  trev_h_->element(i, j) += sign7 * std::conj(civec_i->data(fullpos_i)) * civec_j->data(fullpos_j) * mo_trev->element(Bp, Ep);
                                              }
                                            }
                                          }
                                        }
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }

                    // 4 electron version - Loop over sets of 8 MOs
                    if (zfci.nele() == 4) {
                      for (int Ap = 0; Ap != 2*norb; ++Ap) {
                        const int A = Ap % norb;
                        if ((ia[A] && Ap < norb) || (ib[A] && Ap >= norb)) {

                          for (int Bp = 0; Bp != 2*norb; ++Bp) {
                            if (Ap == Bp) continue;
                            const int B = Bp % norb;
                            if ((ia[B] && Bp < norb) || (ib[B] && Bp >= norb)) {

                              for (int Cp = 0; Cp != 2*norb; ++Cp) {
                                if ((Ap == Cp) || (Bp == Cp)) continue;
                                const int C = Cp % norb;
                                if ((ia[C] && Cp < norb) || (ib[C] && Cp >= norb)) {

                                  for (int Dp = 0; Dp != 2*norb; ++Dp) {
                                    if ((Ap == Dp) || (Bp == Dp) || (Cp == Dp)) continue;
                                    const int D = Dp % norb;
                                    if ((ia[D] && Dp < norb) || (ib[D] && Dp >= norb)) {

                                      for (int Ep = 0; Ep != 2*norb; ++Ep) {
                                        const int E = Ep % norb;
                                        if ((ja[E] && Ep < norb) || (jb[E] && Ep >= norb)) {

                                          for (int Fp = 0; Fp != 2*norb; ++Fp) {
                                            if (Ep == Fp) continue;
                                            const int F = Fp % norb;
                                            if ((ja[F] && Fp < norb) || (jb[F] && Fp >= norb)) {

                                              for (int Gp = 0; Gp != 2*norb; ++Gp) {
                                                if ((Ep == Gp) || (Fp == Gp)) continue;
                                                const int G = Gp % norb;
                                                if ((ja[G] && Gp < norb) || (jb[G] && Gp >= norb)) {

                                                  for (int Hp = 0; Hp != 2*norb; ++Hp) {
                                                    if ((Ep == Hp) || (Fp == Hp) || (Gp == Hp)) continue;
                                                    const int H = Hp % norb;
                                                    if ((ja[H] && Hp < norb) || (jb[H] && Hp >= norb)) {

                                                      const double sign1 = Ap < Bp ? 1.0 : -1.0;
                                                      const double sign2 = Ap < Cp ? 1.0 : -1.0;
                                                      const double sign3 = Ap < Dp ? 1.0 : -1.0;
                                                      const double sign4 = Bp < Cp ? 1.0 : -1.0;
                                                      const double sign5 = Bp < Dp ? 1.0 : -1.0;
                                                      const double sign6 = Cp < Dp ? 1.0 : -1.0;
                                                      const double sign7 = Ep < Fp ? 1.0 : -1.0;
                                                      const double sign8 = Ep < Gp ? 1.0 : -1.0;
                                                      const double sign9 = Ep < Hp ? 1.0 : -1.0;
                                                      const double sign10 = Fp < Gp ? 1.0 : -1.0;
                                                      const double sign11 = Fp < Hp ? 1.0 : -1.0;
                                                      const double sign12 = Gp < Hp ? 1.0 : -1.0;
                                                      const double sign13 = sign1 * sign2 * sign3 * sign4 * sign5 * sign6 * sign7 * sign8 * sign9 * sign10 * sign11 * sign12 / 24.0;

                                                      if (conjop) {
                                                        trev_h_->element(i, j) += sign13 * std::conj(civec_i->data(fullpos_i)) * std::conj(civec_j->data(fullpos_j)) * mo_trev->element(Dp, Hp) * mo_trev->element(Ap, Ep) * mo_trev->element(Bp, Fp) * mo_trev->element(Cp, Gp);
                                                      } else {
                                                        if ((Ap == Ep) && (Bp == Fp) && (Cp == Gp))
                                                          trev_h_->element(i, j) += sign13 * std::conj(civec_i->data(fullpos_i)) * civec_j->data(fullpos_j) * mo_trev->element(Dp, Hp);
                                                        if ((Dp == Hp) && (Bp == Fp) && (Cp == Gp))
                                                          trev_h_->element(i, j) += sign13 * std::conj(civec_i->data(fullpos_i)) * civec_j->data(fullpos_j) * mo_trev->element(Ap, Ep);
                                                        if ((Ap == Ep) && (Dp == Hp) && (Cp == Gp))
                                                          trev_h_->element(i, j) += sign13 * std::conj(civec_i->data(fullpos_i)) * civec_j->data(fullpos_j) * mo_trev->element(Bp, Fp);
                                                        if ((Ap == Ep) && (Bp == Fp) && (Dp == Hp))
                                                          trev_h_->element(i, j) += sign13 * std::conj(civec_i->data(fullpos_i)) * civec_j->data(fullpos_j) * mo_trev->element(Cp, Gp);
                                                      }
                                                    }
                                                  }
                                                }
                                              }
                                            }
                                          }
                                        }
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    trev_h_->print("NEW CALCULATION of time-reversal matrix in ZFCI states", 24);
  }
#endif

# if 1
  { // Attempt to compute time-reversal matrix in Determinant basis
    const static Comb combination;
    const size_t ndet = combination(2*zfci.norb(), zfci.nele());
    auto trev_det = make_shared<ZMatrix>(ndet, ndet);
    trev_det->zero();
    const int maxa = std::min(zfci.nele(), zfci.norb());
    const int mina = std::max(zfci.nele() - zfci.norb(), 0);
    /***********************/
    const bool conjop = !zfci.idata()->get<bool>("aniso_testconj", false);
    /***********************/

    /********************/
    // Redundant, but need this info to be available
    vector<int> aniso_state;
    aniso_state.resize(nspin1_);
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
    }
    /********************/

    vector<array<int,2>> ab = {};
    for (int j = maxa; j >= mina; --j)
      ab.push_back({{j, zfci.nele()-j}});

    int ipos = -1;
    int jpos = -1;

    // Loop over pairs of spin sectors
    for (int k = 0; k != ab.size(); ++k) {
      const int ipos_save = ipos;
      jpos = -1;
      for (int l = 0; l != ab.size(); ++l) {
        auto det_i = zfci.cc()->find(ab[k][0], ab[k][1])->data(0)->det();
        auto det_j = zfci.cc()->find(ab[l][0], ab[l][1])->data(0)->det();
        const int jpos_save = jpos;
        ipos = ipos_save;

        // Loop over pairs of determinants
        for (auto& ia : det_i->string_bits_a()) {
          for (auto& ib : det_i->string_bits_b()) {
            ++ipos;
            jpos = jpos_save;
            if (jpos == -1) cout << " ** Determinant in position ipos = " << std::setw(4) << ipos << " is " << print_bit(ia, ib, det_i->norb()) << endl;

            for (auto& ja : det_j->string_bits_a()) {
              for (auto& jb : det_j->string_bits_b()) {
                ++jpos;
                //cout << " matrix element to be set = (" << std::setw(4) << ipos << ", " << std::setw(4) << jpos << ")" << endl;

                // 1 electron version - Loop over pairs of MOs
                if (zfci.nele() == 1) {
                  for (int Ap = 0; Ap != 2*norb; ++Ap) {
                    const int A = Ap % norb;
                    if ((ia[A] && Ap < norb) || (ib[A] && Ap >= norb)) {
                      for (int Bp = 0; Bp != 2*norb; ++Bp) {
                        const int B = Bp % norb;
                        if ((ja[B] && Bp < norb) || (jb[B] && Bp >= norb)) {
                            trev_det->element(ipos, jpos) += mo_trev->element(Ap, Bp);
                        }
                      }
                    }
                  }
                }

                // 2 electron version - Loop over sets of 4 MOs
                if (zfci.nele() == 2) {
                  for (int Ap = 0; Ap != 2*norb; ++Ap) {
                    const int A = Ap % norb;
                    if ((ia[A] && Ap < norb) || (ib[A] && Ap >= norb)) {
                      for (int Bp = 0; Bp != 2*norb; ++Bp) {
                        if (Ap == Bp) continue;
                        const int B = Bp % norb;
                        if ((ia[B] && Bp < norb) || (ib[B] && Bp >= norb)) {
                          for (int Cp = 0; Cp != 2*norb; ++Cp) {
                            const int C = Cp % norb;
                            if ((ja[C] && Cp < norb) || (jb[C] && Cp >= norb)) {
                              for (int Dp = 0; Dp != 2*norb; ++Dp) {
                                if (Cp == Dp) continue;
                                const int D = Dp % norb;
                                if ((ja[D] && Dp < norb) || (jb[D] && Dp >= norb)) {
                                  const double sign1 = Ap < Bp ? 1.0 : -1.0;
                                  const double sign2 = Cp < Dp ? 1.0 : -1.0;
                                  const double sign3 = sign1 * sign2 / 2.0;
                                  if (conjop) {
                                      trev_det->element(ipos, jpos) += sign3 * mo_trev->element(Bp, Dp) * mo_trev->element(Ap, Cp);
                                  } else {
                                    if (Ap == Cp)
                                      trev_det->element(ipos, jpos) += sign3 * mo_trev->element(Bp, Dp);
                                    if (Bp == Dp)
                                      trev_det->element(ipos, jpos) += sign3 * mo_trev->element(Ap, Cp);
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }

                // 3 electron version - Loop over sets of 6 MOs
                if (zfci.nele() == 3) {
                  for (int Ap = 0; Ap != 2*norb; ++Ap) {
                    const int A = Ap % norb;
                    if ((ia[A] && Ap < norb) || (ib[A] && Ap >= norb)) {

                      for (int Bp = 0; Bp != 2*norb; ++Bp) {
                        if (Ap == Bp) continue;
                        const int B = Bp % norb;
                        if ((ia[B] && Bp < norb) || (ib[B] && Bp >= norb)) {

                          for (int Cp = 0; Cp != 2*norb; ++Cp) {
                            if ((Ap == Cp) || (Bp == Cp)) continue;
                            const int C = Cp % norb;
                            if ((ia[C] && Cp < norb) || (ib[C] && Cp >= norb)) {

                              for (int Dp = 0; Dp != 2*norb; ++Dp) {
                                const int D = Dp % norb;
                                if ((ja[D] && Dp < norb) || (jb[D] && Dp >= norb)) {

                                  for (int Ep = 0; Ep != 2*norb; ++Ep) {
                                    if (Dp == Ep) continue;
                                    const int E = Ep % norb;
                                    if ((ja[E] && Ep < norb) || (jb[E] && Ep >= norb)) {

                                      for (int Fp = 0; Fp != 2*norb; ++Fp) {
                                        if ((Dp == Fp) || (Ep == Fp)) continue;
                                        const int F = Fp % norb;
                                        if ((ja[F] && Fp < norb) || (jb[F] && Fp >= norb)) {

                                          const double sign1 = Ap < Bp ? 1.0 : -1.0;
                                          const double sign2 = Ap < Cp ? 1.0 : -1.0;
                                          const double sign3 = Bp < Cp ? 1.0 : -1.0;
                                          const double sign4 = Dp < Ep ? 1.0 : -1.0;
                                          const double sign5 = Dp < Fp ? 1.0 : -1.0;
                                          const double sign6 = Ep < Fp ? 1.0 : -1.0;
                                          const double sign7 = sign1 * sign2 * sign3 * sign4 * sign5 * sign6 / 6.0;
                                          if (conjop) {
                                            trev_det->element(ipos, jpos) += sign7 * mo_trev->element(Cp, Fp) * mo_trev->element(Ap, Dp) * mo_trev->element(Bp, Ep);
                                          } else {
                                            if ((Ap == Dp) && (Bp == Ep))
                                              trev_det->element(ipos, jpos) += sign7 * mo_trev->element(Cp, Fp);
                                            if ((Bp == Ep) && (Cp == Fp))
                                              trev_det->element(ipos, jpos) += sign7 * mo_trev->element(Ap, Dp);
                                            if ((Ap == Dp) && (Cp == Fp))
                                              trev_det->element(ipos, jpos) += sign7 * mo_trev->element(Bp, Ep);
                                          }
                                        }
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }

                // 4 electron version - Loop over sets of 8 MOs
                if (zfci.nele() == 4) {
                  for (int Ap = 0; Ap != 2*norb; ++Ap) {
                    const int A = Ap % norb;
                    if ((ia[A] && Ap < norb) || (ib[A] && Ap >= norb)) {

                      for (int Bp = 0; Bp != 2*norb; ++Bp) {
                        if (Ap == Bp) continue;
                        const int B = Bp % norb;
                        if ((ia[B] && Bp < norb) || (ib[B] && Bp >= norb)) {

                          for (int Cp = 0; Cp != 2*norb; ++Cp) {
                            if ((Ap == Cp) || (Bp == Cp)) continue;
                            const int C = Cp % norb;
                            if ((ia[C] && Cp < norb) || (ib[C] && Cp >= norb)) {

                              for (int Dp = 0; Dp != 2*norb; ++Dp) {
                                if ((Ap == Dp) || (Bp == Dp) || (Cp == Dp)) continue;
                                const int D = Dp % norb;
                                if ((ia[D] && Dp < norb) || (ib[D] && Dp >= norb)) {

                                  for (int Ep = 0; Ep != 2*norb; ++Ep) {
                                    const int E = Ep % norb;
                                    if ((ja[E] && Ep < norb) || (jb[E] && Ep >= norb)) {

                                      for (int Fp = 0; Fp != 2*norb; ++Fp) {
                                        if (Ep == Fp) continue;
                                        const int F = Fp % norb;
                                        if ((ja[F] && Fp < norb) || (jb[F] && Fp >= norb)) {

                                          for (int Gp = 0; Gp != 2*norb; ++Gp) {
                                            if ((Ep == Gp) || (Fp == Gp)) continue;
                                            const int G = Gp % norb;
                                            if ((ja[G] && Gp < norb) || (jb[G] && Gp >= norb)) {

                                              for (int Hp = 0; Hp != 2*norb; ++Hp) {
                                                if ((Ep == Hp) || (Fp == Hp) || (Gp == Hp)) continue;
                                                const int H = Hp % norb;
                                                if ((ja[H] && Hp < norb) || (jb[H] && Hp >= norb)) {

                                                  const double sign1 = Ap < Bp ? 1.0 : -1.0;
                                                  const double sign2 = Ap < Cp ? 1.0 : -1.0;
                                                  const double sign3 = Ap < Dp ? 1.0 : -1.0;
                                                  const double sign4 = Bp < Cp ? 1.0 : -1.0;
                                                  const double sign5 = Bp < Dp ? 1.0 : -1.0;
                                                  const double sign6 = Cp < Dp ? 1.0 : -1.0;
                                                  const double sign7 = Ep < Fp ? 1.0 : -1.0;
                                                  const double sign8 = Ep < Gp ? 1.0 : -1.0;
                                                  const double sign9 = Ep < Hp ? 1.0 : -1.0;
                                                  const double sign10 = Fp < Gp ? 1.0 : -1.0;
                                                  const double sign11 = Fp < Hp ? 1.0 : -1.0;
                                                  const double sign12 = Gp < Hp ? 1.0 : -1.0;
                                                  const double sign13 = sign1 * sign2 * sign3 * sign4 * sign5 * sign6 * sign7 * sign8 * sign9 * sign10 * sign11 * sign12 / 24.0;

                                                  if (conjop) {
                                                      trev_det->element(ipos, jpos) += sign13 * mo_trev->element(Dp, Hp) * mo_trev->element(Ap, Ep) * mo_trev->element(Bp, Fp) * mo_trev->element(Cp , Gp);
                                                  } else {
                                                    if ((Ap == Ep) && (Bp == Fp) && (Cp == Gp))
                                                      trev_det->element(ipos, jpos) += sign13 * mo_trev->element(Dp, Hp);
                                                    if ((Dp == Hp) && (Bp == Fp) && (Cp == Gp))
                                                      trev_det->element(ipos, jpos) += sign13 * mo_trev->element(Ap, Ep);
                                                    if ((Ap == Ep) && (Dp == Hp) && (Cp == Gp))
                                                      trev_det->element(ipos, jpos) += sign13 * mo_trev->element(Bp, Fp);
                                                    if ((Ap == Ep) && (Bp == Fp) && (Dp == Hp))
                                                      trev_det->element(ipos, jpos) += sign13 * mo_trev->element(Cp, Gp);
                                                  }
                                                }
                                              }
                                            }
                                          }
                                        }
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    ++ipos;
    ++jpos;
    cout << " Final ipos = " << ipos << ", final jpos = " << jpos << ", ndet = " << ndet << endl;
    assert(ipos == jpos);
    assert(ipos == ndet);
    trev_det->print("NEW CALCULATION of time-reversal matrix in a basis of determinants", 50);
  }
#endif

  { // spin angular momentum
    for (int i = 0; i != 3; ++i) {
      magnetic_moment[i] = ao_spin(i)->copy();
      magnetic_moment[i]->scale(g_elec__);
    }
  }

  array<shared_ptr<ZMatrix>,3> ao_orbang;
  { // orbital angular momentum
    // TODO For geometries with only one metal atom, use that atom's position as default mcoord
    const array<double, 3> mcoord = zfci.idata()->get_array<double,3>("aniso_center", array<double, 3>({{0.0, 0.0, 0.0}}));
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
    zfci_spin_[i] = make_shared<ZMatrix>(nspin1_, nspin1_);
    zfci_orbang_[i] = make_shared<ZMatrix>(nspin1_, nspin1_);
  }
  trev_h_ = make_shared<ZMatrix>(nspin1_, nspin1_);
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

      trev_h_->element(i,j) = aodensity.dot_product(*ao_trev);
      for (int k = 0; k != 3; ++k) {
        spinop_h_[k]->element(i,j) = aodenconj.dot_product(*magnetic_moment[k]);
        zfci_spin_[k]->element(i,j) = aodenconj.dot_product(*ao_spin(k));
        zfci_orbang_[k]->element(i,j) = aodenconj.dot_product(*ao_orbang[k]);
      }
    }
  }

#if 0
  for (int k = 0; k != 3; ++k) {
    spinop_h_[k]->print("Magnetic moment in ZFCI basis");
    zfci_spin_[k]->print("Spin angular momentum in ZFCI basis");
    zfci_orbang_[k]->print("Orbital angular momentum in ZFCI basis");
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
            temp->element(i, j) += 0.5 * spinop_h_[i]->element(k, l) * spinop_h_[j]->element(l, k);

    Atensor = temp->get_real_part();
    assert(temp->get_imag_part()->rms() < 1.0e-10);
    Atransform = Atensor->copy();
    Atransform->diagonalize(Aeig);

    Atensor->print("A tensor");
    cout << endl;
    for (int i = 0; i != 3; ++i)
      cout << " *** Atensor eigenvalue " << i << " = " << Aeig[i] << endl;
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
    auto Gtensor = make_shared<Matrix>(*gtensor ^ *gtensor);

    gtensor->print("g-tensor");
    cout << endl;
    Gtensor->print("G-tensor");
    cout << endl;

    assert((*Gtensor - 4.0 * factor * *Atensor).rms() < 1.0e-8);

    cout << "  Main axes of magnetic anisotropy:" << endl;
    for (int i = 0; i != 3; ++i) {
      cout << "   " << i << " |g_" << i << "| = " << setw(12) << gval[i] << ",  axis = ( ";
      cout << setw(12) << Atransform->element(i, 0) << ", ";
      cout << setw(12) << Atransform->element(i, 1) << ", ";
      cout << setw(12) << Atransform->element(i, 2) << " )" << endl;
    }
    cout << endl;
  }
  return Atransform;
}


shared_ptr<const ZMatrix> Pseudospin::compute_spin_eigenvalues(const array<double, 3> rotation, const ZHarrison& zfci) const {

  // Diagonalize S_z to get pseudospin eigenstates as combinations of ZFCI Hamiltonian eigenstates
  auto transform = make_shared<ZMatrix>(nspin1_, nspin1_);
  const complex<double> scale = 1.0 / std::sqrt(rotation[0]*rotation[0] + rotation[1]*rotation[1] + rotation[2]*rotation[2]);
  assert(std::abs(std::imag(scale)) < 1.0e-10);
  for (int i = 0; i != 3; ++i)
    *transform += scale * rotation[i] * *spinop_h_[i];
  VectorB zeig(nspin1_);
#ifndef NDEBUG
  auto spinmat_to_diag = transform->copy();
#endif
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

#if 1
  // For testing arbitrary phase shifts applied to pseudospin eigenfunctions
  const complex<double> dadjust = polar(1.0, zfci.idata()->get<double>("aniso_dadjust", 0.0));
  cout << " ** Phase shift = " << std::arg(dadjust) << endl;
  for (int i = 0; i != nspin1_; ++i) {
    transform->element(i, 0) = dadjust * transform->element(i, 0);
    //transform->element(i, nspin_) = std::conj(dadjust) * transform->element(i, nspin_);
#endif
#if 1
  {
    // Attempt to fix time-reversal symmetry from the Ci vectors
    // TODO Probably this could be made significantly more efficient - see if it's a problematic
    const int maxa = std::min(zfci.nele(), zfci.norb());
    const int mina = std::max(zfci.nele() - zfci.norb(), 0);
    cout << setprecision(8);

    // Loop over pairs of Kramers conjugates
    // km starts with -1/2 (or 0) and goes to -S
    // kp is the Kramers conjugate of km
    for (int km = nspin_; km >= nspin1_ / 2; --km) {
      const int kp = nspin_ - km;
      cout << endl;

      //cout << " nele = " << zfci.nele() << ", norb = " << zfci.norb() << ", number alpha/beta runs from " << maxa << " to " << mina << endl;
      vector<array<int,2>> ab = {};
      for (int j = maxa; j >= mina; --j)
        ab.push_back({{j, zfci.nele()-j}});

      double maxval2 = 0.0;
      bool signflip = false;
      complex<double> maxvalm;
      complex<double> maxvalp;

      // Loop over sectors of different spin values
      for (int j = 0; j != ab.size(); ++j) {

        auto det_m = zfci.cc()->find(ab[j][0], ab[j][1])->data(0)->det();
        auto det_p = zfci.cc()->find(ab[j][1], ab[j][0])->data(0)->det();

        // Loop over different Slater determinants in this sector
        for (auto& ia : det_m->string_bits_a()) {
          for (auto& ib : det_m->string_bits_b()) {
            complex<double> valm = 0.0;
            complex<double> valp = 0.0;

            // Loop over contributions from the different ZFCI eigenstates
            for (int k = 0; k != nspin1_; ++k) {

              auto civec_m = zfci.cc()->find(ab[j][0], ab[j][1])->data(aniso_state[k]);
              auto civec_p = zfci.cc()->find(ab[j][1], ab[j][0])->data(aniso_state[k]);

              const int pos1m = det_m->lexical<0>(ia);
              const int pos2m = det_m->lexical<1>(ib);
              const int pos1p = det_p->lexical<0>(ib);
              const int pos2p = det_p->lexical<1>(ia);
              const int fullpos_m = pos1m*det_m->lenb() + pos2m;
              const int fullpos_p = pos1p*det_p->lenb() + pos2p;
              const complex<double> cival_m = civec_m->data(fullpos_m) * transform->element(k, km);
              const complex<double> cival_p = civec_p->data(fullpos_p) * transform->element(k, kp);
              valm += cival_m;
              valp += cival_p;
            }
            const double mag2 = (std::abs(valm) + std::abs(valp));
            if (mag2 > maxval2) {
              maxval2 = mag2;
              maxvalm = valm;
              maxvalp = valp;
              signflip = (ab[j][1] % 2 == 0) ? false : true;
            }
          }
        }
      }
      assert(std::abs(std::abs(maxvalm) - std::abs(maxvalp)) < 1.0e-6);
      double phase_sum = std::arg(maxvalm) + std::arg(maxvalp);
      if (signflip == true) phase_sum += pi__;

      cout << " km = " << km << ", phase shift = " << sign * phase_sum << endl;
      const complex<double> shift = std::polar(1.0, -phase_sum);
      for (int i = 0; i != nspin1_; ++i)
        transform->element(i, km) = shift * transform->element(i, km);
    }
  }
#endif

#if 1
  // Enforce time-reversal symmetry by looking at the matrices - misses the {1/2, -1/2} element
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

  // Check the spin matrices as well with basically the same procedure, since sometimes we miss a phase due to numerically zero entries in the Hamiltonian
  for (int j = 0; j != 3; ++j) {
    const ZMatrix spinop_x = *transform % *spinop_h_[j] * *transform;
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
#endif

#ifndef NDEBUG
  {
    auto diagspin = transform->clone();
    for (int i = 0; i != nspin1_; ++i)
      diagspin->element(i, i) = zeig[i];
    *diagspin = *transform * *diagspin ^ *transform;
    /*
    spinmat_to_diag->print("Original spin matrix");
    diagspin->print("Recalculated spin matrix");
    (*spinmat_to_diag - *diagspin).print("Error");
    cout << endl;
    */
    assert((*diagspin - *spinmat_to_diag).rms() < 1.0e-6);
  }
#endif

  cout << "    The z-axis is set to (" << rotation[0] << ", " << rotation[1] << ", " << rotation[2] << ")." << endl << endl;
  for (int i = 0; i != nspin1_; ++i)
    cout << "    Pseudospin eigenvalue " << i+1 << " = " << setw(12) << zeig[i] << endl;

  // We can no longer use this option, since I made this function const...
  //if (numerical_eig) {
  //  cout << "  **  By request, we compute the Hamiltonian using eigenvalues of S_z, rather than the canonical m_s values" << endl;
  //  update_spin_matrices(zeig);
  //}

  /***********************/
  const bool conjop = !zfci.idata()->get<bool>("aniso_testconj", false);
  /***********************/

  shared_ptr<ZMatrix> spinham_s = make_shared<ZMatrix>(*transform % *spinham_h_ * *transform);
  array<shared_ptr<ZMatrix>, 3> spinop_s;
  auto trev_s = conjop ? make_shared<ZMatrix>(*transform % *trev_h_ * *transform->get_conjg()) : make_shared<ZMatrix>(*transform % *trev_h_ * *transform);
  for (int i = 0; i != 3; ++i) {
    spinop_s[i] = make_shared<ZMatrix>(*transform % *spinop_h_[i] * *transform);
  }

  spinham_s->print("Spin Hamiltonian");
  cout << endl;
  trev_s->print("Time-reversal operator!  (in pseudospin states)");
  cout << endl;
  for (int i = 0; i != 3; ++i) {
    spinop_s[i]->print("Spin matrix, component " + to_string(i));
    cout << endl;
  }

  if (!is_t_symmetric(*spinham_s, /*hermitian*/true, /*time reversal*/true))
    throw runtime_error("The spin Hamiltonian seems to not have proper time-reversal symmetry.  Check that your spin value and states mapped are reasonable.");

  cout << endl;
  for (int i = 0; i != 3; ++i) {
    //spinop_h_[i]->print("ZFCI basis spin matrix " + to_string(i+1));
#if 0
    spinop_s[i]->print("Spin eigenfunction basis magnetic moment matrix " + to_string(i+1));
    (*transform % *zfci_spin_[i] * *transform).print("Spin eigenfunction basis spin angular momentum matrix " + to_string(i+1));
    (*transform % *zfci_orbang_[i] * *transform).print("Spin eigenfunction basis orbital angular momentum matrix " + to_string(i+1));
    cout << endl;
#endif
    assert(is_t_symmetric(*spinop_s[i], /*hermitian*/true, /*time reversal*/false));
  }
  assert(is_t_symmetric(*spinham_s, /*hermitian*/true, /*time reversal*/true));

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

