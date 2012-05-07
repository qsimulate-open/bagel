//
// Newint - Parallel electron correlation program.
// Filename: fci.h
// Copyright (C) 2011 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//

// Desc :: The implementation closely follows Knowles and Handy 1984 CPL.
//         It is amazing how easy it is to implement FCI !!
//

#ifndef __NEWINT_FCI_FCI_H
#define __NEWINT_FCI_FCI_H

#include <tuple>
#include <src/scf/scf.h>
#ifdef USE_SSE42_INTRINSICS
#include <nmmintrin.h>
#endif
#include <cassert>
#include <iostream>
#include <memory>
#include <src/util/input.h>
#include <src/fci/civec.h>
#include <src/fci/mofile.h>
#include <src/wfn/rdm.h>
#include <src/wfn/reference.h>

class FCI {

  protected:
    // input
    std::multimap<std::string, std::string> idata_; 
    // reference
    std::shared_ptr<Reference> ref_;
    // geometry file
    const std::shared_ptr<Geometry> geom_;
    // number of states
    int nstate_;
    // max #iteration
    int max_iter_;
    // threshold for variants
    double thresh_;

    // numbers of electrons
    int nelea_;
    int neleb_;
    int ncore_;
    int norb_;
    // core
    double core_energy_; // <- usually initialized in create_iiii
    // total energy
    std::vector<double> energy_;

    // CI vector at convergence
    std::shared_ptr<Dvec> cc_;
    // RDMs; should be resized in constructors
    std::vector<std::shared_ptr<RDM<1> > > rdm1_;
    std::vector<std::shared_ptr<RDM<2> > > rdm2_;
    // state averaged RDM
    std::vector<double> weight_;
    std::shared_ptr<RDM<1> > rdm1_av_;
    std::shared_ptr<RDM<2> > rdm2_av_;
    // MO integrals 
    std::shared_ptr<MOFile> jop_;

    //
    // Below here, internal private members
    //

    // Knowles & Handy lexical mapping
    std::vector<unsigned int> zkl_; // contains zkl (Two dimenional array. See the public function).
    // string lists
    std::vector<unsigned int> stringa_;
    std::vector<unsigned int> stringb_;
    // denominator
    std::shared_ptr<Civec> denom_;

    // configuration list
    std::vector<std::vector<std::tuple<unsigned int, int, unsigned int> > > phia_;
    std::vector<std::vector<std::tuple<unsigned int, int, unsigned int> > > phib_;

    // some init functions
    void common_init();
    void create_Jiiii();
    // lexical maps (Zkl)
    void const_lexical_mapping_();
    // alpha and beta string lists
    void const_string_lists_();
    // single displacement vectors Phi's
    template <int>
    void const_phis_(const std::vector<unsigned int>& string,
                     std::vector<std::vector<std::tuple<unsigned int, int, unsigned int> > >& target, bool compress=true);
    // generate spin-adapted guess configurations
    void generate_guess(const int nspin, const int nstate, std::shared_ptr<Dvec>);
    // denominator
    void const_denom();
    // obtain determinants for guess generation
    std::vector<std::pair<int, int> > detseeds(const int ndet);

    int numofbits(unsigned int bits) { //
#ifndef USE_SSE42_INTRINSICS
      bits = (bits & 0x55555555) + (bits >> 1 & 0x55555555); bits = (bits & 0x33333333) + (bits >> 2 & 0x33333333);
      bits = (bits & 0x0f0f0f0f) + (bits >> 4 & 0x0f0f0f0f); bits = (bits & 0x00ff00ff) + (bits >> 8 & 0x00ff00ff);
      return (bits & 0x0000ffff) + (bits >>16 & 0x0000ffff); // can be cheaper, but it is fine for the time being...
#else
      return _mm_popcnt_u32(bits); // not tested. i7 or later - good reason to buy a new laptop.
#endif
    };

    int sign(unsigned int bit, int i, int j) {
      // masking irrelevant bits
      const unsigned int ii = ~((1 << (std::min(i,j)+1)) - 1);
      const unsigned int jj = ((1 << (std::max(i,j))) - 1); 
      bit = (bit & ii) & jj;
      return 1 - ((numofbits(bit) & 1) << 1);
    };

    // maps bit to lexical numbers.
    template <int spin> unsigned int lexical(int bit) {
      unsigned int out = 0, k = 0;
      for (int i = 0; i != norb_; ++i, bit >>= 1) 
        if (bit & 1) { out += zkl(k,i, spin); ++k; }
      return out;
    };

    // this is slow but robust implementation of bit to number converter.
    std::vector<int> bit_to_numbers(unsigned int bit) {
      std::vector<int> out;
      for (int i = 0; bit != 0u; ++i, bit >>= 1) if (bit & 1) out.push_back(i); 
      return out;
    };
  
    // some utility functions
    unsigned int& zkl(int i, int j, int spin) { return zkl_[i*norb_+j+spin*nelea_*norb_]; };

    unsigned int stringa(int i) const { return stringa_[i]; };
    unsigned int stringb(int i) const { return stringb_[i]; };

    // run-time functions
    void form_sigma(std::shared_ptr<Dvec> c, std::shared_ptr<Dvec> sigma, std::shared_ptr<Dvec> d, std::shared_ptr<Dvec> e,
                    const std::vector<int>& conv);
    void sigma_1(std::shared_ptr<Civec> cc, std::shared_ptr<Civec> sigma);
    void sigma_3(std::shared_ptr<Civec> cc, std::shared_ptr<Civec> sigma);
    void sigma_2a1(std::shared_ptr<Civec> cc, std::shared_ptr<Dvec> d);
    void sigma_2a2(std::shared_ptr<Civec> cc, std::shared_ptr<Dvec> d);
    void sigma_2b (std::shared_ptr<Dvec> d, std::shared_ptr<Dvec> e);
    void sigma_2c1(std::shared_ptr<Civec> sigma, std::shared_ptr<Dvec> e);
    void sigma_2c2(std::shared_ptr<Civec> sigma, std::shared_ptr<Dvec> e);

    // functions related to natural orbitals
    void update_rdms(const std::vector<double>& coeff); 

    // print functions
    void print_header() const;
    void print_civectors(const std::vector<std::shared_ptr<Civec> >, const double thr = 0.05) const;
    void print_timing_(const std::string, int& time, std::vector<std::pair<std::string, double> >&) const; 

  public:
    FCI(const std::multimap<std::string, std::string>, const std::shared_ptr<Geometry>, std::shared_ptr<Reference>);
    ~FCI();
    void compute();

    // returns members
    int norb() const { return norb_; };
    int ncore() const { return ncore_; };
    double core_energy() const { return core_energy_; };
    // sets members
    void set_core_energy(const double a) { core_energy_ = a; };
    // rdms
    void compute_rdm12(); // compute all states at once + averaged rdm
    void compute_rdm12(const int istate);
    std::vector<std::shared_ptr<RDM<1> > > rdm1() { return rdm1_; };
    std::vector<std::shared_ptr<RDM<2> > > rdm2() { return rdm2_; };
    std::shared_ptr<RDM<1> > rdm1(const int i) { return rdm1_.at(i); };
    std::shared_ptr<RDM<2> > rdm2(const int i) { return rdm2_.at(i); };
    std::shared_ptr<RDM<1> > rdm1_av() { return rdm1_av_; };
    std::shared_ptr<RDM<2> > rdm2_av() { return rdm2_av_; };
    // move to natural orbitals
    std::pair<std::vector<double>, std::vector<double> > natorb_convert();

    const std::shared_ptr<Geometry> geom() const { return geom_; };

    // returns integral files
    std::shared_ptr<MOFile> jop() { return jop_; };

    // returns total energy
    std::vector<double> energy() const { return energy_; };

    // static constants
    static const int Alpha = 0;
    static const int Beta = 1;

    std::string print_bit(unsigned int bit) const {
      std::string out; 
      for (int i = 0; i != norb_; ++i, bit >>=1) { if (bit&1) { out += "1"; } else { out += "."; } }
      return out;
    };
    std::string print_bit(unsigned int bit1, unsigned int bit2) const {
      std::string out; 
      for (int i = 0; i != norb_; ++i, bit1 >>=1, bit2 >>=1) {
        if (bit1&1 && bit2&1) { out += "2"; }
        else if (bit1&1) { out += "a"; }
        else if (bit2&1) { out += "b"; }
        else { out += "."; }
      }
      return out;
    };
};


// Template function that creates the single-displacement lists (step a and b in Knowles & Handy paper).
template <int spin>
void FCI::const_phis_(const std::vector<unsigned int>& string,
      std::vector<std::vector<std::tuple<unsigned int, int, unsigned int> > >& phi, bool compress) {

  phi.clear();
  phi.resize(compress ? norb_*(norb_+1)/2 : norb_*norb_);
  for (auto iter = phi.begin(); iter != phi.end(); ++iter) {
    iter->reserve(string.size());
  }

  for (auto iter = string.begin(); iter != string.end(); ++iter) {
    for (unsigned int i = 0; i != norb_; ++i) { // annihilation
      const unsigned int ibit = (1 << i);
      if (ibit & *iter && compress) {
        const unsigned int source = lexical<spin>(*iter); 
        const unsigned int nbit = (ibit^*iter); // annihilated.
        for (unsigned int j = 0; j != norb_; ++j) { // creation
          const unsigned int jbit = (1 << j); 
          if (!(jbit & nbit)) {
            const unsigned int mbit = jbit^nbit;
            const int minij = std::min(i,j); 
            const int maxij = std::max(i,j);
            phi[minij+((maxij*(maxij+1))>>1)].push_back(std::make_tuple(lexical<spin>(mbit), sign(mbit, i, j), source));
          }
        }
      } else if (ibit & *iter) {
        const unsigned int source = lexical<spin>(*iter); 
        const unsigned int nbit = (ibit^*iter); // annihilated.
        for (unsigned int j = 0; j != norb_; ++j) { // creation
          const unsigned int jbit = (1 << j); 
          if (!(jbit & nbit)) {
            const unsigned int mbit = jbit^nbit;
            phi[i+j*norb_].push_back(std::make_tuple(lexical<spin>(mbit), sign(mbit, i, j), source));
          }
        }
      }
    }
  }

#if 0
  // sort each vectors
  for (auto iter = phi.begin(); iter != phi.end(); ++iter) {
    std::sort(iter->begin(), iter->end());
  }
#endif
};


#endif

