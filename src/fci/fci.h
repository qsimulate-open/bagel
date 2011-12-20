//
// Author : Toru Shiozaki
// Date   : Dec 2011
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
#include <src/fci/civec.h>
#include <src/fci/mofile.h>

class FCI {
  protected:
    std::shared_ptr<SCF> ref_;
    const std::shared_ptr<Geometry> geom_;

    // Knowles & Handy lexical mapping
    std::vector<unsigned int> zkl_; // contains zkl (Two dimenional array. See the public function).
    // string lists
    std::vector<unsigned int> stringa_;
    std::vector<unsigned int> stringb_;

    // some init functions
    // lexical maps (Zkl)
    void const_lexical_mapping_();
    // alpha and beta string lists
    void const_string_lists_();
    // single displacement vectors Phi's
    template <int>
    void const_phis_(const std::vector<unsigned int>&,
                     std::vector<std::tuple<unsigned int, int, unsigned int, unsigned int> >&);

    // numbers of electrons
    const int nelea_;
    const int neleb_;
    const int ncore_;
    const int norb_;

    // configuration list
    std::vector<std::tuple<unsigned int, int, unsigned int, unsigned int> > phia_;
    std::vector<std::tuple<unsigned int, int, unsigned int, unsigned int> > phib_;

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
      for (int i = 0; i != norb_; ++i, bit = (bit >> 1)) 
        if (bit & 1) { out += zkl(k,i, spin); ++k; }
      return out;
    };
  
    // some utility functions
    unsigned int& zkl(int i, int j, int spin) { return zkl_[i*nelea_+j+spin*nelea_*norb_]; };

    unsigned int stringa(int i) const { return stringa_[i]; };
    unsigned int stringb(int i) const { return stringb_[i]; };

    // run-time functions
    void form_sigma(std::shared_ptr<Dvec>, std::shared_ptr<Dvec>, std::shared_ptr<Dvec>, std::shared_ptr<Dvec>,
                    std::shared_ptr<MOFile>);

    // print functions
    void print_header() const;

  public:
    FCI(const std::shared_ptr<Geometry>);
    ~FCI();
    void compute();

    // static constants
    static const int Alpha = 0;
    static const int Beta = 1;
};


// Template function that creates the single-displacement lists (step a and b in Knowles & Handy paper).
template <int spin>
void FCI::const_phis_(const std::vector<unsigned int>& string,
                      std::vector<std::tuple<unsigned int, int, unsigned int, unsigned int> >& phi) {

  phi.resize(string.size()*norb_*norb_);
  auto piter = phi.begin();

  for (auto iter = string.begin(); iter != string.end(); ++iter) {
    for (unsigned int i = 0; i != norb_; ++i) { // annihilation
      const unsigned int ibit = (1 << i);
      if (ibit & *iter) {
        const unsigned int source = lexical<spin>(*iter); 
        const unsigned int nbit = (ibit^*iter); // annihilated.
        for (unsigned int j = 0; j != norb_; ++j) { // creation
          const unsigned int jbit = (1 << j); 
          if (!(jbit & nbit)) {
            const unsigned int mbit = jbit^nbit;
            *piter = make_tuple(lexical<spin>(mbit), sign(mbit, i, j), j*norb_+i, source); 
#if 0
            std::cout << i << j << " " << (*iter & 1)  << ((*iter >> 1) & 1) << ((*iter >> 2) & 1) << ((*iter >> 3) & 1) 
                 << ((*iter >> 4) & 1) << " " << (mbit & 1) << ((mbit >> 1) & 1) << ((mbit >> 2) & 1) <<
                 ((mbit >> 3) & 1) << ((mbit >> 4) & 1) << " " << mbit << " " << sign(mbit, i, j) << " " << lexical<spin>(mbit) << std::endl;
#endif
            ++piter;
          }
        }
      }
    }
  }
}

#endif

