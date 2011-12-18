//
// Author : Toru Shiozaki
// Date   : Dec 2011
//
// Desc :: The implementation closely follows Knowles and Handy 1984 CPL.
//         It is amazing how easy it is to implement FCI !!
//

#ifndef __NEWINT_FCI_FCI_H
#define __NEWINT_FCI_FCI_H

#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <src/scf/scf.h>
#ifdef USE_SSE42_INTRINSICS
#include <nmmintrin.h>
#endif
#include <cassert>
#include <iostream>

class FCI {
  protected:
    boost::shared_ptr<SCF> ref_;
    const boost::shared_ptr<Geometry> geom_;

    // Knowles & Handy lexical mapping
    std::vector<unsigned int> zkla_; // contains zkl and corresponding bits
    std::vector<unsigned int> zklb_;
    // string list
    std::vector<unsigned int> stringa_;
    std::vector<unsigned int> stringb_;

    // some init functions
    // lexical maps (Zkl)
    void const_lexical_mapping_();
    // alpha and beta string lists
    void const_string_lists_();
    // single displacement vectors Phi's
    void const_phis_();

    // numbers of electrons
    const int nelea_;
    const int neleb_;
    const int ncore_;
    const int norb_;

    // configuration list
    std::vector<boost::tuple<unsigned int, int, unsigned int, unsigned int> > phia_;
    std::vector<boost::tuple<unsigned int, int, unsigned int, unsigned int> > phib_;

    int numofbits(unsigned int bits) { //
#ifndef USE_SSE42_INTRINSICS
      bits = (bits & 0x55555555) + (bits >> 1 & 0x55555555);
      bits = (bits & 0x33333333) + (bits >> 2 & 0x33333333);
      bits = (bits & 0x0f0f0f0f) + (bits >> 4 & 0x0f0f0f0f);
      bits = (bits & 0x00ff00ff) + (bits >> 8 & 0x00ff00ff);
      return (bits & 0x0000ffff) + (bits >>16 & 0x0000ffff);
#else
      return _mm_popcnt_u32(bits); // not tested. i7 or later - good reason to buy a new laptop.
#endif
    };

    int sign(unsigned int bit, int i, int j) {
      const int ii = 32 - std::max(i,j); // TODO 32 is hard wired...
      const int jj = std::min(i,j) + 1;
      bit = (((bit >> jj) << (jj+ii)) >> ii);
      return 1 - ((numofbits(bit) & 1) << 1);
    };

    // maps bit to lexical numbers.
    unsigned int lexicala(int bit) {
      unsigned int out = 0, k = 0;
      for (int i = 0; i != norb_; ++i, bit = (bit >> 1)) 
        if (bit & 1) { out += zkla(k,i); ++k; }
      assert(k == nelea_);
      return out;
    };
    unsigned int lexicalb(int bit) {
      unsigned int out = 0, k = 0;
      for (int i = 0; i != norb_; ++i, bit = (bit >> 1))
        if (bit & 1) { out += zklb(k,i); ++k; }
      assert(k == neleb_);
      return out;
    };
  
    // some utility functions
    unsigned int& zkla(int i, int j) { return zkla_[i*nelea_+j]; };
    unsigned int& zklb(int i, int j) { return zklb_[i*neleb_+j]; };
    unsigned int stringa(int i) const { return stringa_[i]; };
    unsigned int stringb(int i) const { return stringb_[i]; };

  public:
    FCI(const boost::shared_ptr<Geometry>);
    ~FCI();
    void compute();

};


#endif

