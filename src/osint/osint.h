//
// Author: Toru Shiozaki
// Date  : May 2009
//

#ifndef __src_osint_osint_h
#define __src_osint_osint_h

#include <vector>
#include <src/scf/shell.h>
#include <src/rysint/sortlist.h>
#include <boost/shared_ptr.hpp>

class OSInt {
  protected:
    bool spherical_;

    std::vector<boost::shared_ptr<Shell> > basisinfo_;

    double* data_;
    std::vector<double> xp_, xa_, xb_, rho_, p_; 
    std::vector<double> coeffsx_, coeffsy_, coeffsz_;
    std::vector<double> coefftx_, coeffty_, coefftz_;
    double AB_[3];

    int ang0_, ang1_, cont0_, cont1_, prim0_, prim1_;

    int amax_, amax1_, amin_, asize_, asize_final_, asize_intermediate_;
    std::vector<int> amapping_;

    bool swap01_;

    SortList sort_;

    virtual void perform_VRR(double*) {};
    void perform_contraction(const int, const double*, const int, const int, double*, 
                             const std::vector<std::vector<double> >&, const std::vector<std::pair<int, int> >&, const int, 
                             const std::vector<std::vector<double> >&, const std::vector<std::pair<int, int> >&, const int);

  public:
    OSInt(const std::vector<boost::shared_ptr<Shell> >&);
    ~OSInt();

    virtual void compute() {}; 

    const double* data() const { return data_; }; 

};

#endif

