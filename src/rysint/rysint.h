//
// Author: Toru Shiozaki
// Date  : May 2009
//

// Base class for the Rys-type integral evaluator

#ifndef __src_rysint_rysint_h
#define __src_rysint_rysint_h

#include <src/rysint/vrrlist.h>
#include <src/rysint/hrrlist.h>
#include <src/rysint/sortlist.h>
#include <src/scf/shell.h>
#include <vector>
#include <boost/shared_ptr.hpp>

class RysInt {
  protected:
    VRRList vrr_; 
    HRRList hrr_; 
    SortList sort_;

    bool spherical_;

    std::vector<boost::shared_ptr<Shell> > basisinfo_;
    double *data_;
    double *data2_;
    unsigned int size_final_;

    /// info for Rys quadruture
    double *roots_; 
    double *weights_;
    int rank_;

    // for screening
    int* screening_;
    int screening_size_;

  public:
    RysInt(const std::vector<boost::shared_ptr<Shell> >);
    ~RysInt();

    virtual void compute() {};

    /// retrieve a batch of integrals
    const double* data() const { return data_; };
    const double* data2() const { return data2_; };
    const bool data2_exists() const { return data2_ != NULL; };
    const unsigned int data_size() const { return size_final_; };

};

#endif

