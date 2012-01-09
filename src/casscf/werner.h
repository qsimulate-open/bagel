//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#ifndef __NEWINT_CASSCF_WERNER_H
#define __NEWINT_CASSCF_WERNER_H

#include <src/casscf/casscf.h>
#include <src/casscf/rotfile.h>
#include <src/casscf/jvec.h>

class WernerKnowles : public CASSCF {
  protected:
    void common_init() {
      std::cout << "    * Using the two-step Werner-Knowles algorithm (see JCP 1985)" << std::endl << std::endl;
    };

    std::shared_ptr<Matrix1e> compute_bvec(std::shared_ptr<FCI>, std::shared_ptr<Jvec>, std::shared_ptr<Matrix1e>, 
                                           std::shared_ptr<Coeff>);
    std::shared_ptr<Matrix1e> compute_bvec(std::shared_ptr<FCI>, std::shared_ptr<Jvec>, std::shared_ptr<Matrix1e>, 
                                           std::shared_ptr<Matrix1e>,std::shared_ptr<Coeff>);

    double thresh_mmicro_;
    int max_mmicro_iter_;

  public:
    WernerKnowles(const std::multimap<std::string, std::string> idat, const std::shared_ptr<Geometry> geom, std::shared_ptr<Reference> ref)
      : CASSCF(idat, geom, ref) {common_init(); 
      // get thresh (for micro iteration) from the input
      thresh_mmicro_ = read_input<double>(idat, "thresh_mmicro", thresh_micro_);
      max_mmicro_iter_ = read_input<int>(idat, "maxiter_mmicro", 5);
    };
    ~WernerKnowles() {};

    void compute();


}; 

#endif

