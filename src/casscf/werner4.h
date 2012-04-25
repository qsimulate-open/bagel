//
// Author : Toru Shiozaki
// Date   : April 2012
//

// Implements second-order CASSCF (two step)
// Uses Werner-Knowles algorithm, but construcging 4-index MO integrals.
// This is, therefere, only for reference

#ifndef __SRC_CASSCF_WERNER4_H
#define __SRC_CASSCF_WERNER4_H

#include <src/casscf/werner.h>
#include <src/casscf/jkop.h>

class Werner4 : public WernerKnowles {
  protected:
    void common_init() {
      std::cout << "    * Using the two-step 4-index Werner-Knowles algorithm" << std::endl << std::endl;
    };

  public:
    Werner4(const std::multimap<std::string, std::string> idat, const std::shared_ptr<Geometry> geom, std::shared_ptr<Reference> ref)
      : WernerKnowles(idat, geom, ref) { };
    ~Werner4() {};

    void compute();

}; 


#endif
