//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#ifndef __NEWINT_CASSCF_WERNER_H
#define __NEWINT_CASSCF_WERNER_H

#include <src/casscf/casscf.h>

class WernerKnowles : public CASSCF {
  protected:
    void common_init() {
      std::cout << "    * Using the two-step Werner-Knowles algorithm (see JCP 1985)" << std::endl;
    };

  public:
    WernerKnowles(const std::multimap<std::string, std::string> idat, const std::shared_ptr<Geometry> geom, std::shared_ptr<Reference> ref)
      : CASSCF(idat, geom, ref) {common_init(); };
    ~WernerKnowles() {};

    void compute();


}; 

#endif

