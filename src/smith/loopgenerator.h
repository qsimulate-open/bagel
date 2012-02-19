//
// Author : Toru Shiozaki
// Date   : Feb 2012
//

#ifndef __SRC_SMITH_LOOPGENERATOR_H
#define __SRC_SMITH_LOOPGENERATOR_H

#include <list>
#include <src/smith/indexrange.h>

namespace SMITH {

class LoopGenerator {
  protected:
    std::list<IndexRange> loop_;
  public:
    LoopGenerator(const std::list<IndexRange>& o) : loop_(o) {};
    ~LoopGenerator() {};

    std::vector<std::vector<Index> > block_loop() const; 
};

}

#endif
