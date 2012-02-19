//
// Author : Toru Shiozaki
// Date   : Feb 2012
//

#ifndef __SRC_SMITH_TENSOR_H
#define __SRC_SMITH_TENSOR_H

#include <list>
#include <map>
#include <memory>
#include <iostream>
#include <src/smith/storage.h>
#include <src/smith/indexrange.h>
#include <src/smith/loopgenerator.h>

namespace SMITH {

template <typename T>
class Tensor {
  protected:
   std::shared_ptr<T> data_; 

  public:
    Tensor(std::list<IndexRange> in) {
      // make loop 
      LoopGenerator lg(in);
      std::vector<std::vector<Index> > loop = lg.block_loop();
#if 0
      for (auto i = loop.begin(); i != loop.end(); ++i) {
        for (auto j = i->begin(); j != i->end(); ++j)
          std::cout << "(" << j->offset() << "," << j->size() << ")";
        std::cout << std::endl;
      }
#endif
      // first compute hashtags and 
      std::map<size_t, size_t> hashmap;

      std::shared_ptr<T> tmp(new T(hashmap));
      data_ = tmp;
    };
    ~Tensor() {};
};

}

#endif

