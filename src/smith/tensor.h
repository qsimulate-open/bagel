//
// Author : Toru Shiozaki
// Date   : Feb 2012
//

#ifndef __SRC_SMITH_TENSOR_H
#define __SRC_SMITH_TENSOR_H

#include <list>
#include <map>
#include <memory>
#include <src/smith/storage.h>
#include <src/smith/indexrange.h>

namespace SMITH {

template <typename T>
class Tensor {
  protected:
   std::shared_ptr<T> data_; 

  public:
    Tensor(std::list<IndexRange> in) {
      // first compute hashtags and 
      std::map<size_t, size_t> hashmap;
      for (auto i = in.begin(); i != in.end(); ++i) {

      }

      std::shared_ptr<T> tmp(new T(hashmap));
      data_ = tmp;
    };
    ~Tensor() {};
};

}

#endif

