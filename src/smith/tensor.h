//
// Author : Toru Shiozaki
// Date   : Feb 2012
//

#ifndef __SRC_SMITH_TENSOR_H
#define __SRC_SMITH_TENSOR_H

#include <src/smith/storage.h>

namespace SMITH {

class Tensor {
  protected:
   std::shared_ptr<Storage_base> data_; 

  public:
    Tensor();
    ~Tensor() {};


};

}

#endif

