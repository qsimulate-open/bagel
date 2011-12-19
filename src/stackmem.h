//
// Author : Toru Shiozaki
// Date   : Feb 2011
//

#ifndef __STACKMEM_H
#define __STACKMEM_H

// CAUTION last-in-first-out stack to avoid the overhead of new'ing every time

class StackMem {
  protected:
    double* stack_area_;
    size_t pointer_;

  public:
    StackMem(const size_t size) {
      stack_area_ = new double[size];
      pointer_ = 0lu;
    };

    ~StackMem() {
      delete[] stack_area_;
    }

    double* get(const size_t size) {
      double* out = stack_area_ + pointer_;
      pointer_ += size; 
      return out;
    };

    void release(const size_t size) {
      pointer_ -= size; 
    };
    
};


#endif
