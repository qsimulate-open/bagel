//
// Author : Toru Shiozaki
// Date   : Feb 2012
//

#ifndef __SRC_SMITH_TASK_H
#define __SRC_SMITH_TASK_H

namespace SMITH {

// base class for Task objects
// assumes that the operation table is static (not adjustable at runtime). 
template <typename T>
class Task {
  protected:
//  std::list<std::shared_ptr<Task> > depend_;

  public:
    Task() {}; 
    ~Task() { };
    virtual void compute() = 0;

    bool ready() const {
      bool out = true;
//    for (auto i = depend_.begin(); i != depend_.end(); ++i) out &= (*i)->ready();
      return out;
    };
};

}

#endif
