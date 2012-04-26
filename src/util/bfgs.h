//
// Author : Toru Shiozaki
// Date   : April 2012
//

#ifndef __SRC_UTIL_BFGS_H
#define __SRC_UTIL_BFGS_H

// implements BFGS based on Fischer & Almlof JPC 1992.
// T needs clone, ddot and daxpy, along with overloaded operators and a copy constructor

#include <vector>
#include <memory>

template<typename T>
class BFGS {
  protected:
    std::vector<std::shared_ptr<const T> > delta; 
    std::vector<std::shared_ptr<const T> > y;
    std::vector<std::shared_ptr<const T> > D;

    std::shared_ptr<const T> prev_grad;
    std::shared_ptr<const T> prev_value;

    const std::shared_ptr<const T> denom_;

  public:
    BFGS(std::shared_ptr<T> denom) : denom_(denom) {};
    ~BFGS() {};

    std::shared_ptr<T> extrapolate(std::shared_ptr<T> _grad, std::shared_ptr<T> _value) {
      // to make sure, inputs are copied.
      std::shared_ptr<T> grad(new T(*_grad));
      std::shared_ptr<T> value(new T(*_value));

      std::shared_ptr<T> out(new T(*grad));
      assert(grad->size() == grad->ndim() * grad->mdim());
      // (1)
      for (size_t i = 0; i != grad->size(); ++i) out->data(i) /= denom_->data(i);

      if (prev_value) {
        // (3)
        std::shared_ptr<T> yy = grad->clone(); 
        {
          std::shared_ptr<T> DD(new T(*grad - *prev_grad));
          D.push_back(DD);

          for (size_t i = 0; i != grad->size(); ++i) yy->data(i) = DD->data(i) / denom_->data(i);

          std::shared_ptr<T> vv(new T(*value - *prev_value));
          delta.push_back(vv);

          prev_grad = grad;
          prev_value = value;
        }
        const int n = delta.size()-1;
        assert(delta.size() == y.size()+1 && y.size()+1 == D.size());
        
        // (4)
        for (int i = 0; i < n; ++i) {
          const double s1 = 1.0 / delta[i]->ddot(D[i]);
          const double s2 = 1.0 / D[i]->ddot(y[i]);
          const double s3 = delta[i]->ddot(grad);
          const double s4 =     y[i]->ddot(grad);
          const double s5 = delta[i]->ddot(D[n]);
          const double s6 =     y[i]->ddot(D[n]);
          const double t1 = (1.0 + s1/s2) * s1 * s3 - s1 * s4;
          const double t2 = s1 * s3;
          const double t3 = (1.0 + s1/s2) * s1 * s5 - s1 * s6;
          const double t4 = s1 * s5;
          out->daxpy(t1, delta[i]);
          out->daxpy(-t2, y[i]);
          yy->daxpy(t3, delta[i]);
          yy->daxpy(-t4, y[i]);
        }
        { // (5)
          const double s1 = 1.0 / delta[n]->ddot(D[n]);
          const double s2 = 1.0 /     D[n]->ddot(yy);
          const double s3 = delta[n]->ddot(grad);
          const double s4 =       yy->ddot(grad);
          const double t1 = (1.0 + s1/s2) * s1 * s3 - s1 * s4;
          const double t2 = s1 * s3;
          out->daxpy(t1, delta[n]);
          out->daxpy(-t2, yy);
        }
        y.push_back(yy);
      }
      prev_grad = grad; 
      prev_value = value;
      return out;
    };


};

#endif
