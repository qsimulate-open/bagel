//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#ifndef __NEWINT_SRC_UTIL_Linear_H
#define __NEWINT_SRC_UTIL_Linear_H

#include <memory>
#include <list>
#include <stdexcept>
#include <src/util/f77.h>

template<typename T>
class Linear {
  typedef std::shared_ptr<T> RefT;
  typedef std::list<std::pair<RefT, RefT> > Container_type_;
  typedef typename Container_type_::iterator iterator;

  protected:
    std::list<std::shared_ptr<T> > c_;
    std::list<std::shared_ptr<T> > sigma_;

    const int max_;    
    const std::shared_ptr<T> grad_;

    // contains 
    std::vector<double> mat_;
    std::vector<double> prod_;
    // scratch area for diagonalization
    std::vector<double> scr_;
    std::vector<double> vec_; 
    std::vector<int> ipiv_;
    int info;

    int size_;
    // for convenience below
    double& mat(int i, int j) { return mat_[i+j*max_]; };
    double& scr(int i, int j) { return scr_[i+j*max_]; };


  public:
    Linear(const int ndim, std::shared_ptr<T> grad) : max_(ndim), size_(0), grad_(grad) {
      mat_.resize(max_*max_);
      scr_.resize(max_*max_);
      vec_.resize(max_);
      prod_.resize(max_);
      ipiv_.resize(max_*2);
    };
    ~Linear() {};

    std::shared_ptr<T> compute_residual(std::shared_ptr<T> c, std::shared_ptr<T> s) {
      if (size_+2 == max_) throw std::runtime_error("max size reached in Linear");
      // register new vectors
      c_.push_back(c);
      sigma_.push_back(s);
      // first set mat (=x(i)Ax(j)) and prod (= x(i)*y)
      ++size_;
      auto citer = c_.begin();
      for (int i = 0; i != size_; ++i, ++citer) {
        mat(i,size_-1) = mat(size_-1,i) = s->ddot(**citer);
      } 
      prod_[size_-1] = -c->ddot(*grad_); 

      // set to scr_
      std::copy(mat_.begin(), mat_.end(), scr_.begin());
      std::copy(prod_.begin(), prod_.end(), vec_.begin());
      dgesv_(size_, 1, &scr_[0], max_, &ipiv_[0], &vec_[0], size_, info); 
      if (info) throw std::runtime_error("dsyev failed in Linear");

      std::shared_ptr<T> out(new T(*grad_)); 
      int cnt = 0;
      for (auto j = sigma_.begin(); j != sigma_.end(); ++j, ++cnt)
        out->daxpy(vec_[cnt], *j);
      assert(cnt == size_);
      return out;
    };

    std::shared_ptr<T> civec() const {
      std::shared_ptr<T> out = c_.front()->clone();
      int cnt = 0;
      for (auto i = c_.begin(); i != c_.end(); ++i, ++cnt) {
        out->daxpy(vec_[cnt], *i); 
      } 
      return out;
    };

    // make cc orthogonal to cc_ vectors
    double orthog(std::shared_ptr<T> cc) { return cc->orthog(c_); }

};

#endif
