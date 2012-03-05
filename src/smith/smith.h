//
// Author : Toru Shiozaki
// Date   : Feb 2012
//

// compiles some input data for the smith routines.

namespace SMITH {

class SMITH_info {
  protected:
    int maxiter_;
    double thresh_residual_;

  public:
    SMITH_info(); 
    ~SMITH_info();

    int maxiter() const { return maxiter_; };
    double thresh_residual() const { return thresh_residual_; };
};

}
