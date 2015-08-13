// Author: Hai-Anh Le
// Date  : August 2015

#include <src/util/constants.h>
#include <src/util/parallel/resources.h>

using namespace bagel;

class LocalExps {

  private:
    int ws_;
    int lmax_;
    int limit_;
    double thresh_;

    int rank_;
    double beta_;
    double *rvec_, *kvec_;
    double* T_;
    double* Rsq_;
    double* coeff_;
    double *roots_, *weights_;

    std::vector<std::complex<double>> mlm_;
    bool allocated_here_;
    std::shared_ptr<StackMem> stack_;
    size_t size_allocated_;
    double* buff_;

    void allocate_arrays(const size_t ps);
    void compute_mlm();

    double dot(std::array<double, 3> b, std::array<double, 3> c) { return b[0] * c[0] + b[1] * c[1] + b[2] * c[2]; }
    std::array<double, 3> cross(std::array<double, 3> b, std::array<double, 3> c, double s) {
      std::array<double, 3> out;
      out[0] = (b[1] * c[2] - b[2] * c[1]) * s;
      out[1] = (b[2] * c[0] - b[0] * c[2]) * s;
      out[2] = (b[0] * c[1] - b[1] * c[0]) * s;

      return out;
    }

    bool is_in_cff(const int ws, const int n0, const int n1, const int n2);
    void root_weight(const int size);

  public:
    LocalExps(const int ws, const int lmax, const int limit, const double thresh = PRIM_SCREEN_THRESH);
    ~LocalExps() {
      stack_->release(size_allocated_, buff_);
      resources__->release(stack_);
    }

   std::vector<std::complex<double>> mlm() const { return mlm_; }
   std::complex<double> mlm(const int i) const { return mlm_[i]; }
   double mlm_real(const int i) const { return mlm_[i].real(); }
};
