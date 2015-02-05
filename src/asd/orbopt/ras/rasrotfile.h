
#ifndef __BAGEL_ASD_ORBOPT_RASROTFILE_H
#define __BAGEL_ASD_ORBOPT_RASROTFILE_H

#include <src/util/math/matrix.h>

namespace bagel {

template<typename DataType>
class ASDRASRotationMatrix {
  public:
    using data_type = DataType;

  protected:
    const int nclosed_;
    const int nact_;
    const int nvirt_;
    const int nactA_;
    const int nactB_;
    const int size_;
    std::unique_ptr<DataType[]> data_;

  public:
    ASDRASRotationMatrix(const int iclos, const int iact, const int ivirt, const int iactA, const int iactB)
     : nclosed_(iclos), nact_(iact), nvirt_(ivirt), nactA_(iactA), nactB_(iactB), size_(iclos*iact+iclos*ivirt+iact*ivirt+iactA*iactB), data_(new DataType[size_]) {
      zero();
    }
    ASDRASRotationMatrix(const ASDRASRotationMatrix& o) : nclosed_(o.nclosed_), nact_(o.nact_), nvirt_(o.nvirt_), nactA_(o.nactA_), nactB_(o.nactB_),  size_(o.size_), data_(new DataType[o.size_]) {
      *this = o;
    }
    ASDRASRotationMatrix(std::shared_ptr<const ASDRASRotationMatrix> o)
      : nclosed_(o->nclosed_), nact_(o->nact_), nvirt_(o->nvirt_), nactA_(o->nactA_), nactB_(o->nactB_), size_(o->size_), data_(new DataType[o->size_]) {
      *this = *o;
    }
    //Matrix to ASDRotationMatrix conversion
    ASDRASRotationMatrix(std::shared_ptr<const Matrix_base<DataType>> o, const int iclos, const int iact, const int ivirt, const int iactA, const int iactB)
      : nclosed_(iclos), nact_(iact), nvirt_(ivirt), nactA_(iactA), nactB_(iactB), size_(iclos*iact+iclos*ivirt+iact*ivirt+iactA*iactB), data_(new DataType[size_]) {
      const int nocc = nclosed_ + nact_;
      assert(iact == iactA+iactB);
      for (int i = 0; i != nact_; ++i) {
        for (int j = 0; j != nvirt_;   ++j) {
          ele_va(j, i) = o->element(j+nocc, i+nclosed_);
        }
        for (int j = 0; j != nclosed_; ++j) {
          ele_ca(j, i) = o->element(i+nclosed_, j);
        }
      }
      for (int i = 0; i != nclosed_; ++i) {
        for (int j = 0; j != nvirt_;   ++j) {
          ele_vc(j, i) = o->element(j+nocc, i);
        }
      }
      for (int i = 0; i != nactB_; ++i) {
        for (int j = 0; j != nactA_; ++j) {
          ele_aa(j, i) = o->element(j+nclosed_,i+nclosed_+nactA_);
        }
      }
    }

    std::shared_ptr<ASDRASRotationMatrix<DataType>> clone() const {
      return std::make_shared<ASDRASRotationMatrix<DataType>>(nclosed_, nact_, nvirt_, nactA_, nactB_);
    }
    std::shared_ptr<ASDRASRotationMatrix<DataType>> copy() const {
      return std::make_shared<ASDRASRotationMatrix<DataType>>(*this);
    }

    // overloaded operators
    ASDRASRotationMatrix<DataType> operator+(const ASDRASRotationMatrix<DataType>& o) const { ASDRASRotationMatrix<DataType> out(*this); out.ax_plus_y(1.0, o); return out; }
    ASDRASRotationMatrix<DataType> operator-(const ASDRASRotationMatrix<DataType>& o) const { ASDRASRotationMatrix<DataType> out(*this); out.ax_plus_y(-1.0, o); return out; }
    ASDRASRotationMatrix<DataType>& operator+=(const ASDRASRotationMatrix<DataType>& o) { ax_plus_y(1.0, o); return *this; }
    ASDRASRotationMatrix<DataType>& operator-=(const ASDRASRotationMatrix<DataType>& o) { ax_plus_y(-1.0, o); return *this; }
    ASDRASRotationMatrix<DataType>& operator*=(const DataType a) { scale(a); return *this; }
    ASDRASRotationMatrix<DataType>& operator/=(const ASDRASRotationMatrix<DataType>& o) { for (int i = 0; i != size(); ++i) data(i)/= o.data(i); return *this; }
    ASDRASRotationMatrix<DataType>& operator*=(const ASDRASRotationMatrix<DataType>& o) { for (int i = 0; i != size(); ++i) data(i)*= o.data(i); return *this; }
    ASDRASRotationMatrix<DataType> operator/(const ASDRASRotationMatrix<DataType>& o) const { ASDRASRotationMatrix<DataType> out(*this); return out /= o; }
    ASDRASRotationMatrix<DataType> operator*(const ASDRASRotationMatrix<DataType>& o) const { ASDRASRotationMatrix<DataType> out(*this); return out *= o; }
    ASDRASRotationMatrix<DataType>& operator=(const ASDRASRotationMatrix<DataType>& o) { std::copy_n(o.data(), size(), data());  return *this; }

    // size of the file
    int size() const { return size_; }
    // zero out
    void zero() { fill(DataType(0.0)); }
    void fill(const DataType& a) { std::fill_n(data(), size_, a); }
    // returns dot product
    DataType dot_product(const ASDRASRotationMatrix<DataType>& o) const { return blas::dot_product(data(), size_, o.data()); }
    DataType dot_product(const std::shared_ptr<const ASDRASRotationMatrix<DataType>> o) const { return dot_product(*o); }
    // scale function
    void scale(const DataType& a) { std::for_each(data(), data()+size_, [&a](DataType& p) { p *= a; }); }
    // returns norm of the vector
    double norm() const { return std::sqrt(detail::real(dot_product(*this))); }
    // returns a root mean square
    double rms() const { return norm() / std::sqrt(static_cast<double>(size())); }
    // daxpy added to self
    void ax_plus_y(const DataType& a, const ASDRASRotationMatrix& o) { std::transform(o.data(), o.data()+size_, data(), data(), [&a](DataType p, DataType q) { return p*a+q; }); }
    void ax_plus_y(const DataType& a, const std::shared_ptr<const ASDRASRotationMatrix> o) { ax_plus_y(a, *o); }

    // orthogonalize to the liset of ASDRASRotationMatrix's
    double orthog(std::list<std::shared_ptr<const ASDRASRotationMatrix<DataType>>> c) {
      for (auto iter = c.begin(); iter != c.end(); ++iter)
        this->ax_plus_y(- this->dot_product(**iter), **iter);
      return normalize();
    }

    double normalize() {
      const double scal = 1.0/this->norm();
      scale(scal);
      return 1.0/scal;
    }

    void synchronize() {
#ifdef HAVE_MPI_H
      mpi__->broadcast(data(), size(), 0);
#endif
    }

    // return data_
    DataType* data() { return data_.get(); }
    const DataType* data() const { return data_.get(); }
    DataType& data(const size_t i) { return data_[i]; }
    const DataType& data(const size_t i) const { return data_[i]; }
    // return data_
    DataType* begin() { return data(); }
    // return data_
    DataType* end() { return data()+size_; }

    // ROTATION MATRIX
    //   C   A   V        A(A) A(B)   ORDERED AS:
    // C .  ca  cv    A(A)  .  [aa]   ca -> va -> vc -> aa
    // A . [aa] av    A(B)  .    .
    // V .   .   .

    // closed-active block. closed runs first
    DataType* ptr_ca() { return data(); }
    DataType& ele_ca(const int ic, const int ia) { return data_[ic + ia*nclosed_]; }
    // active-virtual block. virtual runs first
    DataType* ptr_va() { return data() + nclosed_*nact_; }
    DataType& ele_va(const int iv, const int ia) { return data_[nclosed_*nact_ + iv + ia*nvirt_]; }
    // closed-virtual block. virtual runs first
    DataType* ptr_vc() { return data() + (nclosed_+nvirt_)*nact_; }
    DataType& ele_vc(const int iv, const int ic) { return data_[(nclosed_+nvirt_)*nact_ + iv + ic*nvirt_]; }
    // active-active block. A index runs first
    DataType* ptr_aa() { return data() + (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_; }
    DataType& ele_aa(const int ia, const int ib) { return data_[(nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + ia + ib*nactA_]; }

    // const references and pointers
    const DataType& ele_ca(const int ic, const int ia) const { return data_[ic + ia*nclosed_]; }
    const DataType& ele_va(const int iv, const int ia) const { return data_[nclosed_*nact_ + iv + ia*nvirt_]; }
    const DataType& ele_vc(const int iv, const int ic) const { return data_[(nclosed_+nvirt_)*nact_ + iv + ic*nvirt_]; }
    const DataType& ele_aa(const int ia, const int ib) const { return data_[(nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + ia + ib*nactA_]; }

    const DataType* ptr_ca() const { return data(); }
    const DataType* ptr_va() const { return data() + nclosed_*nact_; }
    const DataType* ptr_vc() const { return data() + (nclosed_+nvirt_)*nact_; }
    const DataType* ptr_aa() const { return data() + (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_; }

    // unpack to Matrix
    template<class MatType>
    std::shared_ptr<MatType> unpack(const DataType a = 0.0) const {
      const int nocc = nclosed_ + nact_;
      const int nbasis = nclosed_ + nact_ + nvirt_;
      auto out = std::make_shared<MatType>(nbasis, nbasis);
      std::fill_n(out->data(), out->size(), a);

      for (int i = 0; i != nact_; ++i) {
        //Virtual-Active
        for (int j = 0; j != nvirt_;   ++j) {
          out->element(j+nocc, i+nclosed_) = ele_va(j, i);
        }
        //Active-Closed
        for (int j = 0; j != nclosed_; ++j) {
          out->element(i+nclosed_, j) = ele_ca(j, i);
        }
      }
      //Virtual-Closed
      for (int i = 0; i != nclosed_; ++i) {
        for (int j = 0; j != nvirt_;   ++j) {
          out->element(j+nocc, i) = ele_vc(j, i);
        }
      }
      //Active-Active 
      for (int i = 0; i != nactB_; ++i) {
        for (int j = 0; j != nactA_; ++j) {
        //out->element(j+nclosed_, i+nclosed_+nactA_) = ele_aa(j, i);
          out->element(i+nclosed_+nactA_,j+nclosed_) = -ele_aa(j,i);
        //out->element(i+nclosed_+nactA_,j+nclosed_) = ele_aa(j,i);
        }
      }
      //Anti-symmetric
      for (int i = 0; i != nbasis; ++i) {
        for (int j = 0; j <= i; ++j) {
          out->element(j, i) = - detail::conj(out->element(i, j));
        }
      }

      return out;
    }

    // print matrix
    void print(const std::string in = "") const {
      std::cout << "++++ " + in + " ++++" << std::endl;
      if (nact_ && nclosed_) {
        std::cout << " printing closed-active block" << std::endl;
        for (int i = 0; i != nact_; ++i) {
          for (int j = 0; j != nclosed_; ++j) {
            std::cout << std::setw(10) << std::setprecision(6) << ele_ca(j,i);
          }
          std::cout << std::endl;
        }
      }
      if (nact_ && nvirt_) {
        std::cout << " printing virtual-active block" << std::endl;
        for (int i = 0; i != nact_; ++i) {
          for (int j = 0; j != nvirt_; ++j) {
            std::cout << std::setw(10) << std::setprecision(6) << ele_va(j,i);
          }
          std::cout << std::endl;
        }
      }
      if (nclosed_ && nvirt_) {
        std::cout << " printing virtual-closed block" << std::endl;
        for (int i = 0; i != nclosed_; ++i) {
          for (int j = 0; j != nvirt_; ++j) {
            std::cout << std::setw(10) << std::setprecision(6) << ele_vc(j,i);
          }
          std::cout << std::endl;
        }
      }
      if (nactA_ && nactB_) {
        std::cout << " printing active(A)-active(B) block" << std::endl;
        for (int i = 0; i != nactB_; ++i) {
          for (int j = 0; j != nactA_; ++j) {
            std::cout << std::setw(10) << std::setprecision(6) << ele_aa(j,i);
          }
          std::cout << std::endl;
        }
      }

    }
};

using ASDRASRotFile = ASDRASRotationMatrix<double>;

}

#endif
