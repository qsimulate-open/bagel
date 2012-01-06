//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#ifndef __NEWINT_DF_DensityFit_H
#define __NEWINT_DF_DensityFit_H

#include <vector>
#include <memory>
#include <src/scf/atom.h>
#include <src/util/f77.h>

class DF_Half;
class DF_Full;

class DensityFit : public std::enable_shared_from_this<DensityFit> {
  protected:
    std::unique_ptr<double[]> data_;
    std::unique_ptr<double[]> data2_;
    size_t nbasis_;
    size_t naux_;

    double* data() { return data_.get(); };
    double* data2() { return data2_.get(); };
    const double* const data() const { return data_.get(); };
    const double* const data2() const { return data2_.get(); };

  public:
    DensityFit(const int nbas, const int naux,
       const std::vector<std::shared_ptr<Atom> >& atoms,  const std::vector<std::vector<int> >& offsets,
       const std::vector<std::shared_ptr<Atom> >& aux_atoms,  const std::vector<std::vector<int> >& aux_offsets, const double thr);
    ~DensityFit() {};

    const double* const data_3index() const { return data(); };
    const double* const data_2index() const { return data2(); };

    size_t nbasis() const { return nbasis_; };
    size_t naux() const { return naux_; };

    std::shared_ptr<DF_Half> compute_half_transform(const double* c, const size_t nocc);

}; 


class DF_Half {

  protected:
    std::unique_ptr<double[]> data_;
    const std::shared_ptr<DensityFit> df_;
    const int nocc_;

  public:
    DF_Half(const std::shared_ptr<DensityFit> df, const int nocc, std::unique_ptr<double[]>& in)
     : df_(df), nocc_(nocc), data_(std::move(in)) {}; 

    ~DF_Half() {};

    const double* const data() const { return data_.get(); };
    std::unique_ptr<double[]> move_data() { return std::move(data_); };

    std::shared_ptr<DF_Half> apply_J() {
      const int naux = df_->naux();
      const int nbasis = df_->nbasis();
      std::unique_ptr<double[]> out(new double[nocc_*naux*nbasis]);
      dgemm_("N", "N", naux, nocc_*nbasis, naux, 1.0, df_->data_2index(), naux, data_.get(), naux, 0.0, out.get(), naux); 
      std::shared_ptr<DF_Half> tmp(new DF_Half(df_, nocc_, out));
      return tmp; 
    } 

    std::shared_ptr<DF_Full> compute_second_transform(const double* c, const size_t nocc);
};


class DF_Full {

  protected:
    std::unique_ptr<double[]> data_;
    const std::shared_ptr<DensityFit> df_;
    const int nocc1_; // inner
    const int nocc2_; // outer

  public:
    DF_Full(const std::shared_ptr<DensityFit> df, const int nocc1, const int nocc2, std::unique_ptr<double[]>& in)
      : df_(df), nocc1_(nocc1), nocc2_(nocc2), data_(std::move(in)) {};

    ~DF_Full() {};

    std::shared_ptr<DF_Full> apply_J() {
      const int naux = df_->naux();
      const int nbasis = df_->nbasis();
      std::unique_ptr<double[]> out(new double[nocc1_*nocc2_*naux]);
      dgemm_("N", "N", naux, nocc1_*nocc2_, naux, 1.0, df_->data_2index(), naux, data_.get(), naux, 0.0, out.get(), naux); 
      std::shared_ptr<DF_Full> tmp(new DF_Full(df_, nocc1_, nocc2_, out));
      return tmp; 
    };

    void form_4index(std::unique_ptr<double[]>& target) const;

    const double* const data() const { return data_.get(); };
    const std::shared_ptr<DensityFit> df() const { return df_; };

};

#endif

