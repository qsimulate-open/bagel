//
// Author : Toru Shiozaki
// Date   : October 2009
//

#pragma once
#ifndef __src_util_pcompcabsfile_h
#define __src_util_pcompcabsfile_h

#include <src/util/pcompfile.h>

template<class T>
class PCompCABSFile : public PCompFile<T> {
  protected:
    void init_cabs_schwarz();
    std::vector<double> cabs_schwarz_;
    std::vector<boost::shared_ptr<Shell> > cabs_basis_;
    std::vector<int> cabs_offset_;
    void calculate_num_int_each();

  public:
    PCompCABSFile(boost::shared_ptr<PGeometry>, const bool late_init = false, const std::string jobname = "source");

    const double cabs_schwarz(int i) const { return cabs_schwarz_[i]; };
    std::vector<int> cabs_offset() const { return cabs_offset_; };
    int cabs_offset(size_t i) const { return cabs_offset_[i]; };

};


template<class T>
PCompCABSFile<T>::PCompCABSFile(boost::shared_ptr<PGeometry> pg, const bool late_init, const std::string jobname)
 : PCompFile<T>(pg, true, jobname) {

  { // prepare offset and basis
    typedef boost::shared_ptr<Atom> RefAtom;
    typedef boost::shared_ptr<Shell> RefShell;

    const std::vector<RefAtom> atoms = pg->cabs_atoms();
    int cnt = 0;
    for (std::vector<RefAtom>::const_iterator aiter = atoms.begin(); aiter != atoms.end(); ++aiter, ++cnt) {
      const std::vector<RefShell> tmp = (*aiter)->shells();
      cabs_basis_.insert(cabs_basis_.end(), tmp.begin(), tmp.end());  
      const std::vector<int> tmpoff = pg->cabs_offset(cnt); 
      cabs_offset_.insert(cabs_offset_.end(), tmpoff.begin(), tmpoff.end());
    }
  }
  if (!late_init) {
    this->init_schwarz();
    init_cabs_schwarz();
    calculate_num_int_each();
  }

};


template<class T>
void PCompCABSFile<T>::init_cabs_schwarz() {
  typedef boost::shared_ptr<Shell> RefShell;
  typedef boost::shared_ptr<Atom> RefAtom;

  const int size = this->basis_.size(); // the number of shells per unit cell
  const int cabs_size = cabs_basis_.size();

  cabs_schwarz_.resize(size * cabs_size * (2 * this->K_ + 1));

  #pragma omp prallel for
  for (int m = - this->K_; m <= this->K_; ++m) { 
    const double disp[3] = {0.0, 0.0, m * this->A_};
    for (int i0 = 0; i0 != size; ++i0) { // center unit cell
      const RefShell b0 = this->basis_[i0];
      for (int i1 = 0; i1 != cabs_size; ++i1) {
        const RefShell b1 = cabs_basis_[i1]->move_atom(disp);

        std::vector<RefShell> input;
        input.push_back(b0);
        input.push_back(b1);
        input.push_back(b0);
        input.push_back(b1);
        T batch(input, 1.0e100);
        batch.compute();
        const double* data = batch.data();
        const int datasize = batch.data_size();
        double cmax = 0.0;
        for (int xi = 0; xi != datasize; ++xi, ++data) {
          const double absed = (*data) > 0.0 ? *data : -*data;
          if (absed > cmax) cmax = absed;
        }
        cabs_schwarz_[(m + this->K_) * size * cabs_size + i0 * cabs_size + i1] = cmax;
      }
    }
  }
};


template<class T>
void PCompCABSFile<T>::calculate_num_int_each() {

  typedef boost::shared_ptr<Shell> RefShell;

  const int s = this->S_;
  const int l = this->L_;
  const int k = this->K_;
  const double a = this->A_;
  unsigned long data_written = 0ul;
  this->num_int_each_.resize((s + s + 1) * (s + s + 1) * (l + l + 1));

  const int size = this->basis_.size(); // number of shells
  const int cabs_size = cabs_basis_.size(); // number of shells

  #pragma omp parallel for reduction(+:data_written)
  for (int m1 = - s; m1 <= s; ++m1) {
    const double m1disp[3] = {0.0, 0.0, m1 * a}; 
    size_t offset = (m1 + s) * (l + 1) * (s * 2 + 1);
    for (int m2 = - l; m2 <= l; ++m2) { // use bra-ket symmetry!!!
      const double m2disp[3] = {0.0, 0.0, m2 * a}; 
      for (int m3 = m2 - s; m3 <= m2 + s; ++m3, ++offset) {
        const double m3disp[3] = {0.0, 0.0, m3 * a}; 
        size_t thisblock = 0ul; 
        for (int i0 = 0; i0 != size; ++i0) {
          const int b0offset = this->offset_[i0]; 
          const int b0size = this->basis_[i0]->nbasis();

          for (int i1 = 0; i1 != size; ++i1) {
            const int b1offset = this->offset_[i1];
            const int b1size = this->basis_[i1]->nbasis();

            for (int i2 = 0; i2 != size; ++i2) {
              const int b2offset = this->offset_[i2];
              const int b2size = this->basis_[i2]->nbasis();

              for (int i3 = 0; i3 != cabs_size; ++i3) {
                const int b3offset = cabs_offset_[i3]; 
                const int b3size = cabs_basis_[i3]->nbasis();

                const double integral_bound = this->schwarz_[(m1 + k) * size * size + i0 * size + i1]
                                            * cabs_schwarz_[(m3 - m2 + k) * size * cabs_size + i2 * cabs_size + i3];
                const bool skip_schwarz = integral_bound < SCHWARZ_THRESH;
                if (skip_schwarz) continue;
                data_written += b0size * b1size * b2size * b3size;
                thisblock += b0size * b1size * b2size * b3size; 
              }
            }
          }
        }
        this->num_int_each_[offset] = thisblock;
      }
    }
  }

  this->max_num_int_ = *std::max_element(this->num_int_each_.begin(), this->num_int_each_.end());

  std::cout << "  Using ";
  const size_t data_written_byte = data_written * sizeof(double);
  if (data_written_byte > 1.0e9) {
    std::cout << std::setprecision(1) << data_written_byte / 1.0e9 << " GB";
  } else if (data_written_byte > 1.0e6) {
    std::cout << std::setprecision(1) << data_written_byte / 1.0e6 << " MB";
  } else {
    std::cout << std::setprecision(1) << data_written_byte / 1.0e3 << " KB";
  }

  std::cout << " hard disk for storing \"" << this->jobname_ << "\"" << std::endl; 
  std::cout << std::endl;
  assert(data_written < 5.0e9); // 40GB
};

#endif

