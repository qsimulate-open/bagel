//
// BAGEL - Parallel electron correlation program.
// Filename: pcompcabsfile.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#ifndef __src_util_pcompcabsfile_h
#define __src_util_pcompcabsfile_h

#include <vector>
#include <src/pscf/pcompfile.h>

namespace bagel {

template<class T>
class PCompCABSFile : public PCompFile<T> {
  protected:
    void init_schwarz_jb();
    void init_schwarz_ia();
    std::vector<double> schwarz_jb_;
    std::vector<double> schwarz_ia_;
    std::vector<std::shared_ptr<const Shell>> cabs_basis_;
    std::vector<int> aux_offset_;

    const bool i_is_cabs_;
    const bool j_is_cabs_;
    const bool a_is_cabs_;
    const bool b_is_cabs_;

    std::vector<int> offset_i_;
    std::vector<int> offset_j_;
    std::vector<int> offset_a_;
    std::vector<int> offset_b_;

    std::vector<std::shared_ptr<const Shell>> basis_i_;
    std::vector<std::shared_ptr<const Shell>> basis_j_;
    std::vector<std::shared_ptr<const Shell>> basis_a_;
    std::vector<std::shared_ptr<const Shell>> basis_b_;

    int size_i_;
    int size_j_;
    int size_a_;
    int size_b_;

    int nbasis_i_;
    int nbasis_j_;
    int nbasis_a_;
    int nbasis_b_;

  public:
    PCompCABSFile(std::shared_ptr<PGeometry> geom, const double gamma,
        const bool i_is_cabs, const bool j_is_cabs, const bool a_is_cabs, const bool b_is_cabs,
        const bool late_init = false, const std::string jobname = "source");

    int offset_i(const size_t i) const { return offset_i_[i]; };
    int offset_j(const size_t j) const { return offset_j_[j]; };
    int offset_a(const size_t a) const { return offset_a_[a]; };
    int offset_b(const size_t b) const { return offset_b_[b]; };

    const std::shared_ptr<const Shell> basis_i(const size_t i) const { return basis_i_[i]; };
    const std::shared_ptr<const Shell> basis_j(const size_t j) const { return basis_j_[j]; };
    const std::shared_ptr<const Shell> basis_a(const size_t a) const { return basis_a_[a]; };
    const std::shared_ptr<const Shell> basis_b(const size_t b) const { return basis_b_[b]; };

    int size_i() const { return size_i_; };
    int size_j() const { return size_j_; };
    int size_a() const { return size_a_; };
    int size_b() const { return size_b_; };

    int nbasis_i() const { return nbasis_i_; };
    int nbasis_j() const { return nbasis_j_; };
    int nbasis_a() const { return nbasis_a_; };
    int nbasis_b() const { return nbasis_b_; };

    double schwarz_ia(const size_t i) const { return schwarz_ia_[i]; };
    double schwarz_jb(const size_t i) const { return schwarz_jb_[i]; };

    std::vector<int> aux_offset() const { return aux_offset_; };
    int aux_offset(size_t i) const { return aux_offset_[i]; };
    size_t cabs_nbasis(size_t i) const { return cabs_basis_[i]->nbasis(); };

    // virtual functions
    void calculate_num_int_each();
    void store_integrals();
    void eval_new_block(double*, int, int, int);

    std::shared_ptr<PMOFile<std::complex<double>>>
      mo_transform_cabs_aux(std::shared_ptr<PCoeff>,
                            std::shared_ptr<PCoeff>,
                            std::shared_ptr<PCoeff>,
                            std::shared_ptr<PCoeff>,
                            const int istart, const int ifence,
                            const int jstart, const int jfence,
                            const int astart, const int afence,
                            const int bstart, const int bfence,
                            const std::string jobname = "intermediate",
                            const bool direct = true);

};


template<class T>
PCompCABSFile<T>::PCompCABSFile(std::shared_ptr<PGeometry> pg, const double gam,
    const bool i_c, const bool j_c, const bool a_c, const bool b_c,
    const bool late_init, const std::string jobname)
 : PCompFile<T>(pg, gam, true, jobname), i_is_cabs_(i_c), j_is_cabs_(j_c), a_is_cabs_(a_c), b_is_cabs_(b_c) {

  { // prepare offset and basis
    typedef std::shared_ptr<const Atom> RefAtom;
    typedef std::shared_ptr<const Shell> RefShell;

    const std::vector<RefAtom> atoms = pg->aux_atoms();
    int cnt = 0;
    for (std::vector<RefAtom>::const_iterator aiter = atoms.begin(); aiter != atoms.end(); ++aiter, ++cnt) {
      const std::vector<RefShell> tmp = (*aiter)->shells();
      cabs_basis_.insert(cabs_basis_.end(), tmp.begin(), tmp.end());
      const std::vector<int> tmpoff = pg->aux_offset(cnt);
      aux_offset_.insert(aux_offset_.end(), tmpoff.begin(), tmpoff.end());
    }
  }

  {
    const std::vector<int> tmpi = i_is_cabs_ ? aux_offset_ : this->offset_;
    const std::vector<int> tmpj = j_is_cabs_ ? aux_offset_ : this->offset_;
    const std::vector<int> tmpa = a_is_cabs_ ? aux_offset_ : this->offset_;
    const std::vector<int> tmpb = b_is_cabs_ ? aux_offset_ : this->offset_;
    offset_i_.insert(offset_i_.end(), tmpi.begin(), tmpi.end());
    offset_j_.insert(offset_j_.end(), tmpj.begin(), tmpj.end());
    offset_a_.insert(offset_a_.end(), tmpa.begin(), tmpa.end());
    offset_b_.insert(offset_b_.end(), tmpb.begin(), tmpb.end());
  }
  {
    const std::vector<std::shared_ptr<const Shell>> tmpi = i_is_cabs_ ? cabs_basis_ : this->basis_;
    const std::vector<std::shared_ptr<const Shell>> tmpj = j_is_cabs_ ? cabs_basis_ : this->basis_;
    const std::vector<std::shared_ptr<const Shell>> tmpa = a_is_cabs_ ? cabs_basis_ : this->basis_;
    const std::vector<std::shared_ptr<const Shell>> tmpb = b_is_cabs_ ? cabs_basis_ : this->basis_;
    basis_i_.insert(basis_i_.end(), tmpi.begin(), tmpi.end());
    basis_j_.insert(basis_j_.end(), tmpj.begin(), tmpj.end());
    basis_a_.insert(basis_a_.end(), tmpa.begin(), tmpa.end());
    basis_b_.insert(basis_b_.end(), tmpb.begin(), tmpb.end());
    size_i_ = basis_i_.size();
    size_j_ = basis_j_.size();
    size_a_ = basis_a_.size();
    size_b_ = basis_b_.size();
    nbasis_i_ = i_is_cabs_ ? this->geom_->naux() : this->geom_->nbasis();
    nbasis_j_ = j_is_cabs_ ? this->geom_->naux() : this->geom_->nbasis();
    nbasis_a_ = a_is_cabs_ ? this->geom_->naux() : this->geom_->nbasis();
    nbasis_b_ = b_is_cabs_ ? this->geom_->naux() : this->geom_->nbasis();
  }

  if (!late_init) {
    init_schwarz_ia();
    init_schwarz_jb();
    calculate_num_int_each();
  }

};


template<class T>
void PCompCABSFile<T>::init_schwarz_jb() {
  typedef std::shared_ptr<const Shell> RefShell;
  typedef std::shared_ptr<const Atom> RefAtom;

  const int size = this->basis_.size(); // the number of shells per unit cell

  schwarz_jb_.resize(size_b_ * size_j_ * (2 * this->K_ + 1));

//#pragma omp parallel for
  for (int m = - this->K_; m <= this->K_; ++m) {
    const double disp[3] = {0.0, 0.0, m * this->A_};
    for (int i0 = 0; i0 != size_j_; ++i0) { // center unit cell
      const RefShell b0 = basis_j_[i0];
      for (int i1 = 0; i1 != size_b_; ++i1) {
        const RefShell b1 = basis_b_[i1]->move_atom(disp);

        std::array<RefShell,4> input = {{ b0, b1, b0, b1 }};
        T batch(input, 1.0, this->gamma_, false);
        batch.compute();
        const double* data = batch.data();
        const int datasize = batch.data_size();
        double cmax = 0.0;
        for (int xi = 0; xi != datasize; ++xi, ++data) {
          const double absed = std::fabs(*data);
          if (absed > cmax) cmax = absed;
        }
        schwarz_jb_[((m + this->K_) * size_j_ + i0) * size_b_ + i1] = std::sqrt(cmax);
      }
    }
  }
};


template<class T>
void PCompCABSFile<T>::init_schwarz_ia() {
  typedef std::shared_ptr<const Shell> RefShell;
  typedef std::shared_ptr<const Atom> RefAtom;

  schwarz_ia_.resize(size_a_ * size_i_ * (2 * this->K_ + 1));

//#pragma omp parallel for
  for (int m = - this->K_; m <= this->K_; ++m) {
    const double disp[3] = {0.0, 0.0, m * this->A_};
    for (int i0 = 0; i0 != size_i_; ++i0) { // center unit cell
      const RefShell b0 = basis_i_[i0];
      for (int i1 = 0; i1 != size_a_; ++i1) {
        const RefShell b1 = basis_a_[i1]->move_atom(disp);

        std::array<RefShell,4> input = {{ b0, b1, b0, b1}};
        T batch(input, 1.0, this->gamma_, false);
        batch.compute();
        const double* data = batch.data();
        const int datasize = batch.data_size();
        double cmax = 0.0;
        for (int xi = 0; xi != datasize; ++xi, ++data) {
          const double absed = std::fabs(*data);
          if (absed > cmax) cmax = absed;
        }
        schwarz_ia_[((m + this->K_) * size_i_ + i0) * size_a_ + i1] = std::sqrt(cmax);
      }
    }
  }
};


template<class T>
void PCompCABSFile<T>::calculate_num_int_each() {

  typedef std::shared_ptr<const Shell> RefShell;

  const int s = this->S_;
  const int l = this->L_;
  const int k = this->K_;
  const double a = this->A_;
  unsigned long data_written = 0ul;
  this->num_int_each_.resize((s+s+1) * (s+s+1) * (l+l+1));

//#pragma omp parallel for reduction (+:data_written)
  for (int m1 = - s; m1 <= s; ++m1) {
    const double m1disp[3] = {0.0, 0.0, m1*a};
    size_t offset = (m1+s) * (l*2+1) * (s*2+1);
    for (int m2 = - l; m2 <= l; ++m2) { // NO bra-ket symmetry!!!
      const double m2disp[3] = {0.0, 0.0, m2*a};
      for (int m3 = m2 - s; m3 <= m2 + s; ++m3, ++offset) {
        const double m3disp[3] = {0.0, 0.0, m3*a};
        size_t thisblock = 0ul;
        for (int i0 = 0; i0 != size_i_; ++i0) {
          const int b0offset = offset_i_[i0];
          const int b0size = basis_i_[i0]->nbasis();

          for (int i1 = 0; i1 != size_a_; ++i1) {
            const int b1offset = offset_a_[i1];
            const int b1size = basis_a_[i1]->nbasis();

            for (int i2 = 0; i2 != size_j_; ++i2) {
              const int b2offset = offset_j_[i2];
              const int b2size = basis_j_[i2]->nbasis();

              for (int i3 = 0; i3 != size_b_; ++i3) {
                const int b3offset = offset_b_[i3];
                const int b3size = basis_b_[i3]->nbasis();

                const double integral_bound = schwarz_ia_[((m1 + k) * size_i_ + i0) * size_a_ + i1]
                                            * schwarz_jb_[((m3 - m2 + k) * size_j_ + i2) * size_b_ + i3];
                const bool skip_schwarz = integral_bound < schwarz_thresh__;
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

  /*
  if (this->jobname_ != "NULL") {
    std::cout << std::fixed << "  Using ";
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
  }
  */

};


template<class T>
void PCompCABSFile<T>::store_integrals() {
  // TODO control the cache size.
  const size_t cachesize_max = 120000000lu;
  const size_t cachesize = std::max(cachesize_max, this->max_num_int_);

  const int s = this->S_;
  const int l = this->L_;
  const int k = this->K_;

  double* dcache = new double[cachesize];

  size_t remaining = cachesize;
  size_t current = 0lu;
  size_t cnt = 0lu;
  for (int m1 = -s; m1 <= s; ++m1) {
    for (int m2 = -l; m2 <= l; ++m2) { // NO bra-ket symmetry owing to CABS index!!!
      for (int m3 = m2 - s; m3 <= m2 + s; ++m3, ++cnt) {
        const size_t blocksize = this->num_int_each(cnt);;
        if (remaining < blocksize) {
          this->append(current, dcache);
          current = 0lu;
          remaining = cachesize;
        }
        eval_new_block(dcache+current, m1, m2, m3);
        current += blocksize;
        remaining -= blocksize;
      }
    }
  }
  this->append(current, dcache);
  delete[] dcache;
};


template<class T>
void PCompCABSFile<T>::eval_new_block(double* out, int m1, int m2, int m3) {

  typedef std::shared_ptr<const Shell> RefShell;

  const int s = this->S_;
  const int l = this->L_;
  const int k = this->K_;
  const double a = this->A_;

  const double m1disp[3] = {0.0, 0.0, m1 * a};
  const double m2disp[3] = {0.0, 0.0, m2 * a};
  const double m3disp[3] = {0.0, 0.0, m3 * a};

  size_t* blocks = new size_t[size_i_*size_j_*size_a_*size_b_+1];
  blocks[0] = 0;
  int iall = 0;
  for (int i0 = 0; i0 != size_i_; ++i0) {
    const int b0size = basis_i_[i0]->nbasis();
    for (int i1 = 0; i1 != size_a_; ++i1) {
      const int b1size = basis_a_[i1]->nbasis();
      for (int i2 = 0; i2 != size_j_; ++i2) {
        const int b2size = basis_j_[i2]->nbasis();
        for (int i3 = 0; i3 != size_b_; ++i3, ++iall) {
          const int b3size = basis_b_[i3]->nbasis();
          const double integral_bound = schwarz_ia_[((m1 + k) * size_i_ + i0) * size_a_ + i1]
                                      * schwarz_jb_[((m3 - m2 + k) * size_j_ + i2) * size_b_ + i3];
          const bool skip_schwarz = integral_bound < schwarz_thresh__;
          blocks[iall + 1] = blocks[iall] + (skip_schwarz ? 0 : (b0size * b1size * b2size * b3size));
        }
      }
    }
  }
//#pragma omp parallel for
  for (int i01 = 0; i01 < size_i_ * size_a_; ++i01) {
    const int i1 = i01 % size_a_;
    const int i0 = (i01 - i1) / size_a_;
    {
      int offset = i01 * size_j_ * size_b_;
      const RefShell b0 = basis_i_[i0]; // b0 is the center cell
      const RefShell b1 = basis_a_[i1]->move_atom(m1disp);
      for (int i2 = 0; i2 != size_j_; ++i2) {
        const RefShell b2 = basis_j_[i2]->move_atom(m2disp);
        for (int i3 = 0; i3 != size_b_; ++i3, ++offset) {
          const RefShell b3 = basis_b_[i3]->move_atom(m3disp);

          if (blocks[offset] == blocks[offset + 1]) continue;

          std::array<RefShell,4> input = {{ b3, b2, b1, b0 }};

          T batch(input, 1.0, this->gamma_, false);
          batch.compute();
          const double* bdata = batch.data();
          ::memcpy(out + blocks[offset], bdata, (blocks[offset + 1] - blocks[offset]) * sizeof(double));
        }
      }
    }
  }
  delete[] blocks;
};


template<class T>
std::shared_ptr<PMOFile<std::complex<double>>>
  PCompCABSFile<T>::mo_transform_cabs_aux(std::shared_ptr<PCoeff> coeff_i,
                                          std::shared_ptr<PCoeff> coeff_j,
                                          std::shared_ptr<PCoeff> coeff_a,
                                          std::shared_ptr<PCoeff> coeff_b,
                                          const int istart, const int ifence,
                                          const int jstart, const int jfence,
                                          const int astart, const int afence,
                                          const int bstart, const int bfence,
                                          const std::string jobname,
                                          const bool direct) {

  // What is different is that the coefficient of b3 is replaced by cabs_coeff.
  // Other than that, they should be the same as mo_transform.

  // Loading a (2K * 2K * nov) quantity on memory

  assert(coeff_i->ndim() == nbasis_i_);
  assert(coeff_j->ndim() == nbasis_j_);
  assert(coeff_a->ndim() == nbasis_a_);
  assert(coeff_b->ndim() == nbasis_b_);

  assert(bfence <= coeff_b->mdim());

  const int isize = ifence - istart;
  const int jsize = jfence - jstart;
  const int asize = afence - astart;
  const int bsize = bfence - bstart;

  const size_t noovv = static_cast<size_t>(isize) * jsize * asize * bsize;
  assert(noovv > 0);

  const double pi = 3.14159265358979323846;
  const int k = this->K_;
  const int s = this->S_;
  const int l = this->L_;

  const int KK = k + k;

  const int maxK1 = std::max(k, 1);
  const std::complex<double> czero(0.0, 0.0);
  const std::complex<double> cone(1.0, 0.0);

  const int unit = 1;

  const size_t filesize = noovv * std::max(KK, 1) * std::max(KK, 1) * std::max(KK, 1);
  std::cout << "  Creating " << jobname << "  of size ";
  const size_t filesize_byte = filesize * sizeof(std::complex<double>);
  if (filesize_byte > 1.0e9) {
    std::cout << std::setprecision(1) << filesize_byte / 1.0e9 << " GB" << std::endl;
  } else if (filesize_byte > 1.0e6) {
    std::cout << std::setprecision(1) << filesize_byte / 1.0e6 << " MB" << std::endl;
  } else {
    std::cout << std::setprecision(1) << filesize_byte / 1.0e3 << " KB" << std::endl;
  }
  std::shared_ptr<PMOFile<std::complex<double>>>
    mo_int(new PMOFile<std::complex<double>>(this->geom_, filesize, k,
                                              istart, ifence, jstart, jfence,
                                              astart, afence, bstart, bfence, true));

  // we are assuming that the (c.a. two-electron integrals for a unit cell)*K^2 can be
  // held in core. If that is not the case, this must be rewritten.

  const size_t nbasis4 = nbasis_a_ * nbasis_i_ * nbasis_b_ * nbasis_j_;
  const size_t nv = nbasis_i_ * nbasis_a_ * nbasis_j_ * bsize;
  const size_t nov = nbasis_i_ * nbasis_a_ * jsize * bsize;
  const size_t novv = nbasis_i_ * jsize * asize * bsize;

  const size_t nmax = std::max(nbasis4, std::max(std::max(std::max(nv, nov), novv), noovv));

  // allocating a temp array
  const size_t alloc = std::max(nmax, std::max(std::max(nov,novv), noovv) * std::max(KK, 1));
  std::complex<double>* data = new std::complex<double>[nmax];
  std::complex<double>* datas = new std::complex<double>[alloc];
  double* data_read = new double[nbasis4];

  size_t* blocks = new size_t[size_i_*size_j_*size_a_*size_b_ + 1];

  const size_t sizem1 = s*2+1lu;
  const size_t sizem2 = l*2+1lu;
  const size_t num_loops = sizem1 * sizem1 * sizem2;
  size_t loop_counter = 0lu;
  size_t loop_mod10 = 0lu;

  size_t allocsize = *std::max_element(this->num_int_each_.begin(), this->num_int_each_.end());

  PFile<std::complex<double>> intermediate_mmK(std::max(KK*nv, nv), k, false);
  std::complex<double>* intermediate_novv = new std::complex<double>[novv];
  PFile<std::complex<double>> intermediate_mKK(std::max(nov * KK * KK, nov), k, false);
  PFile<std::complex<double>> intermediate_KKK(std::max(novv * KK * KK * KK, novv), k, true);

  for (int q1 = -s; q1 <= s; ++q1) {
    const bool q1_front = (q1 == -s);
    const int m1 = q1;
    intermediate_mKK.clear();
    for (int q2 = -l; q2 <= l; ++q2) {
      const int m2 = q2;
      intermediate_mmK.clear();
      for (int q3 = -s; q3 <= s; ++q3, ++loop_counter) {
        const int m3 = m2 + q3;

        const double* cdata;
        if (!direct) {
          const size_t key = q3 + s + sizem1 * (q2 + l + sizem2 * (q1 + s));
          size_t datasize_acc = 0lu;
          for (size_t i = 0; i != key; ++i) datasize_acc += this->num_int_each_[i];
          this->get_block(datasize_acc, this->num_int_each_[key], data_read);
          cdata = data_read;
        } else {
          this->eval_new_block(data_read, m1, m2, m3);
          cdata = data_read;
        }

        if (loop_mod10 < loop_counter * 10lu / num_loops) {
          loop_mod10++;
          std::cout << "  Loop " << loop_mod10 * 10lu << " percent done" << std::endl;
        }

        blocks[0] = 0;
        int iall = 0;
        for (int i0 = 0; i0 != size_i_; ++i0) {
          const int b0size = basis_i_[i0]->nbasis();
          for (int i1 = 0; i1 != size_a_; ++i1) {
            const int b1size = basis_a_[i1]->nbasis();
            for (int i2 = 0; i2 != size_j_; ++i2) {
              const int b2size = basis_j_[i2]->nbasis();
              for (int i3 = 0; i3 != size_b_; ++i3, ++iall) {
                const int b3size = basis_b_[i3]->nbasis();

                double integral_bound = schwarz_ia_[((q1 + k) * size_i_ + i0) * size_a_ + i1]
                                      * schwarz_jb_[((q3 + k) * size_j_ + i2) * size_b_ + i3];
                const bool skip_schwarz = integral_bound < schwarz_thresh__;
                blocks[iall + 1] = blocks[iall] + (skip_schwarz ? 0 : (b0size * b1size * b2size * b3size));
              }
            }
          }
        }

//      #pragma omp parallel for
        for (int i0 = 0; i0 < size_i_; ++i0) {
          int noffset = i0 * size_j_ * size_a_ * size_b_;
          const int b0offset = offset_i_[i0];
          const int b0size = basis_i_[i0]->nbasis();
          for (int i1 = 0; i1 != size_a_; ++i1) {
            const int b1offset = offset_a_[i1];
            const int b1size = basis_a_[i1]->nbasis();
            for (int i2 = 0; i2 != size_j_; ++i2) {
              const int b2offset = offset_j_[i2];
              const int b2size = basis_j_[i2]->nbasis();
              for (int i3 = 0; i3 != size_b_; ++i3, ++noffset) {
                const int b3offset = offset_b_[i3];
                const int b3size = basis_b_[i3]->nbasis();

                const size_t cn3 = nbasis_a_ * nbasis_j_ * nbasis_b_;
                const size_t cn2 = nbasis_j_ * nbasis_b_;
                const size_t cn1 = nbasis_b_;
                if (blocks[noffset] != blocks[noffset + 1]) {
                  const double* ndata = cdata + blocks[noffset];
                  std::complex<double>* label = data + b3offset + nbasis_b_
                      * (b2offset + nbasis_j_ * (b1offset + nbasis_a_ * b0offset));
                  for (int j0 = 0; j0 != b0size; ++j0, label += cn3) {
                    for (int j1 = 0; j1 != b1size; ++j1, label += cn2) {
                      for (int j2 = 0; j2 != b2size; ++j2, label += cn1) {
                        for (int j3 = 0; j3 != b3size; ++j3, ++label, ++ndata) {
                          *label = static_cast<std::complex<double>>(*ndata);
                        }
                        label -= b3size;
                      }
                      label -= b2size * cn1;
                    }
                    label -= b1size * cn2;
                  }
                } else {
                  std::complex<double>* label = data + b3offset + nbasis_b_
                      * (b2offset + nbasis_j_ * (b1offset + nbasis_a_ * b0offset));
                  for (int j0 = 0; j0 != b0size; ++j0, label += cn3) {
                    for (int j1 = 0; j1 != b1size; ++j1, label += cn2) {
                      for (int j2 = 0; j2 != b2size; ++j2, label += cn1) {
                        std::fill(label, label + b3size, czero);
                      }
                      label -= b2size * cn1;
                    }
                    label -= b1size * cn2;
                  }
                } // end of if skip_schwarz

              }
            }
          }
        } // end of shell loops; now data is ready

        if (true) { // start from j & b contractions -- needs transposition
          const int mn = nbasis_i_ * nbasis_a_;
          const int mnc = nbasis_j_ * nbasis_b_;
          blas::transpose(data, mnc, mn, datas);
          ::memcpy(data, datas, nbasis4 * sizeof(std::complex<double>));
        }

        for (int nkb = -k; nkb < maxK1; ++nkb) {
          const int nkbc = nkb + k;
          const double img = k != 0 ? ((m3 * nkb * pi) / k) : 0.0;
          const std::complex<double> exponent(0.0, img);
          const std::complex<double> prefac = exp(exponent);
          size_t offset1 = 0;
          size_t offset2 = 0;
          const int cn2 = nbasis_i_ * nbasis_a_;
          std::fill(datas, datas+nv, czero);
//        #pragma omp parallel for
          for (int ii = 0; ii < nbasis_j_; ++ii) {
            const int offset1 = cn2 * nbasis_b_ * ii;
            const int offset2 = cn2 * bsize * ii;
            zgemm3m_("N", "N", &cn2, &bsize, &nbasis_b_, &prefac, data + offset1, &cn2,
                                           coeff_b->bp(nkb) + nbasis_b_ * bstart, &nbasis_b_, &cone,
                                           datas + offset2, &cn2);
          }
          intermediate_mmK.add_block(nkbc*nv, nv, datas);
        } // end of contraction b for given m3

      } // end of m3 loop

      // intermediate_mmK is ready

      for (int nkb = -k; nkb < maxK1; ++nkb) {
        const int nkbc = nkb + k;
        int nbj = nkbc * KK;
        intermediate_mmK.get_block(nkbc*nv, nv, data);
//      #pragma omp parallel for
        for (int nkj = -k; nkj < maxK1; ++nkj, ++nbj) {
          fill(datas, datas+nov, czero);
          std::complex<double>* conjc2 = new std::complex<double>[nbasis_j_ * jsize];
          const double img = k != 0 ? ((- m2 * nkj * pi) / k) : 0.0;
          const std::complex<double> exponent(0.0, img);
          const std::complex<double> prefac = exp(exponent);

          const std::complex<double>* jdata = coeff_j->bp(nkj) + nbasis_j_ * jstart;
          for (int ii = 0; ii != nbasis_j_ * jsize; ++ii) conjc2[ii] = conj(jdata[ii]);
          const int nsize = nbasis_a_ * nbasis_i_ * bsize;
          zgemm3m_("N", "N", &nsize, &jsize, &nbasis_j_, &prefac, data, &nsize,
                                                                conjc2, &nbasis_j_, &cone,
                                                                datas, &nsize);
          delete[] conjc2;
          intermediate_mKK.add_block(nbj*nov, nov, datas);
        }
      } // end of contraction j for given m2

    } // end of m2 loop

    // intermediate_mKK is ready
    for (int nkb = -k, nbja = 0, nbj = 0; nkb != maxK1; ++nkb) {
      for (int nkj = -k; nkj != maxK1; ++nkj, ++nbj) {
        intermediate_mKK.get_block(nov*nbj, nov, datas);
        const int m = nbasis_i_ * nbasis_a_;
        const int n = jsize * bsize;
        blas::transpose(datas, m, n, data);

        for (int nka = -k; nka < maxK1; ++nka, ++nbja) {
          const int nkac = nka + k;
          const double img = k != 0 ? ((m1 * nka * pi)/ k) : 0.0;
          const std::complex<double> exponent(0.0, img);
          const std::complex<double> prefac = exp(exponent);
          const int nsize = jsize * bsize;
//        #pragma omp parallel for
          for (int ii = 0; ii < nbasis_i_; ++ii) {
            const size_t offset1 = nsize * nbasis_a_ * ii;
            const size_t offset2 = nsize * asize * ii;
            zgemm3m_("N", "N", &nsize, &asize, &nbasis_a_, &prefac, data + offset1, &nsize,
                                                                coeff_a->bp(nka) + nbasis_a_ * astart, &nbasis_a_, &czero,
                                                                datas + offset2, &nsize);
          }
          if (q1_front)
            intermediate_KKK.append(novv, datas);
          else
            intermediate_KKK.add_block(novv * nbja, novv, datas);
        }
      }
    } // end of contraction a for given m1

    if (q1_front)
      intermediate_KKK.reopen_with_inout();

  } // end of m1 loop

  // now intermediate_KKK is ready.
  for (int nkb = -k; nkb != maxK1; ++nkb) {
    int nbja = (nkb + k) * KK * KK;
    for (int nkj = -k; nkj != maxK1; ++nkj) {
      for (int nka = -k; nka < maxK1; ++nka, ++nbja) {
        intermediate_KKK.get_block(novv * nbja, novv, intermediate_novv);
        std::complex<double>* conjc2 = new std::complex<double>[nbasis_i_ * isize];

        // momentum conservation
        int nki = nka + nkb - nkj;
        if (nki < - k) nki += k * 2;
        else if (nki >=  k) nki -= k * 2;

        const std::complex<double>* idata = coeff_i->bp(nki) + nbasis_i_ * istart;
        for (int ii = 0; ii != nbasis_i_ * isize; ++ii) conjc2[ii] = conj(idata[ii]);
        const int nsize = jsize * asize * bsize;
        zgemm3m_("N", "N", &nsize, &isize, &nbasis_i_, &cone, intermediate_novv, &nsize,
                                                            conjc2, &nbasis_i_, &czero,
                                                            datas, &nsize);
        delete[] conjc2;
        mo_int->append(noovv, datas);
      }
    }
  } // end of contraction i

  delete[] data;
  delete[] datas;
  delete[] data_read;
  delete[] intermediate_novv;
  delete[] blocks;

  std::cout << "  done" << std::endl <<std::endl;
  mo_int->reopen_with_inout();
  mo_int->sort_inside_blocks();
  return mo_int;

};

}

#endif

