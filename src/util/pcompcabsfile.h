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

  public:
    PCompCABSFile(boost::shared_ptr<PGeometry>, const bool late_init = false, const std::string jobname = "source");

    const double cabs_schwarz(int i) const { return cabs_schwarz_[i]; };
    std::vector<int> cabs_offset() const { return cabs_offset_; };
    int cabs_offset(size_t i) const { return cabs_offset_[i]; };
    const size_t cabs_nbasis(size_t i) const { return cabs_basis_[i]->nbasis(); };

    // virtual functions
    void calculate_num_int_each();
    void store_integrals();
    void eval_new_block(double*, int, int, int);

    boost::shared_ptr<PMOFile<std::complex<double> > >
      mo_transform_cabs_aux(boost::shared_ptr<PCoeff>,
                            boost::shared_ptr<PMatrix1e>,
                            const int istart, const int ifence,
                            const int jstart, const int jfence,
                            const int astart, const int afence,
                            const int bstart, const int bfence,
                            const std::string jobname = "intermediate");

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
  this->num_int_each_.resize((s+s+1) * (s+s+1) * (l+l+1));

  const int size = this->basis_.size(); // number of shells
  const int cabs_size = cabs_basis_.size(); // number of shells

  #pragma omp parallel for reduction(+:data_written)
  for (int m1 = - s; m1 <= s; ++m1) {
    const double m1disp[3] = {0.0, 0.0, m1*a};
    size_t offset = (m1 + s) * (l + 1) * (s * 2 + 1);
    for (int m2 = - l; m2 <= l; ++m2) { // use bra-ket symmetry!!!
      const double m2disp[3] = {0.0, 0.0, m2*a};
      for (int m3 = m2 - s; m3 <= m2 + s; ++m3, ++offset) {
        const double m3disp[3] = {0.0, 0.0, m3*a};
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


template<class T>
void PCompCABSFile<T>::store_integrals() {
  const size_t cachesize_max = 100000000lu;
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
        eval_new_block(&dcache[current], m1, m2, m3);
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

  typedef boost::shared_ptr<Shell> RefShell;

  const int s = this->S_;
  const int l = this->L_;
  const int k = this->K_;
  const double a = this->A_;

  const double m1disp[3] = {0.0, 0.0, m1 * a};
  const double m2disp[3] = {0.0, 0.0, m2 * a};
  const double m3disp[3] = {0.0, 0.0, m3 * a};

  const int size = this->basis_.size(); // number of shells
  const int cabs_size = cabs_basis_.size(); // number of shells

  int* blocks = new int[size*size*size*cabs_size+1];
  blocks[0] = 0;
  int iall = 0;
  for (int i0 = 0; i0 != size; ++i0) {
    const int b0size = this->basis_[i0]->nbasis();
    for (int i1 = 0; i1 != size; ++i1) {
      const int b1size = this->basis_[i1]->nbasis();
      for (int i2 = 0; i2 != size; ++i2) {
        const int b2size = this->basis_[i2]->nbasis();
        for (int i3 = 0; i3 != cabs_size; ++i3, ++iall) {
          const int b3size = cabs_basis_[i3]->nbasis();
          const double integral_bound = this->schwarz_[(m1 + k) * size * size + i0 * size + i1]
                                      * cabs_schwarz_[(m3 - m2 + k) * size * size + i2 * size + i3];
          const bool skip_schwarz = integral_bound < SCHWARZ_THRESH;
          blocks[iall + 1] = blocks[iall] + (skip_schwarz ? 0 : (b0size * b1size * b2size * b3size));
        }
      }
    }
  }
  #pragma omp parallel for
  for (int i0 = 0; i0 < size; ++i0) {
    int offset = i0 * size * size * cabs_size;
    const RefShell b0 = this->basis_[i0]; // b0 is the center cell
    for (int i1 = 0; i1 != size; ++i1) {
      const RefShell b1 = this->basis_[i1]->move_atom(m1disp);
      for (int i2 = 0; i2 != size; ++i2) {
        const RefShell b2 = this->basis_[i2]->move_atom(m2disp);
        for (int i3 = 0; i3 != cabs_size; ++i3, ++offset) {
          const RefShell b3 = cabs_basis_[i3]->move_atom(m3disp);

          if (blocks[offset] == blocks[offset + 1]) continue;

          std::vector<RefShell> input;
          input.push_back(b3);
          input.push_back(b2);
          input.push_back(b1);
          input.push_back(b0);

          T batch(input, 1.0);
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
boost::shared_ptr<PMOFile<std::complex<double> > >
  PCompCABSFile<T>::mo_transform_cabs_aux(boost::shared_ptr<PCoeff> coeff,
                                          boost::shared_ptr<PMatrix1e> cabs_coeff,
                                          const int istart, const int ifence,
                                          const int jstart, const int jfence,
                                          const int astart, const int afence,
                                          const int bstart, const int bfence,
                                          const std::string jobname) {

  // What is different is that the coefficient of b3 is replaced by cabs_coeff.
  // Other than that, they should be the same as mo_transform.

  // Loading a (2K * 2K * nov) quantity on memory

  assert(bfence <= cabs_coeff->mdim());

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

  const int nbasis1 = this->geom_->nbasis();
  const int cabs_nbasis1 = cabs_coeff->ndim();
  const int nbasis2 = nbasis1 * nbasis1;
  const int cabs_nbasis2 = nbasis1 * cabs_nbasis1;
  const int nbasis3 = nbasis2 * nbasis1;
  const int cabs_nbasis3 = nbasis2 * cabs_nbasis1;
  const size_t nbasis4 = static_cast<size_t>(nbasis2) * nbasis2;
  const size_t cabs_nbasis4 = static_cast<size_t>(cabs_nbasis1) * nbasis3;

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
  boost::shared_ptr<PMOFile<std::complex<double> > >
    mo_int(new PMOFile<std::complex<double> >(filesize, k,
                                              istart, ifence, jstart, jfence,
                                              astart, afence, bstart, bfence, true));

  // we are assuming that the (c.a. two-electron integrals for a unit cell)*K^2 can be
  // held in core. If that is not the case, this must be rewritten.

  // allocating a temp array
  const int bsizemax = std::max(std::max(nbasis1, bsize), cabs_nbasis1);
  std::complex<double>* data = new std::complex<double>[nbasis3 * bsizemax * std::max(KK, 1)];
  std::complex<double>* datas = new std::complex<double>[nbasis3 * bsizemax * std::max(KK, 1)];
  std::complex<double>* conjc = new std::complex<double>[nbasis1 * std::max(isize, jsize)];
  double* data_read = new double[this->max_num_int_ * (s*2+1)];

  const int size = this->basis_.size();
  const int cabs_size = cabs_basis_.size();
  int* blocks = new int[size * size * size * cabs_size + 1];

  const int nv = nbasis3 * bsize;
  const int nov = nbasis2 * jsize * bsize;
  const int novv = nbasis1 * jsize * asize * bsize;

  const size_t sizem1 = s*2+1lu;
  const size_t sizem2 = l*2+1lu;
  const size_t num_loops = sizem1 * sizem1 * sizem2;
  size_t loop_counter = 0lu;
  size_t loop_mod10 = 0lu;

  size_t allocsize = *std::max_element(this->num_int_each_.begin(), this->num_int_each_.end());

  std::complex<double>* intermediate_mmK = new std::complex<double>[nv * std::max(KK, 1)];
  std::complex<double>* intermediate_mKK = new std::complex<double>[nov * std::max(KK * KK, 1)];
  PFile<std::complex<double> > intermediate_KKK(std::max(novv * KK * KK * KK, novv), k, true);

  for (int q1 = -s; q1 <= s; ++q1) {
    const bool q1_front = (q1 == -s);
    const int m1 = q1;
    fill(intermediate_mKK, intermediate_mKK + nov * std::max(KK * KK, 1), czero);
    for (int q2 = -l; q2 <= l; ++q2) {
      const int m2 = q2;
      fill(intermediate_mmK, intermediate_mmK + nv * std::max(KK, 1), czero);
      for (int q3 = -s; q3 <= s; ++q3, ++loop_counter) {
        const int m3 = m2 + q3;

        const double* cdata;

        if (q3 == -s) {
          const size_t key = sizem1 * ((q2+s) + sizem2 * (q1+s));
          size_t readsize = 0lu;
          for (int i = 0; i <= s*2; ++i) readsize += this->num_int_each_[key + i];
          size_t datasize_acc = 0lu;
          for (int i = 0; i != key; ++i) datasize_acc += this->num_int_each_[i];
          this->get_block(datasize_acc, readsize, data_read);
          cdata = data_read;
        } else {
          cdata = data_read;
          const size_t key = s + sizem1 * ((q2+s) + sizem2 * (q1+s));
          for (int i = key - s; i != key + q3; ++i) cdata += this->num_int_each_[i];
        }

        if (loop_mod10 < loop_counter * 10lu / num_loops) {
          loop_mod10++;
          std::cout << "  Loop " << loop_mod10 * 10lu << " percent done" << std::endl;
        }

        size_t local_counter = 0lu;

        blocks[0] = 0;
        int iall = 0;
        for (int i0 = 0; i0 != size; ++i0) {
          const int b0size = this->basis_[i0]->nbasis();
          for (int i1 = 0; i1 != size; ++i1) {
            const int b1size = this->basis_[i1]->nbasis();
            for (int i2 = 0; i2 != size; ++i2) {
              const int b2size = this->basis_[i2]->nbasis();
              for (int i3 = 0; i3 != cabs_size; ++i3, ++iall) {
                const int b3size = cabs_basis_[i3]->nbasis();

                double integral_bound = this->schwarz_[((q1 + k) * size + i0) * size + i1]
                                      * cabs_schwarz_[((q3 + k) * size + i2) * size + i3];
                const bool skip_schwarz = integral_bound < SCHWARZ_THRESH;
                blocks[iall + 1] = blocks[iall] + (skip_schwarz ? 0 : (b0size * b1size * b2size * b3size));
              }
            }
          }
        }

        #pragma omp parallel for
        for (int i0 = 0; i0 < size; ++i0) {
          int noffset = i0 * size * size * cabs_size;
          const int b0offset = this->offset(i0);
          const int b0size = this->nbasis(i0);
          for (int i1 = 0; i1 != size; ++i1) {
            const int b1offset = this->offset(i1);
            const int b1size = this->nbasis(i1);
            for (int i2 = 0; i2 != size; ++i2) {
              const int b2offset = this->offset(i2);
              const int b2size = this->nbasis(i2);
              for (int i3 = 0; i3 != cabs_size; ++i3, ++noffset) {
                const int b3offset = cabs_offset(i3);
                const int b3size = cabs_nbasis(i3);

                if (blocks[noffset] != blocks[noffset + 1]) {
                  const double* ndata = cdata + blocks[noffset];
                  std::complex<double>* label = data + b3offset + cabs_nbasis1
                      * (b2offset + nbasis1 * (b1offset + nbasis1 * b0offset));
                  for (int j0 = 0; j0 != b0size; ++j0, label += cabs_nbasis3) {
                    for (int j1 = 0; j1 != b1size; ++j1, label += cabs_nbasis2) {
                      for (int j2 = 0; j2 != b2size; ++j2, label += cabs_nbasis1) {
                        for (int j3 = 0; j3 != b3size; ++j3, ++label, ++ndata) {
                          *label = static_cast<std::complex<double> >(*ndata);
                        }
                        label -= b3size;
                      }
                      label -= b2size * cabs_nbasis1;
                    }
                    label -= b1size * cabs_nbasis2;
                  }
                } else {
                  std::complex<double>* label = data + b3offset + cabs_nbasis1
                      * (b2offset + nbasis1 * (b1offset + nbasis1 * b0offset));
                  for (int j0 = 0; j0 != b0size; ++j0, label += cabs_nbasis3) {
                    for (int j1 = 0; j1 != b1size; ++j1, label += cabs_nbasis2) {
                      for (int j2 = 0; j2 != b2size; ++j2, label += cabs_nbasis1) {
                        std::fill(label, label + b3size, czero);
                      }
                      label -= b2size * cabs_nbasis1;
                    }
                    label -= b1size * cabs_nbasis2;
                  }
                } // end of if skip_schwarz

              }
            }
          }
        } // end of shell loops; now data is ready

        if (true) { // start from j & b contractions -- needs transposition
          const int mn = nbasis2;
          const int mnc = cabs_nbasis2;
          mytranspose_complex_(data, &mnc, &mn, datas);
          ::memcpy(data, datas, cabs_nbasis4 * sizeof(std::complex<double>));
        }

        #pragma omp parallel for
        for (int nkb = -k; nkb < maxK1; ++nkb) {
          const int nkbc = nkb + k;
          const std::complex<double> exponent(0.0, k != 0 ? (m3 * nkb * pi) / k : 0.0);
          const std::complex<double> prefac = exp(exponent);
          int offset1 = 0;
          int offset2 = 0;
          for (int ii = 0; ii != nbasis1; ++ii, offset1 += cabs_nbasis3,
                                                offset2 += nbasis2 * bsize) {
            zgemm_("N", "N", &nbasis2, &bsize, &cabs_nbasis1, &prefac, data + offset1, &nbasis2,
                                               cabs_coeff->bp(nkb) + cabs_nbasis1 * bstart, &cabs_nbasis1, &cone,
                                               intermediate_mmK + nv * nkbc + offset2, &nbasis2);
          }
        } // end of contraction b for given m3

      } // end of m3 loop

      // intermediate_mmK is ready
      #pragma omp parallel for
      for (int nkb = -k; nkb < maxK1; ++nkb) {
        std::complex<double>* conjc2 = new std::complex<double>[nbasis1 * isize];
        const int nkbc = nkb + k;
        int nbj = nkbc * KK;
        for (int nkj = -k, nkjc = 0; nkj != maxK1; ++nkj, ++nbj, ++nkjc) {
          const std::complex<double> exponent(0.0, k != 0 ? (- m2 * nkj * pi) / k : 0.0);
          const std::complex<double> prefac = exp(exponent);

          const std::complex<double>* jdata = coeff->bp(nkj) + nbasis1 * jstart;
          for (int ii = 0; ii != nbasis1 * jsize; ++ii) conjc2[ii] = conj(jdata[ii]);
          const int nsize = nbasis2 * bsize;
          zgemm_("N", "N", &nsize, &jsize, &nbasis1, &prefac, intermediate_mmK + nv * nkbc, &nsize,
                                                              conjc2, &nbasis1, &cone,
                                                              intermediate_mKK + nov * nbj, &nsize);
        }
        delete[] conjc2;
      } // end of contraction j for given m2

    } // end of m2 loop

    // intermediate_mKK is ready
    for (int nkb = -k, nbja = 0, nbj = 0; nkb != maxK1; ++nkb) {
      for (int nkj = -k; nkj != maxK1; ++nkj, ++nbj) {
        ::memcpy(datas, intermediate_mKK + nov * nbj, nov * sizeof(std::complex<double>));
        {
          const int m = nbasis2;
          const int n = jsize * bsize;
          mytranspose_complex_(datas, &m, &n, data);
          std::fill(datas, datas + novv, czero);
        }

        #pragma omp parallel for
        for (int nka = -k; nka < maxK1; ++nka) {
          const int nkac = nka + k;
          const std::complex<double> exponent(0.0, k != 0 ? (m1 * nka * pi)/ k : 0.0);
          const std::complex<double> prefac = exp(exponent);
          const int nsize = jsize * bsize;
          int offset1 = 0;
          int offset2 = 0;
          for (int ii = 0; ii != nbasis1; ++ii, offset1 += nsize * nbasis1,
                                                offset2 += nsize * asize) {
            zgemm_("N", "N", &nsize, &asize, &nbasis1, &prefac, data + offset1, &nsize,
                                                                coeff->bp(nka) + nbasis1 * astart, &nbasis1, &czero,
                                                                datas + nkac * novv + offset2, &nsize);
          }
        }
        if (q1_front)
          intermediate_KKK.append(novv * std::max(KK, 1), datas);
        else
          intermediate_KKK.add_block(novv * nbja, novv * std::max(KK, 1), datas);
        nbja += std::max(KK, 1);
      }
    } // end of contraction a for given m1

    if (q1_front)
      intermediate_KKK.reopen_with_inout();

  } // end of m1 loop

  // now intermediate_KKK is ready.
  for (int nkb = -k; nkb != maxK1; ++nkb) {
    int nbja = (nkb + k) * KK * KK;
    for (int nkj = -k; nkj != maxK1; ++nkj) {
      intermediate_KKK.get_block(novv * nbja, novv * std::max(KK, 1), data);
      #pragma omp parallel for
      for (int nka = -k; nka < maxK1; ++nka) {
        std::complex<double>* conjc2 = new std::complex<double>[nbasis1 * isize];

        // momentum conservation
        int nki = nka + nkb - nkj;
        if (nki < - k) nki += k * 2;
        else if (nki >=  k) nki -= k * 2;

        const std::complex<double>* idata = coeff->bp(nki) + nbasis1 * istart;
        for (int ii = 0; ii != nbasis1 * isize; ++ii) conjc2[ii] = conj(idata[ii]);
        const int nsize = jsize * asize * bsize;
        zgemm_("N", "N", &nsize, &isize, &nbasis1, &cone, data + novv * (nka + k), &nsize,
                                                          conjc2, &nbasis1, &czero,
                                                          datas + noovv * (nka + k), &nsize);
        delete[] conjc2;
      }
      mo_int->append(noovv * std::max(KK, 1), datas);
      nbja += std::max(KK, 1);
    }
  } // end of contraction i

  delete[] data;
  delete[] datas;
  delete[] data_read;
  delete[] conjc;
  delete[] intermediate_mmK;
  delete[] intermediate_mKK;
  delete[] blocks;

  std::cout << "  done" << std::endl <<std::endl;
  mo_int->reopen_with_inout();
  return mo_int;

};


#endif

