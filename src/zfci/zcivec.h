//
// BAGEL - Parallel electron correlation program.
// Filename: civec.h
// Copyright (C) 2013 Michael Caldwell
//
// Author: Michael Caldwell <michaelcaldwell2013@u.northwestern.edu>>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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


#ifndef BAGEL_ZFCI_ZCIVEC_H
#define BAGEL_ZFCI_ZCIVEC_H

#include <list>
#include <src/parallel/staticdist.h>
#include <src/util/f77.h>
#include <src/fci/determinants.h>
#include <src/parallel/accrequest.h>
#include <src/parallel/recvrequest.h>

namespace bagel {

class ZDistCivec;

class ZCivec {
  protected:
    // The determinant space in which this Civec object is defined
    mutable std::shared_ptr<const Determinants> det_;

    size_t lena_;
    size_t lenb_;

    // !!CAUTION!!
    // cc is formated so that B runs first.
    // Also, cc_ can be null if this is constructed by Dvec.
    std::unique_ptr<std::complex<double>[]> cc_;

    std::complex<double>* cc_ptr_;

    std::complex<double>& cc(const size_t& i) { return *(cc_ptr_+i); }
    const std::complex<double>& cc(const size_t& i) const { return *(cc_ptr_+i); }
    std::complex<double>* cc() { return cc_ptr_; }
    const std::complex<double>* cc() const { return cc_ptr_; }

  public:
    ZCivec(std::shared_ptr<const Determinants> det);

    // constructor that is called by Dvec.
    ZCivec(std::shared_ptr<const Determinants> det, std::complex<double>* din_);

    // copy constructor
    ZCivec(const ZCivec& o);
    // from a distribtued Civec
    ZCivec(const ZDistCivec& o);

    // this is not a copy constructor.
    ZCivec(std::shared_ptr<ZCivec> o, std::shared_ptr<const Determinants> det);

    std::complex<double>* data() { return cc(); }
    std::complex<double>& element(size_t i, size_t j) { return cc(i+j*lenb_); }
    std::complex<double>* element_ptr(size_t i, size_t j) { return cc()+i+j*lenb_; }

    std::complex<double>& data(const size_t& i) { return cc(i); }
    const std::complex<double>& data(const size_t& i) const { return cc(i); }

    const std::complex<double>* data() const { return cc(); }
    const std::complex<double>* element_ptr(size_t i, size_t j) const { return cc()+i+j*lenb_; }

    std::shared_ptr<const Determinants> det() const { return det_; }
    void set_det(std::shared_ptr<const Determinants> o) const { det_ = o; }

    void zero() { std::fill(cc(), cc()+lena_*lenb_, 0.0); }

    size_t size() const { return lena_*lenb_; }

    std::shared_ptr<ZCivec> transpose(std::shared_ptr<const Determinants> det = std::shared_ptr<Determinants>()) const;

    size_t lena() const { return lena_; }
    size_t lenb() const { return lenb_; }

    // some functions for convenience
    void zaxpy(std::complex<double> a, const ZCivec& other);
    std::complex<double> zdotc(const ZCivec& other) const;
    double norm() const;
    double variance() const;
    void scale(const std::complex<double> a);

    std::complex<double> spin_expectation() const; // returns < S^2 >
    std::shared_ptr<ZCivec> spin() const; // returns S^2 | civec >
    std::shared_ptr<ZCivec> spin_lower(std::shared_ptr<const Determinants> target_det = std::shared_ptr<Determinants>()) const; // S_-
    std::shared_ptr<ZCivec> spin_raise(std::shared_ptr<const Determinants> target_det = std::shared_ptr<Determinants>()) const; // S_+
    void spin_decontaminate(const double thresh = 1.0e-12);

    ZCivec& operator*=(const std::complex<double>& a) { scale(a); return *this; }
    ZCivec& operator+=(const std::complex<double>& a) { zaxpy_(size(),  1.0, &a, 0, data(), 1); return *this; }
    ZCivec& operator-=(const std::complex<double>& a) { zaxpy_(size(), -1.0, &a, 0, data(), 1); return *this; }

    ZCivec& operator=(const ZCivec& o) { assert(size() == o.size()); std::copy(o.cc(), o.cc()+size(), cc()); return *this; }
    ZCivec& operator+=(const ZCivec& o) { zaxpy( 1.0, o); return *this; }
    ZCivec& operator-=(const ZCivec& o) { zaxpy(-1.0, o); return *this; }
    ZCivec& operator/=(const ZCivec& o);
    ZCivec operator/(const ZCivec& o) const;

    // assumes that Civec's in c are already orthogonal with each other.
    // returns scaling factor (see implementation)

    std::complex<double> orthog(std::list<std::shared_ptr<const ZCivec>> c);
    std::complex<double> orthog(std::shared_ptr<const ZCivec> o);
    void project_out(std::shared_ptr<const ZCivec> o) { throw std::logic_error("not yet implemented"); } //zaxpy(-zdotc(*o), *o); }

    void print(double thresh) const;

    std::shared_ptr<ZDistCivec> zdistcivec() const;
};


class ZDistCivec {
  protected:
    mutable std::shared_ptr<const Determinants> det_;

    // global dimension
    size_t lena_;
    size_t lenb_;

    // local storage
    std::unique_ptr<std::complex<double>[]> local_;

    // local alpha strings
    size_t astart_;
    size_t aend_;

    // allocation size
    size_t alloc_;

    // table for alpha string distribution
    const StaticDist dist_;

    // MPI send/receive management
    mutable std::shared_ptr<AccRequest> accum_;
    mutable std::shared_ptr<SendRequest> send_;
    mutable std::shared_ptr<PutRequest> put_;
    mutable std::shared_ptr<RecvRequest> recv_;

    // mutex for write accesses to local_
    mutable std::vector<std::mutex> mutex_;

    // for transpose, buffer can be appended
    mutable std::shared_ptr<ZDistCivec> buf_;
    mutable std::vector<int> transp_;

  public:
    ZDistCivec(std::shared_ptr<const Determinants> det);
    ZDistCivec(const ZDistCivec& o);

    ZDistCivec& operator=(const ZDistCivec& o);

    std::complex<double>& local(const size_t i) { return local_[i]; }
    const std::complex<double>& local(const size_t i) const { return local_[i]; }

    std::complex<double>* local() { return local_.get(); }
    const std::complex<double>* local() const { return local_.get(); }

    size_t size() const { return lenb_*(aend_-astart_); }
    size_t global_size() const { return lena_*lenb_; }
    size_t lena() const { return lena_; }
    size_t lenb() const { return lenb_; }

    size_t astart() const { return astart_; }
    size_t aend() const { return aend_; }
    size_t asize() const { return aend_ - astart_; }

    void zero() { std::fill_n(local_.get(), size(), 0.0); }

    std::shared_ptr<ZCivec> zcivec() const { return std::make_shared<ZCivec>(*this); }
    std::shared_ptr<const Determinants> det() const { return det_; }

    std::shared_ptr<ZDistCivec> clone() const { return std::make_shared<ZDistCivec>(det_); }

    // MPI Isend Irecv
    void init_mpi_accumulate() const;
    void accumulate_bstring_buf(std::unique_ptr<std::complex<double>[]>& buf, const size_t a) const;
    void terminate_mpi_accumulate() const;

    void init_mpi_recv() const;
    int get_bstring_buf(std::complex<double>* buf, const size_t a) const;
    void terminate_mpi_recv() const;

    // utility functions
    std::complex<double> norm() const;
    std::complex<double> variance() const;
    std::complex<double> zdotc(const ZDistCivec& o) const;
    void scale(const std::complex<double> a);
    void zaxpy(const std::complex<double> a, const ZDistCivec& o);

    void project_out(std::shared_ptr<const ZDistCivec> o) { std::logic_error("not yet implemented"); } // zaxpy(-zdotc(*o), *o); }
    std::complex<double> orthog(std::list<std::shared_ptr<const ZDistCivec>> c);
    std::complex<double> orthog(std::shared_ptr<const ZDistCivec> o);

    // mutex
    std::mutex& cimutex(const size_t& i) const { return mutex_[i]; }

// if we use a dedicated server. currently not using this. defined in src/util/serverflush.h
#ifndef USE_SERVER_THREAD
    void flush() const {
      if (accum_) accum_->flush();
      if (send_ ) send_->flush();
      if (put_  ) put_->flush();
    }
#endif

    std::shared_ptr<ZDistCivec> transpose() const;
    void transpose_wait();

};

}

#endif
