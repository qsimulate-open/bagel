//
// BAGEL - Parallel electron correlation program.
// Filename: dfblock_london.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <rreynoldschem@u.northwestern.edu>
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

#include <src/util/taskqueue.h>
#include <src/df/dfblock_london.h>
#include <src/integral/libint/libint.h>
#include <src/integral/rys/eribatch.h>
#include <src/util/constants.h>
#include <src/util/simple.h>

using namespace bagel;
using namespace std;

constexpr complex<double> zone = 1.0;
constexpr complex<double> znil = 0.0;
const complex<double>* zpone = &zone;
const complex<double>* zpnil = &znil;

DFBlock_London::DFBlock_London(std::shared_ptr<const StaticDist> adist_shell, std::shared_ptr<const StaticDist> adist,
                 const size_t a, const size_t b1, const size_t b2, const int as, const int b1s, const int b2s, const bool averaged)
 : adist_shell_(adist_shell), adist_(adist), averaged_(averaged), asize_(a), b1size_(b1), b2size_(b2), astart_(as), b1start_(b1s), b2start_(b2s) {

  assert(asize_ == adist_shell->size(mpi__->rank()) || asize_ == adist_->size(mpi__->rank()) || asize_ == adist_->nele());

  const size_t amax = max(adist_shell_->size(mpi__->rank()), max(adist_->size(mpi__->rank()), asize_));
  data_ = unique_ptr<complex<double>[]>(new complex<double>[amax*b1size_*b2size_]);
}


DFBlock_London::DFBlock_London(const DFBlock_London& o)
 : adist_shell_(o.adist_shell_), adist_(o.adist_), averaged_(o.averaged_), asize_(o.asize_), b1size_(o.b1size_), b2size_(o.b2size_),
   astart_(o.astart_), b1start_(o.b1start_), b2start_(o.b2start_) {

  const size_t amax = max(adist_shell_->size(mpi__->rank()), max(adist_->size(mpi__->rank()), asize_));
  data_ = unique_ptr<complex<double>[]>(new complex<double>[amax*b1size_*b2size_]);
  copy_n(o.data_.get(), size(), data_.get());
}


void DFBlock_London::average() {
  if (averaged_) return;
  averaged_ = true;

  // first make a send and receive buffer
  const size_t o_start = astart_;
  const size_t o_end   = o_start + asize_;
  const int myrank = mpi__->rank();
  size_t t_start, t_end;
  tie(t_start, t_end) = adist_->range(myrank);

  assert(o_end - t_end >= 0);
  assert(o_start - t_start >= 0);

  // TODO so far I am not considering the cases when data must be sent to the next neighbor; CAUTION
  const size_t asendsize = o_end - t_end;
  const size_t arecvsize = o_start - t_start;

  assert(asendsize < t_end-t_start && arecvsize < t_end-t_start);

  unique_ptr<complex<double>[]> sendbuf;
  unique_ptr<complex<double>[]> recvbuf;
  int sendtag = 0;
  int recvtag = 0;

  if (asendsize) {
    TaskQueue<CopyBlockTask_complex> task(b2size_);

    sendbuf = unique_ptr<complex<double>[]>(new complex<double>[asendsize*b1size_*b2size_]);
    const size_t retsize = asize_ - asendsize;
    for (size_t b2 = 0; b2 != b2size_; ++b2)
      task.emplace_back(data_.get()+retsize+asize_*b1size_*b2, asize_, sendbuf.get()+asendsize*b1size_*b2, asendsize, asendsize, b1size_);

    task.compute();

    // send to the next node
    sendtag = mpi__->request_send(sendbuf.get(), asendsize*b1size_*b2size_, myrank+1, myrank);
  }

  if (arecvsize) {
    recvbuf = unique_ptr<complex<double>[]>(new complex<double>[arecvsize*b1size_*b2size_]);
    // recv from the previous node
    recvtag = mpi__->request_recv(recvbuf.get(), arecvsize*b1size_*b2size_, myrank-1, myrank-1);
  }

  // second move local data
  if (arecvsize || asendsize) {
    const size_t t_size = t_end - t_start;
    const size_t retsize = asize_ - asendsize;
    if (t_size <= asize_) {
      for (size_t i = 0; i != b1size_*b2size_; ++i) {
        if (i*asize_ < (i+1)*t_size-retsize) {
          copy_backward(data_.get()+i*asize_, data_.get()+i*asize_+retsize, data_.get()+(i+1)*t_size);
        } else if (i*asize_ > (i+1)*t_size-retsize) {
          copy_n(data_.get()+i*asize_, retsize, data_.get()+(i+1)*t_size-retsize);
        }
      }
    } else {
      for (long long int i = b1size_*b2size_-1; i >= 0; --i) {
        assert(i*asize_ < (i+1)*t_size-retsize);
        copy_backward(data_.get()+i*asize_, data_.get()+i*asize_+retsize, data_.get()+(i+1)*t_size);
      }
    }
  }

  // set new astart_ and asize_
  asize_ = t_end - t_start;
  astart_ = t_start;

  // set received data
  if (arecvsize) {
    // wait for recv communication
    mpi__->wait(recvtag);

    TaskQueue<CopyBlockTask_complex> task(b2size_);
    for (size_t b2 = 0; b2 != b2size_; ++b2)
      task.emplace_back(recvbuf.get()+arecvsize*b1size_*b2, arecvsize, data_.get()+asize_*b1size_*b2, asize_, arecvsize, b1size_);
    task.compute();
  }

  // wait for send communication
  if (asendsize) mpi__->wait(sendtag);

}


// reverse operation of average() function
void DFBlock_London::shell_boundary() {
  if (!averaged_) return;
  averaged_ = false;
  const size_t o_start = astart_;
  const size_t o_end = o_start + asize_;
  const int myrank = mpi__->rank();
  size_t t_start, t_end;
  tie(t_start, t_end) = adist_shell_->range(myrank);

  const size_t asendsize = t_start - o_start;
  const size_t arecvsize = t_end - o_end;
  assert(t_start >= o_start && t_end >= o_end);

  unique_ptr<complex<double>[]> sendbuf, recvbuf;
  int sendtag = 0;
  int recvtag = 0;

  if (asendsize) {
    TaskQueue<CopyBlockTask_complex> task(b2size_);
    sendbuf = unique_ptr<complex<double>[]>(new complex<double>[asendsize*b1size_*b2size_]);
    for (size_t b2 = 0; b2 != b2size_; ++b2)
      task.emplace_back(data_.get()+asize_*b1size_*b2, asize_, sendbuf.get()+asendsize*b1size_*b2, asendsize, asendsize, b1size_);

    task.compute();
    assert(myrank > 0);
    sendtag = mpi__->request_send(sendbuf.get(), asendsize*b1size_*b2size_, myrank-1, myrank);
  }
  if (arecvsize) {
    assert(myrank+1 < mpi__->size());
    recvbuf = unique_ptr<complex<double>[]>(new complex<double>[arecvsize*b1size_*b2size_]);
    recvtag = mpi__->request_recv(recvbuf.get(), arecvsize*b1size_*b2size_, myrank+1, myrank+1);
  }

  if (arecvsize || asendsize) {
    const size_t t_size = t_end - t_start;
    const size_t retsize = asize_ - asendsize;
    assert(t_size >= retsize);
    if (t_size <= asize_) {
      for (size_t i = 0; i != b1size_*b2size_; ++i) {
        assert(i*asize_+asendsize > i*t_size);
        copy_n(data_.get()+i*asize_+asendsize, retsize, data_.get()+i*t_size);
      }
    } else {
      for (long long int i = b1size_*b2size_-1; i >= 0; --i) {
        if (i*asize_+asendsize > i*t_size) {
          copy_n(data_.get()+i*asize_+asendsize, retsize, data_.get()+i*t_size);
        } else if (i*asize_+asendsize < i*t_size) {
          copy_backward(data_.get()+i*asize_+asendsize, data_.get()+(i+1)*asize_, data_.get()+i*t_size+retsize);
        }
      }
    }
  }

  // set new astart_ and asize_
  asize_ = t_end - t_start;
  astart_ = t_start;

  // set received data
  if (arecvsize) {
    // wait for recv communication
    mpi__->wait(recvtag);

    TaskQueue<CopyBlockTask_complex> task(b2size_);
    for (size_t b2 = 0; b2 != b2size_; ++b2)
      task.emplace_back(recvbuf.get()+arecvsize*b1size_*b2, arecvsize, data_.get()+asize_*b1size_*b2+(asize_-arecvsize), asize_, arecvsize, b1size_);
    task.compute();
  }

  // wait for send communication
  if (asendsize) mpi__->wait(sendtag);
}


shared_ptr<DFBlock_London> DFBlock_London::transform_second(std::shared_ptr<const ZMatrix> cmat, const bool trans) const {
  assert(trans ? cmat->mdim() : cmat->ndim() == b1size_);
  const complex<double>* const c = cmat->data();
  const int nocc = trans ? cmat->ndim() : cmat->mdim();

  // so far I only consider the following case
  assert(b1start_ == 0);
  auto out = make_shared<DFBlock_London>(adist_shell_, adist_, asize_, nocc, b2size_, astart_, 0, b2start_, averaged_);

  for (size_t i = 0; i != b2size_; ++i) {
    if (!trans)
      zgemm3m_("N", "N", asize_, nocc, b1size_, 1.0, data_.get()+i*asize_*b1size_, asize_, c, b1size_, 0.0, out->get()+i*asize_*nocc, asize_);
    else
      zgemm3m_("N", "T", asize_, nocc, b1size_, 1.0, data_.get()+i*asize_*b1size_, asize_, c, nocc, 0.0, out->get()+i*asize_*nocc, asize_);
  }
  return out;
}


shared_ptr<DFBlock_London> DFBlock_London::transform_third(std::shared_ptr<const ZMatrix> cmat, const bool trans) const {
  assert(trans ? cmat->mdim() : cmat->ndim() == b2size_);
  const complex<double>* const c = cmat->data();
  const int nocc = trans ? cmat->ndim() : cmat->mdim();

  // so far I only consider the following case
  assert(b2start_ == 0);
  auto out = make_shared<DFBlock_London>(adist_shell_, adist_, asize_, b1size_, nocc, astart_, b1start_, 0, averaged_);

  if (!trans)
    zgemm3m_("N", "N", asize_*b1size_, nocc, b2size_, 1.0, data_.get(), asize_*b1size_, c, b2size_, 0.0, out->get(), asize_*b1size_);
  else  // trans -> back transform
    zgemm3m_("N", "T", asize_*b1size_, nocc, b2size_, 1.0, data_.get(), asize_*b1size_, c, nocc, 0.0, out->get(), asize_*b1size_);

  return out;
}


shared_ptr<DFBlock_London> DFBlock_London::clone() const {
  auto out = make_shared<DFBlock_London>(adist_shell_, adist_, asize_, b1size_, b2size_, astart_, b1start_, b2start_, averaged_);
  out->zero();
  return out;
}


shared_ptr<DFBlock_London> DFBlock_London::copy() const {
  return make_shared<DFBlock_London>(*this);
}


DFBlock_London& DFBlock_London::operator+=(const DFBlock_London& o) { ax_plus_y( 1.0, o); return *this; }
DFBlock_London& DFBlock_London::operator-=(const DFBlock_London& o) { ax_plus_y(-1.0, o); return *this; }


void DFBlock_London::ax_plus_y(const double a, const DFBlock_London& o) {
  if (size() != o.size()) throw logic_error("DFBlock_London::daxpy called illegally");
  blas::ax_plus_y_n(a, o.data_.get(), size(), data_.get());
}


void DFBlock_London::scale(const double a) {
  blas::scale_n(a, data_.get(), size());
}


void DFBlock_London::add_direct_product(const shared_ptr<const ZMatrix> a, const shared_ptr<const ZMatrix> b, const double fac) {
  assert(asize_ == a->ndim() && b1size_*b2size_ == b->size());
  zgeru_(asize_, b1size_*b2size_, fac, a->data(), 1, b->data(), 1, data_.get(), asize_);
}


void DFBlock_London::symmetrize() {
  if (b1size_ != b2size_) throw logic_error("illegal call of DFBlock_London::symmetrize()");
  const int n = b1size_;
  for (int i = 0; i != n; ++i)
    for (int j = i; j != n; ++j) {
      blas::ax_plus_y_n(1.0, data_.get()+asize_*(j+n*i), asize_, data_.get()+asize_*(i+n*j));
      copy_n(data_.get()+asize_*(i+n*j), asize_, data_.get()+asize_*(j+n*i));
    }
}


shared_ptr<DFBlock_London> DFBlock_London::swap() const {
  auto out = make_shared<DFBlock_London>(adist_shell_, adist_, asize_, b2size_, b1size_, astart_, b2start_, b1start_, averaged_);
  for (size_t b2 = b2start_; b2 != b2start_+b2size_; ++b2)
    for (size_t b1 = b1start_; b1 != b1start_+b1size_; ++b1)
      copy_n(data_.get()+asize_*(b1+b1size_*b2), asize_, out->get()+asize_*(b2+b2size_*b1));
  return out;
}


/*
shared_ptr<DFBlock_London> DFBlock_London::apply_rhf_2RDM(const double scale_exch) const {
  assert(b1size_ == b2size_);
  const int nocc = b1size_;
  shared_ptr<DFBlock_London> out = clone();
  out->zero();
  // exchange contributions
  out->ax_plus_y(-2.0*scale_exch, *this);
  // coulomb contributions (diagonal to diagonal)
  unique_ptr<double[]> diagsum(new double[asize_]);
  fill_n(diagsum.get(), asize_, 0.0);
  for (int i = 0; i != nocc; ++i)
    blas::ax_plus_y_n(1.0, data_.get()+asize_*(i+nocc*i), asize_, diagsum.get());
  for (int i = 0; i != nocc; ++i)
    blas::ax_plus_y_n(4.0, diagsum.get(), asize_, out->get()+asize_*(i+nocc*i));
  return out;
}


// Caution
//   o strictly assuming that we are using natural orbitals.
//
shared_ptr<DFBlock_London> DFBlock_London::apply_uhf_2RDM(const double* amat, const double* bmat) const {
  assert(b1size_ == b2size_);
  const int nocc = b1size_;
  shared_ptr<DFBlock_London> out = clone();
  {
    unique_ptr<double[]> d2(new double[size()]);
    // exchange contributions
    zgemm3m_("N", "N", asize_*nocc, nocc, nocc, 1.0, data_.get(), asize_*nocc, amat, nocc, 0.0, d2.get(), asize_*nocc);
    for (int i = 0; i != nocc; ++i)
      zgemm3m_("N", "N", asize_, nocc, nocc, -1.0, d2.get()+asize_*nocc*i, asize_, amat, nocc, 0.0, out->get()+asize_*nocc*i, asize_);
    zgemm3m_("N", "N", asize_*nocc, nocc, nocc, 1.0, data_.get(), asize_*nocc, bmat, nocc, 0.0, d2.get(), asize_*nocc);
    for (int i = 0; i != nocc; ++i)
      zgemm3m_("N", "N", asize_, nocc, nocc, -1.0, d2.get()+asize_*nocc*i, asize_, bmat, nocc, 1.0, out->get()+asize_*nocc*i, asize_);
  }

  unique_ptr<double[]> sum(new double[nocc]);
  for (int i = 0; i != nocc; ++i) sum[i] = amat[i+i*nocc] + bmat[i+i*nocc];
  // coulomb contributions (diagonal to diagonal)
  unique_ptr<double[]> diagsum(new double[asize_]);
  fill_n(diagsum.get(), asize_, 0.0);
  for (int i = 0; i != nocc; ++i)
    blas::ax_plus_y_n(sum[i], data_.get()+asize_*(i+nocc*i), asize_, diagsum.get());
  for (int i = 0; i != nocc; ++i)
    blas::ax_plus_y_n(sum[i], diagsum.get(), asize_, out->get()+asize_*(i+nocc*i));
  return out;
}


shared_ptr<DFBlock_London> DFBlock_London::apply_2RDM(const double* rdm, const double* rdm1, const int nclosed, const int nact) const {
  assert(nclosed+nact == b1size_ && b1size_ == b2size_);
  // checking if natural orbitals...
  bool natural = true;
  {
    const double a = ddot_(nact*nact, rdm1, 1, rdm1, 1);
    double sum = 0.0;
    for (int i = 0; i != nact; ++i) sum += rdm1[i+nact*i]*rdm1[i+nact*i];
    if (fabs(a-sum) > numerical_zero__*100)
      natural = false;
  }
  shared_ptr<DFBlock_London> out = clone();
  out->zero();
  // closed-closed part
  // exchange contribution
  for (int i = 0; i != nclosed; ++i)
    for (int j = 0; j != nclosed; ++j)
      daxpy_(asize_, -2.0, data_.get()+asize_*(j+b1size_*i), 1, out->get()+asize_*(j+b1size_*i), 1);
  // coulomb contribution
  unique_ptr<double[]> diagsum(new double[asize_]);
  fill_n(diagsum.get(), asize_, 0.0);
  for (int i = 0; i != nclosed; ++i)
    daxpy_(asize_, 1.0, data_.get()+asize_*(i+b1size_*i), 1, diagsum.get(), 1);
  for (int i = 0; i != nclosed; ++i)
    daxpy_(asize_, 4.0, diagsum.get(), 1, out->get()+asize_*(i+b1size_*i), 1);

  // act-act part
  // compress
  unique_ptr<double[]> buf(new double[nact*nact*asize_]);
  unique_ptr<double[]> buf2(new double[nact*nact*asize_]);
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      copy_n(data_.get()+asize_*(j+nclosed+b1size_*(i+nclosed)), asize_, buf.get()+asize_*(j+nact*i));
  // multiply
  zgemm3m_("N", "N", asize_, nact*nact, nact*nact, 1.0, buf.get(), asize_, rdm, nact*nact, 0.0, buf2.get(), asize_);
  // slot in
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      copy_n(buf2.get()+asize_*(j+nact*i), asize_, out->get()+asize_*(j+nclosed+b1size_*(i+nclosed)));

  // closed-act part
  // coulomb contribution G^ia_ia = 2*gamma_ab
  // ASSUMING natural orbitals
  if (natural) {
    for (int i = 0; i != nact; ++i)
      daxpy_(asize_, 2.0*rdm1[i+nact*i], diagsum.get(), 1, out->get()+asize_*(i+nclosed+b1size_*(i+nclosed)), 1);
  } else {
    for (int i = 0; i != nact; ++i)
      for (int j = 0; j != nact; ++j)
        daxpy_(asize_, 2.0*rdm1[j+nact*i], diagsum.get(), 1, out->get()+asize_*(j+nclosed+b1size_*(i+nclosed)), 1);
  }
  unique_ptr<double[]> diagsum2(new double[asize_]);
  zgemv_("N", asize_, nact*nact, 1.0, buf.get(), asize_, rdm1, 1, 0.0, diagsum2.get(), 1);
  for (int i = 0; i != nclosed; ++i)
    daxpy_(asize_, 2.0, diagsum2.get(), 1, out->get()+asize_*(i+b1size_*i), 1);
  // exchange contribution
  if (natural) {
    for (int i = 0; i != nact; ++i)
      for (int j = 0; j != nclosed; ++j) {
        daxpy_(asize_, -rdm1[i+nact*i], data_.get()+asize_*(j+b1size_*(i+nclosed)), 1, out->get()+asize_*(j+b1size_*(i+nclosed)), 1);
        daxpy_(asize_, -rdm1[i+nact*i], data_.get()+asize_*(i+nclosed+b1size_*j), 1, out->get()+asize_*(i+nclosed+b1size_*j), 1);
      }
  } else {
    for (int i0 = 0; i0 != nact; ++i0)
      for (int i1 = 0; i1 != nact; ++i1)
        for (int j = 0; j != nclosed; ++j) { // TODO be careful when rdm1 is not symmetric (e.g., transition density matrices)
          daxpy_(asize_, -rdm1[i1+nact*i0], data_.get()+asize_*(j+b1size_*(i0+nclosed)), 1, out->get()+asize_*(j+b1size_*(i1+nclosed)), 1);
          daxpy_(asize_, -rdm1[i1+nact*i0], data_.get()+asize_*(i0+nclosed+b1size_*j), 1, out->get()+asize_*(i1+nclosed+b1size_*j), 1);
        }
  }
  return out;
}


shared_ptr<DFBlock_London> DFBlock_London::apply_2RDM(const double* rdm) const {
  shared_ptr<DFBlock_London> out = clone();
  zgemm3m_("N", "T", asize_, b1size_*b2size_, b1size_*b2size_, 1.0, data_.get(), asize_, rdm, b1size_*b2size_, 0.0, out->get(), asize_);
  return out;
}
*/


shared_ptr<ZMatrix> DFBlock_London::form_2index(const shared_ptr<const DFBlock_London> o, const double a) const {
  if (asize_ != o->asize_ || (b1size_ != o->b1size_ && b2size_ != o->b2size_)) throw logic_error("illegal call of DFBlock_London::form_2index");
  shared_ptr<ZMatrix> target;

  if (b1size_ == o->b1size_) {
    target = make_shared<ZMatrix>(b2size_,o->b2size_);
    zgemm3m_("T", "N", b2size_, o->b2size_, asize_*b1size_, a, data_.get(), asize_*b1size_, o->data_.get(), asize_*b1size_, 0.0, target->data(), b2size_);
  } else {
    assert(b2size_ == o->b2size_);
    target = make_shared<ZMatrix>(b1size_,o->b1size_);
    for (int i = 0; i != b2size_; ++i)
      zgemm3m_("T", "N", b1size_, o->b1size_, asize_, a, data_.get()+i*asize_*b1size_, asize_, o->data_.get()+i*asize_*o->b1size_, asize_, 1.0, target->data(), b1size_);
  }

  return target;
}


shared_ptr<ZMatrix> DFBlock_London::form_4index(const shared_ptr<const DFBlock_London> o, const double a) const {
  if (asize_ != o->asize_) throw logic_error("illegal call of DFBlock_London::form_4index");
  auto target = make_shared<ZMatrix>(b1size_*b2size_, o->b1size_*o->b2size_);
  zgemm3m_("T", "N", b1size_*b2size_, o->b1size_*o->b2size_, asize_, a, data_.get(), asize_, o->data_.get(), asize_, 0.0, target->data(), b1size_*b2size_);
  return target;
}


// slowest index of o is fixed to n
shared_ptr<ZMatrix> DFBlock_London::form_4index_1fixed(const shared_ptr<const DFBlock_London> o, const double a, const size_t n) const {
  if (asize_ != o->asize_) throw logic_error("illegal call of DFBlock_London::form_4index_1fixed");
  auto target = make_shared<ZMatrix>(b2size_*b1size_, o->b1size_);
  zgemm3m_("T", "N", b1size_*b2size_, o->b1size_, asize_, a, data_.get(), asize_, o->data_.get()+n*asize_*o->b1size_, asize_, 0.0, target->data(), b1size_*b2size_);
  return target;
}


shared_ptr<ZMatrix> DFBlock_London::form_aux_2index(const shared_ptr<const DFBlock_London> o, const double a) const {
  if (b1size_ != o->b1size_ || b2size_ != o->b2size_) throw logic_error("illegal call of DFBlock_London::form_aux_2index");
  auto target = make_shared<ZMatrix>(asize_, o->asize_);
  zgemm3m_("N", "T", asize_, o->asize_, b1size_*b2size_, a, data_.get(), asize_, o->data_.get(), o->asize_, 0.0, target->data(), asize_);
  return target;
}


unique_ptr<complex<double>[]> DFBlock_London::form_vec(const shared_ptr<const ZMatrix> den) const {
  unique_ptr<complex<double>[]> out(new complex<double>[asize_]);
  assert(den->ndim() == b1size_ && den->mdim() == b2size_);
  zgemv_("N", asize_, b1size_*b2size_, zone, data_.get(), asize_, den->data(), 1, znil, out.get(), 1);
  return out;
}


shared_ptr<ZMatrix> DFBlock_London::form_mat(const complex<double>* fit) const {
  auto out = make_shared<ZMatrix>(b1size_,b2size_);
  zgemv_("T", asize_, b1size_*b2size_, zone, data_.get(), asize_, fit, 1, znil, out->data(), 1);
  return out;
}


void DFBlock_London::contrib_apply_J(const shared_ptr<const DFBlock_London> o, const shared_ptr<const ZMatrix> d) {
  if (b1size_ != o->b1size_ || b2size_ != o->b2size_) throw logic_error("illegal call of DFBlock_London::contrib_apply_J");
  zgemm3m_("N", "N", asize_, b1size_*b2size_, o->asize_, zone, d->element_ptr(astart_, o->astart_), d->ndim(), o->data_.get(), o->asize_,
                                                        zone, data_.get(), asize_);
}


void DFBlock_London::copy_block(const std::shared_ptr<const ZMatrix> o, const int jdim, const size_t offset) {
  assert(o->size() == asize_*jdim);
  copy_n(o->data(), asize_*jdim, data_.get()+offset);
}


void DFBlock_London::add_block(const std::shared_ptr<const ZMatrix> o, const int jdim, const size_t offset, const double fac) {
  assert(o->size() == asize_*jdim);
  blas::ax_plus_y_n(fac, o->data(), asize_*jdim, data_.get()+offset);
}


shared_ptr<ZMatrix> DFBlock_London::form_Dj(const shared_ptr<const ZMatrix> o, const int jdim) const {
  assert(o->size() == b1size_*b2size_*jdim);
  auto out = make_shared<ZMatrix>(asize_, jdim);
  zgemm3m_("N", "N", asize_, jdim, b1size_*b2size_, 1.0, data_.get(), asize_, o->data(), b1size_*b2size_, 0.0, out->data(), asize_);
  return out;
}


shared_ptr<ZMatrix> DFBlock_London::get_block(const int ist, const int i, const int jst, const int j, const int kst, const int k) const {
  const int ista = ist - astart_;
  const int jsta = jst - b1start_;
  const int ksta = kst - b2start_;
  const int ifen = ist + i - astart_;
  const int jfen = jst + j - b1start_;
  const int kfen = kst + k - b2start_;
  if (ista < 0 || jsta < 0 || ksta < 0 || ifen > asize_ || jfen > b1size_ || kfen > b2size_)
    throw logic_error("illegal call of DFBlock_London::get_block");

  // TODO we need 3-index tensor class here!
  auto out = make_shared<ZMatrix>(i, j*k);
  complex<double>* d = out->data();
  for (int kk = ksta; kk != kfen; ++kk)
    for (int jj = jsta; jj != jfen; ++jj, d += i)
      copy_n(data_.get()+ista+asize_*(jj+b1size_*kk), i, d);

  return out;
}


shared_ptr<ZMatrix> DFBlock_London::get_block_conj(const int ist, const int i, const int jst, const int j, const int kst, const int k) const {
  if (ist != 0 || jst != 0 || kst != 0 || astart_ != 0 || b1start_ != 0 || b2start_ != 0) throw logic_error("DFBlock_London::get_block_conj currently is not designed to work with >1 block");
  if (b1size_ != b2size_) throw logic_error ("DFBLock_London::get_block_conj assumes b1 and b2 contain the same set of basis functions");

  const int ista = ist - astart_;
  const int jsta = jst - b1start_;
  const int ksta = kst - b2start_;
  const int ifen = ist + i - astart_;
  const int jfen = jst + j - b1start_;
  const int kfen = kst + k - b2start_;
  if (ista < 0 || jsta < 0 || ksta < 0 || ifen > asize_ || jfen > b1size_ || kfen > b2size_)
    throw logic_error("illegal call of DFBlock_London::get_block");

  // TODO we need 3-index tensor class here!
  auto out = make_shared<ZMatrix>(i, j*k);
  complex<double>* d = out->data();

  for (int ii=ista; ii!=ifen; ii++) {
    for (int jj=jsta; jj!=jfen; jj++) {
      for (int kk=ksta; kk!=kfen; kk++) d[kk*b1size_*asize_ + jj*asize_ + ii] = conj(data_[kk*asize_ + jj*b1size_*asize_ + ii]);
    }
  }

  return out;
}
