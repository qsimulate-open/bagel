//
// BAGEL - Parallel electron correlation program.
// Filename: dfblock.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
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

#include <numeric>
#include <iomanip>
#include <src/util/taskqueue.h>
#include <src/df/dfblock.h>
#include <src/rysint/libint.h>
#include <src/rysint/eribatch.h>
#include <src/util/constants.h>
#include <src/util/f77.h>
#include <src/util/simple.h>

using namespace bagel;
using namespace std;


DFBlock::DFBlock(std::shared_ptr<const StaticDist> adist, std::shared_ptr<const StaticDist> adist_shell,
                 const size_t a, const size_t b1, const size_t b2, const int as, const int b1s, const int b2s)
 : adist_(adist), adist_shell_(adist_shell), asize_(a), b1size_(b1), b2size_(b2), astart_(as), b1start_(b1s), b2start_(b2s) {

  // for consistency
  assert(adist_shell_->size(mpi__->rank()) == asize_);

  size_alloc_ = max(adist_->size(mpi__->rank()), asize_) * b1size_*b2size_; 
  data_ = unique_ptr<double[]>(new double[size_alloc_]);

}


void DFBlock::average() {
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

  unique_ptr<double[]> sendbuf;
  unique_ptr<double[]> recvbuf;
  int sendtag = 0;
  int recvtag = 0;

  if (asendsize) { 
    vector<CopyBlockTask> task;
    task.reserve(b2size_);

    sendbuf = unique_ptr<double[]>(new double[asendsize*b1size_*b2size_]);
    const size_t retsize = asize_ - asendsize;
    for (size_t b2 = 0, i = 0; b2 != b2size_; ++b2)
      task.push_back(CopyBlockTask(data_.get()+retsize+asize_*b1size_*b2, asize_, sendbuf.get()+asendsize*b1size_*b2, asendsize, asendsize, b1size_));

    TaskQueue<CopyBlockTask> tq(task);
    tq.compute(resources__->max_num_threads());

    // send to the next node
    sendtag = mpi__->request_send(sendbuf.get(), asendsize*b1size_*b2size_, myrank+1, myrank);
  }

  if (arecvsize) {
    recvbuf = unique_ptr<double[]>(new double[arecvsize*b1size_*b2size_]);
    // recv from the previous node
    recvtag = mpi__->request_recv(recvbuf.get(), arecvsize*b1size_*b2size_, myrank-1, myrank-1);
  }

  // second move local data
  if (arecvsize || asendsize) {
    const size_t t_size = t_end - t_start;
    const size_t retsize = asize_ - asendsize;
    if (t_size <= asize_) { 
      for (size_t i = 0; i != b1size_*b2size_; ++i)
        copy_backward(data_.get()+i*asize_, data_.get()+i*asize_+retsize, data_.get()+(i+1)*t_size); 
    } else {
      for (long long int i = b1size_*b2size_-1; i >= 0; --i)
        copy_backward(data_.get()+i*asize_, data_.get()+i*asize_+retsize, data_.get()+(i+1)*t_size); 
    }
  }

  // set new astart_ and asize_
  asize_ = t_end - t_start; 
  astart_ = t_start;

  // set received data
  if (arecvsize) {
    // wait for recv communication 
    mpi__->wait(recvtag);

    vector<CopyBlockTask> task;
    task.reserve(b2size_);
    for (size_t b2 = 0, i = 0; b2 != b2size_; ++b2)
      task.push_back(CopyBlockTask(recvbuf.get()+arecvsize*b1size_*b2, arecvsize, data_.get()+asize_*b1size_*b2, asize_, arecvsize, b1size_));
    TaskQueue<CopyBlockTask> tq(task);
    tq.compute(resources__->max_num_threads());
  }

  // wait for send communication
  if (asendsize) mpi__->wait(sendtag);

}


shared_ptr<DFBlock> DFBlock::transform_second(const double* const c, const int nocc, const bool trans) const {
  // so far I only consider the following case
  assert(b1start_ == 0);
  unique_ptr<double[]> tmp(new double[asize_*nocc*b2size_]);

  for (size_t i = 0; i != b2size_; ++i) {
    if (!trans)
      dgemm_("N", "N", asize_, nocc, b1size_, 1.0, data_.get()+i*asize_*b1size_, asize_, c, b1size_, 0.0, tmp.get()+i*asize_*nocc, asize_);
    else
      dgemm_("N", "T", asize_, nocc, b1size_, 1.0, data_.get()+i*asize_*b1size_, asize_, c, nocc, 0.0, tmp.get()+i*asize_*nocc, asize_);
  }

  return shared_ptr<DFBlock>(new DFBlock(tmp, adist_, adist_shell_, asize_, nocc, b2size_, astart_, 0, b2start_));
}


shared_ptr<DFBlock> DFBlock::transform_third(const double* const c, const int nocc, const bool trans) const {
  // so far I only consider the following case
  assert(b2start_ == 0);
  unique_ptr<double[]> tmp(new double[asize_*b1size_*nocc]);

  if (!trans)
    dgemm_("N", "N", asize_*b1size_, nocc, b2size_, 1.0, data_.get(), asize_*b1size_, c, b2size_, 0.0, tmp.get(), asize_*b1size_);
  else  // trans -> back transform
    dgemm_("N", "T", asize_*b1size_, nocc, b2size_, 1.0, data_.get(), asize_*b1size_, c, nocc, 0.0, tmp.get(), asize_*b1size_);

  return shared_ptr<DFBlock>(new DFBlock(tmp, adist_, adist_shell_, asize_, b1size_, nocc, astart_, b1start_, 0));
}


shared_ptr<DFBlock> DFBlock::clone() const {
  unique_ptr<double[]> tmp(new double[asize_*b1size_*b2size_]);
  return shared_ptr<DFBlock>(new DFBlock(tmp, adist_, adist_shell_, asize_, b1size_, b2size_, astart_, b1start_, b2start_));
}


shared_ptr<DFBlock> DFBlock::copy() const {
  unique_ptr<double[]> tmp(new double[asize_*b1size_*b2size_]);
  copy_n(data_.get(), asize_*b1size_*b2size_, tmp.get());
  return shared_ptr<DFBlock>(new DFBlock(tmp, adist_, adist_shell_, asize_, b1size_, b2size_, astart_, b1start_, b2start_));
}


DFBlock& DFBlock::operator+=(const DFBlock& o) { daxpy( 1.0, o); return *this; }
DFBlock& DFBlock::operator-=(const DFBlock& o) { daxpy(-1.0, o); return *this; }

void DFBlock::daxpy(const double a, const DFBlock& o) {
  if (size() != o.size()) throw logic_error("DFBlock::daxpy called illegally");
  daxpy_(size(), a, o.data_.get(), 1, data_.get(), 1);
}


void DFBlock::scale(const double a) {
  dscal_(size(), a, data_, 1);
}


void DFBlock::add_direct_product(const double* a, const double* b, const double fac) {
  dger_(asize_, b1size_*b2size_, fac, a, 1, b, 1, data_.get(), asize_);
}


// TODO not efficient
void DFBlock::symmetrize() {
  if (b1size_ != b2size_) throw logic_error("illegal call of DFBlock::symmetrize()");
  const int n = b1size_;
  for (int i = 0; i != n; ++i)
    for (int j = i; j != n; ++j)
      for (int k = 0; k != asize_; ++k)
        data_[k+asize_*(j+n*i)] = data_[k+asize_*(i+n*j)] = (data_[k+asize_*(j+n*i)] + data_[k+asize_*(i+n*j)]);
}


shared_ptr<DFBlock> DFBlock::swap() const {
  unique_ptr<double[]> dat(new double[asize_*b1size_*b2size_]);
  for (size_t b2 = b2start_; b2 != b2start_+b2size_; ++b2)
    for (size_t b1 = b1start_; b1 != b1start_+b1size_; ++b1)
      copy_n(data_.get()+asize_*(b1+b1size_*b2), asize_, dat.get()+asize_*(b2+b2size_*b1));

  shared_ptr<DFBlock> out(new DFBlock(dat, adist_, adist_shell_, asize_, b2size_, b1size_, astart_, b2start_, b1start_));
  return out;
}


shared_ptr<DFBlock> DFBlock::apply_rhf_2RDM() const {
  assert(b1size_ == b2size_);
  const int nocc = b1size_;
  shared_ptr<DFBlock> out = clone();
  out->zero();
  // exchange contributions
  out->daxpy(-2.0, *this); 
  // coulomb contributions (diagonal to diagonal)
  unique_ptr<double[]> diagsum(new double[asize_]);
  fill_n(diagsum.get(), asize_, 0.0);
  for (int i = 0; i != nocc; ++i)
    daxpy_(asize_, 1.0, data_.get()+asize_*(i+nocc*i), 1, diagsum.get(), 1);
  for (int i = 0; i != nocc; ++i)
    daxpy_(asize_, 4.0, diagsum.get(), 1, out->get()+asize_*(i+nocc*i), 1);
  return out;
}


// Caution
//   o strictly assuming that we are using natural orbitals.
//
shared_ptr<DFBlock> DFBlock::apply_uhf_2RDM(const double* amat, const double* bmat) const {
  assert(b1size_ == b2size_);
  const int nocc = b1size_;
  shared_ptr<DFBlock> out = clone();
  {
    unique_ptr<double[]> d2(new double[size()]);
    // exchange contributions
    dgemm_("N", "N", asize_*nocc, nocc, nocc, 1.0, data_.get(), asize_*nocc, amat, nocc, 0.0, d2.get(), asize_*nocc);
    for (int i = 0; i != nocc; ++i)
      dgemm_("N", "N", asize_, nocc, nocc, -1.0, d2.get()+asize_*nocc*i, asize_, amat, nocc, 0.0, out->get()+asize_*nocc*i, asize_);
    dgemm_("N", "N", asize_*nocc, nocc, nocc, 1.0, data_.get(), asize_*nocc, bmat, nocc, 0.0, d2.get(), asize_*nocc);
    for (int i = 0; i != nocc; ++i)
      dgemm_("N", "N", asize_, nocc, nocc, -1.0, d2.get()+asize_*nocc*i, asize_, bmat, nocc, 1.0, out->get()+asize_*nocc*i, asize_);
  }

  unique_ptr<double[]> sum(new double[nocc]);
  for (int i = 0; i != nocc; ++i) sum[i] = amat[i+i*nocc] + bmat[i+i*nocc];
  // coulomb contributions (diagonal to diagonal)
  unique_ptr<double[]> diagsum(new double[asize_]);
  fill_n(diagsum.get(), asize_, 0.0);
  for (int i = 0; i != nocc; ++i)
    daxpy_(asize_, sum[i], data_.get()+asize_*(i+nocc*i), 1, diagsum.get(), 1);
  for (int i = 0; i != nocc; ++i)
    daxpy_(asize_, sum[i], diagsum.get(), 1, out->get()+asize_*(i+nocc*i), 1);
  return out;
}



shared_ptr<DFBlock> DFBlock::apply_2RDM(const double* rdm, const double* rdm1, const int nclosed, const int nact) const {
  assert(nclosed+nact == b1size_ && b1size_ == b2size_);
  // checking if natural orbitals...
  {
    const double a = ddot_(nact*nact, rdm1, 1, rdm1, 1);
    double sum = 0.0;
    for (int i = 0; i != nact; ++i) sum += rdm1[i+nact*i]*rdm1[i+nact*i];
    if (fabs(a-sum) > numerical_zero__*100) {
      stringstream ss; ss << "DFFullDist::apply_2rdm should be called with natural orbitals " << scientific << setprecision(3) << fabs(a-sum) - numerical_zero__;
      throw logic_error(ss.str());
    }
  }
  shared_ptr<DFBlock> out = clone();
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
  dgemm_("N", "N", asize_, nact*nact, nact*nact, 1.0, buf.get(), asize_, rdm, nact*nact, 0.0, buf2.get(), asize_);
  // slot in
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      copy_n(buf2.get()+asize_*(j+nact*i), asize_, out->get()+asize_*(j+nclosed+b1size_*(i+nclosed)));

  // closed-act part
  // coulomb contribution G^ia_ia = 2*gamma_ab
  // ASSUMING natural orbitals
  for (int i = 0; i != nact; ++i)
    daxpy_(asize_, 2.0*rdm1[i+nact*i], diagsum.get(), 1, out->get()+asize_*(i+nclosed+b1size_*(i+nclosed)), 1);
  unique_ptr<double[]> diagsum2(new double[asize_]);
  dgemv_("N", asize_, nact*nact, 1.0, buf.get(), asize_, rdm1, 1, 0.0, diagsum2.get(), 1);
  for (int i = 0; i != nclosed; ++i)
    daxpy_(asize_, 2.0, diagsum2.get(), 1, out->get()+asize_*(i+b1size_*i), 1);
  // exchange contribution
  for (int i = 0; i != nact; ++i) {
    for (int j = 0; j != nclosed; ++j) {
      daxpy_(asize_, -rdm1[i+nact*i], data_.get()+asize_*(j+b1size_*(i+nclosed)), 1, out->get()+asize_*(j+b1size_*(i+nclosed)), 1);
      daxpy_(asize_, -rdm1[i+nact*i], data_.get()+asize_*(i+nclosed+b1size_*j), 1, out->get()+asize_*(i+nclosed+b1size_*j), 1);
    }
  }
  return out;
}


shared_ptr<DFBlock> DFBlock::apply_2RDM(const double* rdm) const {
  shared_ptr<DFBlock> out = clone();
  dgemm_("N", "T", asize_, b1size_*b2size_, b1size_*b2size_, 1.0, data_.get(), asize_, rdm, b1size_*b2size_, 0.0, out->get(), asize_);
  return out;
}


shared_ptr<Matrix> DFBlock::form_2index(const shared_ptr<const DFBlock> o, const double a) const {
  if (asize_ != o->asize_ || (b1size_ != o->b1size_ && b2size_ != o->b2size_)) throw logic_error("illegal call of DFBlock::form_2index");
  shared_ptr<Matrix> target;

  if (b1size_ == o->b1size_) {
    target = shared_ptr<Matrix>(new Matrix(b2size_,o->b2size_));
    dgemm_("T", "N", b2size_, o->b2size_, asize_*b1size_, a, data_.get(), asize_*b1size_, o->data_.get(), asize_*b1size_, 0.0, target->data(), b2size_);
  } else {
    assert(b2size_ == o->b2size_);
    target = shared_ptr<Matrix>(new Matrix(b1size_,o->b1size_));
    for (int i = 0; i != b2size_; ++i)
      dgemm_("T", "N", b1size_, o->b1size_, asize_, a, data_.get()+i*asize_*b1size_, asize_, o->data_.get()+i*asize_*o->b1size_, asize_, 1.0, target->data(), b1size_);
  }

  return target;
}


unique_ptr<double[]> DFBlock::form_4index(const shared_ptr<const DFBlock> o, const double a) const {
  if (asize_ != o->asize_) throw logic_error("illegal call of DFBlock::form_4index");
  unique_ptr<double[]> target(new double[b2size_*o->b2size_*b1size_*o->b1size_]);
  dgemm_("T", "N", b1size_*b2size_, o->b1size_*o->b2size_, asize_, a, data_.get(), asize_, o->data_.get(), asize_, 0.0, target.get(), b1size_*b2size_);
  return target;
}


// slowest index of o is fixed to n
unique_ptr<double[]> DFBlock::form_4index_1fixed(const shared_ptr<const DFBlock> o, const double a, const size_t n) const {
  if (asize_ != o->asize_) throw logic_error("illegal call of DFBlock::form_4index");
  unique_ptr<double[]> target(new double[b2size_*b1size_*o->b1size_]);
  dgemm_("T", "N", b1size_*b2size_, o->b1size_, asize_, a, data_.get(), asize_, o->data_.get()+n*asize_*o->b1size_, asize_, 0.0, target.get(), b1size_*b2size_);
  return target;
}


shared_ptr<Matrix> DFBlock::form_aux_2index(const shared_ptr<const DFBlock> o, const double a) const {
  if (b1size_ != o->b1size_ || b2size_ != o->b2size_) throw logic_error("illegal call of DFBlock::form_aux_2index");
  shared_ptr<Matrix> target(new Matrix(asize_, o->asize_));
  dgemm_("N", "T", asize_, o->asize_, b1size_*b2size_, a, data_.get(), asize_, o->data_.get(), o->asize_, 0.0, target->data(), asize_);
  return target;
}


unique_ptr<double[]> DFBlock::form_vec(const shared_ptr<const Matrix> den) const {
  unique_ptr<double[]> out(new double[asize_]);
  assert(den->ndim() == b1size_ && den->mdim() == b2size_);
  dgemv_("N", asize_, b1size_*b2size_, 1.0, data_.get(), asize_, den->data(), 1, 0.0, out.get(), 1);
  return out;
}


shared_ptr<Matrix> DFBlock::form_mat(const double* fit) const {
  shared_ptr<Matrix> out(new Matrix(b1size_,b2size_));
  dgemv_("T", asize_, b1size_*b2size_, 1.0, data_.get(), asize_, fit, 1, 0.0, out->data(), 1);
  return out;
}


void DFBlock::contrib_apply_J(const shared_ptr<const DFBlock> o, const shared_ptr<const Matrix> d) {
  if (b1size_ != o->b1size_ || b2size_ != o->b2size_) throw logic_error("illegal call of DFBlock::contrib_apply_J");
  dgemm_("N", "N", asize_, b1size_*b2size_, o->asize_, 1.0, d->element_ptr(astart_, o->astart_), d->ndim(), o->data_.get(), o->asize_,
                                                        1.0, data_.get(), asize_); 
}


void DFBlock::copy_block(const std::unique_ptr<double[]>& o, const int jdim, const size_t offset) {
  copy_n(o.get(), asize_*jdim, data_.get()+offset);
}


unique_ptr<double[]> DFBlock::form_Dj(const unique_ptr<double[]>& o, const int jdim) const {
  unique_ptr<double[]> out(new double[asize_*jdim]); 
  dgemm_("N", "N", asize_, jdim, b1size_*b2size_, 1.0, data_.get(), asize_, o.get(), b1size_*b2size_, 0.0, out.get(), asize_);
  return out;
}


unique_ptr<double[]> DFBlock::get_block(const int ist, const int i, const int jst, const int j, const int kst, const int k) const {
  const int ista = ist - astart_;
  const int jsta = jst - b1start_;
  const int ksta = kst - b2start_;
  const int ifen = ist + i - astart_;
  const int jfen = jst + j - b1start_;
  const int kfen = kst + k - b2start_;
  if (ista < 0 || jsta < 0 || ksta < 0 || ifen > asize_ || jfen > b1size_ || kfen > b2size_) 
    throw logic_error("illegal call of DFBlock::get_block");

  unique_ptr<double[]> out(new double[i*j*k]);
  double* d = out.get();
  for (int kk = ksta; kk != kfen; ++kk)
    for (int jj = jsta; jj != jfen; ++jj, d += i)
      copy_n(data_.get()+ista+asize_*(jj+b1size_*kk), i, d); 

  return out;
}
