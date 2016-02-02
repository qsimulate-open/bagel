//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: distcivec.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <src/ci/fci/distcivec.h>
#include <src/util/parallel/distqueue.h>
#include <src/util/parallel/mpi_interface.h>

using namespace std;
using namespace bagel;


template<typename DataType>
void GA_Task<DataType>::wait() {
#ifdef HAVE_GA
  NGA_NbWait(&tag);
#else
  assert(false);
#endif
}


template<typename DataType>
bool GA_Task<DataType>::test() {
#ifdef HAVE_GA
  return NGA_NbTest(&tag);
#else
  assert(false);
  return true;
#endif
}


template<typename DataType>
DistCivector<DataType>::DistCivector(shared_ptr<const Determinants> det) : det_(det), lena_(det->lena()), lenb_(det->lenb()), dist_(lena_, mpi__->size()),
                                                                           astart_(dist_.start(mpi__->rank())), aend_(astart_ + dist_.size(mpi__->rank())) {
#ifdef HAVE_GA
  // create GA
  auto type = is_same<double,DataType>::value ? MT_C_DBL : MT_C_DCPL;
  auto atable = dist_.atable();
  int64_t tsize = 0;
  for (auto& i : atable) {
    blocks_.push_back(i.first * lenb_);
    tsize += i.second * lenb_;
  }
  blocks_.push_back(tsize);

  int64_t nblocks = blocks_.size() - 1;
  assert(nblocks <= mpi__->size());
  ga_ = NGA_Create_irreg64(type, 1, &blocks_.back(), const_cast<char*>(""), &nblocks, blocks_.data());

#ifndef NDEBUG
  int64_t offset = astart_*lenb_;
  int64_t last   = astart_*lenb_ + size() - 1;
  assert(NGA_Nodeid() == NGA_Locate64(ga_, &offset));
  assert(NGA_Nodeid() == NGA_Locate64(ga_, &last));
#endif
#else
  assert(false);
#endif
}


template<typename DataType>
DistCivector<DataType>& DistCivector<DataType>::operator=(const DistCivector<DataType>& o) {
#ifdef HAVE_GA
  GA_Copy(o.ga_, ga_);
#endif
  return *this;
}


template<typename DataType>
shared_ptr<Civector<DataType>> DistCivector<DataType>::civec() const {
  auto out = make_shared<Civector<DataType>>(det_);
  const unique_ptr<DataType[]> loc = local();
  copy_n(loc.get(), asize()*lenb_, out->data()+astart()*lenb_);
  mpi__->allreduce(out->data(), out->size());
  return out;
}


template<typename DataType>
unique_ptr<DataType[]> DistCivector<DataType>::local() const {
  unique_ptr<DataType[]> out(new DataType[size()]);
#ifdef HAVE_GA
  int64_t start = astart_*lenb_;
  int64_t last  = start+size()-1;
  int64_t ld  = size();
  assert(is_local(astart_));
  NGA_Get64(ga_, &start, &last, out.get(), &ld);
#endif
  return move(out);
}


template<typename DataType>
bool DistCivector<DataType>::is_local(const size_t a) const {
#ifdef HAVE_GA
  const size_t lo = a*lenb_;
  int64_t nodeid = NGA_Nodeid();
  return (nodeid+1 < blocks_.size()) && (lo >= blocks_[nodeid] && lo < blocks_[nodeid+1]);
#else
  return true;
#endif
}


template<typename DataType>
void DistCivector<DataType>::set_local(const size_t la, const size_t lb, const DataType a) {
#ifdef HAVE_GA
  int64_t element = astart_*lenb_ + lb + la*lenb_;
  int64_t one = 1;
  NGA_Put64(ga_, &element, &element, const_cast<DataType*>(&a), &one);
#endif
}


template<typename DataType>
shared_ptr<GA_Task<DataType>> DistCivector<DataType>::accumulate_bstring_buf(unique_ptr<DataType[]>&& buf, const size_t a) {
#ifdef HAVE_GA
  int64_t offset = a*lenb_;
  int64_t last = offset + lenb_ - 1;
  int64_t size = lenb_;
  DataType one = 1.0;
  ga_nbhdl_t handle;
  // non-blocking accumulate call
  NGA_NbAcc64(ga_, &offset, &last, buf.get(), &size, &one, &handle);
  return make_shared<GA_Task<DataType>>(handle, move(buf));
#else
  return nullptr;
#endif
}


template<typename DataType>
shared_ptr<GA_Task<DataType>> DistCivector<DataType>::get_bstring_buf(DataType* buf, const size_t a) const {
#ifdef HAVE_GA
  int64_t offset = a*lenb_;
  int64_t last = offset + lenb_ - 1;
  int64_t size = lenb_;
  ga_nbhdl_t handle;
  // non-blocking get acall
  NGA_NbGet64(ga_, &offset, &last, buf, &size, &handle);
  return make_shared<GA_Task<DataType>>(handle);
#else
  return nullptr;
#endif
}


template<typename DataType>
void DistCivector<DataType>::local_accumulate(const DataType a, const unique_ptr<DataType[]>& buf) {
#ifdef HAVE_GA
  // asusming that the size is identical to local tile
  int64_t start = astart_ * lenb_;
  int64_t last = start + size() - 1;
  int64_t ld = size();
  NGA_Acc64(ga_, &start, &last, const_cast<DataType*>(buf.get()), &ld, const_cast<DataType*>(&a));
#endif
}


template<typename DataType>
void DistCivector<DataType>::zero() {
#ifdef HAVE_GA
  GA_Zero(ga_);
#endif
}


template<>
double DistCivector<double>::dot_product(const DistCivector<double>& o) const {
#ifdef HAVE_GA
  return GA_Ddot(ga_, o.ga_);
#else
  return 0.0;
#endif
}


template<>
complex<double> DistCivector<complex<double>>::dot_product(const DistCivector<complex<double>>& o) const {
#ifdef HAVE_GA
  assert(size() == o.size());
  // GA_Zdot is zdotu, so we have to compute this manually
  int64_t start = astart_ * lenb_;
  int64_t ld  = size();
  int64_t last  = start + ld - 1;
  unique_ptr<complex<double>[]> a(new complex<double>[size()]);
  unique_ptr<complex<double>[]> b(new complex<double>[size()]);
  NGA_Get64(  ga_, &start, &last, a.get(), &ld);
  NGA_Get64(o.ga_, &start, &last, b.get(), &ld);
  complex<double> sum = blas::dot_product(a.get(), ld, b.get());
  mpi__->allreduce(&sum, 1);
  return sum;
#else
  return 0.0;
#endif
}


template<typename DataType>
void DistCivector<DataType>::scale(const DataType a) {
#ifdef HAVE_GA
  GA_Scale(ga_, const_cast<DataType*>(&a));
#endif
}


template<typename DataType>
void DistCivector<DataType>::ax_plus_y(const DataType a, const DistCivector<DataType>& o) {
#ifdef HAVE_GA
  assert(size() == o.size());
  DataType one = 1.0;
  GA_Add(const_cast<DataType*>(&a), o.ga_, &one, ga_, ga_);
#endif
}


template<typename DataType>
shared_ptr<DistCivector<DataType>> DistCivector<DataType>::transpose() const {
  auto out = make_shared<DistCivector<DataType>>(det_->transpose());
#ifdef HAVE_GA
  // send buffer
  unique_ptr<DataType[]> send(new DataType[max(size(),out->size())]);
  {
    unique_ptr<DataType[]> tmp = local();
    blas::transpose(tmp.get(), lenb_, asize(), send.get());
  }
  // recieve buffer
  unique_ptr<DataType[]> recv = out->local();
  // issue send and recv requests
  vector<int> rqs;
  const int tagoff = 32767 - mpi__->size()*mpi__->size();
  if (tagoff < 0)
    throw runtime_error("Set larger MPI tag offset in DistCivector<DataType>::transpose.");

  for (int i = 0; i != mpi__->size(); ++i) {
    const size_t soffset = out->dist_.start(i) * asize();
    const size_t ssize   = out->dist_.size(i)  * asize();
    const size_t roffset = dist_.start(i) * out->asize();
    const size_t rsize   = dist_.size(i)  * out->asize();
    if (i != mpi__->rank()) {
      rqs.push_back(mpi__->request_send(send.get()+soffset, ssize, i, tagoff+mpi__->rank()+i*mpi__->size()));
      rqs.push_back(mpi__->request_recv(recv.get()+roffset, rsize, i, tagoff+i+mpi__->rank()*mpi__->size()));
    } else {
      assert(rsize == ssize);
      copy_n(send.get()+soffset, ssize, recv.get()+roffset);
    }
  }
  for (auto& i : rqs)
    mpi__->wait(i);

  // rearrange recv buffer
  for (int i = 0; i != mpi__->size(); ++i) {
    const size_t roffset = dist_.start(i) * out->asize();
    const size_t size1   = dist_.size(i);
    for (int j = 0; j != out->asize(); ++j)
      copy_n(recv.get()+roffset+j*size1, size1, send.get()+dist_.start(i)+j*out->lenb_);
  }
  if (det_->nelea()*det_->neleb() & 1)
    out->scale(-1.0);
  out->local_accumulate(1.0, send);
#endif
  return out;
}


template<typename DataType>
shared_ptr<DistCivector<DataType>> DistCivector<DataType>::spin() const {
  auto out = make_shared<DistCivector<DataType>>(*this);
#ifdef HAVE_GA
  // First the easy part, S_z^2 + S_z
  const double sz = 0.5*static_cast<double>(det_->nspin());
  out->scale(sz*sz + sz + det_->neleb());

  const unique_ptr<DataType[]> source = local();
  list<shared_ptr<GA_Task<DataType>>> acc;

  const int norb = det_->norb();
  auto intermediate = make_shared<DistCivector<DataType>>(det_);
  for (int i = 0; i < norb; ++i) {
    for (int j = 0; j < norb; ++j) {
      intermediate->zero();
      for (auto& iter : det_->phia(i,j)) {
        if (is_local(iter.source)) {
          unique_ptr<DataType[]> target(new DataType[lenb_]);
          fill_n(target.get(), lenb_, 0.0);
          blas::ax_plus_y_n(static_cast<double>(iter.sign), source.get()+(iter.source-astart_)*lenb_, lenb_, target.get());
          acc.push_back(intermediate->accumulate_bstring_buf(move(target), iter.target));

          // if done, remove the buffer
          for (auto i = acc.begin(); i != acc.end(); )
            i = (*i)->test() ? acc.erase(i) : ++i;
        }
      }
      for (auto i = acc.begin(); i != acc.end(); ) {
        (*i)->wait();
        i = acc.erase(i);
      }
      GA_Sync();

      unique_ptr<DataType[]> sbuf = intermediate->local();
      unique_ptr<DataType[]> tbuf(new DataType[size()]);
      fill_n(tbuf.get(), size(), 0.0);
      for (int ia = astart_; ia < aend_; ++ia) {
        DataType* target_base = tbuf.get() + (ia-astart_)*lenb_;
        const DataType* source_base = sbuf.get() + (ia-astart_)*lenb_;
        for (auto& iter : det_->phib(j,i)) {
          target_base[iter.target] -= static_cast<double>(iter.sign) * source_base[iter.source];
        }
      }
      out->local_accumulate(1.0, tbuf);
    }
  }
#endif
  return out;
}


template<>
void DistCivector<double>::spin_decontaminate(const double thresh) {
  const int nspin = det_->nspin();
  const int max_spin = det_->nelea() + det_->neleb();
  const double expectation = static_cast<double>(nspin * (nspin + 2)) * 0.25;

  shared_ptr<DistCivec> S2 = spin();

  int k = nspin + 2;
  while(fabs(dot_product(*S2) - expectation) > thresh) {
    if (k > max_spin) throw runtime_error("Spin decontamination failed.");
    const double factor = -4.0/(static_cast<double>(k*(k+2)));
    ax_plus_y(factor, *S2);
    const double norm = this->norm();
    const double rescale = (norm*norm > 1.0e-60) ? 1.0/norm : 0.0;
    scale(rescale);

    S2 = spin();
    k += 2;
  }
}


template<>
void DistCivector<complex<double>>::spin_decontaminate(const double thresh) {
  assert(false);
}


template<typename DataType>
void DistCivector<DataType>::print(const double thresh) const {
#ifdef HAVE_GA
  vector<DataType> data;
  vector<size_t> abits;
  vector<size_t> bbits;

  const unique_ptr<DataType[]> loc = local();
  DataType* d = loc.get();

  for (size_t ia = astart_; ia < aend_; ++ia)
    for (size_t ib = 0; ib < lenb_; ++ib, ++d)
      if (abs(*d) >= thresh) {
        data.push_back(*d);
        abits.push_back(ia);
        bbits.push_back(ib);
      }

  vector<size_t> nelements(mpi__->size(), 0);
  const size_t nn = data.size();
  mpi__->allgather(&nn, 1, nelements.data(), 1);

  const size_t chunk = *max_element(nelements.begin(), nelements.end());
  data.resize(chunk, 0);
  abits.resize(chunk, 0);
  bbits.resize(chunk, 0);

  vector<DataType> alldata(chunk * mpi__->size());
  mpi__->allgather(data.data(), chunk, alldata.data(), chunk);
  vector<size_t> allabits(chunk * mpi__->size());
  mpi__->allgather(abits.data(), chunk, allabits.data(), chunk);
  vector<size_t> allbbits(chunk * mpi__->size());
  mpi__->allgather(bbits.data(), chunk, allbbits.data(), chunk);

  if (mpi__->rank() == 0) {
    multimap<double, tuple<DataType, bitset<nbit__>, bitset<nbit__>>> tmp;
    for (int i = 0; i < chunk * mpi__->size(); ++i)
      if (alldata[i] != 0.0)
        tmp.emplace(-abs(alldata[i]), make_tuple(alldata[i], det_->string_bits_a(allabits[i]), det_->string_bits_b(allbbits[i])));

    for (auto& i : tmp)
      cout << "       " << print_bit(get<1>(i.second), get<2>(i.second), det()->norb())
                << "  " << setprecision(10) << setw(15) << get<0>(i.second) << endl;

  }
#endif
}

template class bagel::GA_Task<double>;
template class bagel::GA_Task<complex<double>>;
template class bagel::DistCivector<double>;
template class bagel::DistCivector<complex<double>>;
