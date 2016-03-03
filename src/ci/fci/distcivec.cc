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

using namespace std;
using namespace bagel;


template<typename DataType>
DistCivector<DataType>::DistCivector(shared_ptr<const Determinants> det)
  : RMAWindow<DataType>(), det_(det), lena_(det->lena()), lenb_(det->lenb()), dist_(lena_, mpi__->size()),
    astart_(dist_.start(mpi__->rank())), aend_(astart_ + dist_.size(mpi__->rank())) {
  if (size() == 0)
    throw runtime_error("Use either Knowles or Harrison for FCI");
  // create an window
  this->initialize();
}


template<typename DataType>
shared_ptr<Civector<DataType>> DistCivector<DataType>::civec() const {
  auto out = make_shared<Civector<DataType>>(det_);
  copy_n(local_data(), localsize(), out->data()+astart()*lenb_);
  mpi__->allreduce(out->data(), out->size());
  return out;
}


template<typename DataType>
bool DistCivector<DataType>::is_local(const size_t a) const {
  return a >= astart_ && a < aend_;
}


template<typename DataType>
tuple<size_t, size_t, size_t> DistCivector<DataType>::locate(const size_t a) const {
  size_t rank, aoff;
  tie(rank, aoff) = dist_.locate(a);
  const size_t off = aoff * lenb_;
  return tie(rank, off, lenb_);
}


template<typename DataType>
void DistCivector<DataType>::set_local(const size_t la, const size_t lb, const DataType a) {
  RMAWindow<DataType>::set_element(mpi__->rank(), lb+la*lenb_, a);
}


template<typename DataType>
shared_ptr<DistCivector<DataType>> DistCivector<DataType>::transpose() const {
  auto out = make_shared<DistCivector<DataType>>(det_->transpose());
#ifdef HAVE_MPI_H
  // send buffer
  unique_ptr<DataType[]> send(new DataType[max(size(),out->size())]);
  blas::transpose(local_data(), lenb_, asize(), send.get());

  MPI_Win win;
  MPI_Win_create(send.get(), size()*sizeof(DataType), sizeof(DataType), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
  MPI_Win_lock_all(MPI_MODE_NOCHECK, win);

  // recieve buffer
  unique_ptr<DataType[]> recv(new DataType[out->size()]);
  // issue send and recv requests
  vector<int> rqs;

  auto type = is_same<double,DataType>::value ? MPI_DOUBLE : MPI_CXX_DOUBLE_COMPLEX;
  for (int i = 0; i != mpi__->size(); ++i) {
    const size_t roffset = dist_.start(i) * out->asize();
    const size_t rsize   = dist_.size(i)  * out->asize();
    MPI_Get(recv.get()+roffset, rsize, type, i, out->dist_.start(mpi__->rank())*dist_.size(i), rsize, type, win);
  }
  MPI_Win_unlock_all(win);
  MPI_Win_free(&win);

  // rearrange recv buffer
  for (int i = 0; i != mpi__->size(); ++i) {
    const size_t roffset = dist_.start(i) * out->asize();
    const size_t size1   = dist_.size(i);
    for (int j = 0; j != out->asize(); ++j)
      copy_n(recv.get()+roffset+j*size1, size1, send.get()+dist_.start(i)+j*out->lenb_);
  }
  out->accumulate_buffer(1.0, send);
  if (det_->nelea()*det_->neleb() & 1)
    out->scale(-1.0);
#endif
  return out;
}


template<typename DataType>
shared_ptr<DistCivector<DataType>> DistCivector<DataType>::spin() const {
  auto out = make_shared<DistCivector<DataType>>(*this);
#ifdef HAVE_MPI_H
  // First the easy part, S_z^2 + S_z
  const double sz = 0.5*static_cast<double>(det_->nspin());
  out->scale(sz*sz + sz + det_->neleb());

  list<shared_ptr<RMATask<DataType>>> acc;

  const int norb = det_->norb();
  auto intermediate = make_shared<DistCivector<DataType>>(det_);
  const DataType* base = local_data();
  for (int i = 0; i < norb; ++i) {
    for (int j = 0; j < norb; ++j) {
      intermediate->zero();
      for (auto& iter : det_->phia(i,j)) {
        if (is_local(iter.source)) {
          unique_ptr<DataType[]> target(new DataType[lenb_]);
          fill_n(target.get(), lenb_, 0.0);
          blas::ax_plus_y_n(static_cast<double>(iter.sign), base+(iter.source-astart_)*lenb_, lenb_, target.get());
          shared_ptr<RMATask<DataType>> acctask = intermediate->rma_radd(move(target), iter.target);
          if (acctask)
            acc.push_back(acctask);

          // if done, remove the buffer
          for (auto i = acc.begin(); i != acc.end(); )
            i = (*i)->test() ? acc.erase(i) : ++i;
        }
      }
      for (auto i = acc.begin(); i != acc.end(); ) {
        (*i)->wait();
        i = acc.erase(i);
      }

      const DataType* sbuf = intermediate->local_data();
      unique_ptr<DataType[]> tbuf(new DataType[size()]);
      fill_n(tbuf.get(), size(), 0.0);
      for (int ia = astart_; ia < aend_; ++ia) {
        DataType* target_base = tbuf.get() + (ia-astart_)*lenb_;
        const DataType* source_base = sbuf + (ia-astart_)*lenb_;
        for (auto& iter : det_->phib(j,i)) {
          target_base[iter.target] -= static_cast<double>(iter.sign) * source_base[iter.source];
        }
      }
      out->accumulate_buffer(1.0, tbuf);
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
#ifdef HAVE_MPI_H
  vector<DataType> data;
  vector<size_t> abits;
  vector<size_t> bbits;

  const DataType* d = local_data();

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

template class bagel::DistCivector<double>;
template class bagel::DistCivector<complex<double>>;
