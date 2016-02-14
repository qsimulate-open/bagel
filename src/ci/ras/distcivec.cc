//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ras/distcivec.cc
// Copyright (C) 2013 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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


#include <iomanip>
#include <unordered_map>
#include <src/ci/ras/distcivector.h>
#include <src/util/parallel/distqueue.h>

using namespace std;
using namespace bagel;


template<typename DataType>
DistRASCivector<DataType>::DistRASCivector(shared_ptr<const RASDeterminants> det) : RASCivector_base<DistCIBlock<DataType>>(det), global_size_(det->size()) {
#if 0
  size_t block_offset = 0;
  for (auto& ipair : det->blockinfo()) {
    if (!ipair->empty())
      blocks_.push_back(make_shared<RBlock>(ipair->stringsa(), ipair->stringsb(), block_offset));
    else
      blocks_.push_back(nullptr);
    ++block_offset;
  }
#endif
}


template<typename DataType>
DistRASCivector<DataType>::DistRASCivector(const DistRASCivector<DataType>& o) : DistRASCivector(o.det_) {
#if 0
  auto j = o.blocks_.begin();
  for (auto i = blocks_.begin(); i != blocks_.end(); ++i, ++j) {
    if (*i) copy_n((*j)->local(), (*i)->size(), (*i)->local());
  }
#endif
}


template<typename DataType>
DistRASCivector<DataType>::DistRASCivector(DistRASCivector<DataType>&& o) : RASCivector_base<DistCIBlock<DataType>>(o.det_), global_size_(det_->size()) {
  for (auto& iblock : o.blocks()) {
    blocks_.push_back(iblock);
  }
}


// Copy assignment
template<typename DataType>
DistRASCivector<DataType>& DistRASCivector<DataType>::operator=(const DistRASCivector<DataType>& o) {
#if 0
  assert(*det_ == *o.det_);
  for (auto i = blocks_.begin(), j = o.blocks_.begin(); i != blocks_.end(); ++i)
    if (*i) copy_n((*j)->local(), (*i)->size(), (*i)->local());
#endif
  return *this;
}

// Move assignment
template<typename DataType>
DistRASCivector<DataType>& DistRASCivector<DataType>::operator=(DistRASCivector<DataType>&& o) {
#if 0
  assert(*det_ == *o.det_);
  for (auto i = blocks_.begin(), j = o.blocks_.begin(); i != blocks_.end(); ++i) { *i = *j; }
#endif
  return *this;
}


template<typename DataType>
int DistRASCivector<DataType>::get_bstring_buf(double* buf, const size_t a) const {
#if 0
  assert(put_ && recv_);
  const size_t mpirank = mpi__->rank();
  shared_ptr<const RASString> aspace = det_->template space<0>(det_->string_bits_a(a));
  size_t rank, off;
  tie(rank, off) = aspace->dist()->locate(a - aspace->offset());

  int out = -1;
  if (mpirank == rank) {
    fill_n(buf, det_->lenb(), 0.0);
    for (auto b : this->template allowed_blocks<0>(aspace))
      copy_n(b->local()+off*b->lenb(), b->lenb(), buf + b->stringsb()->offset());
  } else {
    out = recv_->request_recv(buf, det_->lenb(), rank, a);
  }
  return out;
#else
  return 0.0;
#endif
}


template<typename DataType>
void DistRASCivector<DataType>::zero() {
#if 0
  this->for_each_block( [] (shared_ptr<RBlock> i) { fill_n(i->local(), i->size(), 0.0 ); } );
#endif
}


template<typename DataType>
shared_ptr<DistRASCivector<DataType>> DistRASCivector<DataType>::transpose(shared_ptr<const RASDeterminants> det) const {
  if (!det) det = det_->transpose();
  auto out = make_shared<DistRASCivector<DataType>>(det);
#if 0
  const int myrank = mpi__->rank();

  shared_ptr<DistRASCivector<DataType>> trans = clone();
  for (auto& sblock : blocks_) {
    if (!sblock) continue;
    shared_ptr<RBlock> tblock = out->block(sblock->stringsa(), sblock->stringsb());
    shared_ptr<RBlock> bufblock = trans->block(sblock->stringsb(), sblock->stringsa());
    assert(tblock->global_size() == sblock->global_size() && bufblock->global_size() == sblock->global_size());

    for (int i = 0; i < mpi__->size(); ++i) {
      tuple<size_t, size_t> outrange = tblock->dist().range(i);
      tuple<size_t, size_t> thisrange = sblock->dist().range(i);

      unique_ptr<DataType[]> tmp(new DataType[tblock->dist().size(i)*sblock->asize()]);
      for (size_t j = 0; j != sblock->asize(); ++j)
        copy_n(sblock->local()+get<0>(outrange)+j*sblock->lenb(), tblock->dist().size(i), tmp.get()+j*tblock->dist().size(i));

      const size_t off = get<0>(outrange) * sblock->asize();
      copy_n(tmp.get(), tblock->dist().size(i)*sblock->asize(), bufblock->local()+off);
      if ( i != myrank ) {
        const int tag_offset = sblock->block_offset() * mpi__->size();
        const size_t sendsize = tblock->dist().size(i) * sblock->asize();
        if (sendsize)
          out->transp_.push_back(mpi__->request_send(bufblock->local()+off, sendsize, i, tag_offset + myrank));
        const size_t recvsize = tblock->asize() * sblock->dist().size(i);
        if (recvsize)
          out->transp_.push_back(mpi__->request_recv(tblock->local()+tblock->asize()*get<0>(thisrange), recvsize, i, tag_offset + i));
      }
      else {
        copy_n(bufblock->local() + off, tblock->asize() * sblock->asize(), tblock->local() + sblock->astart() * tblock->asize());
      }
    }
  }
  out->buf_ = trans;
#endif
  return out;
}


template<typename DataType>
DataType DistRASCivector<DataType>::dot_product(const DistRASCivector<DataType>& o) const {
  assert(det_->nelea() == o.det()->nelea() && det_->neleb() == o.det()->neleb() && det_->norb() == o.det()->norb());
  DataType out(0.0);
#if 0
  for (auto& iblock : this->blocks()) {
    if (!iblock) continue;
    shared_ptr<const RBlock> jblock = o.block(iblock->stringsb(), iblock->stringsa());
    if (jblock) out += blas::dot_product(iblock->local(), iblock->size(), jblock->local());
  }
  mpi__->allreduce(&out, 1);
#endif
  return out;
}


template<typename DataType>
void DistRASCivector<DataType>::scale(const DataType a) {
#if 0
  this->for_each_block( [&a] (shared_ptr<RBlock> b) { for_each(b->local(), b->local()+b->size(), [&a] (DataType& p) { p*= a; }); });
#endif
}


template<typename DataType>
void DistRASCivector<DataType>::ax_plus_y(const DataType a, const DistRASCivector<DataType>& o) {
#if 0
  this->for_each_block( [&a, &o] (shared_ptr<RBlock> iblock) {
    shared_ptr<const RBlock> jblock = o.block(iblock->stringsb(), iblock->stringsa());
    assert(jblock);
    blas::ax_plus_y_n(a, jblock->local(), iblock->size(), iblock->local());
  } );
#endif
}


// Returns S^2 | civec >
// S^2 = S_z^2 + S_z + S_-S_+ with S_-S_+ = nbeta - \sum_{ij} j_alpha^dagger i_alpha i_beta^dagger j_beta
template<>
shared_ptr<DistRASCivector<double>> DistRASCivector<double>::spin() const {
#if 0
#if 1
  auto local = civec();
  auto spun = local->spin();
  return spun->distcivec();
#else
  auto out = make_shared<DistRASCivector<double>>(det_);

  unordered_map<bitset<nbit__>, size_t> lexicalmap;
  for (size_t i = 0; i < det_->lenb(); ++i)
    lexicalmap[det_->string_bits_b(i)] = i;

  this->init_mpi_recv();

  DistQueue<RAS::DistSpinTask, const DistRASCivector<double>*> tasks(this);

//This is the code that WOULD be included if you wanted to thread this. Threading doesn't really work yet.
//#define THREADING
#ifdef THREADING
  atomic<bool> still_flushing(true);
  thread flush_thread([this, &still_flushing] () { do { this->flush(); this_thread::sleep_for(sleeptime__); } while(still_flushing); });

  const size_t nthreads = resources__->max_num_threads();
  auto finish_worker = [&]() { tasks.finish_worker(); };

  vector<thread> thread_tasks;
  for (size_t ith = 1; ith < nthreads; ++ith)
    thread_tasks.push_back( thread(finish_worker));
#endif

  // task construction by master
  for (auto& ispace : *this->det()->stringspacea()) {
    size_t astart, aend;
    tie(astart, aend) = ispace->dist().range(mpi__->rank());
    if (astart == aend) continue;

    for (size_t ia = astart; ia < aend; ++ia) {
      tasks.emplace_and_compute(this->det()->string_bits_a(ia + ispace->offset()), this, out, this->det(), &lexicalmap);
    }
  }


  const double sz = static_cast<double>(det_->nspin()) * 0.5;
  const double fac = sz*sz + sz + static_cast<double>(det_->neleb());

  out->ax_plus_y(fac, *this);

  tasks.finish_master();

#ifdef THREADING
  for (auto& th : thread_tasks)
    th.join();

  still_flushing = false;
  flush_thread.join();
#endif

  this->terminate_mpi_recv();

  return out;
#endif
#else
  return nullptr;
#endif
}


template<>
void DistRASCivector<double>::spin_decontaminate(const double thresh) {
#if 0
  const int nspin = det_->nspin();
  const int max_spin = det_->nelea() + det_->neleb();

  const double pure_expectation = static_cast<double>(nspin * (nspin + 2)) * 0.25;

  shared_ptr<DistRASCivec> S2 = spin();
  double actual_expectation = dot_product(*S2);

  int k = nspin + 2;
  while( fabs(actual_expectation - pure_expectation) > thresh ) {
    if ( k > max_spin ) { this->print(0.05); throw runtime_error("Spin decontamination failed."); }

    const double factor = -4.0/(static_cast<double>(k*(k+2)));
    ax_plus_y(factor, *S2);

    const double norm = this->norm();
    const double rescale = (norm*norm > 1.0e-60) ? 1.0/norm : 0.0;
    scale(rescale);

    S2 = spin();
    actual_expectation = dot_product(*S2);

    k += 2;
  }
#endif
}


template<typename DataType>
void DistRASCivector<DataType>::print(const double thr) const {
#if 0
  vector<DataType> data;
  vector<size_t> abits;
  vector<size_t> bbits;
  // multimap sorts elements so that they will be shown in the descending order in magnitude
  multimap<double, tuple<DataType, bitset<nbit__>, bitset<nbit__>>> tmp;
  for (auto& iblock : blocks_) {
    if (!iblock) continue;
    double* i = iblock->local();
    for (size_t ia = iblock->astart(); ia < iblock->aend(); ++ia) {
      for (size_t ib = 0; ib < iblock->lenb(); ++ib) {
        if (abs(*i) >= thr) {
          data.push_back(*i);
          abits.push_back(ia + iblock->stringsa()->offset());
          bbits.push_back(ib + iblock->stringsb()->offset());
        }
        ++i;
      }
    }
  }
  vector<size_t> nelements(mpi__->size(), 0);
  const size_t nn = data.size();
  mpi__->allgather(&nn, 1, nelements.data(), 1);

  const size_t chunk = *max_element(nelements.begin(), nelements.end());
  data.resize(chunk, 0);
  abits.resize(chunk, 0);
  bbits.resize(chunk, 0);

  vector<double> alldata(chunk * mpi__->size());
  mpi__->allgather(data.data(), chunk, alldata.data(), chunk);
  vector<size_t> allabits(chunk * mpi__->size());
  mpi__->allgather(abits.data(), chunk, allabits.data(), chunk);
  vector<size_t> allbbits(chunk * mpi__->size());
  mpi__->allgather(bbits.data(), chunk, allbbits.data(), chunk);

  if (mpi__->rank() == 0) {
    multimap<double, tuple<double, bitset<nbit__>, bitset<nbit__>>> tmp;
    for (int i = 0; i < chunk * mpi__->size(); ++i) {
      if (alldata[i] != 0.0)
        tmp.emplace(-abs(alldata[i]), make_tuple(alldata[i], det_->string_bits_a(allabits[i]), det_->string_bits_b(allbbits[i])));
    }

    for (auto& i : tmp) {
      cout << "       " << print_bit(get<1>(i.second), get<2>(i.second), det_->ras(0))
                << "-" << print_bit(get<1>(i.second), get<2>(i.second), det_->ras(0), det_->ras(0)+det_->ras(1))
                << "-" << print_bit(get<1>(i.second), get<2>(i.second), det_->ras(0)+det_->ras(1), det_->norb())
                << "  " << setprecision(10) << setw(15) << get<0>(i.second) << endl;

    }
  }
#endif
}


template class bagel::DistRASCivector<double>;
