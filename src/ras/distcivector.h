//
// BAGEL - Parallel electron correlation program.
// Filename: ras/distcivector.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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


#ifndef __BAGEL_RAS_DISTRASCIVECTOR_H
#define __BAGEL_RAS_DISTRASCIVECTOR_H

#include <list>
#include <map>
#include <cassert>
#include <iostream>
#include <iomanip>

#include <src/parallel/staticdist.h>
#include <src/parallel/recvrequest.h>
#include <src/ras/dvector_base.h>
#include <src/ras/determinants.h>
#include <src/math/algo.h>

namespace bagel {

template <typename DataType> class DistRASCivector;

// Contains and owns all the data and information for a sub block of the CI coefficient matrix
template <typename DataType>
class DistRASBlock {
  friend class DistRASCivector<DataType>;
  protected:
    std::shared_ptr<const StringSpace> astrings_;
    std::shared_ptr<const StringSpace> bstrings_;

    const StaticDist dist_;

    // allocation size
    size_t alloc_;
    size_t astart_;
    size_t aend_;

    std::unique_ptr<DataType[]> local_;

    // Used during MPI routines
    DistRASCivector<DataType>* cc_; // parent
    const size_t block_offset_;

  public:
    DistRASBlock(std::shared_ptr<const StringSpace> astrings, std::shared_ptr<const StringSpace> bstrings, DistRASCivector<DataType>* cc, const size_t o) :
      astrings_(astrings), bstrings_(bstrings), dist_(astrings->size(), mpi__->size()), cc_(cc), block_offset_(o)
    {
      std::tie(astart_, aend_) = dist_.range(mpi__->rank());
      alloc_ = size();
      local_ = std::unique_ptr<DataType[]>(new DataType[alloc_]);
      std::fill_n(local_.get(), alloc_, 0.0);
      mutex_ = std::vector<std::mutex>(asize());
    }

    DistRASBlock(const DistRASBlock<DataType>& o, DistRASCivector<DataType>* cc) : DistRASBlock<DataType>(o.stringa(), o.stringb(), cc, o.block_offset_) {
      std::copy_n(o.local(), alloc_, local_.get());
    }

    DistRASBlock(DistRASBlock<DataType>&& o, DistRASCivector<DataType>* cc) : astrings_(o.astrings_), bstrings_(o.bstrings_), dist_(astrings_->size(), mpi__->size()), cc_(cc), block_offset_(o.block_offset_) {
      std::tie(astart_, aend_) = dist_.range(mpi__->rank());
      alloc_ = size();
      local_ = std::move(o.local_);
      mutex_ = std::vector<std::mutex>(asize());
    }

    // mutex for write accesses to local_
    mutable std::vector<std::mutex> mutex_;

    const size_t asize() const { return aend_ - astart_; }
    const size_t astart() const { return astart_; }
    const size_t aend() const { return aend_; }

    const size_t size() const { return (aend_ - astart_) * lenb(); }
    const size_t global_size() const { return lena() * lenb(); }

    const size_t lena() const { return astrings_->size(); }
    const size_t lenb() const { return bstrings_->size(); }

    DataType* local() { return local_.get(); }
    const DataType* local() const { return local_.get(); }
    void set_local(std::unique_ptr<DataType[]>&& d) { local_ = std::move(d); }

    std::shared_ptr<const StringSpace> stringa() const { return astrings_; }
    std::shared_ptr<const StringSpace> stringb() const { return bstrings_; }
};

template <typename DataType>
class DistRASCivector {
  friend class DistRASBlock<DataType>;
  public: using DetType = RASDeterminants;
  public: using RBlock = DistRASBlock<DataType>;
  protected:
    std::vector<std::shared_ptr<RBlock>> blocks_;

    std::shared_ptr<const RASDeterminants> det_;

    mutable std::shared_ptr<RecvRequest> recv_;
    mutable std::shared_ptr<BufferPutRequest> put_;

    const size_t global_size_;

    const int hpaddress(const int na, const int nb) const {
      const int N = na + nb;
      return ( (N*(N+1))/2 + nb );
    }

    // for transpose, buffer can be appended
    mutable std::shared_ptr<DistRASCivector<DataType>> buf_;
    mutable std::vector<int> transp_;

    mutable std::mutex mutex_;

  public:
    DistRASCivector(std::shared_ptr<const RASDeterminants> det) : det_(det), global_size_(det->size()) {
      size_t block_offset = 0;
      for (auto& ipair : det->stringpairs()) {
        if (ipair.first && ipair.second)
          blocks_.push_back(std::make_shared<RBlock>(ipair.first, ipair.second, this, block_offset));
        else
          blocks_.push_back(std::shared_ptr<RBlock>());
        ++block_offset;
      }
    }

    DistRASCivector(const DistRASCivector<DataType>& o) : det_(o.det_), global_size_(det_->size()) {
      for (auto& iblock : o.blocks()) {
        if (iblock)
          blocks_.push_back(std::make_shared<RBlock>(*iblock, this));
        else
          blocks_.push_back(std::shared_ptr<RBlock>());
      }
    }

    DistRASCivector(DistRASCivector<DataType>&& o) : det_(o.det_), global_size_(det_->size()) {
      for (auto& iblock : o.blocks()) {
        if (iblock)
          blocks_.push_back(std::make_shared<RBlock>(std::move(*iblock), this));
        else
          blocks_.push_back(std::shared_ptr<RBlock>());
      }
    }

    // Copy assignment
    DistRASCivector<DataType>& operator=(const DistRASCivector<DataType>& o) {
      assert(*det_ == *o.det_);
      for (auto i = blocks_.begin(), j = o.blocks_.begin(); i != blocks_.end(); ++i)
        if (*i) std::copy_n((*j)->local(), (*i)->size(), (*i)->local());
      return *this;
    }

    // Move assignment
    DistRASCivector<DataType>& operator=(DistRASCivector<DataType>&& o) {
      assert(*det_ == *o.det_);
      for (auto i = blocks_.begin(), j = o.blocks_.begin(); i != blocks_.end(); ++i) {
        if (*i) {
          *i = *j;
          (*i)->cc_ = this;
        }
      }
      return *this;
    }

    // Access to individual blocks
    const std::vector<std::shared_ptr<RBlock>>& blocks() const { return blocks_; }
    std::vector<std::shared_ptr<RBlock>>& blocks() { return blocks_; }

    std::shared_ptr<RBlock> block(const int nha, const int nhb, const int npa, const int npb) {
      if ( det_->allowed(nha, nhb, npa, npb) ) {
        const int lp = det_->lenparts();
        return blocks_[ hpaddress(npa, npb) + lp * hpaddress(nha, nhb) ];
      }
      else return std::shared_ptr<RBlock>();
    }
    std::shared_ptr<RBlock> block(const std::bitset<nbit__> bstring, const std::bitset<nbit__> astring) {
      return block( det_->nholes(astring), det_->nholes(bstring), det_->nparticles(astring), det_->nparticles(bstring) );
    }
    std::shared_ptr<RBlock> block(std::shared_ptr<const StringSpace> beta, std::shared_ptr<const StringSpace> alpha) {
      return block( alpha->nholes(), beta->nholes(), alpha->nparticles(), beta->nparticles() );
    }

    std::shared_ptr<const RBlock> block(const int nha, const int nhb, const int npa, const int npb) const {
      if ( det_->allowed(nha, nhb, npa, npb) ) {
        const int lp = det_->lenparts();
        return blocks_[ hpaddress(npa, npb) + lp * hpaddress(nha, nhb) ];
      }
      else return std::shared_ptr<const RBlock>();
    }
    std::shared_ptr<const RBlock> block(const std::bitset<nbit__> bstring, const std::bitset<nbit__> astring) const {
      return block( det_->nholes(astring), det_->nholes(bstring), det_->nparticles(astring), det_->nparticles(bstring) );
    }
    std::shared_ptr<const RBlock> block(std::shared_ptr<const StringSpace> beta, std::shared_ptr<const StringSpace> alpha) const {
      return block( alpha->nholes(), beta->nholes(), alpha->nparticles(), beta->nparticles() );
    }

    template <int spin>
    const std::vector<std::shared_ptr<RBlock>> allowed_blocks(const std::bitset<nbit__> bit) { return allowed_blocks<spin>(det_->nholes(bit), det_->nparticles(bit)); }

    template <int spin>
    const std::vector<std::shared_ptr<RBlock>> allowed_blocks(const int nh, const int np) {
      std::vector<std::shared_ptr<RBlock>> out;
      for (int jp = 0; jp + np <= det_->max_particles(); ++jp) {
        for (int ih = 0; ih + nh <= det_->max_holes(); ++ih) {
          std::shared_ptr<RBlock> blk;
          if (spin == 0) blk = block(nh, ih, np, jp);
          else           blk = block(ih, nh, jp, np);

          if (blk) out.push_back(blk);
        }
      }
      return out;
    }

    template <int spin>
    const std::vector<std::shared_ptr<const RBlock>> allowed_blocks(const std::bitset<nbit__> bit) const { return allowed_blocks<spin>(det_->nholes(bit), det_->nparticles(bit)); }
    template <int spin>
    const std::vector<std::shared_ptr<const RBlock>> allowed_blocks(const std::shared_ptr<const StringSpace> space) const { return allowed_blocks<spin>(space->nholes(), space->nparticles()); }

    template <int spin>
    const std::vector<std::shared_ptr<const RBlock>> allowed_blocks(const int nh, const int np) const {
      std::vector<std::shared_ptr<const RBlock>> out;
      for (int jp = 0; jp + np <= det_->max_particles(); ++jp) {
        for (int ih = 0; ih + nh <= det_->max_holes(); ++ih) {
          std::shared_ptr<const RBlock> blk;
          if (spin == 0) blk = block(nh, ih, np, jp);
          else           blk = block(ih, nh, jp, np);

          if (blk) out.push_back(blk);
        }
      }
      return out;
    }

    // MPI routines
    // Never call concurrently
    void init_mpi_recv() const {
      std::lock_guard<std::mutex> lock(mutex_);
      put_ = std::make_shared<BufferPutRequest>();
      recv_ = std::make_shared<RecvRequest>();
    }

    // Never call concurrently
    void terminate_mpi_recv() const {
      std::lock_guard<std::mutex> lock(mutex_);
      assert( put_ && recv_);
      bool done;
      do {
        done = recv_->test();
#ifndef USE_SERVER_THREAD
        // in case no thread is running behind, we need to cycle this to flush
        size_t d = done ? 0 : 1;
        mpi__->soft_allreduce(&d, 1);
        done = d == 0;
#endif
        if (!done) this->flush();
        if (!done) std::this_thread::sleep_for(sleeptime__);
      } while (!done);
      // cancel all MPI calls
      recv_.reset();
      put_.reset();
    }

    void flush() const {
      std::lock_guard<std::mutex> lock(mutex_);
      for (auto i : put_->get_calls()) {
        // off is interpreted as lexical number of the alpha string
        const size_t tag = i[1];
        const size_t dest = i[2];
        const size_t astring = i[3];
        std::unique_ptr<double[]> buf(new double[det_->lenb()]);
        std::fill_n(buf.get(), det_->lenb(), 0.0);
        // locate astring
        std::shared_ptr<const StringSpace> aspace = det_->space<0>(det_->stringa(astring));
        size_t rank, off;
        std::tie(rank, off) = aspace->dist().locate(astring - aspace->offset());
        assert(rank == mpi__->rank());
        for (auto b : allowed_blocks<0>(aspace))
          std::copy_n(b->local() + off * b->lenb(), b->lenb(), buf.get() + b->stringb()->offset());
        put_->request_send(std::move(buf), det_->lenb(), dest, tag);
      }
#ifndef USE_SERVER_THREAD
      put_->flush();
#endif
    }

    int get_bstring_buf(double* buf, const size_t a) const {
      assert(put_ && recv_);
      const size_t mpirank = mpi__->rank();
      std::shared_ptr<const StringSpace> aspace = det_->space<0>(det_->stringa(a));
      size_t rank, off;
      std::tie(rank, off) = aspace->dist().locate(a - aspace->offset());

      int out = -1;
      if (mpirank == rank) {
        std::fill_n(buf, det_->lenb(), 0.0);
        for (auto b : allowed_blocks<0>(aspace))
          std::copy_n(b->local()+off*b->lenb(), b->lenb(), buf + b->stringb()->offset());
      } else {
        out = recv_->request_recv(buf, det_->lenb(), rank, a);
      }
      return out;
    }

    void zero() { for (auto& i : blocks_) if (i) std::fill_n(i->local(), i->size(), 0.0); }

    std::shared_ptr<const RASDeterminants> det() const { return det_; }
    std::shared_ptr<DistRASCivector<DataType>> clone() const { return std::make_shared<DistRASCivector<DataType>>(det_); }
    std::shared_ptr<DistRASCivector<DataType>> transpose(std::shared_ptr<const RASDeterminants> det = std::shared_ptr<const RASDeterminants>()) const {
      if (!det) det = det_->transpose();
      auto out = std::make_shared<DistRASCivector<DataType>>(det);
      const int myrank = mpi__->rank();

      std::shared_ptr<DistRASCivector<DataType>> trans = clone();
      for (auto& sblock : blocks_) {
        if (!sblock) continue;
        std::shared_ptr<RBlock> tblock = out->block(sblock->stringa(), sblock->stringb());
        std::shared_ptr<RBlock> bufblock = trans->block(sblock->stringb(), sblock->stringa());
        assert(tblock->global_size() == sblock->global_size() && bufblock->global_size() == sblock->global_size());

        for (int i = 0; i < mpi__->size(); ++i) {
          std::tuple<size_t, size_t> outrange = tblock->dist_.range(i);
          std::tuple<size_t, size_t> thisrange = sblock->dist_.range(i);

          std::unique_ptr<DataType[]> tmp(new DataType[tblock->dist_.size(i)*sblock->asize()]);
          for (size_t j = 0; j != sblock->asize(); ++j)
            std::copy_n(sblock->local()+std::get<0>(outrange)+j*sblock->lenb(), tblock->dist_.size(i), tmp.get()+j*tblock->dist_.size(i));

          const size_t off = std::get<0>(outrange) * sblock->asize();
          std::copy_n(tmp.get(), tblock->dist_.size(i)*sblock->asize(), bufblock->local()+off);
          if ( i != myrank ) {
            const int tag_offset = sblock->block_offset_ * mpi__->size();
            const size_t sendsize = tblock->dist_.size(i) * sblock->asize();
            if (sendsize)
              out->transp_.push_back(mpi__->request_send(bufblock->local()+off, sendsize, i, tag_offset + myrank));
            const size_t recvsize = tblock->asize() * sblock->dist_.size(i);
            if (recvsize)
              out->transp_.push_back(mpi__->request_recv(tblock->local()+tblock->asize()*std::get<0>(thisrange), recvsize, i, tag_offset + i));
          }
          else {
            std::copy_n(bufblock->local() + off, tblock->asize() * sblock->asize(), tblock->local() + sblock->astart() * tblock->asize());
          }
        }
      }

      out->buf_ = trans;
      return out;
    }

    void transpose_wait() {
      for (auto& i: transp_)
        mpi__->wait(i);
      buf_ = clone();
      for (auto i = blocks_.begin(), j = buf_->blocks().begin(); i != blocks_.end(); ++i, ++j) {
        if (!(*i)) continue;
        const size_t asize = (*i)->asize();
        const size_t lb = (*i)->lenb();
        if ( asize * lb == 0 ) continue;
        mytranspose_((*i)->local(), asize, lb, (*j)->local());
        std::copy_n((*j)->local(), asize * lb, (*i)->local());
      }
      buf_.reset();
    }

    // Safe for any structure of blocks.
    DataType dot_product(const DistRASCivector<DataType>& o) const {
      assert( det_->nelea() == o.det()->nelea() && det_->neleb() == o.det()->neleb() && det_->norb() == o.det()->norb() );
      DataType out(0.0);
      for (auto& iblock : this->blocks()) {
        if (!iblock) continue;
        std::shared_ptr<const RBlock> jblock = o.block(iblock->stringb(), iblock->stringa());
        if (!jblock) continue;

        out += blas::dot_product(iblock->local(), iblock->size(), jblock->local());
      }

      mpi__->allreduce(&out, 1);
      return out;
    }

    double norm() const { return std::sqrt(dot_product(*this)); }
    double variance() const { return dot_product(*this) / global_size_; }

    void set_det(std::shared_ptr<const RASDeterminants> det) { det_ = det; }
    void scale(const DataType a) {
      for (auto& i : blocks_)
        if (i) std::for_each( i->local(), i->local() + i->size(), [&a] (DataType& p) { p *= a; } );
    }
    void ax_plus_y(const DataType a, const DistRASCivector<DataType>& o) {
      for (auto& iblock : this->blocks()) {
        if (!iblock) continue;
        std::shared_ptr<const RBlock> jblock = o.block(iblock->stringb(), iblock->stringa());
        if (!jblock) continue;
        std::transform(iblock->local(), iblock->local() + iblock->size(), jblock->local(), iblock->local(), [&a] (const DataType& p, const DataType& q) { return p + a * q;});
      }
    }
    void ax_plus_y(const DataType a, std::shared_ptr<const DistRASCivector<DataType>> o) { ax_plus_y(a, *o); }

    // Spin functions are only implememted as specialized functions for double (see civec.cc)
    // returns < S^2 >
    DataType spin_expectation() const {
      std::shared_ptr<const DistRASCivector<DataType>> S2 = spin();
      return this->dot_product(*S2);
    }
    std::shared_ptr<DistRASCivector<DataType>> spin() const { assert(false); return std::shared_ptr<DistRASCivector<DataType>>();} // returns S^2 | civec >
    std::shared_ptr<DistRASCivector<DataType>> spin_lower(std::shared_ptr<const RASDeterminants> target_det = std::shared_ptr<RASDeterminants>()) const
      { assert(false); return std::shared_ptr<DistRASCivector<DataType>>(); } // S_-
    std::shared_ptr<DistRASCivector<DataType>> spin_raise(std::shared_ptr<const RASDeterminants> target_det = std::shared_ptr<RASDeterminants>()) const
      { assert(false); return std::shared_ptr<DistRASCivector<DataType>>(); } // S_+
    void spin_decontaminate(const double thresh = 1.0e-8) { assert(false); }

    std::shared_ptr<DistRASCivector<DataType>> apply(const int orbital, const bool action, const bool spin) const {
      // action: true -> create; false -> annihilate
      // spin: true -> alpha; false -> beta
      auto out = std::shared_ptr<DistRASCivector<DataType>>();

      return out;
    }

    void project_out(std::shared_ptr<const DistRASCivector<DataType>> o) { ax_plus_y(-dot_product(*o), *o); }

    double orthog(std::list<std::shared_ptr<const DistRASCivector<DataType>>> c) {
      for (auto& iter : c)
        project_out(iter);
      const double norm = this->norm();
      const double scal = (norm*norm<1.0e-60 ? 0.0 : 1.0/norm);
      scale(DataType(scal));
      return norm;
    }

    double orthog(std::shared_ptr<const DistRASCivector<DataType>> o) {
      return orthog(std::list<std::shared_ptr<const DistRASCivector<DataType>>>{o});
    }

    void print(const double thr = 0.05) const {
      std::vector<DataType> data;
      std::vector<size_t> abits;
      std::vector<size_t> bbits;
      // multimap sorts elements so that they will be shown in the descending order in magnitude
      std::multimap<double, std::tuple<DataType, std::bitset<nbit__>, std::bitset<nbit__>>> tmp;
      for (auto& iblock : blocks_) {
        if (!iblock) continue;
        double* i = iblock->local();
        for (size_t ia = iblock->astart(); ia < iblock->aend(); ++ia) {
          for (size_t ib = 0; ib < iblock->lenb(); ++ib) {
            if (std::abs(*i) >= thr) {
              data.push_back(*i);
              abits.push_back(ia + iblock->stringa()->offset());
              bbits.push_back(ib + iblock->stringb()->offset());
            }
            ++i;
          }
        }
      }
      std::vector<size_t> nelements(mpi__->size(), 0);
      const size_t nn = data.size();
      mpi__->allgather(&nn, 1, nelements.data(), 1);

      const size_t chunk = *std::max_element(nelements.begin(), nelements.end());
      data.resize(chunk, 0);
      abits.resize(chunk, 0);
      bbits.resize(chunk, 0);

      std::vector<double> alldata(chunk * mpi__->size());
      mpi__->allgather(data.data(), chunk, alldata.data(), chunk);
      std::vector<size_t> allabits(chunk * mpi__->size());
      mpi__->allgather(abits.data(), chunk, allabits.data(), chunk);
      std::vector<size_t> allbbits(chunk * mpi__->size());
      mpi__->allgather(bbits.data(), chunk, allbbits.data(), chunk);

      if (mpi__->rank() == 0) {
        std::multimap<double, std::tuple<double, std::bitset<nbit__>, std::bitset<nbit__>>> tmp;
        for (int i = 0; i < chunk * mpi__->size(); ++i) {
          if (alldata[i] != 0.0)
            tmp.emplace(-std::abs(alldata[i]), std::make_tuple(alldata[i], det_->stringa(allabits[i]), det_->stringb(allbbits[i])));
        }

        for (auto& i : tmp) {
          std::cout << "       " << det_->print_bit(std::get<1>(i.second), std::get<2>(i.second))
                    << "  " << std::setprecision(10) << std::setw(15) << std::get<0>(i.second) << std::endl;

        }
      }
    }
};

template<> std::shared_ptr<DistRASCivector<double>> DistRASCivector<double>::spin() const; // returns S^2 | civec >
template<> std::shared_ptr<DistRASCivector<double>> DistRASCivector<double>::spin_lower(std::shared_ptr<const RASDeterminants>) const; // S_-
template<> std::shared_ptr<DistRASCivector<double>> DistRASCivector<double>::spin_raise(std::shared_ptr<const RASDeterminants>) const; // S_+
template<> void DistRASCivector<double>::spin_decontaminate(const double thresh);

using DistRASCivec = DistRASCivector<double>;
// RASZCivec may come at some later time

using DistRASDvec = Dvector_base<DistRASCivec>;

}

#endif
