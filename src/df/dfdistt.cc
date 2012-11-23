//
// BAGEL - Parallel electron correlation program.
// Filename: dfblockt.cc
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

#include <src/util/f77.h>
#include <src/parallel/mpi_interface.h>
#include <src/df/dfdistt.h>

using namespace std;
using namespace bagel;

DFDistT::DFDistT(const DFDist& in, const size_t naux, const vector<int> start, const vector<int> size)
 : data_(new double[naux*size[mpi__->rank()]]), naux_(naux), nindex1_(in.nindex2()), nindex2_(in.nindex2()), start_(start[mpi__->rank()]), size_(size[mpi__->rank()]), df_(in.df()) {

#ifndef HAVE_MPI_H
  assert(false); // this class should be used with MPI
#endif

  // first make index mapping
  for (size_t i = 0lu; i != nindex2_; ++i)
    for (size_t j = 0lu; j != nindex1_; ++j)
      index_.insert(make_pair(j+nindex1_*i, make_pair(j,i)));

  // second form a matrix
  unique_ptr<double[]> buf(new double[naux_*size_]);

  // source block
  shared_ptr<const DFBlock> source = in.block();

  // information on the data layout
  vector<pair<int, int> > atab = df_->atable();

  vector<int> srequest, rrequest;
  // first issue all the send and receive requests
  for (int i = 0; i != mpi__->size(); ++i) {
    if (i != mpi__->rank()) {
      srequest.push_back(mpi__->request_send(source->get()+source->asize()*start[i], source->asize()*size[i], i));
      rrequest.push_back(mpi__->request_recv(buf.get()+atab[i].first*size_, atab[i].second*size_, i));
    } else {
      assert(source->asize()*size[i] == atab[i].second*size_);
      copy_n(source->get()+source->asize()*start[i], source->asize()*size[i], buf.get()+atab[i].first*size_); 
    }
  }
  for (auto& i : rrequest) mpi__->wait(i);

  // second transpose each block
  for (int i = 0; i != mpi__->size(); ++i) {
    const int n = atab[i].second;
    const int m = size_;
    const int o = atab[i].first*m;
    mytranspose_(buf.get()+o, &n, &m, data_.get()+o);
  }

  for (auto& i : srequest) mpi__->wait(i);

}
