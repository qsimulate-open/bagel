//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelCASPT2_tasks5.cc
// Copyright (C) 2014 Toru Shiozaki
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

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#include <src/smith/relcaspt2/RelCASPT2_tasks5.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelCASPT2;

void Task200::Task_local::compute() {
  const Index x1 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x1, a2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(x1, a2, x0, a1), 0.0);
  {
    // tensor label: I259
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a1, a2, x0, x1);
    sort_indices<3,1,2,0,1,1,1,1>(i0data, odata, a1.size(), a2.size(), x0.size(), x1.size());
  }
  out()->add_block(odata, x1, a2, x0, a1);
}

void Task201::Task_local::compute() {
  const Index a1 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I259
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, a2, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, a2, x0, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, a2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, a2, x0, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: Gamma90
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x4, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x0, x4, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x1.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x5, a1, x4, a2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x5, a1, x4, a2)]);
      sort_indices<0,2,1,3,0,1,2,1>(i1data, i1data_sorted, x5.size(), a1.size(), x4.size(), a2.size());
      zgemm3m_("T", "N", x0.size()*x1.size(), a1.size()*a2.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), a2.size());
  out()->add_block(odata, a1, a2, x0, x1);
}

void Task202::Task_local::compute() {
  const Index a1 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I259
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, a2, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, a2, x0, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, a2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, a2, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma60
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, a1, x2, a2);
      zscal_(x3.size()*a1.size()*x2.size()*a2.size(), e0_, i1data.get(), 1);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, a1, x2, a2)]);
      sort_indices<0,2,1,3,0,1,-2,1>(i1data, i1data_sorted, x3.size(), a1.size(), x2.size(), a2.size());
      zgemm3m_("T", "N", x0.size()*x1.size(), a1.size()*a2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), a2.size());
  out()->add_block(odata, a1, a2, x0, x1);
}

void Task204::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index x0 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, x1, c1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, c1, x0), 0.0);
  {
    // tensor label: I288
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, c2, x0, x1);
    sort_indices<1,3,0,2,1,1,1,1>(i0data, odata, c1.size(), c2.size(), x0.size(), x1.size());
  }
  out()->add_block(odata, c2, x1, c1, x0);
}

void Task205::Task_local::compute() {
  const Index c1 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I288
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, c2, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, c2, x0, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c2, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma92
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x3, x1, x2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, x3, x1, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size(), x1.size(), x2.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, x3, c2, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, x3, c2, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), x3.size(), c2.size(), x2.size());
      zgemm3m_("T", "N", x0.size()*x1.size(), c1.size()*c2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c1.size(), c2.size());
  out()->add_block(odata, c1, c2, x0, x1);
}

void Task206::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x2, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x0, x1), 0.0);
  {
    // tensor label: I290
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, x2, x1, x0);
    sort_indices<0,1,3,2,1,1,1,1>(i0data, odata, c1.size(), x2.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, c1, x2, x0, x1);
}

void Task207::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I290
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma105
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x2, x5, x4, x3, x1, x0);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x2, x5, x4, x3, x1, x0)]);
        sort_indices<1,2,3,0,4,5,0,1,1,1>(i0data, i0data_sorted, x2.size(), x5.size(), x4.size(), x3.size(), x1.size(), x0.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, x5, x4, x3);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, x5, x4, x3)]);
        sort_indices<1,2,3,0,0,1,1,2>(i1data, i1data_sorted, c1.size(), x5.size(), x4.size(), x3.size());
        zgemm3m_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task208::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I290
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma6
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x4, x2, x3, x1, x0);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x4, x2, x3, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x5, x4, c1, x3);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x5, x4, c1, x3)]);
        sort_indices<0,1,3,2,0,1,1,2>(i1data, i1data_sorted, x5.size(), x4.size(), c1.size(), x3.size());
        zgemm3m_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task209::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I290
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: Gamma7
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x2, x3, x1, x0);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x2, x3, x1, x0)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, x3);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, x3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), x3.size());
    zgemm3m_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, x2.size()*x1.size()*x0.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task210::Task_local::compute() {
  const Index c3 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, x0, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(c3, x0, c1, a2), 0.0);
  {
    // tensor label: I294
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, c1, a2, x0);
    sort_indices<0,3,1,2,1,1,1,1>(i0data, odata, c3.size(), c1.size(), a2.size(), x0.size());
  }
  out()->add_block(odata, c3, x0, c1, a2);
}

void Task211::Task_local::compute() {
  const Index c3 = b(0);
  const Index c1 = b(1);
  const Index a2 = b(2);
  const Index x0 = b(3);
  // tensor label: I294
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, c1, a2, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, c1, a2, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, c1, a2, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c1, a2, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma16
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I295
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, x1, c1, a2);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, x1, c1, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), x1.size(), c1.size(), a2.size());
    zgemm3m_("T", "N", x0.size(), c3.size()*c1.size()*a2.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x0.size(), c3.size(), c1.size(), a2.size());
  out()->add_block(odata, c3, c1, a2, x0);
}

void Task212::Task_local::compute() {
  const Index c3 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I295
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, x1, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(c3, x1, c1, a2), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, x1, c1, a2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c3.size(), x1.size(), c1.size(), a2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(c1, x1, c3, a2);
    sort_indices<2,1,0,3,1,1,-1,1>(i1data, odata, c1.size(), x1.size(), c3.size(), a2.size());
  }
  out()->add_block(odata, c3, x1, c1, a2);
}

void Task213::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  {
    // tensor label: I298
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c2, a1, x1, x0);
    sort_indices<0,2,3,1,1,1,1,1>(i0data, odata, c2.size(), a1.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, c2, x1, x0, a1);
}

void Task214::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I298
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, a1, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c2, a1, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, a1, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma35
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I299
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, a1, x3, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, a1, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size(), x3.size(), x2.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), c2.size()*a1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a1.size());
  out()->add_block(odata, c2, a1, x1, x0);
}

void Task215::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I299
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, a1, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c2, a1, x3, x2), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c2, a1, x3, x2);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, c2.size(), a1.size(), x3.size(), x2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(x3, x2, c2, a1);
    sort_indices<2,3,0,1,1,1,-1,2>(i1data, odata, x3.size(), x2.size(), c2.size(), a1.size());
  }
  out()->add_block(odata, c2, a1, x3, x2);
}

void Task216::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I298
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, a1, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c2, a1, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, a1, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma29
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x3, x2, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x3, x2, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), x3.size(), x2.size(), x0.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, x3, x2, a1);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, x3, x2, a1)]);
      sort_indices<1,2,0,3,0,1,-1,2>(i1data, i1data_sorted, c2.size(), x3.size(), x2.size(), a1.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), c2.size()*a1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a1.size());
  out()->add_block(odata, c2, a1, x1, x0);
}

void Task217::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I298
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, a1, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c2, a1, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, a1, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma32
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0, x1, x2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x0, x1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, a1, c2, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, a1, c2, x2)]);
      sort_indices<0,3,1,2,0,1,1,2>(i1data, i1data_sorted, x3.size(), a1.size(), c2.size(), x2.size());
      zgemm3m_("T", "N", x0.size()*x1.size(), a1.size()*c2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), c2.size());
  out()->add_block(odata, c2, a1, x1, x0);
}

void Task218::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I298
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, a1, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c2, a1, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, a1, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x1, x0), 0.0);
  // tensor label: Gamma38
  std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
  std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
  // tensor label: h1
  std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, a1);
  std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, a1)]);
  sort_indices<0,1,0,1,-1,1>(i1data, i1data_sorted, c2.size(), a1.size());
  zgemm3m_("T", "N", x1.size()*x0.size(), c2.size()*a1.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, x1.size()*x0.size());
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a1.size());
  out()->add_block(odata, c2, a1, x1, x0);
}

void Task219::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x0, x1, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1, c1, a2), 0.0);
  {
    // tensor label: I306
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a2, x1, x0);
    sort_indices<3,2,0,1,1,1,1,1>(i0data, odata, c1.size(), a2.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, x0, x1, c1, a2);
}

void Task220::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I306
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma35
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I307
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a2, x3, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size(), x3.size(), x2.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), c1.size()*a2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), a2.size());
  out()->add_block(odata, c1, a2, x1, x0);
}

void Task221::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I307
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, x3, x2), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a2, x3, x2);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, c1.size(), a2.size(), x3.size(), x2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(x3, a2, c1, x2);
    sort_indices<2,1,0,3,1,1,-1,2>(i1data, odata, x3.size(), a2.size(), c1.size(), x2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i2data = in(0)->get_block(x3, x2, c1, a2);
    sort_indices<2,3,0,1,1,1,1,2>(i2data, odata, x3.size(), x2.size(), c1.size(), a2.size());
  }
  out()->add_block(odata, c1, a2, x3, x2);
}

void Task222::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I306
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, x1, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: Gamma7
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x2, x3, x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x2, x3, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, x3, x2, a2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, x3, x2, a2)]);
      sort_indices<2,1,0,3,0,1,1,2>(i1data, i1data_sorted, c1.size(), x3.size(), x2.size(), a2.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), c1.size()*a2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), a2.size());
  out()->add_block(odata, c1, a2, x1, x0);
}

void Task223::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I306
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, x1, x0), 0.0);
  // tensor label: Gamma38
  std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
  std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
  // tensor label: h1
  std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a2);
  std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a2)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size());
  zgemm3m_("T", "N", x1.size()*x0.size(), c1.size()*a2.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, x1.size()*x0.size());
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), a2.size());
  out()->add_block(odata, c1, a2, x1, x0);
}

void Task224::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x1, x2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(x1, x2, x0, a1), 0.0);
  {
    // tensor label: I314
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a1, x0, x2, x1);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, a1.size(), x0.size(), x2.size(), x1.size());
  }
  out()->add_block(odata, x1, x2, x0, a1);
}

void Task225::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I314
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma59
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x5, a1, x4, x3);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x5, a1, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,2>(i1data, i1data_sorted, x5.size(), a1.size(), x4.size(), x3.size());
        zgemm3m_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task226::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I314
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma57
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x4, x3, x0, x2, x1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x4, x3, x0, x2, x1)]);
        sort_indices<0,1,2,3,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x5, x4, x3, a1);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x5, x4, x3, a1)]);
        sort_indices<0,1,2,3,0,1,1,2>(i1data, i1data_sorted, x5.size(), x4.size(), x3.size(), a1.size());
        zgemm3m_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task227::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I314
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: Gamma60
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x0, x2, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, a1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, a1)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size());
    zgemm3m_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, x0.size()*x2.size()*x1.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task228::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, a4, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(c3, a4, c1, a2), 0.0);
  {
    // tensor label: I318
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a4, c3, a2);
    sort_indices<2,1,0,3,1,1,1,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
  }
  out()->add_block(odata, c3, a4, c1, a2);
}

void Task229::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I318
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a4, c3, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a4, c3, a2), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a4, c3, a2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<0,3,2,1,1,1,1,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->add_block(odata, c1, a4, c3, a2);
}

void Task230::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  {
    // tensor label: I320
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a3, c2, a1, x0);
    sort_indices<1,0,3,2,1,1,1,1>(i0data, odata, a3.size(), c2.size(), a1.size(), x0.size());
  }
  out()->add_block(odata, c2, a3, x0, a1);
}

void Task231::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I320
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a3, c2, a1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, a1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a3, c2, a1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, a1, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma38
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I321
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x1, a3, c2, a1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x1, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), a3.size(), c2.size(), a1.size());
    zgemm3m_("T", "N", x0.size(), a3.size()*c2.size()*a1.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x0.size(), a3.size(), c2.size(), a1.size());
  out()->add_block(odata, a3, c2, a1, x0);
}

void Task232::Task_local::compute() {
  const Index x1 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I321
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x1, a3, c2, a1)]);
  std::fill_n(odata.get(), out()->get_size(x1, a3, c2, a1), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, a3, c2, a1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x1.size(), a3.size(), c2.size(), a1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(x1, a1, c2, a3);
    sort_indices<0,3,2,1,1,1,1,1>(i1data, odata, x1.size(), a1.size(), c2.size(), a3.size());
  }
  out()->add_block(odata, x1, a3, c2, a1);
}

void Task233::Task_local::compute() {
  const Index x1 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x1, a2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(x1, a2, x0, a1), 0.0);
  {
    // tensor label: I324
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a1, a2, x0, x1);
    sort_indices<3,1,2,0,1,1,1,1>(i0data, odata, a1.size(), a2.size(), x0.size(), x1.size());
  }
  out()->add_block(odata, x1, a2, x0, a1);
}

void Task234::Task_local::compute() {
  const Index a1 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I324
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, a2, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, a2, x0, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, a2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, a2, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma60
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, a1, x2, a2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, a1, x2, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), x2.size(), a2.size());
      zgemm3m_("T", "N", x0.size()*x1.size(), a1.size()*a2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), a2.size());
  out()->add_block(odata, a1, a2, x0, x1);
}

void Task236::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index x0 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, x1, c1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, c1, x0), 0.0);
  {
    // tensor label: I334
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, c2, x0, x1);
    sort_indices<1,3,0,2,1,1,1,1>(i0data, odata, c1.size(), c2.size(), x0.size(), x1.size());
  }
  out()->add_block(odata, c2, x1, c1, x0);
}

void Task237::Task_local::compute() {
  const Index c1 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I334
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, c2, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, c2, x0, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c2, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma92
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x3, x1, x2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, x3, x1, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size(), x1.size(), x2.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, x3, c2, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, x3, c2, x2)]);
      sort_indices<1,3,0,2,0,1,2,1>(i1data, i1data_sorted, c1.size(), x3.size(), c2.size(), x2.size());
      zgemm3m_("T", "N", x0.size()*x1.size(), c1.size()*c2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c1.size(), c2.size());
  out()->add_block(odata, c1, c2, x0, x1);
}

void Task238::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x2, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x0, x1), 0.0);
  {
    // tensor label: I336
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, x2, x1, x0);
    sort_indices<0,1,3,2,1,1,1,1>(i0data, odata, c1.size(), x2.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, c1, x2, x0, x1);
}

void Task239::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I336
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma6
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x4, x2, x3, x1, x0);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x4, x2, x3, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x5, x4, c1, x3);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x5, x4, c1, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), c1.size(), x3.size());
        zgemm3m_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task240::Task_local::compute() {
  const Index c3 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, x0, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(c3, x0, c1, a2), 0.0);
  {
    // tensor label: I338
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a2, c1, x0);
    sort_indices<0,3,2,1,1,1,1,1>(i0data, odata, c3.size(), a2.size(), c1.size(), x0.size());
  }
  out()->add_block(odata, c3, x0, c1, a2);
}

void Task241::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x0 = b(3);
  // tensor label: I338
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, a2, c1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, c1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, a2, c1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a2, c1, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma16
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I339
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, a2, c1, x1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, a2, c1, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), c1.size(), x1.size());
    zgemm3m_("T", "N", x0.size(), c3.size()*a2.size()*c1.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x0.size(), c3.size(), a2.size(), c1.size());
  out()->add_block(odata, c3, a2, c1, x0);
}

void Task242::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x1 = b(3);
  // tensor label: I339
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, a2, c1, x1)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, c1, x1), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a2, c1, x1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c3.size(), a2.size(), c1.size(), x1.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(c1, a2, c3, x1);
    sort_indices<2,1,0,3,1,1,1,1>(i1data, odata, c1.size(), a2.size(), c3.size(), x1.size());
  }
  out()->add_block(odata, c3, a2, c1, x1);
}

void Task243::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  {
    // tensor label: I342
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a1, c2, x0, x1);
    sort_indices<1,3,2,0,1,1,1,1>(i0data, odata, a1.size(), c2.size(), x0.size(), x1.size());
  }
  out()->add_block(odata, c2, x1, x0, a1);
}

void Task244::Task_local::compute() {
  const Index a1 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I342
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, c2, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, c2, x0, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, c2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, c2, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma32
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0, x1, x2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x0, x1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, a1, c2, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, a1, c2, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), c2.size(), x2.size());
      zgemm3m_("T", "N", x0.size()*x1.size(), a1.size()*c2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), c2.size());
  out()->add_block(odata, a1, c2, x0, x1);
}

void Task245::Task_local::compute() {
  const Index a1 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I342
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, c2, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, c2, x0, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, c2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, c2, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma35
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, a1, x3, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, a1, x3, x2)]);
      sort_indices<2,3,0,1,0,1,-1,1>(i1data, i1data_sorted, c2.size(), a1.size(), x3.size(), x2.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), c2.size()*a1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a1.size());
  out()->add_block(odata, a1, c2, x0, x1);
}

void Task246::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x0, x1, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1, c1, a2), 0.0);
  {
    // tensor label: I346
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a2, c1, x1, x0);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, a2.size(), c1.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, x0, x1, c1, a2);
}

void Task247::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I346
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c1, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c1, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma35
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I347
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, a2, c1, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, a2, c1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a2.size(), c1.size(), x2.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), a2.size()*c1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a2.size(), c1.size());
  out()->add_block(odata, a2, c1, x1, x0);
}

void Task248::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x2 = b(3);
  // tensor label: I347
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x3, a2, c1, x2)]);
  std::fill_n(odata.get(), out()->get_size(x3, a2, c1, x2), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, a2, c1, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), a2.size(), c1.size(), x2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(c1, a2, x3, x2);
    sort_indices<2,1,0,3,1,1,1,1>(i1data, odata, c1.size(), a2.size(), x3.size(), x2.size());
  }
  out()->add_block(odata, x3, a2, c1, x2);
}

void Task249::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x1, x2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(x1, x2, x0, a1), 0.0);
  {
    // tensor label: I350
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a1, x0, x2, x1);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, a1.size(), x0.size(), x2.size(), x1.size());
  }
  out()->add_block(odata, x1, x2, x0, a1);
}

#endif
