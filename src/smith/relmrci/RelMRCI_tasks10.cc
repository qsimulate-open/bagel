//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_tasks10.cc
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

#include <src/smith/relmrci/RelMRCI_tasks10.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelMRCI;

void Task450::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& a5 : *range_[2]) {
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a5, a4);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(a5, a4)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a5.size(), a4.size());
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a2, c3, a5);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a2, c3, a5)]);
    sort_indices<3,0,1,2,0,1,2,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), a5.size());
    zgemm3m_("T", "N", a4.size(), c1.size()*a2.size()*c3.size(), a5.size(),
           1.0, i0data_sorted, a5.size(), i1data_sorted, a5.size(),
           1.0, odata_sorted, a4.size());
  }
  sort_indices<2,1,0,3,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), a2.size(), c3.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task451::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, a4, c1, a2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, a4, c1, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c1.size(), a2.size());
    // tensor label: I129
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, x1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), x1.size());
    zgemm3m_("T", "N", a4.size()*c1.size()*a2.size(), c3.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a4.size()*c1.size()*a2.size());
  }
  sort_indices<2,1,0,3,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), a2.size(), c3.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task452::Task_local::compute() {
  const Index c3 = b(0);
  const Index x1 = b(1);
  // tensor label: I129
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c3, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma27
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, x0);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, x0)]);
    sort_indices<1,0,0,1,-1,1>(i1data, i1data_sorted, c3.size(), x0.size());
    zgemm3m_("T", "N", x1.size(), c3.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), c3.size());
  out()->add_block(odata, c3, x1);
}

void Task453::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, a2, c1, a4);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, a2, c1, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c1.size(), a4.size());
    // tensor label: I132
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, x1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), x1.size());
    zgemm3m_("T", "N", a2.size()*c1.size()*a4.size(), c3.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a2.size()*c1.size()*a4.size());
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), a4.size(), c3.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task454::Task_local::compute() {
  const Index c3 = b(0);
  const Index x1 = b(1);
  // tensor label: I132
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c3, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma27
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, x0);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, x0)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), x0.size());
    zgemm3m_("T", "N", x1.size(), c3.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), c3.size());
  out()->add_block(odata, c3, x1);
}

void Task455::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, a4, c3, a2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, a4, c3, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a4.size(), c3.size(), a2.size());
    // tensor label: I777
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, x0);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, x0)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), x0.size());
    zgemm3m_("T", "N", a4.size()*c3.size()*a2.size(), c1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a4.size()*c3.size()*a2.size());
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a4.size(), c3.size(), a2.size(), c1.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task456::Task_local::compute() {
  const Index c1 = b(0);
  const Index x0 = b(1);
  // tensor label: I777
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma9
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x0, x1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x0, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, x2, c1, x1);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, x2, c1, x1)]);
        sort_indices<0,1,3,2,0,1,-1,1>(i1data, i1data_sorted, x3.size(), x2.size(), c1.size(), x1.size());
        zgemm3m_("T", "N", x0.size(), c1.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c1.size());
  out()->add_block(odata, c1, x0);
}

void Task457::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, a2, c3, a4);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, a2, c3, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a2.size(), c3.size(), a4.size());
    // tensor label: I780
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, x0);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, x0)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), x0.size());
    zgemm3m_("T", "N", a2.size()*c3.size()*a4.size(), c1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a2.size()*c3.size()*a4.size());
  }
  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, a2.size(), c3.size(), a4.size(), c1.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task458::Task_local::compute() {
  const Index c1 = b(0);
  const Index x0 = b(1);
  // tensor label: I780
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma9
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x0, x1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x0, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, x2, c1, x1);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, x2, c1, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), c1.size(), x1.size());
        zgemm3m_("T", "N", x0.size(), c1.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c1.size());
  out()->add_block(odata, c1, x0);
}

void Task459::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a4, c3, x3);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a4, c3, x3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x3.size());
    // tensor label: I783
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a2, x3);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a2, x3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a2.size(), x3.size());
    zgemm3m_("T", "N", c1.size()*a4.size()*c3.size(), a2.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, c1.size()*a4.size()*c3.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task460::Task_local::compute() {
  const Index a2 = b(0);
  const Index x3 = b(1);
  // tensor label: I783
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, x3)]);
  std::fill_n(odata.get(), out()->get_size(a2, x3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: Gamma5
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x2, x3, x1, x0);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x2, x3, x1, x0)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x2, a2, x1, x0);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x2, a2, x1, x0)]);
        sort_indices<0,2,3,1,0,1,-1,2>(i1data, i1data_sorted, x2.size(), a2.size(), x1.size(), x0.size());
        zgemm3m_("T", "N", x3.size(), a2.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), a2.size());
  out()->add_block(odata, a2, x3);
}

void Task461::Task_local::compute() {
  const Index a2 = b(0);
  const Index x3 = b(1);
  // tensor label: I783
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, x3)]);
  std::fill_n(odata.get(), out()->get_size(a2, x3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x3), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma160
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x3, x2, x1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, x3, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x2, x1, x0, a2);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x2, x1, x0, a2)]);
        sort_indices<2,0,1,3,0,1,-1,2>(i1data, i1data_sorted, x2.size(), x1.size(), x0.size(), a2.size());
        zgemm3m_("T", "N", x3.size(), a2.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), a2.size());
  out()->add_block(odata, a2, x3);
}

void Task462::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a2, c3, x3);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a2, c3, x3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x3.size());
    // tensor label: I786
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, x3);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, x3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a4.size(), x3.size());
    zgemm3m_("T", "N", c1.size()*a2.size()*c3.size(), a4.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c3.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task463::Task_local::compute() {
  const Index a4 = b(0);
  const Index x3 = b(1);
  // tensor label: I786
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a4, x3)]);
  std::fill_n(odata.get(), out()->get_size(a4, x3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: Gamma5
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x2, x3, x1, x0);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x2, x3, x1, x0)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x2, a4, x1, x0);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x2, a4, x1, x0)]);
        sort_indices<0,2,3,1,0,1,1,2>(i1data, i1data_sorted, x2.size(), a4.size(), x1.size(), x0.size());
        zgemm3m_("T", "N", x3.size(), a4.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), a4.size());
  out()->add_block(odata, a4, x3);
}

void Task464::Task_local::compute() {
  const Index a4 = b(0);
  const Index x3 = b(1);
  // tensor label: I786
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a4, x3)]);
  std::fill_n(odata.get(), out()->get_size(a4, x3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x3), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma160
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x3, x2, x1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, x3, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x2, x1, x0, a4);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x2, x1, x0, a4)]);
        sort_indices<2,0,1,3,0,1,1,2>(i1data, i1data_sorted, x2.size(), x1.size(), x0.size(), a4.size());
        zgemm3m_("T", "N", x3.size(), a4.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), a4.size());
  out()->add_block(odata, a4, x3);
}

void Task465::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& c5 : *range_[0]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c5, a4, c1, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c5, a4, c1, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, c5.size(), a4.size(), c1.size(), x1.size());
      // tensor label: I795
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c5, c3, a2, x1);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c5, c3, a2, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, c5.size(), c3.size(), a2.size(), x1.size());
      zgemm3m_("T", "N", a4.size()*c1.size(), c3.size()*a2.size(), c5.size()*x1.size(),
             1.0, i0data_sorted, c5.size()*x1.size(), i1data_sorted, c5.size()*x1.size(),
             1.0, odata_sorted, a4.size()*c1.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), c3.size(), a2.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task466::Task_local::compute() {
  const Index c5 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index x1 = b(3);
  // tensor label: I795
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c5, c3, a2, x1)]);
  std::fill_n(odata.get(), out()->get_size(c5, c3, a2, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c5, c3, a2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c5, c3, a2, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma11
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I796
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x0, c5, c3, a2);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x0, c5, c3, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c5.size(), c3.size(), a2.size());
    zgemm3m_("T", "N", x1.size(), c5.size()*c3.size()*a2.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), c5.size(), c3.size(), a2.size());
  out()->add_block(odata, c5, c3, a2, x1);
}

void Task467::Task_local::compute() {
  const Index x0 = b(0);
  const Index c5 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I796
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x0, c5, c3, a2)]);
  std::fill_n(odata.get(), out()->get_size(x0, c5, c3, a2), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, c5, c3, a2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x0.size(), c5.size(), c3.size(), a2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(x0, a2, c3, c5);
    sort_indices<0,3,2,1,1,1,-1,1>(i1data, odata, x0.size(), a2.size(), c3.size(), c5.size());
  }
  out()->add_block(odata, x0, c5, c3, a2);
}

void Task468::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& c5 : *range_[0]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c5, a2, c1, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c5, a2, c1, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, c5.size(), a2.size(), c1.size(), x1.size());
      // tensor label: I798
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c5, c3, a4, x1);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c5, c3, a4, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, c5.size(), c3.size(), a4.size(), x1.size());
      zgemm3m_("T", "N", a2.size()*c1.size(), c3.size()*a4.size(), c5.size()*x1.size(),
             1.0, i0data_sorted, c5.size()*x1.size(), i1data_sorted, c5.size()*x1.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), c3.size(), a4.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task469::Task_local::compute() {
  const Index c5 = b(0);
  const Index c3 = b(1);
  const Index a4 = b(2);
  const Index x1 = b(3);
  // tensor label: I798
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c5, c3, a4, x1)]);
  std::fill_n(odata.get(), out()->get_size(c5, c3, a4, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c5, c3, a4, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c5, c3, a4, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma11
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I799
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x0, c5, c3, a4);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x0, c5, c3, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c5.size(), c3.size(), a4.size());
    zgemm3m_("T", "N", x1.size(), c5.size()*c3.size()*a4.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), c5.size(), c3.size(), a4.size());
  out()->add_block(odata, c5, c3, a4, x1);
}

void Task470::Task_local::compute() {
  const Index x0 = b(0);
  const Index c5 = b(1);
  const Index c3 = b(2);
  const Index a4 = b(3);
  // tensor label: I799
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x0, c5, c3, a4)]);
  std::fill_n(odata.get(), out()->get_size(x0, c5, c3, a4), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, c5, c3, a4);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x0.size(), c5.size(), c3.size(), a4.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(x0, a4, c3, c5);
    sort_indices<0,3,2,1,1,1,1,1>(i1data, odata, x0.size(), a4.size(), c3.size(), c5.size());
  }
  out()->add_block(odata, x0, c5, c3, a4);
}

void Task471::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& c5 : *range_[0]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a4, c5, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a4, c5, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c5.size(), x1.size());
      // tensor label: I807
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c5, c3, a2, x1);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c5, c3, a2, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, c5.size(), c3.size(), a2.size(), x1.size());
      zgemm3m_("T", "N", c1.size()*a4.size(), c3.size()*a2.size(), c5.size()*x1.size(),
             1.0, i0data_sorted, c5.size()*x1.size(), i1data_sorted, c5.size()*x1.size(),
             1.0, odata_sorted, c1.size()*a4.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task472::Task_local::compute() {
  const Index c5 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index x1 = b(3);
  // tensor label: I807
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c5, c3, a2, x1)]);
  std::fill_n(odata.get(), out()->get_size(c5, c3, a2, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c5, c3, a2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c5, c3, a2, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma11
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I808
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x0, c5, c3, a2);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x0, c5, c3, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c5.size(), c3.size(), a2.size());
    zgemm3m_("T", "N", x1.size(), c5.size()*c3.size()*a2.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), c5.size(), c3.size(), a2.size());
  out()->add_block(odata, c5, c3, a2, x1);
}

void Task473::Task_local::compute() {
  const Index x0 = b(0);
  const Index c5 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I808
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x0, c5, c3, a2)]);
  std::fill_n(odata.get(), out()->get_size(x0, c5, c3, a2), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, c5, c3, a2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x0.size(), c5.size(), c3.size(), a2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(x0, a2, c3, c5);
    sort_indices<0,3,2,1,1,1,1,1>(i1data, odata, x0.size(), a2.size(), c3.size(), c5.size());
  }
  out()->add_block(odata, x0, c5, c3, a2);
}

void Task474::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& c5 : *range_[0]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a2, c5, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a2, c5, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c5.size(), x1.size());
      // tensor label: I810
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c5, c3, a4, x1);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c5, c3, a4, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, c5.size(), c3.size(), a4.size(), x1.size());
      zgemm3m_("T", "N", c1.size()*a2.size(), c3.size()*a4.size(), c5.size()*x1.size(),
             1.0, i0data_sorted, c5.size()*x1.size(), i1data_sorted, c5.size()*x1.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task475::Task_local::compute() {
  const Index c5 = b(0);
  const Index c3 = b(1);
  const Index a4 = b(2);
  const Index x1 = b(3);
  // tensor label: I810
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c5, c3, a4, x1)]);
  std::fill_n(odata.get(), out()->get_size(c5, c3, a4, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c5, c3, a4, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c5, c3, a4, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma11
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I811
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x0, c5, c3, a4);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x0, c5, c3, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c5.size(), c3.size(), a4.size());
    zgemm3m_("T", "N", x1.size(), c5.size()*c3.size()*a4.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), c5.size(), c3.size(), a4.size());
  out()->add_block(odata, c5, c3, a4, x1);
}

void Task476::Task_local::compute() {
  const Index x0 = b(0);
  const Index c5 = b(1);
  const Index c3 = b(2);
  const Index a4 = b(3);
  // tensor label: I811
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x0, c5, c3, a4)]);
  std::fill_n(odata.get(), out()->get_size(x0, c5, c3, a4), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, c5, c3, a4);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x0.size(), c5.size(), c3.size(), a4.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(x0, a4, c3, c5);
    sort_indices<0,3,2,1,1,1,-1,1>(i1data, odata, x0.size(), a4.size(), c3.size(), c5.size());
  }
  out()->add_block(odata, x0, c5, c3, a4);
}

void Task477::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, a4, a5, a2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, a4, a5, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a4.size(), a5.size(), a2.size());
      // tensor label: I819
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a5, c3, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a5, c3, x0)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), a5.size(), c3.size(), x0.size());
      zgemm3m_("T", "N", a4.size()*a2.size(), c1.size()*c3.size(), a5.size()*x0.size(),
             1.0, i0data_sorted, a5.size()*x0.size(), i1data_sorted, a5.size()*x0.size(),
             1.0, odata_sorted, a4.size()*a2.size());
    }
  }
  sort_indices<1,2,0,3,1,1,1,1>(odata_sorted, odata, a4.size(), a2.size(), c1.size(), c3.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task478::Task_local::compute() {
  const Index c1 = b(0);
  const Index a5 = b(1);
  const Index c3 = b(2);
  const Index x0 = b(3);
  // tensor label: I819
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a5, c3, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, a5, c3, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a5, c3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a5, c3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma11
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a5, c3, x1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a5, c3, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), a5.size(), c3.size(), x1.size());
    zgemm3m_("T", "N", x0.size(), c1.size()*a5.size()*c3.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x0.size(), c1.size(), a5.size(), c3.size());
  out()->add_block(odata, c1, a5, c3, x0);
}

void Task479::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, a2, a5, a4);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, a2, a5, a4)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a2.size(), a5.size(), a4.size());
      // tensor label: I822
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a5, c3, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a5, c3, x0)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), a5.size(), c3.size(), x0.size());
      zgemm3m_("T", "N", a2.size()*a4.size(), c1.size()*c3.size(), a5.size()*x0.size(),
             1.0, i0data_sorted, a5.size()*x0.size(), i1data_sorted, a5.size()*x0.size(),
             1.0, odata_sorted, a2.size()*a4.size());
    }
  }
  sort_indices<0,2,1,3,1,1,1,1>(odata_sorted, odata, a2.size(), a4.size(), c1.size(), c3.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task480::Task_local::compute() {
  const Index c1 = b(0);
  const Index a5 = b(1);
  const Index c3 = b(2);
  const Index x0 = b(3);
  // tensor label: I822
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a5, c3, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, a5, c3, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a5, c3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a5, c3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma11
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a5, c3, x1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a5, c3, x1)]);
    sort_indices<3,0,1,2,0,1,-1,1>(i1data, i1data_sorted, c1.size(), a5.size(), c3.size(), x1.size());
    zgemm3m_("T", "N", x0.size(), c1.size()*a5.size()*c3.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x0.size(), c1.size(), a5.size(), c3.size());
  out()->add_block(odata, c1, a5, c3, x0);
}

void Task481::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a4, x3, x2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a4, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), x3.size(), x2.size());
      // tensor label: I825
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, a2, x3, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, a2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), x3.size(), x2.size());
      zgemm3m_("T", "N", c1.size()*a4.size(), c3.size()*a2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c1.size()*a4.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task482::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I825
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, a2, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, x3, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, a2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a2, x3, x2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma24
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I826
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, a2, x1, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, a2, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), x1.size(), x0.size());
      zgemm3m_("T", "N", x3.size()*x2.size(), c3.size()*a2.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c3.size(), a2.size());
  out()->add_block(odata, c3, a2, x3, x2);
}

void Task483::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I826
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, a2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, x1, x0), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a2, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, c3.size(), a2.size(), x1.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(x1, a2, c3, x0);
    sort_indices<2,1,0,3,1,1,1,2>(i1data, odata, x1.size(), a2.size(), c3.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i2data = in(0)->get_block(x1, x0, c3, a2);
    sort_indices<2,3,0,1,1,1,-1,2>(i2data, odata, x1.size(), x0.size(), c3.size(), a2.size());
  }
  out()->add_block(odata, c3, a2, x1, x0);
}

void Task484::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I825
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, a2, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, x3, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, a2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a2, x3, x2), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma9
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, x1, x0, a2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, x1, x0, a2)]);
      sort_indices<2,1,0,3,0,1,-1,2>(i1data, i1data_sorted, c3.size(), x1.size(), x0.size(), a2.size());
      zgemm3m_("T", "N", x3.size()*x2.size(), c3.size()*a2.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c3.size(), a2.size());
  out()->add_block(odata, c3, a2, x3, x2);
}

void Task485::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a2, x3, x2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x3.size(), x2.size());
      // tensor label: I828
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, a4, x3, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, a4, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), a4.size(), x3.size(), x2.size());
      zgemm3m_("T", "N", c1.size()*a2.size(), c3.size()*a4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task486::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I828
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, a4, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c3, a4, x3, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, a4, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a4, x3, x2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma24
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I829
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, a4, x1, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, a4, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), a4.size(), x1.size(), x0.size());
      zgemm3m_("T", "N", x3.size()*x2.size(), c3.size()*a4.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c3.size(), a4.size());
  out()->add_block(odata, c3, a4, x3, x2);
}

void Task487::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I829
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, a4, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, a4, x1, x0), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a4, x1, x0);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, c3.size(), a4.size(), x1.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(x1, a4, c3, x0);
    sort_indices<2,1,0,3,1,1,-1,2>(i1data, odata, x1.size(), a4.size(), c3.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i2data = in(0)->get_block(x1, x0, c3, a4);
    sort_indices<2,3,0,1,1,1,1,2>(i2data, odata, x1.size(), x0.size(), c3.size(), a4.size());
  }
  out()->add_block(odata, c3, a4, x1, x0);
}

void Task488::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I828
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, a4, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c3, a4, x3, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, a4, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a4, x3, x2), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma9
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, x1, x0, a4);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, x1, x0, a4)]);
      sort_indices<2,1,0,3,0,1,1,2>(i1data, i1data_sorted, c3.size(), x1.size(), x0.size(), a4.size());
      zgemm3m_("T", "N", x3.size()*x2.size(), c3.size()*a4.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c3.size(), a4.size());
  out()->add_block(odata, c3, a4, x3, x2);
}

void Task489::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& c5 : *range_[0]) {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a2, c3, c5);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a2, c3, c5)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), c5.size());
    // tensor label: I849
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c5, a4);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c5, a4)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c5.size(), a4.size());
    zgemm3m_("T", "N", c1.size()*a2.size()*c3.size(), a4.size(), c5.size(),
           1.0, i0data_sorted, c5.size(), i1data_sorted, c5.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c3.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task490::Task_local::compute() {
  const Index c5 = b(0);
  const Index a4 = b(1);
  // tensor label: I849
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c5, a4)]);
  std::fill_n(odata.get(), out()->get_size(c5, a4), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c5, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c5, a4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma27
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c5, a4, x1, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c5, a4, x1, x0)]);
      sort_indices<2,3,0,1,0,1,-1,1>(i1data, i1data_sorted, c5.size(), a4.size(), x1.size(), x0.size());
      zgemm3m_("T", "N", 1, c5.size()*a4.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c5.size(), a4.size());
  out()->add_block(odata, c5, a4);
}

void Task491::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& c5 : *range_[0]) {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a4, c3, c5);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a4, c3, c5)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), c5.size());
    // tensor label: I852
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c5, a2);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c5, a2)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c5.size(), a2.size());
    zgemm3m_("T", "N", c1.size()*a4.size()*c3.size(), a2.size(), c5.size(),
           1.0, i0data_sorted, c5.size(), i1data_sorted, c5.size(),
           1.0, odata_sorted, c1.size()*a4.size()*c3.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task492::Task_local::compute() {
  const Index c5 = b(0);
  const Index a2 = b(1);
  // tensor label: I852
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c5, a2)]);
  std::fill_n(odata.get(), out()->get_size(c5, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c5, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c5, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma27
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c5, a2, x1, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c5, a2, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c5.size(), a2.size(), x1.size(), x0.size());
      zgemm3m_("T", "N", 1, c5.size()*a2.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c5.size(), a2.size());
  out()->add_block(odata, c5, a2);
}

void Task493::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& a5 : *range_[2]) {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a4, a5, a2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c3, a4, a5, a2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), a5.size(), a2.size());
    // tensor label: I855
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a5);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a5)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), a5.size());
    zgemm3m_("T", "N", c3.size()*a4.size()*a2.size(), c1.size(), a5.size(),
           1.0, i0data_sorted, a5.size(), i1data_sorted, a5.size(),
           1.0, odata_sorted, c3.size()*a4.size()*a2.size());
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, c3.size(), a4.size(), a2.size(), c1.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task494::Task_local::compute() {
  const Index c1 = b(0);
  const Index a5 = b(1);
  // tensor label: I855
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a5)]);
  std::fill_n(odata.get(), out()->get_size(c1, a5), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a5), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma27
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a5, x1, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a5, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), a5.size(), x1.size(), x0.size());
      zgemm3m_("T", "N", 1, c1.size()*a5.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a5.size());
  out()->add_block(odata, c1, a5);
}

void Task495::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& a5 : *range_[2]) {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a2, a5, a4);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c3, a2, a5, a4)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), a5.size(), a4.size());
    // tensor label: I858
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a5);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a5)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), a5.size());
    zgemm3m_("T", "N", c3.size()*a2.size()*a4.size(), c1.size(), a5.size(),
           1.0, i0data_sorted, a5.size(), i1data_sorted, a5.size(),
           1.0, odata_sorted, c3.size()*a2.size()*a4.size());
  }
  sort_indices<1,3,2,0,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), a4.size(), c1.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task496::Task_local::compute() {
  const Index c1 = b(0);
  const Index a5 = b(1);
  // tensor label: I858
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a5)]);
  std::fill_n(odata.get(), out()->get_size(c1, a5), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a5), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma27
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a5, x1, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a5, x1, x0)]);
      sort_indices<2,3,0,1,0,1,-1,1>(i1data, i1data_sorted, c1.size(), a5.size(), x1.size(), x0.size());
      zgemm3m_("T", "N", 1, c1.size()*a5.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a5.size());
  out()->add_block(odata, c1, a5);
}

void Task497::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, x0, c3, a2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, x0, c3, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x0.size(), c3.size(), a2.size());
    // tensor label: I861
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, x0);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, x0)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a4.size(), x0.size());
    zgemm3m_("T", "N", c1.size()*c3.size()*a2.size(), a4.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, c1.size()*c3.size()*a2.size());
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), a2.size(), a4.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task498::Task_local::compute() {
  const Index a4 = b(0);
  const Index x0 = b(1);
  // tensor label: I861
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a4, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma33
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0, x2, x1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x0, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, a4, x2, x1);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, a4, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), a4.size(), x2.size(), x1.size());
        zgemm3m_("T", "N", x0.size(), a4.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), a4.size());
  out()->add_block(odata, a4, x0);
}

void Task499::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, x0, c3, a4);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, x0, c3, a4)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x0.size(), c3.size(), a4.size());
    // tensor label: I864
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a2, x0);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a2, x0)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a2.size(), x0.size());
    zgemm3m_("T", "N", c1.size()*c3.size()*a4.size(), a2.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, c1.size()*c3.size()*a4.size());
  }
  sort_indices<3,0,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), a4.size(), a2.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

#endif
