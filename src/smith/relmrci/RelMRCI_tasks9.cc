//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_tasks9.cc
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

#include <src/smith/relmrci/RelMRCI_tasks9.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelMRCI;

void Task400::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  // tensor label: I100
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x3, a1)]);
  std::fill_n(odata.get(), out()->get_size(x3, a1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& c4 : *range_[0]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c2, a1, c4, a3);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c2, a1, c4, a3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c4.size(), a3.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, c4, a3, c2);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, c4, a3, c2)]);
        sort_indices<3,1,2,0,0,1,2,1>(i1data, i1data_sorted, x3.size(), c4.size(), a3.size(), c2.size());
        zgemm3m_("T", "N", a1.size(), x3.size(), c4.size()*a3.size()*c2.size(),
               1.0, i0data_sorted, c4.size()*a3.size()*c2.size(), i1data_sorted, c4.size()*a3.size()*c2.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a1.size(), x3.size());
  out()->add_block(odata, x3, a1);
}

void Task401::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  // tensor label: I100
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x3, a1)]);
  std::fill_n(odata.get(), out()->get_size(x3, a1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, a4, c2, a3);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, a4, c2, a3)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x3.size(), a4.size(), c2.size(), a3.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, a1, a3, c2);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, a1, a3, c2)]);
        sort_indices<0,3,2,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a1.size(), a3.size(), c2.size());
        zgemm3m_("T", "N", x3.size(), a1.size(), a4.size()*a3.size()*c2.size(),
               1.0, i0data_sorted, a4.size()*a3.size()*c2.size(), i1data_sorted, a4.size()*a3.size()*c2.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size());
  out()->add_block(odata, x3, a1);
}

void Task402::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  // tensor label: I100
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x3, a1)]);
  std::fill_n(odata.get(), out()->get_size(x3, a1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a4 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, a3, c2, a4);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, a3, c2, a4)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a4.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, a1, a3, c2);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, a1, a3, c2)]);
        sort_indices<2,3,0,1,0,1,-1,1>(i1data, i1data_sorted, a4.size(), a1.size(), a3.size(), c2.size());
        zgemm3m_("T", "N", x3.size(), a1.size(), a4.size()*a3.size()*c2.size(),
               1.0, i0data_sorted, a4.size()*a3.size()*c2.size(), i1data_sorted, a4.size()*a3.size()*c2.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size());
  out()->add_block(odata, x3, a1);
}

void Task403::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I93
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& c2 : *range_[0]) {
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x4, a1, x3, c2);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x4, a1, x3, c2)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x4.size(), a1.size(), x3.size(), c2.size());
        // tensor label: I699
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, x4, x0, x2, x1, c2);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, x4, x0, x2, x1, c2)]);
        sort_indices<1,0,5,2,3,4,0,1,1,1>(i1data, i1data_sorted, x3.size(), x4.size(), x0.size(), x2.size(), x1.size(), c2.size());
        zgemm3m_("T", "N", a1.size(), x0.size()*x2.size()*x1.size(), x3.size()*x4.size()*c2.size(),
               1.0, i0data_sorted, x3.size()*x4.size()*c2.size(), i1data_sorted, x3.size()*x4.size()*c2.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x2.size(), x1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task404::Task_local::compute() {
  const Index x3 = b(0);
  const Index x4 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  const Index c2 = b(5);
  // tensor label: I699
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x3, x4, x0, x2, x1, c2)]);
  std::fill_n(odata.get(), out()->get_size(x3, x4, x0, x2, x1, c2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x3, x4, x0, x2, x1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x4, x0, x2, x1, c2), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x6, c2, x5);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x7, x6, c2, x5)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), c2.size(), x5.size());
        // tensor label: Gamma230
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x7, x6, x3, x5, x4, x0, x2, x1);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x7, x6, x3, x5, x4, x0, x2, x1)]);
        sort_indices<0,1,3,2,4,5,6,7,0,1,-1,1>(i1data, i1data_sorted, x7.size(), x6.size(), x3.size(), x5.size(), x4.size(), x0.size(), x2.size(), x1.size());
        zgemm3m_("T", "N", c2.size(), x3.size()*x4.size()*x0.size()*x2.size()*x1.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, c2.size());
      }
    }
  }
  sort_indices<1,2,3,4,5,0,1,1,1,1>(odata_sorted, odata, c2.size(), x3.size(), x4.size(), x0.size(), x2.size(), x1.size());
  out()->add_block(odata, x3, x4, x0, x2, x1, c2);
}

void Task405::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I93
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma231
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x4, x5, x3, x0, x2, x1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x4, x5, x3, x0, x2, x1)]);
        sort_indices<0,1,2,3,4,5,0,1,1,1>(i0data, i0data_sorted, x4.size(), x5.size(), x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: I702
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x4, x3, a1, x5);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x4, x3, a1, x5)]);
        sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x4.size(), x3.size(), a1.size(), x5.size());
        zgemm3m_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x4.size()*x3.size()*x5.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x5.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task406::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index x5 = b(3);
  // tensor label: I702
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x4, x3, a1, x5)]);
  std::fill_n(odata.get(), out()->get_size(x4, x3, a1, x5), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x4, x3, a1, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x3, a1, x5), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c2, a1, c3, x5);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c2, a1, c3, x5)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c3.size(), x5.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x4, c3, x3, c2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x4, c3, x3, c2)]);
      sort_indices<3,1,0,2,0,1,-1,1>(i1data, i1data_sorted, x4.size(), c3.size(), x3.size(), c2.size());
      zgemm3m_("T", "N", a1.size()*x5.size(), x4.size()*x3.size(), c3.size()*c2.size(),
             1.0, i0data_sorted, c3.size()*c2.size(), i1data_sorted, c3.size()*c2.size(),
             1.0, odata_sorted, a1.size()*x5.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, x4, x3, a1, x5);
}

void Task407::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I93
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& x7 : *range_[1]) {
      for (auto& x6 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c2, a1, x7, x6);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c2, a1, x7, x6)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x7.size(), x6.size());
        // tensor label: I705
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, x7, x6, x0, x2, x1);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, x7, x6, x0, x2, x1)]);
        sort_indices<0,1,2,3,4,5,0,1,1,1>(i1data, i1data_sorted, c2.size(), x7.size(), x6.size(), x0.size(), x2.size(), x1.size());
        zgemm3m_("T", "N", a1.size(), x0.size()*x2.size()*x1.size(), c2.size()*x7.size()*x6.size(),
               1.0, i0data_sorted, c2.size()*x7.size()*x6.size(), i1data_sorted, c2.size()*x7.size()*x6.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x2.size(), x1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task408::Task_local::compute() {
  const Index c2 = b(0);
  const Index x7 = b(1);
  const Index x6 = b(2);
  const Index x0 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);
  // tensor label: I705
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, x7, x6, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x7, x6, x0, x2, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, x7, x6, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x7, x6, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma232
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x6, x5, x0, x4, x3, x2, x1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x7, x6, x5, x0, x4, x3, x2, x1)]);
        sort_indices<2,4,5,0,1,3,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x5, c2, x4, x3);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x5, c2, x4, x3)]);
        sort_indices<0,2,3,1,0,1,-1,2>(i1data, i1data_sorted, x5.size(), c2.size(), x4.size(), x3.size());
        zgemm3m_("T", "N", x7.size()*x6.size()*x0.size()*x2.size()*x1.size(), c2.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x7.size()*x6.size()*x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x7.size(), x6.size(), x0.size(), x2.size(), x1.size(), c2.size());
  out()->add_block(odata, c2, x7, x6, x0, x2, x1);
}

void Task409::Task_local::compute() {
  const Index c2 = b(0);
  const Index x7 = b(1);
  const Index x6 = b(2);
  const Index x0 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);
  // tensor label: I705
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, x7, x6, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x7, x6, x0, x2, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, x7, x6, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x7, x6, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma233
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x6, x5, x4, x3, x0, x2, x1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x7, x6, x5, x4, x3, x0, x2, x1)]);
        sort_indices<2,3,4,0,1,5,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x5, x4, x3, c2);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x5, x4, x3, c2)]);
        sort_indices<0,1,2,3,0,1,-1,2>(i1data, i1data_sorted, x5.size(), x4.size(), x3.size(), c2.size());
        zgemm3m_("T", "N", x7.size()*x6.size()*x0.size()*x2.size()*x1.size(), c2.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x7.size()*x6.size()*x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x7.size(), x6.size(), x0.size(), x2.size(), x1.size(), c2.size());
  out()->add_block(odata, c2, x7, x6, x0, x2, x1);
}

void Task410::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I93
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        for (auto& x4 : *range_[1]) {
          for (auto& x3 : *range_[1]) {
            // tensor label: Gamma236
            std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x6, x5, x4, x3, x2, x1);
            std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x7, x0, x6, x5, x4, x3, x2, x1)]);
            sort_indices<0,2,3,4,5,1,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x6.size(), x5.size(), x4.size(), x3.size(), x2.size(), x1.size());
            // tensor label: I717
            std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a1, x4, x3, x7, x6, x5);
            std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a1, x4, x3, x7, x6, x5)]);
            sort_indices<3,4,5,1,2,0,0,1,1,1>(i1data, i1data_sorted, a1.size(), x4.size(), x3.size(), x7.size(), x6.size(), x5.size());
            zgemm3m_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, i0data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, odata_sorted, x0.size()*x2.size()*x1.size());
          }
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task411::Task_local::compute() {
  const Index a1 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x7 = b(3);
  const Index x6 = b(4);
  const Index x5 = b(5);
  // tensor label: I717
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x4, x3, x7, x6, x5)]);
  std::fill_n(odata.get(), out()->get_size(a1, x4, x3, x7, x6, x5), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x4, x3, x7, x6, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x4, x3, x7, x6, x5), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, a2, x6, x5);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x7, a2, x6, x5)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), a2.size(), x6.size(), x5.size());
    // tensor label: I718
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a2, a1, x4, x3);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a2, a1, x4, x3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), a1.size(), x4.size(), x3.size());
    zgemm3m_("T", "N", x7.size()*x6.size()*x5.size(), a1.size()*x4.size()*x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x7.size()*x6.size()*x5.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x7.size(), x6.size(), x5.size(), a1.size(), x4.size(), x3.size());
  out()->add_block(odata, a1, x4, x3, x7, x6, x5);
}

void Task412::Task_local::compute() {
  const Index a2 = b(0);
  const Index a1 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I718
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, a1, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(a2, a1, x4, x3), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a2, a1, x4, x3);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, a2.size(), a1.size(), x4.size(), x3.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(x4, x3, a2, a1);
    sort_indices<2,3,0,1,1,1,1,2>(i1data, odata, x4.size(), x3.size(), a2.size(), a1.size());
  }
  out()->add_block(odata, a2, a1, x4, x3);
}

void Task413::Task_local::compute() {
  const Index a1 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x7 = b(3);
  const Index x6 = b(4);
  const Index x5 = b(5);
  // tensor label: I717
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x4, x3, x7, x6, x5)]);
  std::fill_n(odata.get(), out()->get_size(a1, x4, x3, x7, x6, x5), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x4, x3, x7, x6, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x4, x3, x7, x6, x5), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, a1, x6, a2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x7, a1, x6, a2)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x7.size(), a1.size(), x6.size(), a2.size());
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a2, x5, x4, x3);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a2, x5, x4, x3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), x5.size(), x4.size(), x3.size());
    zgemm3m_("T", "N", x7.size()*a1.size()*x6.size(), x5.size()*x4.size()*x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x7.size()*a1.size()*x6.size());
  }
  sort_indices<1,4,5,0,2,3,1,1,1,1>(odata_sorted, odata, x7.size(), a1.size(), x6.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, a1, x4, x3, x7, x6, x5);
}

void Task414::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I93
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x6 : *range_[1]) {
        for (auto& x5 : *range_[1]) {
          for (auto& x3 : *range_[1]) {
            // tensor label: Gamma237
            std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x4, x6, x5, x3, x0, x2, x1);
            std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x7, x4, x6, x5, x3, x0, x2, x1)]);
            sort_indices<0,1,2,3,4,5,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x4.size(), x6.size(), x5.size(), x3.size(), x0.size(), x2.size(), x1.size());
            // tensor label: I720
            std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x4, x3, a1, x7, x6, x5);
            std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x4, x3, a1, x7, x6, x5)]);
            sort_indices<3,0,4,5,1,2,0,1,1,1>(i1data, i1data_sorted, x4.size(), x3.size(), a1.size(), x7.size(), x6.size(), x5.size());
            zgemm3m_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, i0data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, odata_sorted, x0.size()*x2.size()*x1.size());
          }
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task415::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index x7 = b(3);
  const Index x6 = b(4);
  const Index x5 = b(5);
  // tensor label: I720
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x4, x3, a1, x7, x6, x5)]);
  std::fill_n(odata.get(), out()->get_size(x4, x3, a1, x7, x6, x5), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x4, x3, a1, x7, x6, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x3, a1, x7, x6, x5), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, a2, x6, x5);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x7, a2, x6, x5)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), a2.size(), x6.size(), x5.size());
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a2, x4, x3, a1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a2, x4, x3, a1)]);
    sort_indices<0,1,2,3,0,1,1,2>(i1data, i1data_sorted, a2.size(), x4.size(), x3.size(), a1.size());
    zgemm3m_("T", "N", x7.size()*x6.size()*x5.size(), x4.size()*x3.size()*a1.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x7.size()*x6.size()*x5.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x7.size(), x6.size(), x5.size(), x4.size(), x3.size(), a1.size());
  out()->add_block(odata, x4, x3, a1, x7, x6, x5);
}

void Task416::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I93
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x6 : *range_[1]) {
        for (auto& x5 : *range_[1]) {
          for (auto& x4 : *range_[1]) {
            // tensor label: Gamma238
            std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x3, x6, x5, x4, x0, x2, x1);
            std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x7, x3, x6, x5, x4, x0, x2, x1)]);
            sort_indices<0,1,2,3,4,5,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x3.size(), x6.size(), x5.size(), x4.size(), x0.size(), x2.size(), x1.size());
            // tensor label: I723
            std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x4, a1, x3, x7, x6, x5);
            std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x4, a1, x3, x7, x6, x5)]);
            sort_indices<3,2,4,5,0,1,0,1,1,1>(i1data, i1data_sorted, x4.size(), a1.size(), x3.size(), x7.size(), x6.size(), x5.size());
            zgemm3m_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, i0data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, odata_sorted, x0.size()*x2.size()*x1.size());
          }
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task417::Task_local::compute() {
  const Index x4 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x7 = b(3);
  const Index x6 = b(4);
  const Index x5 = b(5);
  // tensor label: I723
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x4, a1, x3, x7, x6, x5)]);
  std::fill_n(odata.get(), out()->get_size(x4, a1, x3, x7, x6, x5), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x4, a1, x3, x7, x6, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, a1, x3, x7, x6, x5), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, a2, x6, x5);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x7, a2, x6, x5)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), a2.size(), x6.size(), x5.size());
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x4, a1, a2, x3);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x4, a1, a2, x3)]);
    sort_indices<2,0,1,3,0,1,-1,2>(i1data, i1data_sorted, x4.size(), a1.size(), a2.size(), x3.size());
    zgemm3m_("T", "N", x7.size()*x6.size()*x5.size(), x4.size()*a1.size()*x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x7.size()*x6.size()*x5.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x7.size(), x6.size(), x5.size(), x4.size(), a1.size(), x3.size());
  out()->add_block(odata, x4, a1, x3, x7, x6, x5);
}

void Task418::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I93
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: Gamma245
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x3, x4, x2, x1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x0, x3, x4, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x3.size(), x4.size(), x2.size(), x1.size());
        // tensor label: I744
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x4, x3, x5, a1);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x4, x3, x5, a1)]);
        sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x4.size(), x3.size(), x5.size(), a1.size());
        zgemm3m_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x4.size()*x3.size()*x5.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x5.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task419::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index x5 = b(2);
  const Index a1 = b(3);
  // tensor label: I744
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x4, x3, x5, a1)]);
  std::fill_n(odata.get(), out()->get_size(x4, x3, x5, a1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x4, x3, x5, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x3, x5, a1), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, a1, c2, a3);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, a1, c2, a3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c2.size(), a3.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a3, x4, x3, c2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a3, x4, x3, c2)]);
      sort_indices<3,0,1,2,0,1,1,2>(i1data, i1data_sorted, a3.size(), x4.size(), x3.size(), c2.size());
      zgemm3m_("T", "N", x5.size()*a1.size(), x4.size()*x3.size(), a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, x5.size()*a1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), x3.size());
  out()->add_block(odata, x4, x3, x5, a1);
}

void Task420::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I93
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: Gamma246
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x3, x4, x0, x2, x1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x3, x4, x0, x2, x1)]);
        sort_indices<0,1,2,3,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x3.size(), x4.size(), x0.size(), x2.size(), x1.size());
        // tensor label: I747
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x4, x3, x5, a1);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x4, x3, x5, a1)]);
        sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x4.size(), x3.size(), x5.size(), a1.size());
        zgemm3m_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x4.size()*x3.size()*x5.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x5.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task421::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index x5 = b(2);
  const Index a1 = b(3);
  // tensor label: I747
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x4, x3, x5, a1)]);
  std::fill_n(odata.get(), out()->get_size(x4, x3, x5, a1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x4, x3, x5, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x3, x5, a1), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, a2, c3, a1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, a2, c3, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), c3.size(), a1.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x4, c3, a2, x3);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x4, c3, a2, x3)]);
      sort_indices<2,1,0,3,0,1,1,2>(i1data, i1data_sorted, x4.size(), c3.size(), a2.size(), x3.size());
      zgemm3m_("T", "N", x5.size()*a1.size(), x4.size()*x3.size(), c3.size()*a2.size(),
             1.0, i0data_sorted, c3.size()*a2.size(), i1data_sorted, c3.size()*a2.size(),
             1.0, odata_sorted, x5.size()*a1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), x3.size());
  out()->add_block(odata, x4, x3, x5, a1);
}

void Task422::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I93
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        for (auto& x5 : *range_[1]) {
          for (auto& x4 : *range_[1]) {
            // tensor label: Gamma253
            std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x6, x3, x5, x4, x2, x1);
            std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x7, x0, x6, x3, x5, x4, x2, x1)]);
            sort_indices<0,2,3,4,5,1,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x6.size(), x3.size(), x5.size(), x4.size(), x2.size(), x1.size());
            // tensor label: I768
            std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x5, x4, x3, x7, a1, x6);
            std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x5, x4, x3, x7, a1, x6)]);
            sort_indices<3,5,2,0,1,4,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x3.size(), x7.size(), a1.size(), x6.size());
            zgemm3m_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x5.size()*x4.size()*x3.size()*x7.size()*x6.size(),
                   1.0, i0data_sorted, x5.size()*x4.size()*x3.size()*x7.size()*x6.size(), i1data_sorted, x5.size()*x4.size()*x3.size()*x7.size()*x6.size(),
                   1.0, odata_sorted, x0.size()*x2.size()*x1.size());
          }
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task423::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x7 = b(3);
  const Index a1 = b(4);
  const Index x6 = b(5);
  // tensor label: I768
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x5, x4, x3, x7, a1, x6)]);
  std::fill_n(odata.get(), out()->get_size(x5, x4, x3, x7, a1, x6), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x5, x4, x3, x7, a1, x6)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x5, x4, x3, x7, a1, x6), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, a1, x6, a2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x7, a1, x6, a2)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x7.size(), a1.size(), x6.size(), a2.size());
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x5, x4, a2, x3);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x5, x4, a2, x3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), a2.size(), x3.size());
    zgemm3m_("T", "N", x7.size()*a1.size()*x6.size(), x5.size()*x4.size()*x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x7.size()*a1.size()*x6.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x7.size(), a1.size(), x6.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, x5, x4, x3, x7, a1, x6);
}

void Task424::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I93
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: Gamma422
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x6, x5, x2, x1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x7, x0, x6, x5, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x6.size(), x5.size(), x2.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x7, a1, x6, x5);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x7, a1, x6, x5)]);
        sort_indices<0,2,3,1,0,1,-1,1>(i1data, i1data_sorted, x7.size(), a1.size(), x6.size(), x5.size());
        zgemm3m_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task425::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I93
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x9 : *range_[1]) {
    for (auto& x8 : *range_[1]) {
      for (auto& x7 : *range_[1]) {
        // tensor label: Gamma423
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x9, x0, x8, x7, x2, x1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x9, x0, x8, x7, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x9.size(), x0.size(), x8.size(), x7.size(), x2.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x9, a1, x8, x7);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x9, a1, x8, x7)]);
        sort_indices<0,2,3,1,0,1,-1,2>(i1data, i1data_sorted, x9.size(), a1.size(), x8.size(), x7.size());
        zgemm3m_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x9.size()*x8.size()*x7.size(),
               1.0, i0data_sorted, x9.size()*x8.size()*x7.size(), i1data_sorted, x9.size()*x8.size()*x7.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task426::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, a4, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(c3, a4, c1, a2), 0.0);
  {
    // tensor label: I108
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a2, c1, a4, c3);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, a2.size(), c1.size(), a4.size(), c3.size());
  }
  {
    // tensor label: I108
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a4, c3, a2, c1);
    sort_indices<1,0,3,2,1,1,1,1>(i0data, odata, a4.size(), c3.size(), a2.size(), c1.size());
  }
  out()->add_block(odata, c3, a4, c1, a2);
}

void Task427::Task_local::compute() {
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
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a4, c3, x1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a4, c3, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x1.size());
    // tensor label: I109
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a2, x1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a2, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a2.size(), x1.size());
    zgemm3m_("T", "N", c1.size()*a4.size()*c3.size(), a2.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c1.size()*a4.size()*c3.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task428::Task_local::compute() {
  const Index a2 = b(0);
  const Index x1 = b(1);
  // tensor label: I109
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a2, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma11
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x0, a2);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x0, a2)]);
    sort_indices<0,1,0,1,-1,1>(i1data, i1data_sorted, x0.size(), a2.size());
    zgemm3m_("T", "N", x1.size(), a2.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), a2.size());
  out()->add_block(odata, a2, x1);
}

void Task429::Task_local::compute() {
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
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a2, c3, x1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a2, c3, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x1.size());
    // tensor label: I112
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, x1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a4.size(), x1.size());
    zgemm3m_("T", "N", c1.size()*a2.size()*c3.size(), a4.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c3.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task430::Task_local::compute() {
  const Index a4 = b(0);
  const Index x1 = b(1);
  // tensor label: I112
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a4, x1)]);
  std::fill_n(odata.get(), out()->get_size(a4, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a4, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma11
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x0, a4);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x0, a4)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size());
    zgemm3m_("T", "N", x1.size(), a4.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), a4.size());
  out()->add_block(odata, a4, x1);
}

void Task431::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  // tensor label: h1
  std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a2);
  std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c3, a2)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size());
  // tensor label: I115
  std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a4);
  std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a4)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), a4.size());
  zgemm3m_("T", "N", c3.size()*a2.size(), c1.size()*a4.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c3.size()*a2.size());
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), a4.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task432::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  // tensor label: I115
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a4)]);
  std::fill_n(odata.get(), out()->get_size(c1, a4), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma27
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a4, x1, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a4, x1, x0)]);
      sort_indices<2,3,0,1,0,1,-1,1>(i1data, i1data_sorted, c1.size(), a4.size(), x1.size(), x0.size());
      zgemm3m_("T", "N", 1, c1.size()*a4.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size());
  out()->add_block(odata, c1, a4);
}

void Task433::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  // tensor label: h1
  std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a4);
  std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c3, a4)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size());
  // tensor label: I118
  std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a2);
  std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a2)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size());
  zgemm3m_("T", "N", c3.size()*a4.size(), c1.size()*a2.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c3.size()*a4.size());
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, c3.size(), a4.size(), c1.size(), a2.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task434::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  // tensor label: I118
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma27
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a2, x1, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a2, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size(), x1.size(), x0.size());
      zgemm3m_("T", "N", 1, c1.size()*a2.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size());
  out()->add_block(odata, c1, a2);
}

void Task435::Task_local::compute() {
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
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a4, c5, a2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a4, c5, a2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c5.size(), a2.size());
    // tensor label: I121
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, c5);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, c5)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), c5.size());
    zgemm3m_("T", "N", c1.size()*a4.size()*a2.size(), c3.size(), c5.size(),
           1.0, i0data_sorted, c5.size(), i1data_sorted, c5.size(),
           1.0, odata_sorted, c1.size()*a4.size()*a2.size());
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), a2.size(), c3.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task436::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  // tensor label: I121
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, c5)]);
  std::fill_n(odata.get(), out()->get_size(c3, c5), 0.0);
  {
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, c5);
    sort_indices<0,1,1,1,2,1>(i0data, odata, c3.size(), c5.size());
  }
  out()->add_block(odata, c3, c5);
}

void Task437::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  // tensor label: I121
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, c5)]);
  std::fill_n(odata.get(), out()->get_size(c3, c5), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, c5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c5), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma27
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: I868
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, c5, x1, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, c5, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), c5.size(), x1.size(), x0.size());
      zgemm3m_("T", "N", 1, c3.size()*c5.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c3.size(), c5.size());
  out()->add_block(odata, c3, c5);
}

void Task438::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I868
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, c5, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, c5, x1, x0), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, c5, x1, x0);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c3.size(), c5.size(), x1.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(x1, c5, c3, x0);
    sort_indices<2,1,0,3,1,1,-1,1>(i1data, odata, x1.size(), c5.size(), c3.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i2data = in(0)->get_block(x1, x0, c3, c5);
    sort_indices<2,3,0,1,1,1,1,1>(i2data, odata, x1.size(), x0.size(), c3.size(), c5.size());
  }
  out()->add_block(odata, c3, c5, x1, x0);
}

void Task439::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  // tensor label: I121
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, c5)]);
  std::fill_n(odata.get(), out()->get_size(c3, c5), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, c5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c5), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma11
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, x1)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, x1, x0, c5);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, x1, x0, c5)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), x1.size(), x0.size(), c5.size());
      zgemm3m_("T", "N", 1, c3.size()*c5.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c3.size(), c5.size());
  out()->add_block(odata, c3, c5);
}

void Task440::Task_local::compute() {
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
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a2, c5, a4);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a2, c5, a4)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c5.size(), a4.size());
    // tensor label: I123
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, c5);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, c5)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), c5.size());
    zgemm3m_("T", "N", c1.size()*a2.size()*a4.size(), c3.size(), c5.size(),
           1.0, i0data_sorted, c5.size(), i1data_sorted, c5.size(),
           1.0, odata_sorted, c1.size()*a2.size()*a4.size());
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), a4.size(), c3.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task441::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  // tensor label: I123
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, c5)]);
  std::fill_n(odata.get(), out()->get_size(c3, c5), 0.0);
  {
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, c5);
    sort_indices<0,1,1,1,-2,1>(i0data, odata, c3.size(), c5.size());
  }
  out()->add_block(odata, c3, c5);
}

void Task442::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  // tensor label: I123
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, c5)]);
  std::fill_n(odata.get(), out()->get_size(c3, c5), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, c5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c5), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma27
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: I871
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, c5, x1, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, c5, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), c5.size(), x1.size(), x0.size());
      zgemm3m_("T", "N", 1, c3.size()*c5.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c3.size(), c5.size());
  out()->add_block(odata, c3, c5);
}

void Task443::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I871
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, c5, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, c5, x1, x0), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, c5, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c3.size(), c5.size(), x1.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(x1, c5, c3, x0);
    sort_indices<2,1,0,3,1,1,1,1>(i1data, odata, x1.size(), c5.size(), c3.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i2data = in(0)->get_block(x1, x0, c3, c5);
    sort_indices<2,3,0,1,1,1,-1,1>(i2data, odata, x1.size(), x0.size(), c3.size(), c5.size());
  }
  out()->add_block(odata, c3, c5, x1, x0);
}

void Task444::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  // tensor label: I123
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, c5)]);
  std::fill_n(odata.get(), out()->get_size(c3, c5), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, c5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c5), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma11
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, x1)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, x1, x0, c5);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, x1, x0, c5)]);
      sort_indices<2,1,0,3,0,1,-1,1>(i1data, i1data_sorted, c3.size(), x1.size(), x0.size(), c5.size());
      zgemm3m_("T", "N", 1, c3.size()*c5.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c3.size(), c5.size());
  out()->add_block(odata, c3, c5);
}

void Task445::Task_local::compute() {
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
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a5, c3, a2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a5, c3, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a5.size(), c3.size(), a2.size());
    // tensor label: I125
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a5, a4);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a5, a4)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a5.size(), a4.size());
    zgemm3m_("T", "N", c1.size()*c3.size()*a2.size(), a4.size(), a5.size(),
           1.0, i0data_sorted, a5.size(), i1data_sorted, a5.size(),
           1.0, odata_sorted, c1.size()*c3.size()*a2.size());
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), a2.size(), a4.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task446::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  // tensor label: I125
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a5, a4)]);
  std::fill_n(odata.get(), out()->get_size(a5, a4), 0.0);
  {
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a5, a4);
    sort_indices<0,1,1,1,-2,1>(i0data, odata, a5.size(), a4.size());
  }
  out()->add_block(odata, a5, a4);
}

void Task447::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  // tensor label: I125
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a5, a4)]);
  std::fill_n(odata.get(), out()->get_size(a5, a4), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a5, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a5, a4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma27
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: I874
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a5, a4, x1, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a5, a4, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a5.size(), a4.size(), x1.size(), x0.size());
      zgemm3m_("T", "N", 1, a5.size()*a4.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a5.size(), a4.size());
  out()->add_block(odata, a5, a4);
}

void Task448::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I874
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a5, a4, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a5, a4, x1, x0), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a5, a4, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, a5.size(), a4.size(), x1.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(x1, a4, a5, x0);
    sort_indices<2,1,0,3,1,1,1,1>(i1data, odata, x1.size(), a4.size(), a5.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i2data = in(0)->get_block(x1, x0, a5, a4);
    sort_indices<2,3,0,1,1,1,-1,1>(i2data, odata, x1.size(), x0.size(), a5.size(), a4.size());
  }
  out()->add_block(odata, a5, a4, x1, x0);
}

void Task449::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  // tensor label: I125
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a5, a4)]);
  std::fill_n(odata.get(), out()->get_size(a5, a4), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a5, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a5, a4), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma11
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, x1)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a5, x1, x0, a4);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a5, x1, x0, a4)]);
      sort_indices<2,1,0,3,0,1,-1,1>(i1data, i1data_sorted, a5.size(), x1.size(), x0.size(), a4.size());
      zgemm3m_("T", "N", 1, a5.size()*a4.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a5.size(), a4.size());
  out()->add_block(odata, a5, a4);
}

#endif
