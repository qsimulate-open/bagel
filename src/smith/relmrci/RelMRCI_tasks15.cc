//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_tasks15.cc
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

#include <src/smith/relmrci/RelMRCI_tasks15.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelMRCI;

void Task700::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index a2 = b(2);
  const Index x5 = b(3);
  const Index a1 = b(4);
  const Index x4 = b(5);
  // tensor label: I1253
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x3, x2, a2, x5, a1, x4)]);
  std::fill_n(odata.get(), out()->get_size(x3, x2, a2, x5, a1, x4), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x3, x2, a2, x5, a1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, a2, x5, a1, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, a1, x4, a3);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, a1, x4, a3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a3.size());
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a3, x3, x2, a2);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a3, x3, x2, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), x3.size(), x2.size(), a2.size());
    zgemm3m_("T", "N", x5.size()*a1.size()*x4.size(), x3.size()*x2.size()*a2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x5.size()*a1.size()*x4.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), x3.size(), x2.size(), a2.size());
  out()->add_block(odata, x3, x2, a2, x5, a1, x4);
}

void Task701::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);
  // tensor label: I173
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        for (auto& x3 : *range_[1]) {
          // tensor label: Gamma409
          std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x4, x2, x3, x1);
          std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x0, x4, x2, x3, x1)]);
          sort_indices<0,2,3,4,1,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x2.size(), x3.size(), x1.size());
          // tensor label: I1256
          std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, a2, x2, x5, a1, x4);
          std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, a2, x2, x5, a1, x4)]);
          sort_indices<3,5,2,0,1,4,0,1,1,1>(i1data, i1data_sorted, x3.size(), a2.size(), x2.size(), x5.size(), a1.size(), x4.size());
          zgemm3m_("T", "N", x0.size()*x1.size(), a2.size()*a1.size(), x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, x0.size()*x1.size());
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a2.size(), a1.size());
  out()->add_block(odata, a1, x0, x1, a2);
}

void Task702::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  const Index x2 = b(2);
  const Index x5 = b(3);
  const Index a1 = b(4);
  const Index x4 = b(5);
  // tensor label: I1256
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x3, a2, x2, x5, a1, x4)]);
  std::fill_n(odata.get(), out()->get_size(x3, a2, x2, x5, a1, x4), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x3, a2, x2, x5, a1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a2, x2, x5, a1, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, a1, x4, a3);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, a1, x4, a3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a3.size());
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, a2, a3, x2);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, a2, a3, x2)]);
    sort_indices<2,0,1,3,0,1,-1,1>(i1data, i1data_sorted, x3.size(), a2.size(), a3.size(), x2.size());
    zgemm3m_("T", "N", x5.size()*a1.size()*x4.size(), x3.size()*a2.size()*x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x5.size()*a1.size()*x4.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), x3.size(), a2.size(), x2.size());
  out()->add_block(odata, x3, a2, x2, x5, a1, x4);
}

void Task703::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index x0 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, x1, c1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, c1, x0), 0.0);
  {
    // tensor label: I194
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, c2, x1, x0);
    sort_indices<1,2,0,3,1,1,1,1>(i0data, odata, c1.size(), c2.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, c2, x1, c1, x0);
}

void Task704::Task_local::compute() {
  const Index c1 = b(0);
  const Index c2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I194
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, c2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, c2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c2, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x3, x0, x2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x3, x0, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), x3.size(), x0.size(), x2.size());
      // tensor label: I195
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, c2, x3, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, c2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), c2.size(), x3.size(), x2.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), c1.size()*c2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), c2.size());
  out()->add_block(odata, c1, c2, x1, x0);
}

void Task705::Task_local::compute() {
  const Index c1 = b(0);
  const Index c2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I195
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, c2, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c2, x3, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c2, x3, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, x3, c4, x2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c3, x3, c4, x2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), x3.size(), c4.size(), x2.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, c4, c2, c3);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, c4, c2, c3)]);
      sort_indices<3,1,0,2,0,1,-2,1>(i1data, i1data_sorted, c1.size(), c4.size(), c2.size(), c3.size());
      zgemm3m_("T", "N", x3.size()*x2.size(), c1.size()*c2.size(), c4.size()*c3.size(),
             1.0, i0data_sorted, c4.size()*c3.size(), i1data_sorted, c4.size()*c3.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c1.size(), c2.size());
  out()->add_block(odata, c1, c2, x3, x2);
}

void Task706::Task_local::compute() {
  const Index c1 = b(0);
  const Index c2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I195
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, c2, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c2, x3, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c2, x3, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a3, c2, a4);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a3, c2, a4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), c2.size(), a4.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, x3, a3, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, x3, a3, x2)]);
      sort_indices<2,0,1,3,0,1,-2,1>(i1data, i1data_sorted, a4.size(), x3.size(), a3.size(), x2.size());
      zgemm3m_("T", "N", c1.size()*c2.size(), x3.size()*x2.size(), a4.size()*a3.size(),
             1.0, i0data_sorted, a4.size()*a3.size(), i1data_sorted, a4.size()*a3.size(),
             1.0, odata_sorted, c1.size()*c2.size());
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), c2.size(), x3.size(), x2.size());
  out()->add_block(odata, c1, c2, x3, x2);
}

void Task707::Task_local::compute() {
  const Index c1 = b(0);
  const Index c2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I194
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, c2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, c2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: Gamma412
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x5, x1, x4);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, x5, x1, x4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), x5.size(), x1.size(), x4.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, x5, c2, x4);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, x5, c2, x4)]);
      sort_indices<1,3,0,2,0,1,-2,1>(i1data, i1data_sorted, c1.size(), x5.size(), c2.size(), x4.size());
      zgemm3m_("T", "N", x0.size()*x1.size(), c1.size()*c2.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c1.size(), c2.size());
  out()->add_block(odata, c1, c2, x1, x0);
}

void Task708::Task_local::compute() {
  const Index c1 = b(0);
  const Index c2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I194
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, c2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, c2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c2, x1, x0), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      // tensor label: Gamma413
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x7, x1, x6);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, x7, x1, x6)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), x7.size(), x1.size(), x6.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, x7, c2, x6);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, x7, c2, x6)]);
      sort_indices<1,3,0,2,0,1,-1,1>(i1data, i1data_sorted, c1.size(), x7.size(), c2.size(), x6.size());
      zgemm3m_("T", "N", x0.size()*x1.size(), c1.size()*c2.size(), x7.size()*x6.size(),
             1.0, i0data_sorted, x7.size()*x6.size(), i1data_sorted, x7.size()*x6.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c1.size(), c2.size());
  out()->add_block(odata, c1, c2, x1, x0);
}

void Task709::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, a4, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(c3, a4, c1, a2), 0.0);
  {
    // tensor label: I773
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, c3, a2, a4);
    sort_indices<1,3,0,2,1,1,1,1>(i0data, odata, c1.size(), c3.size(), a2.size(), a4.size());
  }
  out()->add_block(odata, c3, a4, c1, a2);
}

void Task710::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index a4 = b(3);
  // tensor label: I773
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, c3, a2, a4)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, a2, a4), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c3, a2, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, a2, a4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, a2, x0, a4);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, a2, x0, a4)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), x0.size(), a4.size());
      // tensor label: I774
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, c3, x1, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, c3, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), c3.size(), x1.size(), x0.size());
      zgemm3m_("T", "N", a2.size()*a4.size(), c1.size()*c3.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, a2.size()*a4.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a2.size(), a4.size(), c1.size(), c3.size());
  out()->add_block(odata, c1, c3, a2, a4);
}

void Task711::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I774
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, c3, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x3, x0, x2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x3, x0, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), x3.size(), x0.size(), x2.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, x3, c3, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, x3, c3, x2)]);
      sort_indices<1,3,0,2,0,1,-2,1>(i1data, i1data_sorted, c1.size(), x3.size(), c3.size(), x2.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), c1.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), c3.size());
  out()->add_block(odata, c1, c3, x1, x0);
}

void Task712::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index a4 = b(3);
  // tensor label: I773
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, c3, a2, a4)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, a2, a4), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c3, a2, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, a2, a4), 0.0);
  for (auto& c5 : *range_[0]) {
    for (auto& c6 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c5, a4, c6, a2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c5, a4, c6, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c5.size(), a4.size(), c6.size(), a2.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, c6, c3, c5);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, c6, c3, c5)]);
      sort_indices<3,1,0,2,0,1,2,1>(i1data, i1data_sorted, c1.size(), c6.size(), c3.size(), c5.size());
      zgemm3m_("T", "N", a4.size()*a2.size(), c1.size()*c3.size(), c6.size()*c5.size(),
             1.0, i0data_sorted, c6.size()*c5.size(), i1data_sorted, c6.size()*c5.size(),
             1.0, odata_sorted, a4.size()*a2.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, a4.size(), a2.size(), c1.size(), c3.size());
  out()->add_block(odata, c1, c3, a2, a4);
}

void Task713::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index a4 = b(3);
  // tensor label: I773
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, c3, a2, a4)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, a2, a4), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c3, a2, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, a2, a4), 0.0);
  for (auto& c5 : *range_[0]) {
    for (auto& c6 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c5, a2, c6, a4);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c5, a2, c6, a4)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c5.size(), a2.size(), c6.size(), a4.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, c6, c3, c5);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, c6, c3, c5)]);
      sort_indices<3,1,0,2,0,1,-2,1>(i1data, i1data_sorted, c1.size(), c6.size(), c3.size(), c5.size());
      zgemm3m_("T", "N", a2.size()*a4.size(), c1.size()*c3.size(), c6.size()*c5.size(),
             1.0, i0data_sorted, c6.size()*c5.size(), i1data_sorted, c6.size()*c5.size(),
             1.0, odata_sorted, a2.size()*a4.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a2.size(), a4.size(), c1.size(), c3.size());
  out()->add_block(odata, c1, c3, a2, a4);
}

void Task714::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index a4 = b(3);
  // tensor label: I773
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, c3, a2, a4)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, a2, a4), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c3, a2, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, a2, a4), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, a2, x2, a4);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, a2, x2, a4)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), x2.size(), a4.size());
      // tensor label: I977
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, c3, x3, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, c3, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), c3.size(), x3.size(), x2.size());
      zgemm3m_("T", "N", a2.size()*a4.size(), c1.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a2.size()*a4.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a2.size(), a4.size(), c1.size(), c3.size());
  out()->add_block(odata, c1, c3, a2, a4);
}

void Task715::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I977
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, c3, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x3, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c3, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x3, x2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma368
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x1, x2, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x1, x2, x0)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x2.size(), x0.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, x1, c3, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, x1, c3, x0)]);
      sort_indices<1,3,0,2,0,1,-2,1>(i1data, i1data_sorted, c1.size(), x1.size(), c3.size(), x0.size());
      zgemm3m_("T", "N", x3.size()*x2.size(), c1.size()*c3.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c1.size(), c3.size());
  out()->add_block(odata, c1, c3, x3, x2);
}

void Task716::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index a4 = b(3);
  // tensor label: I773
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, c3, a2, a4)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, a2, a4), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c3, a2, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, a2, a4), 0.0);
  // scalar
  // tensor label: Gamma424
  std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block();
  std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size()]);
  sort_indices<0,1,1,1>(i0data, i0data_sorted);
  // tensor label: I1289
  std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a4, c3, a2);
  std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a4, c3, a2)]);
  sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
  zgemm3m_("T", "N", 1, c1.size()*a4.size()*c3.size()*a2.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, 1);
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->add_block(odata, c1, c3, a2, a4);
}

void Task717::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I1289
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a4, c3, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a4, c3, a2), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a4, c3, a2);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<0,3,2,1,1,1,-2,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->add_block(odata, c1, a4, c3, a2);
}

void Task718::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index a4 = b(3);
  // tensor label: I773
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, c3, a2, a4)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, a2, a4), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c3, a2, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, a2, a4), 0.0);
  // scalar
  // tensor label: Gamma426
  std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block();
  std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size()]);
  sort_indices<0,1,1,1>(i0data, i0data_sorted);
  // tensor label: I1293
  std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a4, c3, a2);
  std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a4, c3, a2)]);
  sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
  zgemm3m_("T", "N", 1, c1.size()*a4.size()*c3.size()*a2.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, 1);
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->add_block(odata, c1, c3, a2, a4);
}

void Task719::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I1293
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a4, c3, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a4, c3, a2), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a4, c3, a2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<0,3,2,1,1,1,-1,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->add_block(odata, c1, a4, c3, a2);
}

void Task720::Task_local::compute() {
  const Index x1 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x1, a2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(x1, a2, x0, a1), 0.0);
  {
    // tensor label: I1228
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0, a1, a2);
    sort_indices<0,3,1,2,1,1,1,1>(i0data, odata, x1.size(), x0.size(), a1.size(), a2.size());
  }
  out()->add_block(odata, x1, a2, x0, a1);
}

void Task721::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);
  // tensor label: I1228
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x1, x0, a1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, a1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x1, x0, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, a1, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a1, c4, a2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c3, a1, c4, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a1.size(), c4.size(), a2.size());
      // tensor label: I1229
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c4, c3, x1, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c4, c3, x1, x0)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c4.size(), c3.size(), x1.size(), x0.size());
      zgemm3m_("T", "N", a1.size()*a2.size(), x1.size()*x0.size(), c4.size()*c3.size(),
             1.0, i0data_sorted, c4.size()*c3.size(), i1data_sorted, c4.size()*c3.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x1.size(), x0.size());
  out()->add_block(odata, x1, x0, a1, a2);
}

void Task722::Task_local::compute() {
  const Index c4 = b(0);
  const Index c3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1229
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c4, c3, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c4, c3, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c4, c3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c4, c3, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma368
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x1, x2, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x1, x2, x0)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x2.size(), x0.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, c4, x2, c3);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, c4, x2, c3)]);
      sort_indices<0,2,1,3,0,1,-2,1>(i1data, i1data_sorted, x3.size(), c4.size(), x2.size(), c3.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), c4.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c4.size(), c3.size());
  out()->add_block(odata, c4, c3, x1, x0);
}

void Task723::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);
  // tensor label: I1228
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x1, x0, a1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, a1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x1, x0, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, a1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma368
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x1, x2, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x1, x2, x0)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x2.size(), x0.size());
      // tensor label: I1262
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a1, a2, x3, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a1, a2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a1.size(), a2.size(), x3.size(), x2.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), a1.size()*a2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a1.size(), a2.size());
  out()->add_block(odata, x1, x0, a1, a2);
}

void Task724::Task_local::compute() {
  const Index a1 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1262
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, a2, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(a1, a2, x3, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, a2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, a2, x3, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, a3, x2, a4);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, a3, x2, a4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), x2.size(), a4.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, a1, a3, a2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, a1, a3, a2)]);
      sort_indices<2,0,1,3,0,1,-2,1>(i1data, i1data_sorted, a4.size(), a1.size(), a3.size(), a2.size());
      zgemm3m_("T", "N", x3.size()*x2.size(), a1.size()*a2.size(), a4.size()*a3.size(),
             1.0, i0data_sorted, a4.size()*a3.size(), i1data_sorted, a4.size()*a3.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), a1.size(), a2.size());
  out()->add_block(odata, a1, a2, x3, x2);
}

void Task725::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);
  // tensor label: I1228
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x1, x0, a1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, a1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x1, x0, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, a1, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: Gamma432
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x4, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x0, x4, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x1.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x5, a1, x4, a2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x5, a1, x4, a2)]);
      sort_indices<0,2,1,3,0,1,-2,1>(i1data, i1data_sorted, x5.size(), a1.size(), x4.size(), a2.size());
      zgemm3m_("T", "N", x0.size()*x1.size(), a1.size()*a2.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), a2.size());
  out()->add_block(odata, x1, x0, a1, a2);
}

void Task726::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);
  // tensor label: I1228
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x1, x0, a1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, a1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x1, x0, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, a1, a2), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      // tensor label: Gamma433
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x6, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x7, x0, x6, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x6.size(), x1.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x7, a1, x6, a2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x7, a1, x6, a2)]);
      sort_indices<0,2,1,3,0,1,-1,1>(i1data, i1data_sorted, x7.size(), a1.size(), x6.size(), a2.size());
      zgemm3m_("T", "N", x0.size()*x1.size(), a1.size()*a2.size(), x7.size()*x6.size(),
             1.0, i0data_sorted, x7.size()*x6.size(), i1data_sorted, x7.size()*x6.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), a2.size());
  out()->add_block(odata, x1, x0, a1, a2);
}

void Task728::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x2, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x0, x1), 0.0);
  {
    // tensor label: I1308
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, x2, x1, x0);
    sort_indices<0,1,3,2,1,1,1,1>(i0data, odata, c1.size(), x2.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, c1, x2, x0, x1);
}

void Task729::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1308
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: Gamma5
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

void Task730::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1308
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, x5, x4, x3);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, x5, x4, x3)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), x5.size(), x4.size(), x3.size());
        // tensor label: Gamma81
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x2, x5, x4, x3, x1, x0);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x2, x5, x4, x3, x1, x0)]);
        sort_indices<1,2,3,0,4,5,0,1,1,2>(i1data, i1data_sorted, x2.size(), x5.size(), x4.size(), x3.size(), x1.size(), x0.size());
        zgemm3m_("T", "N", c1.size(), x2.size()*x1.size()*x0.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), x2.size(), x1.size(), x0.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task731::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1308
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x4, c1, x3);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x4, c1, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), c1.size(), x3.size());
        // tensor label: Gamma4
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x5, x4, x2, x3, x1, x0);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x5, x4, x2, x3, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,2>(i1data, i1data_sorted, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
        zgemm3m_("T", "N", c1.size(), x2.size()*x1.size()*x0.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), x2.size(), x1.size(), x0.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task732::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x0, x1, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1, c1, a2), 0.0);
  {
    // tensor label: I1310
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0, c1, a2);
    sort_indices<1,0,2,3,1,1,1,1>(i0data, odata, x1.size(), x0.size(), c1.size(), a2.size());
  }
  out()->add_block(odata, x0, x1, c1, a2);
}

void Task733::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I1310
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  // tensor label: h1
  std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a2);
  std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a2)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size());
  // tensor label: Gamma27
  std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x1, x0);
  std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
  zgemm3m_("T", "N", c1.size()*a2.size(), x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a2.size());
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), x0.size());
  out()->add_block(odata, x1, x0, c1, a2);
}

void Task734::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I1310
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma24
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a2, x3, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,2>(i1data, i1data_sorted, c1.size(), a2.size(), x3.size(), x2.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), c1.size()*a2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), a2.size());
  out()->add_block(odata, x1, x0, c1, a2);
}

void Task735::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I1310
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, x3, x2, a2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, x3, x2, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x3.size(), x2.size(), a2.size());
      // tensor label: Gamma5
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x2, x3, x1, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x2, x3, x1, x0)]);
      sort_indices<1,0,2,3,0,1,1,2>(i1data, i1data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
      zgemm3m_("T", "N", c1.size()*a2.size(), x1.size()*x0.size(), x2.size()*x3.size(),
             1.0, i0data_sorted, x2.size()*x3.size(), i1data_sorted, x2.size()*x3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), x0.size());
  out()->add_block(odata, x1, x0, c1, a2);
}

void Task736::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I1310
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, a2, c1, x2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, a2, c1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c1.size(), x2.size());
      // tensor label: Gamma24
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, x2, x1, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,-1,2>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      zgemm3m_("T", "N", a2.size()*c1.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x1.size(), x0.size());
  out()->add_block(odata, x1, x0, c1, a2);
}

void Task737::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I1310
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, c1, a2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, c1, a2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), c1.size(), a2.size());
      // tensor label: Gamma24
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, x2, x1, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,2>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      zgemm3m_("T", "N", c1.size()*a2.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), x0.size());
  out()->add_block(odata, x1, x0, c1, a2);
}

void Task738::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x1, x2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(x1, x2, x0, a1), 0.0);
  {
    // tensor label: I1312
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x2, x1, a1);
    sort_indices<2,1,0,3,1,1,1,1>(i0data, odata, x0.size(), x2.size(), x1.size(), a1.size());
  }
  out()->add_block(odata, x1, x2, x0, a1);
}

void Task739::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  // tensor label: I1312
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x0, x2, x1, a1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x2, x1, a1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, a1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size());
    // tensor label: Gamma33
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, x0, x2, x1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, x0, x2, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
    zgemm3m_("T", "N", a1.size(), x0.size()*x2.size()*x1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x2.size(), x1.size());
  out()->add_block(odata, x0, x2, x1, a1);
}

void Task740::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  // tensor label: I1312
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x0, x2, x1, a1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x2, x1, a1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, a1, x4, x3);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, a1, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), x3.size());
        // tensor label: Gamma32
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,2>(i1data, i1data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        zgemm3m_("T", "N", a1.size(), x0.size()*x2.size()*x1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x2.size(), x1.size());
  out()->add_block(odata, x0, x2, x1, a1);
}

void Task741::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  // tensor label: I1312
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x0, x2, x1, a1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x2, x1, a1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma31
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
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, x0, x2, x1, a1);
}

void Task742::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index x0 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, x1, c1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, c1, x0), 0.0);
  {
    // tensor label: I1314
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, c2, x0, x1);
    sort_indices<1,3,0,2,1,1,1,1>(i0data, odata, c1.size(), c2.size(), x0.size(), x1.size());
  }
  out()->add_block(odata, c2, x1, c1, x0);
}

void Task743::Task_local::compute() {
  const Index c1 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1314
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, c2, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, c2, x0, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c2, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma0
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

void Task744::Task_local::compute() {
  const Index c3 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, x0, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(c3, x0, c1, a2), 0.0);
  {
    // tensor label: I1320
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, c3, c1, a2);
    sort_indices<1,0,2,3,1,1,1,1>(i0data, odata, x0.size(), c3.size(), c1.size(), a2.size());
  }
  out()->add_block(odata, c3, x0, c1, a2);
}

void Task745::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I1320
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x0, c3, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x0, c3, c1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x0, c3, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c3, c1, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, x1, c1, a2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c3, x1, c1, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), x1.size(), c1.size(), a2.size());
    // tensor label: Gamma11
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x0, x1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size());
    zgemm3m_("T", "N", c3.size()*c1.size()*a2.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c3.size()*c1.size()*a2.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c3.size(), c1.size(), a2.size(), x0.size());
  out()->add_block(odata, x0, c3, c1, a2);
}

void Task746::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I1320
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x0, c3, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x0, c3, c1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x0, c3, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c3, c1, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, x1, c3, a2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, x1, c3, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x1.size(), c3.size(), a2.size());
    // tensor label: Gamma11
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x0, x1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,-1,1>(i1data, i1data_sorted, x0.size(), x1.size());
    zgemm3m_("T", "N", c1.size()*c3.size()*a2.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c1.size()*c3.size()*a2.size());
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), a2.size(), x0.size());
  out()->add_block(odata, x0, c3, c1, a2);
}

void Task747::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, a4, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(c3, a4, c1, a2), 0.0);
  {
    // tensor label: I1336
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a4, c3, a2);
    sort_indices<2,1,0,3,1,1,1,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
  }
  out()->add_block(odata, c3, a4, c1, a2);
}

void Task748::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I1336
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

void Task749::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  {
    // tensor label: I1338
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, a3, c2, a1);
    sort_indices<2,1,0,3,1,1,1,1>(i0data, odata, x0.size(), a3.size(), c2.size(), a1.size());
  }
  out()->add_block(odata, c2, a3, x0, a1);
}

#endif
