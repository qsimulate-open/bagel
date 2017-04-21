//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelCASA_tasks3.cc
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

#include <src/smith/relcasa/RelCASA_tasks3.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelCASA;

void Task100::Task_local::compute() {
  const Index a1 = b(0);
  const Index c2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I79
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, c2, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(a1, c2, x3, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, c2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, c2, x3, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x2, c3);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x2, c3)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x2.size(), c3.size());
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, a1, c2, x3);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, a1, c2, x3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), a1.size(), c2.size(), x3.size());
    zgemm3m_("T", "N", x2.size(), a1.size()*c2.size()*x3.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, x2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x2.size(), a1.size(), c2.size(), x3.size());
  out()->add_block(odata, a1, c2, x3, x2);
}

void Task101::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index c2 = b(3);
  // tensor label: I72
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x1, c2)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x1, c2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, c2), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: Gamma5
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x2, x3, x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x2, x3, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
      // tensor label: I82
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, a1, x3, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, a1, x3, x2)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size(), x3.size(), x2.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), c2.size()*a1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a1.size());
  out()->add_block(odata, a1, x0, x1, c2);
}

void Task102::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I82
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, a1, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c2, a1, x3, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, a1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x3, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x2, c3);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x2, c3)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x2.size(), c3.size());
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, a1, c3, x3);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, a1, c3, x3)]);
    sort_indices<2,0,1,3,0,1,-1,1>(i1data, i1data_sorted, c2.size(), a1.size(), c3.size(), x3.size());
    zgemm3m_("T", "N", x2.size(), c2.size()*a1.size()*x3.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, x2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x2.size(), c2.size(), a1.size(), x3.size());
  out()->add_block(odata, c2, a1, x3, x2);
}

void Task103::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index c2 = b(3);
  // tensor label: I72
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x1, c2)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x1, c2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, c2), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c2, x2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c2, x2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), x2.size());
    // tensor label: I85
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a1, x0, x1, x2);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a1, x0, x1, x2)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a1.size(), x0.size(), x1.size(), x2.size());
    zgemm3m_("T", "N", c2.size(), a1.size()*x0.size()*x1.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, c2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x0.size(), x1.size());
  out()->add_block(odata, a1, x0, x1, c2);
}

void Task104::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: I85
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x1, x2)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x1, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x1, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, x2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma28
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x4, x3, x1, x2);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x0, x4, x3, x1, x2)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x1.size(), x2.size());
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x5, a1, x4, x3);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x5, a1, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x5.size(), a1.size(), x4.size(), x3.size());
        zgemm3m_("T", "N", x0.size()*x1.size()*x2.size(), a1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x0.size()*x1.size()*x2.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), a1.size());
  out()->add_block(odata, a1, x0, x1, x2);
}

void Task105::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index c2 = b(3);
  // tensor label: I72
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x1, c2)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x1, c2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, c2), 0.0);
  // tensor label: Gamma29
  std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
  std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
  // tensor label: I88
  std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, a1);
  std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, a1)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size());
  zgemm3m_("T", "N", x1.size()*x0.size(), c2.size()*a1.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, x1.size()*x0.size());
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a1.size());
  out()->add_block(odata, a1, x0, x1, c2);
}

void Task106::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  // tensor label: I88
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, a1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: f1
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a4, c3);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(a4, c3)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), c3.size());
      // tensor label: I89
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, a4, c3, a1);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, a4, c3, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), a4.size(), c3.size(), a1.size());
      zgemm3m_("T", "N", 1, c2.size()*a1.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size());
  out()->add_block(odata, c2, a1);
}

void Task107::Task_local::compute() {
  const Index c2 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a1 = b(3);
  // tensor label: I89
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, a4, c3, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, a4, c3, a1), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c2, a4, c3, a1);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, c2.size(), a4.size(), c3.size(), a1.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(c2, a1, c3, a4);
    sort_indices<0,3,2,1,1,1,-2,1>(i1data, odata, c2.size(), a1.size(), c3.size(), a4.size());
  }
  out()->add_block(odata, c2, a4, c3, a1);
}

void Task108::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index c2 = b(3);
  // tensor label: I72
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x1, c2)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x1, c2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, c2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma25
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I94
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, c2, a1, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, c2, a1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), c2.size(), a1.size(), x2.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), c2.size()*a1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a1.size());
  out()->add_block(odata, a1, x0, x1, c2);
}

void Task109::Task_local::compute() {
  const Index x3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x2 = b(3);
  // tensor label: I94
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x3, c2, a1, x2)]);
  std::fill_n(odata.get(), out()->get_size(x3, c2, a1, x2), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c2, a1, x3, x2);
    zscal_(c2.size()*a1.size()*x3.size()*x2.size(), e0_, i0data.get(), 1);
    sort_indices<2,0,1,3,1,1,1,1>(i0data, odata, c2.size(), a1.size(), x3.size(), x2.size());
  }
  out()->add_block(odata, x3, c2, a1, x2);
}

void Task110::Task_local::compute() {
  const Index x3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x2 = b(3);
  // tensor label: I94
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x3, c2, a1, x2)]);
  std::fill_n(odata.get(), out()->get_size(x3, c2, a1, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x3, c2, a1, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, c2, a1, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a3, x2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(a3, x2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), x2.size());
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, a3, c2, a1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, a3, c2, a1)]);
    sort_indices<1,0,2,3,0,1,-1,1>(i1data, i1data_sorted, x3.size(), a3.size(), c2.size(), a1.size());
    zgemm3m_("T", "N", x2.size(), x3.size()*c2.size()*a1.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x2.size(), x3.size(), c2.size(), a1.size());
  out()->add_block(odata, x3, c2, a1, x2);
}

void Task111::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index c2 = b(3);
  // tensor label: I72
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x1, c2)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x1, c2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, c2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma24
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0, x1, x2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x0, x1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
      // tensor label: I97
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, a1, c2, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, a1, c2, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), c2.size(), x2.size());
      zgemm3m_("T", "N", x0.size()*x1.size(), a1.size()*c2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), c2.size());
  out()->add_block(odata, a1, x0, x1, c2);
}

void Task112::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index c2 = b(2);
  const Index x2 = b(3);
  // tensor label: I97
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x3, a1, c2, x2)]);
  std::fill_n(odata.get(), out()->get_size(x3, a1, c2, x2), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, a1, c2, x2);
    zscal_(x3.size()*a1.size()*c2.size()*x2.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), a1.size(), c2.size(), x2.size());
  }
  out()->add_block(odata, x3, a1, c2, x2);
}

void Task113::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index c2 = b(2);
  const Index x2 = b(3);
  // tensor label: I97
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x3, a1, c2, x2)]);
  std::fill_n(odata.get(), out()->get_size(x3, a1, c2, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x3, a1, c2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1, c2, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a3, x2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(a3, x2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), x2.size());
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, a1, c2, a3);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, a1, c2, a3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), c2.size(), a3.size());
    zgemm3m_("T", "N", x2.size(), x3.size()*a1.size()*c2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x2.size(), x3.size(), a1.size(), c2.size());
  out()->add_block(odata, x3, a1, c2, x2);
}

void Task114::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index c2 = b(3);
  // tensor label: I72
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x1, c2)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x1, c2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, c2), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x2, a1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x2, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), a1.size());
    // tensor label: I100
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, x1, x2, x0);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, x1, x2, x0)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), x1.size(), x2.size(), x0.size());
    zgemm3m_("T", "N", a1.size(), c2.size()*x1.size()*x0.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, a1.size());
  }
  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), x1.size(), x0.size());
  out()->add_block(odata, a1, x0, x1, c2);
}

void Task115::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index x0 = b(3);
  // tensor label: I100
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, x1, x2, x0)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x2, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, x1, x2, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x2, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma33
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x4, x1, x3, x2, x0);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x4, x1, x3, x2, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x1.size(), x3.size(), x2.size(), x0.size());
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x5, x4, c2, x3);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x5, x4, c2, x3)]);
        sort_indices<0,1,3,2,0,1,-1,1>(i1data, i1data_sorted, x5.size(), x4.size(), c2.size(), x3.size());
        zgemm3m_("T", "N", x1.size()*x2.size()*x0.size(), c2.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x1.size()*x2.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x2.size(), x0.size(), c2.size());
  out()->add_block(odata, c2, x1, x2, x0);
}

void Task116::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index c2 = b(3);
  // tensor label: I72
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x1, c2)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x1, c2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, c2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a3, a1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(a3, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), a1.size());
    // tensor label: I103
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a3, c2, x0, x1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a3, c2, x0, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), x0.size(), x1.size());
    zgemm3m_("T", "N", a1.size(), c2.size()*x0.size()*x1.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, a1.size());
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), x0.size(), x1.size());
  out()->add_block(odata, a1, x0, x1, c2);
}

void Task117::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I103
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a3, c2, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a3, c2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma24
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0, x1, x2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x0, x1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, a3, c2, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, a3, c2, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a3.size(), c2.size(), x2.size());
      zgemm3m_("T", "N", x0.size()*x1.size(), a3.size()*c2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a3.size(), c2.size());
  out()->add_block(odata, a3, c2, x0, x1);
}

void Task118::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I103
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a3, c2, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a3, c2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma25
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, a3, x3, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, a3, x3, x2)]);
      sort_indices<2,3,0,1,0,1,-1,1>(i1data, i1data_sorted, c2.size(), a3.size(), x3.size(), x2.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), c2.size()*a3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a3.size());
  out()->add_block(odata, a3, c2, x0, x1);
}

void Task119::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index c2 = b(3);
  // tensor label: I72
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x1, c2)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x1, c2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, c2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: Gamma112
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x1, x4);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x0, x1, x4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x1.size(), x4.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x5, a1, c2, x4);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x5, a1, c2, x4)]);
      sort_indices<0,3,1,2,0,1,-1,1>(i1data, i1data_sorted, x5.size(), a1.size(), c2.size(), x4.size());
      zgemm3m_("T", "N", x0.size()*x1.size(), a1.size()*c2.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), c2.size());
  out()->add_block(odata, a1, x0, x1, c2);
}

void Task120::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index c2 = b(3);
  // tensor label: I72
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x1, c2)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x1, c2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, c2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: Gamma113
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x4, x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x4, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x1.size(), x0.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, a1, x5, x4);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, a1, x5, x4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size(), x5.size(), x4.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), c2.size()*a1.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a1.size());
  out()->add_block(odata, a1, x0, x1, c2);
}

void Task121::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index c2 = b(3);
  // tensor label: I72
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x1, c2)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x1, c2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, c2), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      // tensor label: Gamma116
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x1, x6);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x7, x0, x1, x6)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x1.size(), x6.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x7, a1, c2, x6);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x7, a1, c2, x6)]);
      sort_indices<0,3,1,2,0,1,-1,2>(i1data, i1data_sorted, x7.size(), a1.size(), c2.size(), x6.size());
      zgemm3m_("T", "N", x0.size()*x1.size(), a1.size()*c2.size(), x7.size()*x6.size(),
             1.0, i0data_sorted, x7.size()*x6.size(), i1data_sorted, x7.size()*x6.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), c2.size());
  out()->add_block(odata, a1, x0, x1, c2);
}

void Task122::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index c2 = b(3);
  // tensor label: I72
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x1, c2)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x1, c2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, c2), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      // tensor label: Gamma117
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x6, x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x7, x6, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x1.size(), x0.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, a1, x7, x6);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, a1, x7, x6)]);
      sort_indices<2,3,0,1,0,1,1,2>(i1data, i1data_sorted, c2.size(), a1.size(), x7.size(), x6.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), c2.size()*a1.size(), x7.size()*x6.size(),
             1.0, i0data_sorted, x7.size()*x6.size(), i1data_sorted, x7.size()*x6.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a1.size());
  out()->add_block(odata, a1, x0, x1, c2);
}

void Task123::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x0, x1, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1, c1, a2), 0.0);
  {
    // tensor label: I108
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a2, x1, x0, c1);
    sort_indices<2,1,3,0,1,1,1,1>(i0data, odata, a2.size(), x1.size(), x0.size(), c1.size());
  }
  out()->add_block(odata, x0, x1, c1, a2);
}

void Task124::Task_local::compute() {
  const Index a2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, x1, x0, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, x1, x0, c1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, x1, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x1, x0, c1), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, c3);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, c3)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), c3.size());
    // tensor label: I109
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a2, c3, x1, x0);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a2, c3, x1, x0)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), c3.size(), x1.size(), x0.size());
    zgemm3m_("T", "N", c1.size(), a2.size()*x1.size()*x0.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, c1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), x0.size());
  out()->add_block(odata, a2, x1, x0, c1);
}

void Task125::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I109
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c3, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a2, c3, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma25
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I110
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, a2, c3, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, a2, c3, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a2.size(), c3.size(), x2.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), a2.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a2.size(), c3.size());
  out()->add_block(odata, a2, c3, x1, x0);
}

void Task126::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x2 = b(3);
  // tensor label: I110
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x3, a2, c3, x2)]);
  std::fill_n(odata.get(), out()->get_size(x3, a2, c3, x2), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, a2, c3, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), a2.size(), c3.size(), x2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(c3, a2, x3, x2);
    sort_indices<2,1,0,3,1,1,-1,1>(i1data, odata, c3.size(), a2.size(), x3.size(), x2.size());
  }
  out()->add_block(odata, x3, a2, c3, x2);
}

void Task127::Task_local::compute() {
  const Index a2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, x1, x0, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, x1, x0, c1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, x1, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x1, x0, c1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: Gamma5
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x2, x3, x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x2, x3, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
      // tensor label: I115
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a2, c1, x3, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a2, c1, x3, x2)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, a2.size(), c1.size(), x3.size(), x2.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), a2.size()*c1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a2.size(), c1.size());
  out()->add_block(odata, a2, x1, x0, c1);
}

void Task128::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I115
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c1, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, x3, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, x3, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x2, c3);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x2, c3)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x2.size(), c3.size());
    // tensor label: I116
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, a2, c1, x3);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, a2, c1, x3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), c1.size(), x3.size());
    zgemm3m_("T", "N", x2.size(), a2.size()*c1.size()*x3.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, x2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x2.size(), a2.size(), c1.size(), x3.size());
  out()->add_block(odata, a2, c1, x3, x2);
}

void Task129::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I116
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, a2, c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, c1, x3), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a2, c1, x3);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c3.size(), a2.size(), c1.size(), x3.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(c1, a2, c3, x3);
    sort_indices<2,1,0,3,1,1,1,1>(i1data, odata, c1.size(), a2.size(), c3.size(), x3.size());
  }
  out()->add_block(odata, c3, a2, c1, x3);
}

void Task130::Task_local::compute() {
  const Index a2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, x1, x0, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, x1, x0, c1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, x1, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x1, x0, c1), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, x2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, x2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), x2.size());
    // tensor label: I121
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a2, x2, x1, x0);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a2, x2, x1, x0)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), x2.size(), x1.size(), x0.size());
    zgemm3m_("T", "N", c1.size(), a2.size()*x1.size()*x0.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, c1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), x0.size());
  out()->add_block(odata, a2, x1, x0, c1);
}

void Task131::Task_local::compute() {
  const Index a2 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I121
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a2, x2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma40
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x2, x4, x3, x1, x0);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x2, x4, x3, x1, x0)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x5, a2, x4, x3);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x5, a2, x4, x3)]);
        sort_indices<0,2,3,1,0,1,-1,1>(i1data, i1data_sorted, x5.size(), a2.size(), x4.size(), x3.size());
        zgemm3m_("T", "N", x2.size()*x1.size()*x0.size(), a2.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), a2.size());
  out()->add_block(odata, a2, x2, x1, x0);
}

void Task132::Task_local::compute() {
  const Index a2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, x1, x0, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, x1, x0, c1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, x1, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x1, x0, c1), 0.0);
  // tensor label: Gamma29
  std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
  std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
  // tensor label: I124
  std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a2);
  std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a2)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size());
  zgemm3m_("T", "N", x1.size()*x0.size(), c1.size()*a2.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, x1.size()*x0.size());
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), a2.size());
  out()->add_block(odata, a2, x1, x0, c1);
}

void Task133::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  // tensor label: I124
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: f1
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a4, c3);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(a4, c3)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), c3.size());
      // tensor label: I125
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a4, c3, a2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a4, c3, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
      zgemm3m_("T", "N", 1, c1.size()*a2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size());
  out()->add_block(odata, c1, a2);
}

void Task134::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I125
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a4, c3, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a4, c3, a2), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a4, c3, a2);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<0,3,2,1,1,1,2,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->add_block(odata, c1, a4, c3, a2);
}

void Task135::Task_local::compute() {
  const Index a2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, x1, x0, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, x1, x0, c1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, x1, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x1, x0, c1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma25
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I130
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, c1, a2, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, c1, a2, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), c1.size(), a2.size(), x2.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), c1.size()*a2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), a2.size());
  out()->add_block(odata, a2, x1, x0, c1);
}

void Task136::Task_local::compute() {
  const Index x3 = b(0);
  const Index c1 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I130
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x3, c1, a2, x2)]);
  std::fill_n(odata.get(), out()->get_size(x3, c1, a2, x2), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, a2, c1, x2);
    zscal_(x3.size()*a2.size()*c1.size()*x2.size(), e0_, i0data.get(), 1);
    sort_indices<0,2,1,3,1,1,1,1>(i0data, odata, x3.size(), a2.size(), c1.size(), x2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(c1, a2, x3, x2);
    zscal_(c1.size()*a2.size()*x3.size()*x2.size(), e0_, i1data.get(), 1);
    sort_indices<2,0,1,3,1,1,-1,1>(i1data, odata, c1.size(), a2.size(), x3.size(), x2.size());
  }
  out()->add_block(odata, x3, c1, a2, x2);
}

void Task137::Task_local::compute() {
  const Index x3 = b(0);
  const Index c1 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I130
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x3, c1, a2, x2)]);
  std::fill_n(odata.get(), out()->get_size(x3, c1, a2, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x3, c1, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, c1, a2, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a3, x2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(a3, x2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), x2.size());
    // tensor label: I131
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, a3, c1, a2);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, a3, c1, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a3.size(), c1.size(), a2.size());
    zgemm3m_("T", "N", x2.size(), x3.size()*c1.size()*a2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x2.size(), x3.size(), c1.size(), a2.size());
  out()->add_block(odata, x3, c1, a2, x2);
}

void Task138::Task_local::compute() {
  const Index x3 = b(0);
  const Index a3 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I131
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x3, a3, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x3, a3, c1, a2), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, a3, c1, a2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), a3.size(), c1.size(), a2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(x3, a2, c1, a3);
    sort_indices<0,3,2,1,1,1,-1,1>(i1data, odata, x3.size(), a2.size(), c1.size(), a3.size());
  }
  out()->add_block(odata, x3, a3, c1, a2);
}

void Task139::Task_local::compute() {
  const Index a2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, x1, x0, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, x1, x0, c1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, x1, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x1, x0, c1), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x2, a2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x2, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), a2.size());
    // tensor label: I136
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, x2, x1, x0);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, x2, x1, x0)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), x2.size(), x1.size(), x0.size());
    zgemm3m_("T", "N", a2.size(), c1.size()*x1.size()*x0.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x1.size(), x0.size());
  out()->add_block(odata, a2, x1, x0, c1);
}

void Task140::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I136
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma3
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

void Task141::Task_local::compute() {
  const Index a2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, x1, x0, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, x1, x0, c1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, x1, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x1, x0, c1), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a3, a2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(a3, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), a2.size());
    // tensor label: I139
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a3, c1, x1, x0);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a3, c1, x1, x0)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c1.size(), x1.size(), x0.size());
    zgemm3m_("T", "N", a2.size(), c1.size()*x1.size()*x0.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x1.size(), x0.size());
  out()->add_block(odata, a2, x1, x0, c1);
}

void Task142::Task_local::compute() {
  const Index a3 = b(0);
  const Index c1 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I139
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a3, c1, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a3, c1, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a3, c1, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c1, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma25
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I140
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, a3, c1, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, a3, c1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a3.size(), c1.size(), x2.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), a3.size()*c1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a3.size(), c1.size());
  out()->add_block(odata, a3, c1, x1, x0);
}

void Task143::Task_local::compute() {
  const Index x3 = b(0);
  const Index a3 = b(1);
  const Index c1 = b(2);
  const Index x2 = b(3);
  // tensor label: I140
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x3, a3, c1, x2)]);
  std::fill_n(odata.get(), out()->get_size(x3, a3, c1, x2), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, a3, c1, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), a3.size(), c1.size(), x2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(c1, a3, x3, x2);
    sort_indices<2,1,0,3,1,1,1,1>(i1data, odata, c1.size(), a3.size(), x3.size(), x2.size());
  }
  out()->add_block(odata, x3, a3, c1, x2);
}

void Task144::Task_local::compute() {
  const Index a2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, x1, x0, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, x1, x0, c1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, x1, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x1, x0, c1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: Gamma113
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x4, x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x4, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x1.size(), x0.size());
      // tensor label: I319
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x5, a2, c1, x4);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x5, a2, c1, x4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x5.size(), a2.size(), c1.size(), x4.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), a2.size()*c1.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a2.size(), c1.size());
  out()->add_block(odata, a2, x1, x0, c1);
}

void Task145::Task_local::compute() {
  const Index x5 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x4 = b(3);
  // tensor label: I319
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x5, a2, c1, x4)]);
  std::fill_n(odata.get(), out()->get_size(x5, a2, c1, x4), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, a2, c1, x4);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x5.size(), a2.size(), c1.size(), x4.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(c1, a2, x5, x4);
    sort_indices<2,1,0,3,1,1,-1,1>(i1data, odata, c1.size(), a2.size(), x5.size(), x4.size());
  }
  out()->add_block(odata, x5, a2, c1, x4);
}

void Task146::Task_local::compute() {
  const Index a2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, x1, x0, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, x1, x0, c1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, x1, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x1, x0, c1), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      // tensor label: Gamma117
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x6, x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x7, x6, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x1.size(), x0.size());
      // tensor label: I327
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x7, a2, c1, x6);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x7, a2, c1, x6)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x7.size(), a2.size(), c1.size(), x6.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), a2.size()*c1.size(), x7.size()*x6.size(),
             1.0, i0data_sorted, x7.size()*x6.size(), i1data_sorted, x7.size()*x6.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a2.size(), c1.size());
  out()->add_block(odata, a2, x1, x0, c1);
}

void Task147::Task_local::compute() {
  const Index x7 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x6 = b(3);
  // tensor label: I327
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x7, a2, c1, x6)]);
  std::fill_n(odata.get(), out()->get_size(x7, a2, c1, x6), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, a2, c1, x6);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x7.size(), a2.size(), c1.size(), x6.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(c1, a2, x7, x6);
    sort_indices<2,1,0,3,1,1,-1,2>(i1data, odata, c1.size(), a2.size(), x7.size(), x6.size());
  }
  out()->add_block(odata, x7, a2, c1, x6);
}

void Task148::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x1, x2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(x1, x2, x0, a1), 0.0);
  {
    // tensor label: I144
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a1, x0, x2, x1);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, a1.size(), x0.size(), x2.size(), x1.size());
  }
  out()->add_block(odata, x1, x2, x0, a1);
}

void Task149::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I144
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: Gamma48
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x3, x4, x2, x1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x0, x3, x4, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x3.size(), x4.size(), x2.size(), x1.size());
        // tensor label: I145
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x5, a1, x4, x3);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x5, a1, x4, x3)]);
        sort_indices<0,3,2,1,0,1,1,1>(i1data, i1data_sorted, x5.size(), a1.size(), x4.size(), x3.size());
        zgemm3m_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

#endif
