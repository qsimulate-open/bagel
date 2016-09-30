//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelCASA_tasks2.cc
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

#include <src/smith/relcasa/RelCASA_tasks2.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelCASA;

void Task50::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  // tensor label: I10
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(c1, x3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: f1
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a3, c2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(a3, c2)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), c2.size());
      // tensor label: I11
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, a3, c1, x3);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, a3, c1, x3)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), a3.size(), c1.size(), x3.size());
      zgemm3m_("T", "N", 1, c1.size()*x3.size(), c2.size()*a3.size(),
             1.0, i0data_sorted, c2.size()*a3.size(), i1data_sorted, c2.size()*a3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c1.size(), x3.size());
  out()->add_block(odata, c1, x3);
}

void Task51::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I11
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, a3, c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(c2, a3, c1, x3), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c2, a3, c1, x3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c2.size(), a3.size(), c1.size(), x3.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(c1, a3, c2, x3);
    sort_indices<2,1,0,3,1,1,-1,1>(i1data, odata, c1.size(), a3.size(), c2.size(), x3.size());
  }
  out()->add_block(odata, c2, a3, c1, x3);
}

void Task52::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I6
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: Gamma5
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x3, x2, x4, x1, x0);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x3, x2, x4, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x3.size(), x2.size(), x4.size(), x1.size(), x0.size());
        // tensor label: I16
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x5, c1, x4, x3);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x5, c1, x4, x3)]);
        sort_indices<0,3,2,1,0,1,1,1>(i1data, i1data_sorted, x5.size(), c1.size(), x4.size(), x3.size());
        zgemm3m_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task53::Task_local::compute() {
  const Index x5 = b(0);
  const Index c1 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I16
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x5, c1, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(x5, c1, x4, x3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x5, c1, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x5, c1, x4, x3), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a2, x3);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(a2, x3)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a2.size(), x3.size());
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x5, a2, c1, x4);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x5, a2, c1, x4)]);
    sort_indices<1,0,2,3,0,1,-1,1>(i1data, i1data_sorted, x5.size(), a2.size(), c1.size(), x4.size());
    zgemm3m_("T", "N", x3.size(), x5.size()*c1.size()*x4.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x3.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x3.size(), x5.size(), c1.size(), x4.size());
  out()->add_block(odata, x5, c1, x4, x3);
}

void Task54::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I6
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
        // tensor label: I19
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, x5, x4, x3);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, x5, x4, x3)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), x5.size(), x4.size(), x3.size());
        zgemm3m_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task55::Task_local::compute() {
  const Index c1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I19
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x5, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(c1, x5, x4, x3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x5, x4, x3), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a2, x3);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(a2, x3)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a2.size(), x3.size());
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a2, x5, x4);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a2, x5, x4)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size(), x5.size(), x4.size());
    zgemm3m_("T", "N", x3.size(), c1.size()*x5.size()*x4.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x3.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x3.size(), c1.size(), x5.size(), x4.size());
  out()->add_block(odata, c1, x5, x4, x3);
}

void Task56::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I6
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: Gamma70
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x6, x2, x5, x1, x0);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x7, x6, x2, x5, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x2.size(), x5.size(), x1.size(), x0.size());
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x7, x6, c1, x5);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x7, x6, c1, x5)]);
        sort_indices<0,1,3,2,0,1,-1,1>(i1data, i1data_sorted, x7.size(), x6.size(), c1.size(), x5.size());
        zgemm3m_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task57::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I6
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x9 : *range_[1]) {
    for (auto& x8 : *range_[1]) {
      for (auto& x7 : *range_[1]) {
        // tensor label: Gamma71
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x9, x8, x2, x7, x1, x0);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x9, x8, x2, x7, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x9.size(), x8.size(), x2.size(), x7.size(), x1.size(), x0.size());
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x9, x8, c1, x7);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x9, x8, c1, x7)]);
        sort_indices<0,1,3,2,0,1,-1,2>(i1data, i1data_sorted, x9.size(), x8.size(), c1.size(), x7.size());
        zgemm3m_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x9.size()*x8.size()*x7.size(),
               1.0, i0data_sorted, x9.size()*x8.size()*x7.size(), i1data_sorted, x9.size()*x8.size()*x7.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task58::Task_local::compute() {
  const Index c3 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, x0, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(c3, x0, c1, a2), 0.0);
  {
    // tensor label: I21
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a2, c3, x0, c1);
    sort_indices<1,2,3,0,1,1,1,1>(i0data, odata, a2.size(), c3.size(), x0.size(), c1.size());
  }
  out()->add_block(odata, c3, x0, c1, a2);
}

void Task59::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I21
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c3, x0, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c3, x0, c1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c3, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x0, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, x1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), x1.size());
    // tensor label: I22
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a2, c3, x1, x0);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a2, c3, x1, x0)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), c3.size(), x1.size(), x0.size());
    zgemm3m_("T", "N", c1.size(), a2.size()*c3.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), x0.size());
  out()->add_block(odata, a2, c3, x0, c1);
}

void Task60::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I22
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c3, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a2, c3, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma7
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x1, x0, x2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x1, x0, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x0.size(), x2.size());
      // tensor label: t2
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

void Task61::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I22
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c3, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a2, c3, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma9
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, a2, x3, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, a2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,-1,1>(i1data, i1data_sorted, c3.size(), a2.size(), x3.size(), x2.size());
      zgemm3m_("T", "N", x0.size()*x1.size(), c3.size()*a2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c3.size(), a2.size());
  out()->add_block(odata, a2, c3, x1, x0);
}

void Task62::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I21
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c3, x0, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c3, x0, c1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c3, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x0, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, x1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c3, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c3.size(), x1.size());
    // tensor label: I25
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a2, c1, x0, x1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a2, c1, x0, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a2.size(), c1.size(), x0.size(), x1.size());
    zgemm3m_("T", "N", c3.size(), a2.size()*c1.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c3.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), x0.size());
  out()->add_block(odata, a2, c3, x0, c1);
}

void Task63::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I25
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c1, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, x0, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma9
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: I26
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, a2, c1, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, a2, c1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a2.size(), c1.size(), x2.size());
      zgemm3m_("T", "N", x0.size()*x1.size(), a2.size()*c1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a2.size(), c1.size());
  out()->add_block(odata, a2, c1, x0, x1);
}

void Task64::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x2 = b(3);
  // tensor label: I26
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

void Task65::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I21
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c3, x0, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c3, x0, c1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c3, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x0, c1), 0.0);
  // tensor label: f1
  std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a2);
  std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a2)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size());
  // tensor label: I34
  std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, x0);
  std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), x0.size());
  zgemm3m_("T", "N", c1.size()*a2.size(), c3.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a2.size());
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), x0.size());
  out()->add_block(odata, a2, c3, x0, c1);
}

void Task66::Task_local::compute() {
  const Index c3 = b(0);
  const Index x0 = b(1);
  // tensor label: I34
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma9
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x0, x1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x0, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, x2, c3, x1);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, x2, c3, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), c3.size(), x1.size());
        zgemm3m_("T", "N", x0.size(), c3.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c3.size());
  out()->add_block(odata, c3, x0);
}

void Task67::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I21
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c3, x0, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c3, x0, c1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c3, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x0, c1), 0.0);
  // tensor label: f1
  std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a2);
  std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c3, a2)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size());
  // tensor label: I37
  std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, x0);
  std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), x0.size());
  zgemm3m_("T", "N", c3.size()*a2.size(), c1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c3.size()*a2.size());
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), x0.size());
  out()->add_block(odata, a2, c3, x0, c1);
}

void Task68::Task_local::compute() {
  const Index c1 = b(0);
  const Index x0 = b(1);
  // tensor label: I37
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

void Task69::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I21
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c3, x0, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c3, x0, c1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c3, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x0, c1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a4, c3, a2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a4, c3, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
    // tensor label: I40
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, x0);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), x0.size());
    zgemm3m_("T", "N", c1.size()*c3.size()*a2.size(), x0.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, c1.size()*c3.size()*a2.size());
  }
  sort_indices<2,1,3,0,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), a2.size(), x0.size());
  out()->add_block(odata, a2, c3, x0, c1);
}

void Task70::Task_local::compute() {
  const Index a4 = b(0);
  const Index x0 = b(1);
  // tensor label: I40
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a4, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma13
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: f1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, x1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, x1)]);
    sort_indices<1,0,0,1,-2,1>(i1data, i1data_sorted, a4.size(), x1.size());
    zgemm3m_("T", "N", x0.size(), a4.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), a4.size());
  out()->add_block(odata, a4, x0);
}

void Task71::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I21
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c3, x0, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c3, x0, c1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c3, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x0, c1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a2, c3, a4);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a2, c3, a4)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
    // tensor label: I43
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, x0);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), x0.size());
    zgemm3m_("T", "N", c1.size()*a2.size()*c3.size(), x0.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c3.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), x0.size());
  out()->add_block(odata, a2, c3, x0, c1);
}

void Task72::Task_local::compute() {
  const Index a4 = b(0);
  const Index x0 = b(1);
  // tensor label: I43
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a4, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma13
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: f1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, x1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, x1)]);
    sort_indices<1,0,0,1,2,1>(i1data, i1data_sorted, a4.size(), x1.size());
    zgemm3m_("T", "N", x0.size(), a4.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), a4.size());
  out()->add_block(odata, a4, x0);
}

void Task73::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I21
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c3, x0, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c3, x0, c1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c3, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x0, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, a2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size());
    // tensor label: I46
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, c3, x1, x0);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, c3, x1, x0)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), c3.size(), x1.size(), x0.size());
    zgemm3m_("T", "N", a2.size(), c1.size()*c3.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), c3.size(), x0.size());
  out()->add_block(odata, a2, c3, x0, c1);
}

void Task74::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I46
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, c3, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma1
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

void Task75::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I21
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c3, x0, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c3, x0, c1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c3, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x0, c1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a4, a2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(a4, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a2.size());
    // tensor label: I49
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, a4, c1, x0);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, a4, c1, x0)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), a4.size(), c1.size(), x0.size());
    zgemm3m_("T", "N", a2.size(), c3.size()*c1.size()*x0.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, a2.size(), c3.size(), c1.size(), x0.size());
  out()->add_block(odata, a2, c3, x0, c1);
}

void Task76::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index x0 = b(3);
  // tensor label: I49
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, a4, c1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, a4, c1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, a4, c1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a4, c1, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma13
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I50
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, a4, c1, x1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, a4, c1, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, c3.size(), a4.size(), c1.size(), x1.size());
    zgemm3m_("T", "N", x0.size(), c3.size()*a4.size()*c1.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x0.size(), c3.size(), a4.size(), c1.size());
  out()->add_block(odata, c3, a4, c1, x0);
}

void Task77::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index x1 = b(3);
  // tensor label: I50
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, a4, c1, x1)]);
  std::fill_n(odata.get(), out()->get_size(c3, a4, c1, x1), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a4, c1, x1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c3.size(), a4.size(), c1.size(), x1.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(c1, a4, c3, x1);
    sort_indices<2,1,0,3,1,1,1,1>(i1data, odata, c1.size(), a4.size(), c3.size(), x1.size());
  }
  out()->add_block(odata, c3, a4, c1, x1);
}

void Task78::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I21
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c3, x0, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c3, x0, c1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c3, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x0, c1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: Gamma72
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x3);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, x3)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size());
    // tensor label: I217
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, a2, c1, x3);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, a2, c1, x3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), c1.size(), x3.size());
    zgemm3m_("T", "N", x0.size(), c3.size()*a2.size()*c1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<2,1,0,3,1,1,1,1>(odata_sorted, odata, x0.size(), c3.size(), a2.size(), c1.size());
  out()->add_block(odata, a2, c3, x0, c1);
}

void Task79::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I217
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, a2, c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, c1, x3), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a2, c1, x3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c3.size(), a2.size(), c1.size(), x3.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(c1, a2, c3, x3);
    sort_indices<2,1,0,3,1,1,-1,1>(i1data, odata, c1.size(), a2.size(), c3.size(), x3.size());
  }
  out()->add_block(odata, c3, a2, c1, x3);
}

void Task80::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I21
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c3, x0, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c3, x0, c1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c3, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x0, c1), 0.0);
  for (auto& x5 : *range_[1]) {
    // tensor label: Gamma74
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x5);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, x5)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), x5.size());
    // tensor label: I221
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, a2, c1, x5);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, a2, c1, x5)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), c1.size(), x5.size());
    zgemm3m_("T", "N", x0.size(), c3.size()*a2.size()*c1.size(), x5.size(),
           1.0, i0data_sorted, x5.size(), i1data_sorted, x5.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<2,1,0,3,1,1,1,1>(odata_sorted, odata, x0.size(), c3.size(), a2.size(), c1.size());
  out()->add_block(odata, a2, c3, x0, c1);
}

void Task81::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x5 = b(3);
  // tensor label: I221
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, a2, c1, x5)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, c1, x5), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a2, c1, x5);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, c3.size(), a2.size(), c1.size(), x5.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(c1, a2, c3, x5);
    sort_indices<2,1,0,3,1,1,-1,2>(i1data, odata, c1.size(), a2.size(), c3.size(), x5.size());
  }
  out()->add_block(odata, c3, a2, c1, x5);
}

void Task82::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  {
    // tensor label: I54
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a1, c2, x1, x0);
    sort_indices<1,2,3,0,1,1,1,1>(i0data, odata, a1.size(), c2.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, c2, x1, x0, a1);
}

void Task83::Task_local::compute() {
  const Index a1 = b(0);
  const Index c2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I54
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, c2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a1, c2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, c2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, c2, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma18
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x3, x2, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x3, x2, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), x3.size(), x2.size(), x0.size());
      // tensor label: I55
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a1, c2, x3, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a1, c2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a1.size(), c2.size(), x3.size(), x2.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), a1.size()*c2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a1.size(), c2.size());
  out()->add_block(odata, a1, c2, x1, x0);
}

void Task84::Task_local::compute() {
  const Index a1 = b(0);
  const Index c2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I55
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

void Task85::Task_local::compute() {
  const Index a1 = b(0);
  const Index c2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I54
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, c2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a1, c2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, c2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, c2, x1, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: Gamma3
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x2, x3, x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x2, x3, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
      // tensor label: I58
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, a1, x3, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, a1, x3, x2)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size(), x3.size(), x2.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), c2.size()*a1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,2,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a1.size());
  out()->add_block(odata, a1, c2, x1, x0);
}

void Task86::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I58
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

void Task87::Task_local::compute() {
  const Index a1 = b(0);
  const Index c2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I54
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, c2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a1, c2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, c2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, c2, x1, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c2, x2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c2, x2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), x2.size());
    // tensor label: I61
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a1, x0, x1, x2);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a1, x0, x1, x2)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a1.size(), x0.size(), x1.size(), x2.size());
    zgemm3m_("T", "N", c2.size(), a1.size()*x0.size()*x1.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, c2.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x0.size(), x1.size());
  out()->add_block(odata, a1, c2, x1, x0);
}

void Task88::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: I61
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x1, x2)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x1, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x1, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, x2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma20
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

void Task89::Task_local::compute() {
  const Index a1 = b(0);
  const Index c2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I54
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, c2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a1, c2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, c2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, c2, x1, x0), 0.0);
  // tensor label: Gamma21
  std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
  std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
  // tensor label: I64
  std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, a1);
  std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, a1)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size());
  zgemm3m_("T", "N", x1.size()*x0.size(), c2.size()*a1.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, x1.size()*x0.size());
  sort_indices<3,2,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a1.size());
  out()->add_block(odata, a1, c2, x1, x0);
}

void Task90::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  // tensor label: I64
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
      // tensor label: I65
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

void Task91::Task_local::compute() {
  const Index c2 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a1 = b(3);
  // tensor label: I65
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

void Task92::Task_local::compute() {
  const Index a1 = b(0);
  const Index c2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I54
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, c2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a1, c2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, c2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, c2, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma23
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I70
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, c2, a1, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, c2, a1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), c2.size(), a1.size(), x2.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), c2.size()*a1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,2,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a1.size());
  out()->add_block(odata, a1, c2, x1, x0);
}

void Task93::Task_local::compute() {
  const Index x3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x2 = b(3);
  // tensor label: I70
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

void Task94::Task_local::compute() {
  const Index a1 = b(0);
  const Index c2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I54
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, c2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a1, c2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, c2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, c2, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma24
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0, x1, x2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x0, x1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
      // tensor label: I73
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, a1, c2, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, a1, c2, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), c2.size(), x2.size());
      zgemm3m_("T", "N", x0.size()*x1.size(), a1.size()*c2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), c2.size());
  out()->add_block(odata, a1, c2, x1, x0);
}

void Task95::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index c2 = b(2);
  const Index x2 = b(3);
  // tensor label: I73
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

void Task96::Task_local::compute() {
  const Index a1 = b(0);
  const Index c2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I54
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, c2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a1, c2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, c2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, c2, x1, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x2, a1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x2, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), a1.size());
    // tensor label: I76
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, x1, x2, x0);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, x1, x2, x0)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), x1.size(), x2.size(), x0.size());
    zgemm3m_("T", "N", a1.size(), c2.size()*x1.size()*x0.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, a1.size());
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), x1.size(), x0.size());
  out()->add_block(odata, a1, c2, x1, x0);
}

void Task97::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index x0 = b(3);
  // tensor label: I76
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, x1, x2, x0)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x2, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, x1, x2, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x2, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma25
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

void Task98::Task_local::compute() {
  const Index a1 = b(0);
  const Index c2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I54
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, c2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a1, c2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, c2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, c2, x1, x0), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a3, a1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(a3, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), a1.size());
    // tensor label: I79
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a3, c2, x0, x1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a3, c2, x0, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), x0.size(), x1.size());
    zgemm3m_("T", "N", a1.size(), c2.size()*x0.size()*x1.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, a1.size());
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), x0.size(), x1.size());
  out()->add_block(odata, a1, c2, x1, x0);
}

void Task99::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I79
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

#endif
