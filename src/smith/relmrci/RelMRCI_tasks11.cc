//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_tasks11.cc
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

#include <src/smith/relmrci/RelMRCI_tasks11.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelMRCI;

void Task500::Task_local::compute() {
  const Index a2 = b(0);
  const Index x0 = b(1);
  // tensor label: I864
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, x0)]);
  std::fill_n(odata.get(), out()->get_size(a2, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma33
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0, x2, x1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x0, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, a2, x2, x1);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, a2, x2, x1)]);
        sort_indices<0,2,3,1,0,1,-1,1>(i1data, i1data_sorted, x3.size(), a2.size(), x2.size(), x1.size());
        zgemm3m_("T", "N", x0.size(), a2.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), a2.size());
  out()->add_block(odata, a2, x0);
}

void Task501::Task_local::compute() {
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
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a2, c3, a5);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a2, c3, a5)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a5.size());
    // tensor label: I876
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a5, a4);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a5, a4)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a5.size(), a4.size());
    zgemm3m_("T", "N", c1.size()*a2.size()*c3.size(), a4.size(), a5.size(),
           1.0, i0data_sorted, a5.size(), i1data_sorted, a5.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c3.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task502::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  // tensor label: I876
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
      // tensor label: I877
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

void Task503::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I877
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a5, a4, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a5, a4, x1, x0), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a5, a4, x1, x0);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, a5.size(), a4.size(), x1.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(x1, a4, a5, x0);
    sort_indices<2,1,0,3,1,1,-1,1>(i1data, odata, x1.size(), a4.size(), a5.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i2data = in(0)->get_block(x1, x0, a5, a4);
    sort_indices<2,3,0,1,1,1,1,1>(i2data, odata, x1.size(), x0.size(), a5.size(), a4.size());
  }
  out()->add_block(odata, a5, a4, x1, x0);
}

void Task504::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  // tensor label: I876
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
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, a5.size(), x1.size(), x0.size(), a4.size());
      zgemm3m_("T", "N", 1, a5.size()*a4.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a5.size(), a4.size());
  out()->add_block(odata, a5, a4);
}

void Task505::Task_local::compute() {
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
    for (auto& c6 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a5, c6, a4);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a5, c6, a4)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a5.size(), c6.size(), a4.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, c6, a5, a2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, c6, a5, a2)]);
      sort_indices<2,1,0,3,0,1,-2,1>(i1data, i1data_sorted, c3.size(), c6.size(), a5.size(), a2.size());
      zgemm3m_("T", "N", c1.size()*a4.size(), c3.size()*a2.size(), c6.size()*a5.size(),
             1.0, i0data_sorted, c6.size()*a5.size(), i1data_sorted, c6.size()*a5.size(),
             1.0, odata_sorted, c1.size()*a4.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task506::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& c6 : *range_[0]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a4, c6, a5);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a4, c6, a5)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c6.size(), a5.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, c6, a5, a2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, c6, a5, a2)]);
      sort_indices<1,2,0,3,0,1,2,1>(i1data, i1data_sorted, c3.size(), c6.size(), a5.size(), a2.size());
      zgemm3m_("T", "N", c1.size()*a4.size(), c3.size()*a2.size(), c6.size()*a5.size(),
             1.0, i0data_sorted, c6.size()*a5.size(), i1data_sorted, c6.size()*a5.size(),
             1.0, odata_sorted, c1.size()*a4.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task507::Task_local::compute() {
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
    for (auto& c6 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a5, c6, a2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a5, c6, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a5.size(), c6.size(), a2.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, c6, a5, a4);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, c6, a5, a4)]);
      sort_indices<2,1,0,3,0,1,2,1>(i1data, i1data_sorted, c3.size(), c6.size(), a5.size(), a4.size());
      zgemm3m_("T", "N", c1.size()*a2.size(), c3.size()*a4.size(), c6.size()*a5.size(),
             1.0, i0data_sorted, c6.size()*a5.size(), i1data_sorted, c6.size()*a5.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task508::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& c6 : *range_[0]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a2, c6, a5);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a2, c6, a5)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c6.size(), a5.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, c6, a5, a4);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, c6, a5, a4)]);
      sort_indices<1,2,0,3,0,1,-2,1>(i1data, i1data_sorted, c3.size(), c6.size(), a5.size(), a4.size());
      zgemm3m_("T", "N", c1.size()*a2.size(), c3.size()*a4.size(), c6.size()*a5.size(),
             1.0, i0data_sorted, c6.size()*a5.size(), i1data_sorted, c6.size()*a5.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task509::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& a6 : *range_[2]) {
    for (auto& c5 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a6, c5, a4);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a6, c5, a4)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a6.size(), c5.size(), a4.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, a2, a6, c5);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, a2, a6, c5)]);
      sort_indices<2,3,0,1,0,1,2,1>(i1data, i1data_sorted, c3.size(), a2.size(), a6.size(), c5.size());
      zgemm3m_("T", "N", c1.size()*a4.size(), c3.size()*a2.size(), a6.size()*c5.size(),
             1.0, i0data_sorted, a6.size()*c5.size(), i1data_sorted, a6.size()*c5.size(),
             1.0, odata_sorted, c1.size()*a4.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task510::Task_local::compute() {
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
    for (auto& a6 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a4, c5, a6);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a4, c5, a6)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c5.size(), a6.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, a2, a6, c5);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, a2, a6, c5)]);
      sort_indices<3,2,0,1,0,1,-2,1>(i1data, i1data_sorted, c3.size(), a2.size(), a6.size(), c5.size());
      zgemm3m_("T", "N", c1.size()*a4.size(), c3.size()*a2.size(), a6.size()*c5.size(),
             1.0, i0data_sorted, a6.size()*c5.size(), i1data_sorted, a6.size()*c5.size(),
             1.0, odata_sorted, c1.size()*a4.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task511::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I108
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& a6 : *range_[2]) {
    for (auto& c5 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a6, c5, a2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a6, c5, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a6.size(), c5.size(), a2.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, a4, a6, c5);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, a4, a6, c5)]);
      sort_indices<2,3,0,1,0,1,-2,1>(i1data, i1data_sorted, c3.size(), a4.size(), a6.size(), c5.size());
      zgemm3m_("T", "N", c1.size()*a2.size(), c3.size()*a4.size(), a6.size()*c5.size(),
             1.0, i0data_sorted, a6.size()*c5.size(), i1data_sorted, a6.size()*c5.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task512::Task_local::compute() {
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
    for (auto& a6 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a2, c5, a6);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a2, c5, a6)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c5.size(), a6.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, a4, a6, c5);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, a4, a6, c5)]);
      sort_indices<3,2,0,1,0,1,2,1>(i1data, i1data_sorted, c3.size(), a4.size(), a6.size(), c5.size());
      zgemm3m_("T", "N", c1.size()*a2.size(), c3.size()*a4.size(), a6.size()*c5.size(),
             1.0, i0data_sorted, a6.size()*c5.size(), i1data_sorted, a6.size()*c5.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task513::Task_local::compute() {
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
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, a4, c1, a2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, a4, c1, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a4.size(), c1.size(), a2.size());
    // tensor label: I935
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, x3);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, x3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), x3.size());
    zgemm3m_("T", "N", a4.size()*c1.size()*a2.size(), c3.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a4.size()*c1.size()*a2.size());
  }
  sort_indices<2,1,0,3,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), a2.size(), c3.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task514::Task_local::compute() {
  const Index c3 = b(0);
  const Index x3 = b(1);
  // tensor label: I935
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, x3)]);
  std::fill_n(odata.get(), out()->get_size(c3, x3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: Gamma24
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x1, x0);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x1, x0)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, x2, x1, x0);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, x2, x1, x0)]);
        sort_indices<1,2,3,0,0,1,-1,2>(i1data, i1data_sorted, c3.size(), x2.size(), x1.size(), x0.size());
        zgemm3m_("T", "N", x3.size(), c3.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), c3.size());
  out()->add_block(odata, c3, x3);
}

void Task515::Task_local::compute() {
  const Index c3 = b(0);
  const Index x3 = b(1);
  // tensor label: I935
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, x3)]);
  std::fill_n(odata.get(), out()->get_size(c3, x3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x3), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma33
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0, x2, x1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x0, x2, x1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x2, x1, c3, x0);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x2, x1, c3, x0)]);
        sort_indices<3,0,1,2,0,1,-1,2>(i1data, i1data_sorted, x2.size(), x1.size(), c3.size(), x0.size());
        zgemm3m_("T", "N", x3.size(), c3.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), c3.size());
  out()->add_block(odata, c3, x3);
}

void Task516::Task_local::compute() {
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
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, a2, c1, a4);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, a2, c1, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c1.size(), a4.size());
    // tensor label: I938
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, x3);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, x3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), x3.size());
    zgemm3m_("T", "N", a2.size()*c1.size()*a4.size(), c3.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a2.size()*c1.size()*a4.size());
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), a4.size(), c3.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task517::Task_local::compute() {
  const Index c3 = b(0);
  const Index x3 = b(1);
  // tensor label: I938
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, x3)]);
  std::fill_n(odata.get(), out()->get_size(c3, x3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: Gamma24
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x1, x0);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x1, x0)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, x2, x1, x0);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, x2, x1, x0)]);
        sort_indices<1,2,3,0,0,1,1,2>(i1data, i1data_sorted, c3.size(), x2.size(), x1.size(), x0.size());
        zgemm3m_("T", "N", x3.size(), c3.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), c3.size());
  out()->add_block(odata, c3, x3);
}

void Task518::Task_local::compute() {
  const Index c3 = b(0);
  const Index x3 = b(1);
  // tensor label: I938
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, x3)]);
  std::fill_n(odata.get(), out()->get_size(c3, x3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x3), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma33
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0, x2, x1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x0, x2, x1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x2, x1, c3, x0);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x2, x1, c3, x0)]);
        sort_indices<3,0,1,2,0,1,1,2>(i1data, i1data_sorted, x2.size(), x1.size(), c3.size(), x0.size());
        zgemm3m_("T", "N", x3.size(), c3.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), c3.size());
  out()->add_block(odata, c3, x3);
}

void Task519::Task_local::compute() {
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
    for (auto& c5 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, a4, c5, a2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, a4, c5, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c5.size(), a2.size());
      // tensor label: I947
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, c3, c5, x1);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, c3, c5, x1)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), c3.size(), c5.size(), x1.size());
      zgemm3m_("T", "N", a4.size()*a2.size(), c1.size()*c3.size(), c5.size()*x1.size(),
             1.0, i0data_sorted, c5.size()*x1.size(), i1data_sorted, c5.size()*x1.size(),
             1.0, odata_sorted, a4.size()*a2.size());
    }
  }
  sort_indices<1,2,0,3,1,1,1,1>(odata_sorted, odata, a4.size(), a2.size(), c1.size(), c3.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task520::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index c5 = b(2);
  const Index x1 = b(3);
  // tensor label: I947
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, c3, c5, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, c5, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c3, c5, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, c5, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma27
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, x0, c3, c5);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, x0, c3, c5)]);
    sort_indices<1,0,2,3,0,1,-1,1>(i1data, i1data_sorted, c1.size(), x0.size(), c3.size(), c5.size());
    zgemm3m_("T", "N", x1.size(), c1.size()*c3.size()*c5.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), c1.size(), c3.size(), c5.size());
  out()->add_block(odata, c1, c3, c5, x1);
}

void Task521::Task_local::compute() {
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
    for (auto& c5 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, a2, c5, a4);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, a2, c5, a4)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c5.size(), a4.size());
      // tensor label: I950
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, c3, c5, x1);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, c3, c5, x1)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), c3.size(), c5.size(), x1.size());
      zgemm3m_("T", "N", a2.size()*a4.size(), c1.size()*c3.size(), c5.size()*x1.size(),
             1.0, i0data_sorted, c5.size()*x1.size(), i1data_sorted, c5.size()*x1.size(),
             1.0, odata_sorted, a2.size()*a4.size());
    }
  }
  sort_indices<0,2,1,3,1,1,1,1>(odata_sorted, odata, a2.size(), a4.size(), c1.size(), c3.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task522::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index c5 = b(2);
  const Index x1 = b(3);
  // tensor label: I950
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, c3, c5, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, c5, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c3, c5, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, c5, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma27
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, x0, c3, c5);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, x0, c3, c5)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), x0.size(), c3.size(), c5.size());
    zgemm3m_("T", "N", x1.size(), c1.size()*c3.size()*c5.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), c1.size(), c3.size(), c5.size());
  out()->add_block(odata, c1, c3, c5, x1);
}

void Task523::Task_local::compute() {
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
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, a5, c1, a4);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, a5, c1, a4)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a5.size(), c1.size(), a4.size());
      // tensor label: I953
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a5, c3, a2, x1);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a5, c3, a2, x1)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a5.size(), c3.size(), a2.size(), x1.size());
      zgemm3m_("T", "N", c1.size()*a4.size(), c3.size()*a2.size(), a5.size()*x1.size(),
             1.0, i0data_sorted, a5.size()*x1.size(), i1data_sorted, a5.size()*x1.size(),
             1.0, odata_sorted, c1.size()*a4.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task524::Task_local::compute() {
  const Index a5 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index x1 = b(3);
  // tensor label: I953
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a5, c3, a2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a5, c3, a2, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a5, c3, a2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a5, c3, a2, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma27
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I954
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a5, x0, c3, a2);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a5, x0, c3, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a5.size(), x0.size(), c3.size(), a2.size());
    zgemm3m_("T", "N", x1.size(), a5.size()*c3.size()*a2.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), a5.size(), c3.size(), a2.size());
  out()->add_block(odata, a5, c3, a2, x1);
}

void Task525::Task_local::compute() {
  const Index a5 = b(0);
  const Index x0 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I954
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a5, x0, c3, a2)]);
  std::fill_n(odata.get(), out()->get_size(a5, x0, c3, a2), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a5, x0, c3, a2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, a5.size(), x0.size(), c3.size(), a2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(c3, x0, a5, a2);
    sort_indices<2,1,0,3,1,1,1,1>(i1data, odata, c3.size(), x0.size(), a5.size(), a2.size());
  }
  out()->add_block(odata, a5, x0, c3, a2);
}

void Task526::Task_local::compute() {
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
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, a4, c1, a5);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, a4, c1, a5)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c1.size(), a5.size());
      // tensor label: I956
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a5, c3, a2, x1);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a5, c3, a2, x1)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a5.size(), c3.size(), a2.size(), x1.size());
      zgemm3m_("T", "N", a4.size()*c1.size(), c3.size()*a2.size(), a5.size()*x1.size(),
             1.0, i0data_sorted, a5.size()*x1.size(), i1data_sorted, a5.size()*x1.size(),
             1.0, odata_sorted, a4.size()*c1.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), c3.size(), a2.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task527::Task_local::compute() {
  const Index a5 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index x1 = b(3);
  // tensor label: I956
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a5, c3, a2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a5, c3, a2, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a5, c3, a2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a5, c3, a2, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma27
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I957
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a5, x0, c3, a2);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a5, x0, c3, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a5.size(), x0.size(), c3.size(), a2.size());
    zgemm3m_("T", "N", x1.size(), a5.size()*c3.size()*a2.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), a5.size(), c3.size(), a2.size());
  out()->add_block(odata, a5, c3, a2, x1);
}

void Task528::Task_local::compute() {
  const Index a5 = b(0);
  const Index x0 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I957
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a5, x0, c3, a2)]);
  std::fill_n(odata.get(), out()->get_size(a5, x0, c3, a2), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a5, x0, c3, a2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, a5.size(), x0.size(), c3.size(), a2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(c3, x0, a5, a2);
    sort_indices<2,1,0,3,1,1,-1,1>(i1data, odata, c3.size(), x0.size(), a5.size(), a2.size());
  }
  out()->add_block(odata, a5, x0, c3, a2);
}

void Task529::Task_local::compute() {
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
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, a5, c1, a2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, a5, c1, a2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a5.size(), c1.size(), a2.size());
      // tensor label: I959
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a5, c3, a4, x1);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a5, c3, a4, x1)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a5.size(), c3.size(), a4.size(), x1.size());
      zgemm3m_("T", "N", c1.size()*a2.size(), c3.size()*a4.size(), a5.size()*x1.size(),
             1.0, i0data_sorted, a5.size()*x1.size(), i1data_sorted, a5.size()*x1.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task530::Task_local::compute() {
  const Index a5 = b(0);
  const Index c3 = b(1);
  const Index a4 = b(2);
  const Index x1 = b(3);
  // tensor label: I959
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a5, c3, a4, x1)]);
  std::fill_n(odata.get(), out()->get_size(a5, c3, a4, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a5, c3, a4, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a5, c3, a4, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma27
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I960
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a5, x0, c3, a4);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a5, x0, c3, a4)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a5.size(), x0.size(), c3.size(), a4.size());
    zgemm3m_("T", "N", x1.size(), a5.size()*c3.size()*a4.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), a5.size(), c3.size(), a4.size());
  out()->add_block(odata, a5, c3, a4, x1);
}

void Task531::Task_local::compute() {
  const Index a5 = b(0);
  const Index x0 = b(1);
  const Index c3 = b(2);
  const Index a4 = b(3);
  // tensor label: I960
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a5, x0, c3, a4)]);
  std::fill_n(odata.get(), out()->get_size(a5, x0, c3, a4), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a5, x0, c3, a4);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, a5.size(), x0.size(), c3.size(), a4.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(c3, x0, a5, a4);
    sort_indices<2,1,0,3,1,1,-1,1>(i1data, odata, c3.size(), x0.size(), a5.size(), a4.size());
  }
  out()->add_block(odata, a5, x0, c3, a4);
}

void Task532::Task_local::compute() {
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
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, a2, c1, a5);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, a2, c1, a5)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c1.size(), a5.size());
      // tensor label: I962
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a5, c3, a4, x1);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a5, c3, a4, x1)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a5.size(), c3.size(), a4.size(), x1.size());
      zgemm3m_("T", "N", a2.size()*c1.size(), c3.size()*a4.size(), a5.size()*x1.size(),
             1.0, i0data_sorted, a5.size()*x1.size(), i1data_sorted, a5.size()*x1.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), c3.size(), a4.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task533::Task_local::compute() {
  const Index a5 = b(0);
  const Index c3 = b(1);
  const Index a4 = b(2);
  const Index x1 = b(3);
  // tensor label: I962
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a5, c3, a4, x1)]);
  std::fill_n(odata.get(), out()->get_size(a5, c3, a4, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a5, c3, a4, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a5, c3, a4, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma27
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I963
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a5, x0, c3, a4);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a5, x0, c3, a4)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a5.size(), x0.size(), c3.size(), a4.size());
    zgemm3m_("T", "N", x1.size(), a5.size()*c3.size()*a4.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), a5.size(), c3.size(), a4.size());
  out()->add_block(odata, a5, c3, a4, x1);
}

void Task534::Task_local::compute() {
  const Index a5 = b(0);
  const Index x0 = b(1);
  const Index c3 = b(2);
  const Index a4 = b(3);
  // tensor label: I963
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a5, x0, c3, a4)]);
  std::fill_n(odata.get(), out()->get_size(a5, x0, c3, a4), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a5, x0, c3, a4);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, a5.size(), x0.size(), c3.size(), a4.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(c3, x0, a5, a4);
    sort_indices<2,1,0,3,1,1,1,1>(i1data, odata, c3.size(), x0.size(), a5.size(), a4.size());
  }
  out()->add_block(odata, a5, x0, c3, a4);
}

void Task535::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  {
    // tensor label: I134
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c2, a3, x0, a1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c2.size(), a3.size(), x0.size(), a1.size());
  }
  out()->add_block(odata, c2, a3, x0, a1);
}

void Task536::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I134
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, a1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size());
    // tensor label: I135
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, a3, x1, x0);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, a3, x1, x0)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), a3.size(), x1.size(), x0.size());
    zgemm3m_("T", "N", a1.size(), c2.size()*a3.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size());
  out()->add_block(odata, c2, a3, x0, a1);
}

void Task537::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I135
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, a3, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c2, a3, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, a3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a3, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma24
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, a3, x3, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, a3, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), a3.size(), x3.size(), x2.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), c2.size()*a3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a3.size());
  out()->add_block(odata, c2, a3, x1, x0);
}

void Task538::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I134
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, a3);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, a3)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size());
    // tensor label: I138
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, a1, x1, x0);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, a1, x1, x0)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size(), x1.size(), x0.size());
    zgemm3m_("T", "N", a3.size(), c2.size()*a1.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a3.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size());
  out()->add_block(odata, c2, a3, x0, a1);
}

void Task539::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I138
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, a1, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c2, a1, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, a1, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma24
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
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a1.size());
  out()->add_block(odata, c2, a1, x1, x0);
}

void Task540::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I134
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  // tensor label: h1
  std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c2, a1);
  std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c2, a1)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size());
  // tensor label: I141
  std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a3, x0);
  std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a3, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), x0.size());
  zgemm3m_("T", "N", c2.size()*a1.size(), a3.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a1.size());
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), a3.size(), x0.size());
  out()->add_block(odata, c2, a3, x0, a1);
}

void Task541::Task_local::compute() {
  const Index a3 = b(0);
  const Index x0 = b(1);
  // tensor label: I141
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a3, x0)]);
  std::fill_n(odata.get(), out()->get_size(a3, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma33
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0, x2, x1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x0, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, a3, x2, x1);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, a3, x2, x1)]);
        sort_indices<0,2,3,1,0,1,-1,1>(i1data, i1data_sorted, x3.size(), a3.size(), x2.size(), x1.size());
        zgemm3m_("T", "N", x0.size(), a3.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), a3.size());
  out()->add_block(odata, a3, x0);
}

void Task542::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I134
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  // tensor label: h1
  std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c2, a3);
  std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c2, a3)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size());
  // tensor label: I144
  std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a1, x0);
  std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a1.size(), x0.size());
  zgemm3m_("T", "N", c2.size()*a3.size(), a1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a3.size());
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), a1.size(), x0.size());
  out()->add_block(odata, c2, a3, x0, a1);
}

void Task543::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  // tensor label: I144
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma33
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0, x2, x1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x0, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, a1, x2, x1);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, a1, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), x2.size(), x1.size());
        zgemm3m_("T", "N", x0.size(), a1.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), a1.size());
  out()->add_block(odata, a1, x0);
}

void Task544::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I134
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c2, a3, c4, a1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c2, a3, c4, a1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c4.size(), a1.size());
    // tensor label: I147
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c4, x0);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c4, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c4.size(), x0.size());
    zgemm3m_("T", "N", c2.size()*a3.size()*a1.size(), x0.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c2.size()*a3.size()*a1.size());
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), a1.size(), x0.size());
  out()->add_block(odata, c2, a3, x0, a1);
}

void Task545::Task_local::compute() {
  const Index c4 = b(0);
  const Index x0 = b(1);
  // tensor label: I147
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c4, x0)]);
  std::fill_n(odata.get(), out()->get_size(c4, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c4, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma27
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x1, c4);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x1, c4)]);
    sort_indices<0,1,0,1,-2,1>(i1data, i1data_sorted, x1.size(), c4.size());
    zgemm3m_("T", "N", x0.size(), c4.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c4.size());
  out()->add_block(odata, c4, x0);
}

void Task546::Task_local::compute() {
  const Index c4 = b(0);
  const Index x0 = b(1);
  // tensor label: I147
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c4, x0)]);
  std::fill_n(odata.get(), out()->get_size(c4, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma33
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0, x2, x1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x0, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, c4, x2, x1);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, c4, x2, x1)]);
        sort_indices<0,2,3,1,0,1,-1,1>(i1data, i1data_sorted, x3.size(), c4.size(), x2.size(), x1.size());
        zgemm3m_("T", "N", x0.size(), c4.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c4.size());
  out()->add_block(odata, c4, x0);
}

void Task547::Task_local::compute() {
  const Index c4 = b(0);
  const Index x0 = b(1);
  // tensor label: I147
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c4, x0)]);
  std::fill_n(odata.get(), out()->get_size(c4, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma24
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x1, x0);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x1, x0)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, x2, x1, c4);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, x2, x1, c4)]);
        sort_indices<0,1,2,3,0,1,-1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), c4.size());
        zgemm3m_("T", "N", x0.size(), c4.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c4.size());
  out()->add_block(odata, c4, x0);
}

void Task548::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I134
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c2, a1, c4, a3);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c2, a1, c4, a3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c4.size(), a3.size());
    // tensor label: I150
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c4, x0);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c4, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c4.size(), x0.size());
    zgemm3m_("T", "N", c2.size()*a1.size()*a3.size(), x0.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c2.size()*a1.size()*a3.size());
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), a3.size(), x0.size());
  out()->add_block(odata, c2, a3, x0, a1);
}

void Task549::Task_local::compute() {
  const Index c4 = b(0);
  const Index x0 = b(1);
  // tensor label: I150
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c4, x0)]);
  std::fill_n(odata.get(), out()->get_size(c4, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c4, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma27
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x1, c4);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x1, c4)]);
    sort_indices<0,1,0,1,2,1>(i1data, i1data_sorted, x1.size(), c4.size());
    zgemm3m_("T", "N", x0.size(), c4.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c4.size());
  out()->add_block(odata, c4, x0);
}

#endif
