//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_tasks14.cc
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

#include <src/smith/relmrci/RelMRCI_tasks14.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelMRCI;

void Task650::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I134
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, a4, c2, a1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, a4, c2, a1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a4.size(), c2.size(), a1.size());
      // tensor label: I1109
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, a3, x3, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, a3, x3, x0)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), a3.size(), x3.size(), x0.size());
      zgemm3m_("T", "N", c2.size()*a1.size(), a3.size()*x0.size(), a4.size()*x3.size(),
             1.0, i0data_sorted, a4.size()*x3.size(), i1data_sorted, a4.size()*x3.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), a3.size(), x0.size());
  out()->add_block(odata, c2, a3, x0, a1);
}

void Task651::Task_local::compute() {
  const Index a4 = b(0);
  const Index a3 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1109
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a4, a3, x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, a3, x3, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a4, a3, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a3, x3, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma33
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: I1110
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, a3, x2, x1);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, a3, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a3.size(), x2.size(), x1.size());
      zgemm3m_("T", "N", x3.size()*x0.size(), a4.size()*a3.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a4.size(), a3.size());
  out()->add_block(odata, a4, a3, x3, x0);
}

void Task652::Task_local::compute() {
  const Index a4 = b(0);
  const Index a3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1110
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a4, a3, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a4, a3, x2, x1), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a4, a3, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, a4.size(), a3.size(), x2.size(), x1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(x2, x1, a4, a3);
    sort_indices<2,3,0,1,1,1,-1,2>(i1data, odata, x2.size(), x1.size(), a4.size(), a3.size());
  }
  out()->add_block(odata, a4, a3, x2, x1);
}

void Task653::Task_local::compute() {
  const Index a4 = b(0);
  const Index a3 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1109
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a4, a3, x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, a3, x3, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a4, a3, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a3, x3, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma24
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, x2, x1, a3);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, x2, x1, a3)]);
      sort_indices<1,2,0,3,0,1,-1,2>(i1data, i1data_sorted, a4.size(), x2.size(), x1.size(), a3.size());
      zgemm3m_("T", "N", x3.size()*x0.size(), a4.size()*a3.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a4.size(), a3.size());
  out()->add_block(odata, a4, a3, x3, x0);
}

void Task654::Task_local::compute() {
  const Index a4 = b(0);
  const Index a3 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1109
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a4, a3, x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, a3, x3, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a4, a3, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a3, x3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma368
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x1, x2, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x1, x2, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x2.size(), x0.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x2, a3, a4, x1);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x2, a3, a4, x1)]);
      sort_indices<3,0,1,2,0,1,1,2>(i1data, i1data_sorted, x2.size(), a3.size(), a4.size(), x1.size());
      zgemm3m_("T", "N", x3.size()*x0.size(), a3.size()*a4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<3,2,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a3.size(), a4.size());
  out()->add_block(odata, a4, a3, x3, x0);
}

void Task655::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I134
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, a1, c2, a4);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, a1, c2, a4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a4.size());
      // tensor label: I1112
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, a3, x3, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, a3, x3, x0)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), a3.size(), x3.size(), x0.size());
      zgemm3m_("T", "N", a1.size()*c2.size(), a3.size()*x0.size(), a4.size()*x3.size(),
             1.0, i0data_sorted, a4.size()*x3.size(), i1data_sorted, a4.size()*x3.size(),
             1.0, odata_sorted, a1.size()*c2.size());
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size());
  out()->add_block(odata, c2, a3, x0, a1);
}

void Task656::Task_local::compute() {
  const Index a4 = b(0);
  const Index a3 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1112
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a4, a3, x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, a3, x3, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a4, a3, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a3, x3, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma33
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: I1113
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, a3, x2, x1);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, a3, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a3.size(), x2.size(), x1.size());
      zgemm3m_("T", "N", x3.size()*x0.size(), a4.size()*a3.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a4.size(), a3.size());
  out()->add_block(odata, a4, a3, x3, x0);
}

void Task657::Task_local::compute() {
  const Index a4 = b(0);
  const Index a3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1113
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a4, a3, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a4, a3, x2, x1), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a4, a3, x2, x1);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, a4.size(), a3.size(), x2.size(), x1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(x2, x1, a4, a3);
    sort_indices<2,3,0,1,1,1,1,2>(i1data, odata, x2.size(), x1.size(), a4.size(), a3.size());
  }
  out()->add_block(odata, a4, a3, x2, x1);
}

void Task658::Task_local::compute() {
  const Index a4 = b(0);
  const Index a3 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1112
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a4, a3, x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, a3, x3, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a4, a3, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a3, x3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma363
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0, x1, x2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x0, x1, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, x2, x1, a3);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, x2, x1, a3)]);
      sort_indices<2,1,0,3,0,1,1,2>(i1data, i1data_sorted, a4.size(), x2.size(), x1.size(), a3.size());
      zgemm3m_("T", "N", x3.size()*x0.size(), a4.size()*a3.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a4.size(), a3.size());
  out()->add_block(odata, a4, a3, x3, x0);
}

void Task659::Task_local::compute() {
  const Index a4 = b(0);
  const Index a3 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1112
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a4, a3, x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, a3, x3, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a4, a3, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a3, x3, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x2, a3, a4, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x2, a3, a4, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x2.size(), a3.size(), a4.size(), x1.size());
      // tensor label: Gamma33
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, x0, x2, x1);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, x0, x2, x1)]);
      sort_indices<2,3,0,1,0,1,-1,2>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      zgemm3m_("T", "N", a3.size()*a4.size(), x3.size()*x0.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, a3.size()*a4.size());
    }
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, a3.size(), a4.size(), x3.size(), x0.size());
  out()->add_block(odata, a4, a3, x3, x0);
}

void Task660::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I134
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, a1, x4, a3);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, a1, x4, a3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a3.size());
      // tensor label: I1199
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, x5, x0, x4);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, x5, x0, x4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, c2.size(), x5.size(), x0.size(), x4.size());
      zgemm3m_("T", "N", a1.size()*a3.size(), c2.size()*x0.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, a1.size()*a3.size());
    }
  }
  sort_indices<2,1,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), c2.size(), x0.size());
  out()->add_block(odata, c2, a3, x0, a1);
}

void Task661::Task_local::compute() {
  const Index c2 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  // tensor label: I1199
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, x5, x0, x4)]);
  std::fill_n(odata.get(), out()->get_size(c2, x5, x0, x4), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, x5, x0, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x5, x0, x4), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma32
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<3,4,5,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, x3, x2, x1);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, x3, x2, x1)]);
        sort_indices<1,2,3,0,0,1,-1,1>(i1data, i1data_sorted, c2.size(), x3.size(), x2.size(), x1.size());
        zgemm3m_("T", "N", x5.size()*x0.size()*x4.size(), c2.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x0.size()*x4.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x0.size(), x4.size(), c2.size());
  out()->add_block(odata, c2, x5, x0, x4);
}

void Task662::Task_local::compute() {
  const Index c2 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  // tensor label: I1199
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, x5, x0, x4)]);
  std::fill_n(odata.get(), out()->get_size(c2, x5, x0, x4), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, x5, x0, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x5, x0, x4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma391
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x4, x1, x3, x2);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x0, x4, x1, x3, x2)]);
        sort_indices<3,4,5,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x1.size(), x3.size(), x2.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, x2, c2, x1);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, x2, c2, x1)]);
        sort_indices<3,0,1,2,0,1,-1,1>(i1data, i1data_sorted, x3.size(), x2.size(), c2.size(), x1.size());
        zgemm3m_("T", "N", x5.size()*x0.size()*x4.size(), c2.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x0.size()*x4.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x0.size(), x4.size(), c2.size());
  out()->add_block(odata, c2, x5, x0, x4);
}

void Task663::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I134
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a4, x1, c2, a1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(a4, x1, c2, a1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, a4.size(), x1.size(), c2.size(), a1.size());
      // tensor label: I1205
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a3, a4, x0, x1);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a3, a4, x0, x1)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, a3.size(), a4.size(), x0.size(), x1.size());
      zgemm3m_("T", "N", c2.size()*a1.size(), a3.size()*x0.size(), a4.size()*x1.size(),
             1.0, i0data_sorted, a4.size()*x1.size(), i1data_sorted, a4.size()*x1.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), a3.size(), x0.size());
  out()->add_block(odata, c2, a3, x0, a1);
}

void Task664::Task_local::compute() {
  const Index a3 = b(0);
  const Index a4 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1205
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a3, a4, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a3, a4, x0, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a3, a4, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, a4, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma33
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, a3, x2, a4);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, a3, x2, a4)]);
      sort_indices<0,2,1,3,0,1,-2,1>(i1data, i1data_sorted, x3.size(), a3.size(), x2.size(), a4.size());
      zgemm3m_("T", "N", x0.size()*x1.size(), a3.size()*a4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a3.size(), a4.size());
  out()->add_block(odata, a3, a4, x0, x1);
}

void Task665::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I134
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a4, x1, c2, a3);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(a4, x1, c2, a3)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, a4.size(), x1.size(), c2.size(), a3.size());
      // tensor label: I1208
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a1, a4, x0, x1);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a1, a4, x0, x1)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, a1.size(), a4.size(), x0.size(), x1.size());
      zgemm3m_("T", "N", c2.size()*a3.size(), a1.size()*x0.size(), a4.size()*x1.size(),
             1.0, i0data_sorted, a4.size()*x1.size(), i1data_sorted, a4.size()*x1.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), a1.size(), x0.size());
  out()->add_block(odata, c2, a3, x0, a1);
}

void Task666::Task_local::compute() {
  const Index a1 = b(0);
  const Index a4 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1208
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, a4, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, a4, x0, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, a4, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, a4, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma33
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, a1, x2, a4);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, a1, x2, a4)]);
      sort_indices<0,2,1,3,0,1,2,1>(i1data, i1data_sorted, x3.size(), a1.size(), x2.size(), a4.size());
      zgemm3m_("T", "N", x0.size()*x1.size(), a1.size()*a4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), a4.size());
  out()->add_block(odata, a1, a4, x0, x1);
}

void Task667::Task_local::compute() {
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
    for (auto& a4 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c2, x1, a4, a1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c2, x1, a4, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), x1.size(), a4.size(), a1.size());
      // tensor label: I1211
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a3, a4, x1, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a3, a4, x1, x0)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), a4.size(), x1.size(), x0.size());
      zgemm3m_("T", "N", c2.size()*a1.size(), a3.size()*x0.size(), a4.size()*x1.size(),
             1.0, i0data_sorted, a4.size()*x1.size(), i1data_sorted, a4.size()*x1.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), a3.size(), x0.size());
  out()->add_block(odata, c2, a3, x0, a1);
}

void Task668::Task_local::compute() {
  const Index a3 = b(0);
  const Index a4 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1211
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a3, a4, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a3, a4, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a3, a4, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, a4, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma368
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x1, x2, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x1, x2, x0)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x2.size(), x0.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, a3, x2, a4);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, a3, x2, a4)]);
      sort_indices<0,2,1,3,0,1,2,1>(i1data, i1data_sorted, x3.size(), a3.size(), x2.size(), a4.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), a3.size()*a4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a3.size(), a4.size());
  out()->add_block(odata, a3, a4, x1, x0);
}

void Task669::Task_local::compute() {
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
    for (auto& a4 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c2, x1, a4, a3);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c2, x1, a4, a3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), x1.size(), a4.size(), a3.size());
      // tensor label: I1214
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a1, a4, x0, x1);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a1, a4, x0, x1)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, a1.size(), a4.size(), x0.size(), x1.size());
      zgemm3m_("T", "N", c2.size()*a3.size(), a1.size()*x0.size(), a4.size()*x1.size(),
             1.0, i0data_sorted, a4.size()*x1.size(), i1data_sorted, a4.size()*x1.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), a1.size(), x0.size());
  out()->add_block(odata, c2, a3, x0, a1);
}

void Task670::Task_local::compute() {
  const Index a1 = b(0);
  const Index a4 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1214
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, a4, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, a4, x0, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, a4, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, a4, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma33
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, a1, x2, a4);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, a1, x2, a4)]);
      sort_indices<0,2,1,3,0,1,-2,1>(i1data, i1data_sorted, x3.size(), a1.size(), x2.size(), a4.size());
      zgemm3m_("T", "N", x0.size()*x1.size(), a1.size()*a4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), a4.size());
  out()->add_block(odata, a1, a4, x0, x1);
}

void Task671::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I134
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: Gamma428
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size());
    // tensor label: I1297
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, a3, c2, a1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a3.size(), c2.size(), a1.size());
    zgemm3m_("T", "N", x0.size(), a3.size()*c2.size()*a1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<2,1,0,3,1,1,1,1>(odata_sorted, odata, x0.size(), a3.size(), c2.size(), a1.size());
  out()->add_block(odata, c2, a3, x0, a1);
}

void Task672::Task_local::compute() {
  const Index x3 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I1297
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x3, a3, c2, a1)]);
  std::fill_n(odata.get(), out()->get_size(x3, a3, c2, a1), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, a3, c2, a1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), a3.size(), c2.size(), a1.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(x3, a1, c2, a3);
    sort_indices<0,3,2,1,1,1,-1,1>(i1data, odata, x3.size(), a1.size(), c2.size(), a3.size());
  }
  out()->add_block(odata, x3, a3, c2, a1);
}

void Task673::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I134
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    // tensor label: Gamma430
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size());
    // tensor label: I1301
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x5, a3, c2, a1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x5, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), a3.size(), c2.size(), a1.size());
    zgemm3m_("T", "N", x0.size(), a3.size()*c2.size()*a1.size(), x5.size(),
           1.0, i0data_sorted, x5.size(), i1data_sorted, x5.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<2,1,0,3,1,1,1,1>(odata_sorted, odata, x0.size(), a3.size(), c2.size(), a1.size());
  out()->add_block(odata, c2, a3, x0, a1);
}

void Task674::Task_local::compute() {
  const Index x5 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I1301
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x5, a3, c2, a1)]);
  std::fill_n(odata.get(), out()->get_size(x5, a3, c2, a1), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, a3, c2, a1);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x5.size(), a3.size(), c2.size(), a1.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(x5, a1, c2, a3);
    sort_indices<0,3,2,1,1,1,-1,2>(i1data, odata, x5.size(), a1.size(), c2.size(), a3.size());
  }
  out()->add_block(odata, x5, a3, c2, a1);
}

void Task675::Task_local::compute() {
  const Index x1 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x1, a2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(x1, a2, x0, a1), 0.0);
  {
    // tensor label: I173
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a1, x0, x1, a2);
    sort_indices<2,3,1,0,1,1,1,1>(i0data, odata, a1.size(), x0.size(), x1.size(), a2.size());
  }
  {
    // tensor label: I173
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a2, x1, x0, a1);
    sort_indices<1,0,2,3,1,1,1,1>(i0data, odata, a2.size(), x1.size(), x0.size(), a1.size());
  }
  out()->add_block(odata, x1, a2, x0, a1);
}

void Task676::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);
  // tensor label: I173
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x2, a2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x2, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), a2.size());
    // tensor label: I174
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a1, x0, x2, x1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a1, x0, x2, x1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a1.size(), x0.size(), x2.size(), x1.size());
    zgemm3m_("T", "N", a2.size(), a1.size()*x0.size()*x1.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), a1.size(), x0.size(), x1.size());
  out()->add_block(odata, a1, x0, x1, a2);
}

void Task677::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I174
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma32
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x5, a1, x4, x3);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x5, a1, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x5.size(), a1.size(), x4.size(), x3.size());
        zgemm3m_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task678::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);
  // tensor label: I173
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma33
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: I177
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x2, x3, a1, a2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x2, x3, a1, a2)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), x3.size(), a1.size(), a2.size());
      zgemm3m_("T", "N", x0.size()*x1.size(), a1.size()*a2.size(), x2.size()*x3.size(),
             1.0, i0data_sorted, x2.size()*x3.size(), i1data_sorted, x2.size()*x3.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), a2.size());
  out()->add_block(odata, a1, x0, x1, a2);
}

void Task679::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);
  // tensor label: I177
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x2, x3, a1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x2, x3, a1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x2, x3, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, a1, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, a1, c3, a2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, a1, c3, a2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c3.size(), a2.size());
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x2, c3);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x2, c3)]);
    sort_indices<1,0,0,1,-1,1>(i1data, i1data_sorted, x2.size(), c3.size());
    zgemm3m_("T", "N", x3.size()*a1.size()*a2.size(), x2.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, x3.size()*a1.size()*a2.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size(), a2.size(), x2.size());
  out()->add_block(odata, x2, x3, a1, a2);
}

void Task680::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);
  // tensor label: I177
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x2, x3, a1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x2, x3, a1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x2, x3, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, a1, a2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, a1, x2, a3);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, a1, x2, a3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a3.size());
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a3, a2);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a3, a2)]);
    sort_indices<0,1,0,1,2,1>(i1data, i1data_sorted, a3.size(), a2.size());
    zgemm3m_("T", "N", x3.size()*a1.size()*x2.size(), a2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x3.size()*a1.size()*x2.size());
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size(), x2.size(), a2.size());
  out()->add_block(odata, x2, x3, a1, a2);
}

void Task681::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);
  // tensor label: I177
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x2, x3, a1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x2, x3, a1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x2, x3, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, a1, a2), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, a1, c4, a3);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, a1, c4, a3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c4.size(), a3.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x2, c4, a3, a2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x2, c4, a3, a2)]);
      sort_indices<1,2,0,3,0,1,-1,1>(i1data, i1data_sorted, x2.size(), c4.size(), a3.size(), a2.size());
      zgemm3m_("T", "N", x3.size()*a1.size(), x2.size()*a2.size(), c4.size()*a3.size(),
             1.0, i0data_sorted, c4.size()*a3.size(), i1data_sorted, c4.size()*a3.size(),
             1.0, odata_sorted, x3.size()*a1.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size(), x2.size(), a2.size());
  out()->add_block(odata, x2, x3, a1, a2);
}

void Task682::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);
  // tensor label: I177
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x2, x3, a1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x2, x3, a1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x2, x3, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, a1, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, a4, c3, a1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, a4, c3, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a4.size(), c3.size(), a1.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x2, a2, a4, c3);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x2, a2, a4, c3)]);
      sort_indices<2,3,0,1,0,1,-1,1>(i1data, i1data_sorted, x2.size(), a2.size(), a4.size(), c3.size());
      zgemm3m_("T", "N", x3.size()*a1.size(), x2.size()*a2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, x3.size()*a1.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size(), x2.size(), a2.size());
  out()->add_block(odata, x2, x3, a1, a2);
}

void Task683::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);
  // tensor label: I177
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x2, x3, a1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x2, x3, a1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x2, x3, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, a1, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, a1, c3, a4);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, a1, c3, a4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c3.size(), a4.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x2, a2, a4, c3);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x2, a2, a4, c3)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, x2.size(), a2.size(), a4.size(), c3.size());
      zgemm3m_("T", "N", x3.size()*a1.size(), x2.size()*a2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, x3.size()*a1.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size(), x2.size(), a2.size());
  out()->add_block(odata, x2, x3, a1, a2);
}

void Task684::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);
  // tensor label: I173
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& x5 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a1, x5, x4);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c3, a1, x5, x4)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a1.size(), x5.size(), x4.size());
        // tensor label: I1217
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a2, c3, x5, x4, x1, x0);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a2, c3, x5, x4, x1, x0)]);
        sort_indices<1,2,3,0,4,5,0,1,1,1>(i1data, i1data_sorted, a2.size(), c3.size(), x5.size(), x4.size(), x1.size(), x0.size());
        zgemm3m_("T", "N", a1.size(), a2.size()*x1.size()*x0.size(), c3.size()*x5.size()*x4.size(),
               1.0, i0data_sorted, c3.size()*x5.size()*x4.size(), i1data_sorted, c3.size()*x5.size()*x4.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x1.size(), x0.size());
  out()->add_block(odata, a1, x0, x1, a2);
}

void Task685::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I1217
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, c3, x5, x4, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a2, c3, x5, x4, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, c3, x5, x4, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x5, x4, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma396
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x4, x3, x1, x2, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x4, x3, x1, x2, x0)]);
      sort_indices<2,4,0,1,3,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x1.size(), x2.size(), x0.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, a2, x2, c3);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, a2, x2, c3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a2.size(), x2.size(), c3.size());
      zgemm3m_("T", "N", x5.size()*x4.size()*x1.size()*x0.size(), a2.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x5.size()*x4.size()*x1.size()*x0.size());
    }
  }
  sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x1.size(), x0.size(), a2.size(), c3.size());
  out()->add_block(odata, a2, c3, x5, x4, x1, x0);
}

void Task686::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);
  // tensor label: I173
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, a1, x6, x5);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x7, a1, x6, x5)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x7.size(), a1.size(), x6.size(), x5.size());
        // tensor label: I1220
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a2, x7, x0, x6, x5, x1);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a2, x7, x0, x6, x5, x1)]);
        sort_indices<1,3,4,0,2,5,0,1,1,1>(i1data, i1data_sorted, a2.size(), x7.size(), x0.size(), x6.size(), x5.size(), x1.size());
        zgemm3m_("T", "N", a1.size(), a2.size()*x0.size()*x1.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
  out()->add_block(odata, a1, x0, x1, a2);
}

void Task687::Task_local::compute() {
  const Index a2 = b(0);
  const Index x7 = b(1);
  const Index x0 = b(2);
  const Index x6 = b(3);
  const Index x5 = b(4);
  const Index x1 = b(5);
  // tensor label: I1220
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, x7, x0, x6, x5, x1)]);
  std::fill_n(odata.get(), out()->get_size(a2, x7, x0, x6, x5, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, x7, x0, x6, x5, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x7, x0, x6, x5, x1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma397
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x6, x5, x4, x1, x3, x2);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x7, x0, x6, x5, x4, x1, x3, x2)]);
        sort_indices<4,6,7,0,1,2,3,5,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x6.size(), x5.size(), x4.size(), x1.size(), x3.size(), x2.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x4, a2, x3, x2);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x4, a2, x3, x2)]);
        sort_indices<0,2,3,1,0,1,1,2>(i1data, i1data_sorted, x4.size(), a2.size(), x3.size(), x2.size());
        zgemm3m_("T", "N", x7.size()*x0.size()*x6.size()*x5.size()*x1.size(), a2.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x7.size()*x0.size()*x6.size()*x5.size()*x1.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x7.size(), x0.size(), x6.size(), x5.size(), x1.size(), a2.size());
  out()->add_block(odata, a2, x7, x0, x6, x5, x1);
}

void Task688::Task_local::compute() {
  const Index a2 = b(0);
  const Index x7 = b(1);
  const Index x0 = b(2);
  const Index x6 = b(3);
  const Index x5 = b(4);
  const Index x1 = b(5);
  // tensor label: I1220
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, x7, x0, x6, x5, x1)]);
  std::fill_n(odata.get(), out()->get_size(a2, x7, x0, x6, x5, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, x7, x0, x6, x5, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x7, x0, x6, x5, x1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma236
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x6, x5, x4, x3, x2, x1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x7, x0, x6, x5, x4, x3, x2, x1)]);
        sort_indices<4,5,6,0,1,2,3,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x6.size(), x5.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x4, x3, x2, a2);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x4, x3, x2, a2)]);
        sort_indices<0,1,2,3,0,1,1,2>(i1data, i1data_sorted, x4.size(), x3.size(), x2.size(), a2.size());
        zgemm3m_("T", "N", x7.size()*x0.size()*x6.size()*x5.size()*x1.size(), a2.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x7.size()*x0.size()*x6.size()*x5.size()*x1.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x7.size(), x0.size(), x6.size(), x5.size(), x1.size(), a2.size());
  out()->add_block(odata, a2, x7, x0, x6, x5, x1);
}

void Task689::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);
  // tensor label: I173
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x2, a1, a3, a2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x2, a1, a3, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x2.size(), a1.size(), a3.size(), a2.size());
      // tensor label: I1226
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a3, x1, x2, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a3, x1, x2, x0)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), x1.size(), x2.size(), x0.size());
      zgemm3m_("T", "N", a1.size()*a2.size(), x1.size()*x0.size(), a3.size()*x2.size(),
             1.0, i0data_sorted, a3.size()*x2.size(), i1data_sorted, a3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }
  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x1.size(), x0.size());
  out()->add_block(odata, a1, x0, x1, a2);
}

void Task690::Task_local::compute() {
  const Index a3 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index x0 = b(3);
  // tensor label: I1226
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a3, x1, x2, x0)]);
  std::fill_n(odata.get(), out()->get_size(a3, x1, x2, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a3, x1, x2, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x1, x2, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma336
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x1, x4, x3, x2, x0);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x1, x4, x3, x2, x0)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x1.size(), x4.size(), x3.size(), x2.size(), x0.size());
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x5, a3, x4, x3);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x5, a3, x4, x3)]);
        sort_indices<0,2,3,1,0,1,-1,1>(i1data, i1data_sorted, x5.size(), a3.size(), x4.size(), x3.size());
        zgemm3m_("T", "N", x1.size()*x2.size()*x0.size(), a3.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x1.size()*x2.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x2.size(), x0.size(), a3.size());
  out()->add_block(odata, a3, x1, x2, x0);
}

void Task691::Task_local::compute() {
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
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, a1, c3, a2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, a1, c3, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c3.size(), a2.size());
      // tensor label: I1232
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, x5, x0, x1);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, x5, x0, x1)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), x5.size(), x0.size(), x1.size());
      zgemm3m_("T", "N", a1.size()*a2.size(), x0.size()*x1.size(), c3.size()*x5.size(),
             1.0, i0data_sorted, c3.size()*x5.size(), i1data_sorted, c3.size()*x5.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
  out()->add_block(odata, a1, x0, x1, a2);
}

void Task692::Task_local::compute() {
  const Index c3 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1232
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, x5, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c3, x5, x0, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, x5, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x5, x0, x1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma391
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x4, x1, x3, x2);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x0, x4, x1, x3, x2)]);
        sort_indices<2,4,5,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x1.size(), x3.size(), x2.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x4, c3, x3, x2);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x4, c3, x3, x2)]);
        sort_indices<0,2,3,1,0,1,-1,2>(i1data, i1data_sorted, x4.size(), c3.size(), x3.size(), x2.size());
        zgemm3m_("T", "N", x5.size()*x0.size()*x1.size(), c3.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x5.size()*x0.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x0.size(), x1.size(), c3.size());
  out()->add_block(odata, c3, x5, x0, x1);
}

void Task693::Task_local::compute() {
  const Index c3 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1232
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, x5, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c3, x5, x0, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, x5, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x5, x0, x1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma32
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<2,3,4,0,1,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x4, x3, x2, c3);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x4, x3, x2, c3)]);
        sort_indices<0,1,2,3,0,1,-1,2>(i1data, i1data_sorted, x4.size(), x3.size(), x2.size(), c3.size());
        zgemm3m_("T", "N", x5.size()*x0.size()*x1.size(), c3.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x5.size()*x0.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x0.size(), x1.size(), c3.size());
  out()->add_block(odata, c3, x5, x0, x1);
}

void Task694::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);
  // tensor label: I173
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma368
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x1, x2, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x1, x2, x0)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x2.size(), x0.size());
      // tensor label: I1238
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x2, a2, x3, a1);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x2, a2, x3, a1)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), a2.size(), x3.size(), a1.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), a2.size()*a1.size(), x2.size()*x3.size(),
             1.0, i0data_sorted, x2.size()*x3.size(), i1data_sorted, x2.size()*x3.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a2.size(), a1.size());
  out()->add_block(odata, a1, x0, x1, a2);
}

void Task695::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index a1 = b(3);
  // tensor label: I1238
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x2, a2, x3, a1)]);
  std::fill_n(odata.get(), out()->get_size(x2, a2, x3, a1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x2, a2, x3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, x3, a1), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, a3, c4, a1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, a3, c4, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c4.size(), a1.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x2, c4, a3, a2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x2, c4, a3, a2)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), c4.size(), a3.size(), a2.size());
      zgemm3m_("T", "N", x3.size()*a1.size(), x2.size()*a2.size(), c4.size()*a3.size(),
             1.0, i0data_sorted, c4.size()*a3.size(), i1data_sorted, c4.size()*a3.size(),
             1.0, odata_sorted, x3.size()*a1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size(), x2.size(), a2.size());
  out()->add_block(odata, x2, a2, x3, a1);
}

void Task696::Task_local::compute() {
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
      for (auto& x3 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: Gamma391
          std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x4, x1, x3, x2);
          std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x0, x4, x1, x3, x2)]);
          sort_indices<0,2,4,5,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x1.size(), x3.size(), x2.size());
          // tensor label: I1250
          std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a2, x3, x2, x5, a1, x4);
          std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a2, x3, x2, x5, a1, x4)]);
          sort_indices<3,5,1,2,0,4,0,1,1,1>(i1data, i1data_sorted, a2.size(), x3.size(), x2.size(), x5.size(), a1.size(), x4.size());
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

void Task697::Task_local::compute() {
  const Index a2 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x5 = b(3);
  const Index a1 = b(4);
  const Index x4 = b(5);
  // tensor label: I1250
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, x3, x2, x5, a1, x4)]);
  std::fill_n(odata.get(), out()->get_size(a2, x3, x2, x5, a1, x4), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, x3, x2, x5, a1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x3, x2, x5, a1, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, a1, x4, a3);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, a1, x4, a3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a3.size());
    // tensor label: I1251
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a3, a2, x3, x2);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a3, a2, x3, x2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), a2.size(), x3.size(), x2.size());
    zgemm3m_("T", "N", x5.size()*a1.size()*x4.size(), a2.size()*x3.size()*x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x5.size()*a1.size()*x4.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), a2.size(), x3.size(), x2.size());
  out()->add_block(odata, a2, x3, x2, x5, a1, x4);
}

void Task698::Task_local::compute() {
  const Index a3 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1251
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a3, a2, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(a3, a2, x3, x2), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a3, a2, x3, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, a3.size(), a2.size(), x3.size(), x2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(x3, x2, a3, a2);
    sort_indices<2,3,0,1,1,1,1,1>(i1data, odata, x3.size(), x2.size(), a3.size(), a2.size());
  }
  out()->add_block(odata, a3, a2, x3, x2);
}

void Task699::Task_local::compute() {
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
      for (auto& x3 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: Gamma32
          std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
          std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x0, x4, x3, x2, x1)]);
          sort_indices<0,2,3,4,1,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
          // tensor label: I1253
          std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, x2, a2, x5, a1, x4);
          std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, x2, a2, x5, a1, x4)]);
          sort_indices<3,5,0,1,2,4,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), a2.size(), x5.size(), a1.size(), x4.size());
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

#endif
