//
// BAGEL - Parallel electron correlation program.
// Filename: RelMRCI_tasks1.cc
// Copyright (C) 2014 Shiozaki group
//
// Author: Shiozaki group <shiozaki@northwestern.edu>
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

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#include <src/smith/RelMRCI_tasks1.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelMRCI;

void Task0::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: Gamma0
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x1, x0);
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task1::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: Gamma4
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x0, x1);
  {
    // rdm0 non-merged case
    if (x0 == x1) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block();
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        odata[i1+x0.size()*(i1)]  += 1.0 * i0data[0];
      }
    }
  }
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x0, x1);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, x0.size(), x1.size());
  }
  out()->put_block(odata, x0, x1);
}

void Task2::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // scalar
  // tensor label: Gamma16
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block();
  // associated with merged
  std::unique_ptr<std::complex<double>[]> fdata = in(1)->get_block(x1, x0);
  out()->put_block(odata);
}

void Task3::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index x3 = b(3);
  // scalar
  // tensor label: Gamma18
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block();
  // associated with merged
  std::unique_ptr<std::complex<double>[]> fdata = in(2)->get_block(x3, x2, x1, x0);
  out()->put_block(odata);
}

void Task5::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c3, a4, c1, a2);
  {
    // tensor label: I0
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, c1, a4, a2);
    sort_indices<0,2,1,3,1,1,1,1>(i0data, odata, c3.size(), c1.size(), a4.size(), a2.size());
  }
  {
    // tensor label: I0
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, c3, a2, a4);
    sort_indices<1,3,0,2,1,1,1,1>(i0data, odata, c1.size(), c3.size(), a2.size(), a4.size());
  }
  out()->put_block(odata, c3, a4, c1, a2);
}

void Task6::Task_local::compute() {
  const Index c3 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index a2 = b(3);
  // tensor label: I0
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c3, c1, a4, a2);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, c1, a4, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c1, a4, a2), 0.0);
  for (auto& c5 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a4, c5, a2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a4, c5, a2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c5.size(), a2.size());
    // tensor label: I1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, c5);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, c5)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), c5.size());
    zgemm3m_("T", "N", c1.size()*a4.size()*a2.size(), c3.size(), c5.size(),
           1.0, i0data_sorted, c5.size(), i1data_sorted, c5.size(),
           1.0, odata_sorted, c1.size()*a4.size()*a2.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), a2.size(), c3.size());
  out()->put_block(odata, c3, c1, a4, a2);
}

void Task7::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  // tensor label: I1
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c3, c5);
  {
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, c5);
    sort_indices<0,1,1,1,2,1>(i0data, odata, c3.size(), c5.size());
  }
  out()->put_block(odata, c3, c5);
}

void Task8::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  // tensor label: I1
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c3, c5);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, c5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c5), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, c5, x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c3, c5, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), c5.size(), x1.size(), x0.size());
      // tensor label: Gamma0
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x1, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
      zgemm3m_("T", "N", c3.size()*c5.size(), 1, x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, c3.size()*c5.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c3.size(), c5.size());
  out()->put_block(odata, c3, c5);
}

void Task9::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  // tensor label: I1
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c3, c5);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, c5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c5), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma4
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
  out()->put_block(odata, c3, c5);
}

void Task10::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  // tensor label: I1
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c3, c5);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, c5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c5), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, c5, c3, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, c5, c3, x0)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), c5.size(), c3.size(), x0.size());
      // tensor label: Gamma0
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x1, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,-1,1>(i1data, i1data_sorted, x1.size(), x0.size());
      zgemm3m_("T", "N", c5.size()*c3.size(), 1, x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, c5.size()*c3.size());
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c5.size(), c3.size());
  out()->put_block(odata, c3, c5);
}

void Task11::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  // tensor label: I1
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c3, c5);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, c5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c5), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0, c3, c5);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0, c3, c5)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size(), c3.size(), c5.size());
      // tensor label: Gamma0
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x1, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
      zgemm3m_("T", "N", c3.size()*c5.size(), 1, x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, c3.size()*c5.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c3.size(), c5.size());
  out()->put_block(odata, c3, c5);
}

void Task12::Task_local::compute() {
  const Index c3 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index a2 = b(3);
  // tensor label: I0
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c3, c1, a4, a2);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, c1, a4, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c1, a4, a2), 0.0);
  for (auto& c5 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a2, c5, a4);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a2, c5, a4)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c5.size(), a4.size());
    // tensor label: I3
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, c5);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, c5)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), c5.size());
    zgemm3m_("T", "N", c1.size()*a2.size()*a4.size(), c3.size(), c5.size(),
           1.0, i0data_sorted, c5.size(), i1data_sorted, c5.size(),
           1.0, odata_sorted, c1.size()*a2.size()*a4.size());
  }
  sort_indices<3,0,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), a4.size(), c3.size());
  out()->put_block(odata, c3, c1, a4, a2);
}

void Task13::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  // tensor label: I3
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c3, c5);
  {
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, c5);
    sort_indices<0,1,1,1,-2,1>(i0data, odata, c3.size(), c5.size());
  }
  out()->put_block(odata, c3, c5);
}

void Task14::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  // tensor label: I3
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c3, c5);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, c5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c5), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, c5, x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c3, c5, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), c5.size(), x1.size(), x0.size());
      // tensor label: Gamma0
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x1, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,-1,1>(i1data, i1data_sorted, x1.size(), x0.size());
      zgemm3m_("T", "N", c3.size()*c5.size(), 1, x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, c3.size()*c5.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c3.size(), c5.size());
  out()->put_block(odata, c3, c5);
}

void Task15::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  // tensor label: I3
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c3, c5);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, c5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c5), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma4
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
  out()->put_block(odata, c3, c5);
}

void Task16::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  // tensor label: I3
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c3, c5);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, c5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c5), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma0
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: I37
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x1, c5, c3, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x1, c5, c3, x0)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), c5.size(), c3.size(), x0.size());
      zgemm3m_("T", "N", 1, c5.size()*c3.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c5.size(), c3.size());
  out()->put_block(odata, c3, c5);
}

void Task17::Task_local::compute() {
  const Index x1 = b(0);
  const Index c5 = b(1);
  const Index c3 = b(2);
  const Index x0 = b(3);
  // tensor label: I37
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x1, c5, c3, x0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, c5, c3, x0);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x1.size(), c5.size(), c3.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(x1, x0, c3, c5);
    sort_indices<0,3,2,1,1,1,-1,1>(i1data, odata, x1.size(), x0.size(), c3.size(), c5.size());
  }
  out()->put_block(odata, x1, c5, c3, x0);
}

void Task18::Task_local::compute() {
  const Index c3 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index a2 = b(3);
  // tensor label: I0
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c3, c1, a4, a2);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, c1, a4, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c1, a4, a2), 0.0);
  for (auto& a5 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a5, c3, a2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a5, c3, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a5.size(), c3.size(), a2.size());
    // tensor label: I5
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a5, a4);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a5, a4)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a5.size(), a4.size());
    zgemm3m_("T", "N", c1.size()*c3.size()*a2.size(), a4.size(), a5.size(),
           1.0, i0data_sorted, a5.size(), i1data_sorted, a5.size(),
           1.0, odata_sorted, c1.size()*c3.size()*a2.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), a2.size(), a4.size());
  out()->put_block(odata, c3, c1, a4, a2);
}

void Task19::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  // tensor label: I5
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(a5, a4);
  {
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a5, a4);
    sort_indices<0,1,1,1,-2,1>(i0data, odata, a5.size(), a4.size());
  }
  out()->put_block(odata, a5, a4);
}

void Task20::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  // tensor label: I5
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(a5, a4);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a5, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a5, a4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a5, a4, x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(a5, a4, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, a5.size(), a4.size(), x1.size(), x0.size());
      // tensor label: Gamma0
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x1, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,-1,1>(i1data, i1data_sorted, x1.size(), x0.size());
      zgemm3m_("T", "N", a5.size()*a4.size(), 1, x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, a5.size()*a4.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a5.size(), a4.size());
  out()->put_block(odata, a5, a4);
}

void Task21::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  // tensor label: I5
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(a5, a4);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a5, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a5, a4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a5, x1, x0, a4);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(a5, x1, x0, a4)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, a5.size(), x1.size(), x0.size(), a4.size());
      // tensor label: Gamma4
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x0, x1);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x0, x1)]);
      sort_indices<1,0,0,1,-1,1>(i1data, i1data_sorted, x0.size(), x1.size());
      zgemm3m_("T", "N", a5.size()*a4.size(), 1, x0.size()*x1.size(),
             1.0, i0data_sorted, x0.size()*x1.size(), i1data_sorted, x0.size()*x1.size(),
             1.0, odata_sorted, a5.size()*a4.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a5.size(), a4.size());
  out()->put_block(odata, a5, a4);
}

void Task22::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  // tensor label: I5
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(a5, a4);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a5, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a5, a4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma0
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: I40
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x1, a4, a5, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x1, a4, a5, x0)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), a4.size(), a5.size(), x0.size());
      zgemm3m_("T", "N", 1, a4.size()*a5.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a4.size(), a5.size());
  out()->put_block(odata, a5, a4);
}

void Task23::Task_local::compute() {
  const Index x1 = b(0);
  const Index a4 = b(1);
  const Index a5 = b(2);
  const Index x0 = b(3);
  // tensor label: I40
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x1, a4, a5, x0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, a4, a5, x0);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x1.size(), a4.size(), a5.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(x1, x0, a5, a4);
    sort_indices<0,3,2,1,1,1,-1,1>(i1data, odata, x1.size(), x0.size(), a5.size(), a4.size());
  }
  out()->put_block(odata, x1, a4, a5, x0);
}

void Task24::Task_local::compute() {
  const Index c3 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index a2 = b(3);
  // tensor label: I0
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c3, c1, a4, a2);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, c1, a4, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c1, a4, a2), 0.0);
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
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), a2.size(), c3.size());
  out()->put_block(odata, c3, c1, a4, a2);
}

void Task25::Task_local::compute() {
  const Index c3 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index a2 = b(3);
  // tensor label: I0
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c3, c1, a4, a2);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, c1, a4, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c1, a4, a2), 0.0);
  for (auto& a5 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a2, c3, a5);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a2, c3, a5)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a5.size());
    // tensor label: I18
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a5, a4);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a5, a4)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a5.size(), a4.size());
    zgemm3m_("T", "N", c1.size()*a2.size()*c3.size(), a4.size(), a5.size(),
           1.0, i0data_sorted, a5.size(), i1data_sorted, a5.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c3.size());
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->put_block(odata, c3, c1, a4, a2);
}

void Task26::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  // tensor label: I18
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(a5, a4);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a5, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a5, a4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma0
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: I19
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a5, a4, x1, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a5, a4, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a5.size(), a4.size(), x1.size(), x0.size());
      zgemm3m_("T", "N", 1, a5.size()*a4.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a5.size(), a4.size());
  out()->put_block(odata, a5, a4);
}

void Task27::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I19
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(a5, a4, x1, x0);
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
  out()->put_block(odata, a5, a4, x1, x0);
}

void Task28::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  // tensor label: I18
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(a5, a4);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a5, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a5, a4), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma4
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
  out()->put_block(odata, a5, a4);
}

void Task29::Task_local::compute() {
  const Index c3 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index a2 = b(3);
  // tensor label: I0
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c3, c1, a4, a2);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, c1, a4, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c1, a4, a2), 0.0);
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
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->put_block(odata, c3, c1, a4, a2);
}

void Task30::Task_local::compute() {
  const Index c3 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index a2 = b(3);
  // tensor label: I0
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c3, c1, a4, a2);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, c1, a4, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c1, a4, a2), 0.0);
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
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->put_block(odata, c3, c1, a4, a2);
}

void Task31::Task_local::compute() {
  const Index c3 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index a2 = b(3);
  // tensor label: I0
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c3, c1, a4, a2);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, c1, a4, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c1, a4, a2), 0.0);
  for (auto& c6 : *range_[0]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, c6, a5, a4);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c3, c6, a5, a4)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), c6.size(), a5.size(), a4.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a5, c6, a2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a5, c6, a2)]);
      sort_indices<2,1,0,3,0,1,2,1>(i1data, i1data_sorted, c1.size(), a5.size(), c6.size(), a2.size());
      zgemm3m_("T", "N", c3.size()*a4.size(), c1.size()*a2.size(), a5.size()*c6.size(),
             1.0, i0data_sorted, a5.size()*c6.size(), i1data_sorted, a5.size()*c6.size(),
             1.0, odata_sorted, c3.size()*a4.size());
    }
  }
  sort_indices<0,2,1,3,1,1,1,1>(odata_sorted, odata, c3.size(), a4.size(), c1.size(), a2.size());
  out()->put_block(odata, c3, c1, a4, a2);
}

void Task32::Task_local::compute() {
  const Index c3 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index a2 = b(3);
  // tensor label: I0
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c3, c1, a4, a2);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, c1, a4, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c1, a4, a2), 0.0);
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
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->put_block(odata, c3, c1, a4, a2);
}

void Task33::Task_local::compute() {
  const Index c3 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index a2 = b(3);
  // tensor label: I0
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c3, c1, a4, a2);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, c1, a4, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c1, a4, a2), 0.0);
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
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->put_block(odata, c3, c1, a4, a2);
}

void Task34::Task_local::compute() {
  const Index c3 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index a2 = b(3);
  // tensor label: I0
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c3, c1, a4, a2);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, c1, a4, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c1, a4, a2), 0.0);
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
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->put_block(odata, c3, c1, a4, a2);
}

void Task35::Task_local::compute() {
  const Index c3 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index a2 = b(3);
  // tensor label: I0
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c3, c1, a4, a2);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, c1, a4, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c1, a4, a2), 0.0);
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
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->put_block(odata, c3, c1, a4, a2);
}

void Task36::Task_local::compute() {
  const Index c3 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index a2 = b(3);
  // tensor label: I0
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c3, c1, a4, a2);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, c1, a4, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c1, a4, a2), 0.0);
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
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->put_block(odata, c3, c1, a4, a2);
}

void Task37::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c3, a4, c1, a2);
  {
    // tensor label: I56
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, c3, a4, a2);
    sort_indices<1,2,0,3,1,1,1,1>(i0data, odata, c1.size(), c3.size(), a4.size(), a2.size());
  }
  out()->put_block(odata, c3, a4, c1, a2);
}

void Task38::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index a4 = b(2);
  const Index a2 = b(3);
  // tensor label: I56
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c1, c3, a4, a2);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c3, a4, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, a4, a2), 0.0);
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
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a4.size(), a2.size(), c1.size(), c3.size());
  out()->put_block(odata, c1, c3, a4, a2);
}

void Task39::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index a4 = b(2);
  const Index a2 = b(3);
  // tensor label: I56
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c1, c3, a4, a2);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c3, a4, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, a4, a2), 0.0);
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
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), a4.size(), c1.size(), c3.size());
  out()->put_block(odata, c1, c3, a4, a2);
}

void Task40::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index a4 = b(2);
  const Index a2 = b(3);
  // tensor label: I56
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c1, c3, a4, a2);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c3, a4, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, a4, a2), 0.0);
  for (auto& a6 : *range_[2]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a6, c3, a5);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a6, c3, a5)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a6.size(), c3.size(), a5.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a6, a2, a5, a4);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a6, a2, a5, a4)]);
      sort_indices<0,2,1,3,0,1,2,1>(i1data, i1data_sorted, a6.size(), a2.size(), a5.size(), a4.size());
      zgemm3m_("T", "N", c1.size()*c3.size(), a2.size()*a4.size(), a6.size()*a5.size(),
             1.0, i0data_sorted, a6.size()*a5.size(), i1data_sorted, a6.size()*a5.size(),
             1.0, odata_sorted, c1.size()*c3.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), a2.size(), a4.size());
  out()->put_block(odata, c1, c3, a4, a2);
}

void Task41::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index a4 = b(2);
  const Index a2 = b(3);
  // tensor label: I56
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c1, c3, a4, a2);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c3, a4, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, a4, a2), 0.0);
  for (auto& a5 : *range_[2]) {
    for (auto& a6 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a5, c3, a6);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a5, c3, a6)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a5.size(), c3.size(), a6.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a6, a2, a5, a4);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a6, a2, a5, a4)]);
      sort_indices<2,0,1,3,0,1,-2,1>(i1data, i1data_sorted, a6.size(), a2.size(), a5.size(), a4.size());
      zgemm3m_("T", "N", c1.size()*c3.size(), a2.size()*a4.size(), a6.size()*a5.size(),
             1.0, i0data_sorted, a6.size()*a5.size(), i1data_sorted, a6.size()*a5.size(),
             1.0, odata_sorted, c1.size()*c3.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), a2.size(), a4.size());
  out()->put_block(odata, c1, c3, a4, a2);
}

void Task42::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index a4 = b(2);
  const Index a2 = b(3);
  // tensor label: I56
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c1, c3, a4, a2);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c3, a4, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, a4, a2), 0.0);
  // scalar
  // tensor label: Gamma16
  std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block();
  std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size()]);
  sort_indices<0,1,1,1>(i0data, i0data_sorted);
  // tensor label: I81
  std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a4, c3, a2);
  std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a4, c3, a2)]);
  sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
  zgemm3m_("T", "N", 1, c1.size()*a4.size()*c3.size()*a2.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, 1);
  sort_indices<0,2,1,3,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->put_block(odata, c1, c3, a4, a2);
}

void Task43::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I81
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c1, a4, c3, a2);
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
  out()->put_block(odata, c1, a4, c3, a2);
}

void Task44::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index a4 = b(2);
  const Index a2 = b(3);
  // tensor label: I56
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c1, c3, a4, a2);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c3, a4, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, a4, a2), 0.0);
  // scalar
  // tensor label: Gamma18
  std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block();
  std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size()]);
  sort_indices<0,1,1,1>(i0data, i0data_sorted);
  // tensor label: I85
  std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a4, c3, a2);
  std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a4, c3, a2)]);
  sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
  zgemm3m_("T", "N", 1, c1.size()*a4.size()*c3.size()*a2.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, 1);
  sort_indices<0,2,1,3,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->put_block(odata, c1, c3, a4, a2);
}

void Task45::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I85
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c1, a4, c3, a2);
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
  out()->put_block(odata, c1, a4, c3, a2);
}

void Task47::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c3, a4, c1, a2);
  {
    // tensor label: I88
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a4, c3, a2);
    sort_indices<2,1,0,3,1,1,1,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
  }
  out()->put_block(odata, c3, a4, c1, a2);
}

void Task48::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I88
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c1, a4, c3, a2);
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
  out()->put_block(odata, c1, a4, c3, a2);
}

#endif
