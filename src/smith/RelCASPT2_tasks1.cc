//
// BAGEL - Parallel electron correlation program.
// Filename: RelCASPT2_tasks1.cc
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

#include <src/smith/RelCASPT2_tasks1.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelCASPT2;

void Task0::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: Gamma0
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x3, x0);
  // associated with merged
  std::unique_ptr<std::complex<double>[]> fdata = in(1)->get_block(x2, x1);
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            odata[i3+x3.size()*(i0)]
              += (1.0) * i0data[i3+x3.size()*(i0+x0.size()*(i2+x2.size()*(i1)))] * fdata[i2+x2.size()*(i1)];
          }
        }
      }
    }
  }
  out()->put_block(odata, x3, x0);
}

void Task1::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: Gamma2
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x1, x0);
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task3::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c2, a3, x0, a1);
  {
    // tensor label: I0
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a3, c2, a1, x0);
    sort_indices<1,0,3,2,1,1,1,1>(i0data, odata, a3.size(), c2.size(), a1.size(), x0.size());
  }
  out()->put_block(odata, c2, a3, x0, a1);
}

void Task4::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I0
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(a3, c2, a1, x0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a3, c2, a1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, a1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: Gamma0
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size());
    // tensor label: I1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, a3, c2, a1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a3.size(), c2.size(), a1.size());
    zgemm3m_("T", "N", x0.size(), a3.size()*c2.size()*a1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x0.size(), a3.size(), c2.size(), a1.size());
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task5::Task_local::compute() {
  const Index x3 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I1
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x3, a3, c2, a1);
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, a3, c2, a1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), a3.size(), c2.size(), a1.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(x3, a1, c2, a3);
    sort_indices<0,3,2,1,1,1,1,1>(i1data, odata, x3.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, x3, a3, c2, a1);
}

void Task6::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I0
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(a3, c2, a1, x0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a3, c2, a1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, a1, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I5
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, x1, a3, a1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, x1, a3, a1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), x1.size(), a3.size(), a1.size());
    zgemm3m_("T", "N", x0.size(), c2.size()*a3.size()*a1.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<2,1,3,0,1,1,1,1>(odata_sorted, odata, x0.size(), c2.size(), a3.size(), a1.size());
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task7::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I5
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c2, x1, a3, a1);
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, a3, c2, a1);
    zscal_(x1.size()*a3.size()*c2.size()*a1.size(), e0_, i0data.get(), 1);
    sort_indices<2,0,1,3,1,1,1,1>(i0data, odata, x1.size(), a3.size(), c2.size(), a1.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(x1, a1, c2, a3);
    zscal_(x1.size()*a1.size()*c2.size()*a3.size(), e0_, i1data.get(), 1);
    sort_indices<2,0,3,1,1,1,-1,1>(i1data, odata, x1.size(), a1.size(), c2.size(), a3.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i2data = in(1)->get_block(x1, a3, c2, a1);
    sort_indices<2,0,1,3,1,1,-1,1>(i2data, odata, x1.size(), a3.size(), c2.size(), a1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i3data = in(1)->get_block(x1, a1, c2, a3);
    sort_indices<2,0,3,1,1,1,1,1>(i3data, odata, x1.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task8::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I5
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, a3, c4, a1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, a3, c4, a1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c4.size(), a1.size());
    // tensor label: f1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, c4);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, c4)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c2.size(), c4.size());
    zgemm3m_("T", "N", x1.size()*a3.size()*a1.size(), c2.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, x1.size()*a3.size()*a1.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), a3.size(), a1.size(), c2.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task9::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I5
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, a1, c4, a3);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, a1, c4, a3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c4.size(), a3.size());
    // tensor label: f1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, c4);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, c4)]);
    sort_indices<1,0,0,1,-1,1>(i1data, i1data_sorted, c2.size(), c4.size());
    zgemm3m_("T", "N", x1.size()*a1.size()*a3.size(), c2.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, x1.size()*a1.size()*a3.size());
  }
  sort_indices<3,0,2,1,1,1,1,1>(odata_sorted, odata, x1.size(), a1.size(), a3.size(), c2.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task10::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I5
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, a4, c2, a3);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, a4, c2, a3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a3.size());
    // tensor label: f1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, a1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, a1)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a1.size());
    zgemm3m_("T", "N", x1.size()*c2.size()*a3.size(), a1.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, x1.size()*c2.size()*a3.size());
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, x1.size(), c2.size(), a3.size(), a1.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task11::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I5
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, a3, c2, a4);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, a3, c2, a4)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a4.size());
    // tensor label: f1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, a1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, a1)]);
    sort_indices<0,1,0,1,-1,1>(i1data, i1data_sorted, a4.size(), a1.size());
    zgemm3m_("T", "N", x1.size()*a3.size()*c2.size(), a1.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, x1.size()*a3.size()*c2.size());
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x1.size(), a3.size(), c2.size(), a1.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task12::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I5
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, a4, c2, a1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, a4, c2, a1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a1.size());
    // tensor label: f1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, a3);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, a3)]);
    sort_indices<0,1,0,1,-1,1>(i1data, i1data_sorted, a4.size(), a3.size());
    zgemm3m_("T", "N", x1.size()*c2.size()*a1.size(), a3.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, x1.size()*c2.size()*a1.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, x1.size(), c2.size(), a1.size(), a3.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task13::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I5
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, a1, c2, a4);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, a1, c2, a4)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a4.size());
    // tensor label: f1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, a3);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, a3)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a3.size());
    zgemm3m_("T", "N", x1.size()*a1.size()*c2.size(), a3.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, x1.size()*a1.size()*c2.size());
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, x1.size(), a1.size(), c2.size(), a3.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task14::Task_local::compute() {
  target_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: Gamma2
  std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
  std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
  // tensor label: I31
  std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x0, x1);
  std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x0, x1)]);
  sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size());
  target_ += zdotu_(x0.size()*x1.size(), i0data_sorted, 1, i1data_sorted, 1);
}

void Task15::Task_local::compute() {
  target_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I31
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x0, x1);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a1 : *range_[2]) {
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, a3, c2, a1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, a3, c2, a1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a1.size());
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x0, a1, c2, a3);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x0, a1, c2, a3)]);
        sort_indices<3,2,1,0,0,1,-1,1>(i1data, i1data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
        zgemm3m_("T", "N", x1.size(), x0.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x0, x1);
}

void Task16::Task_local::compute() {
  target_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I31
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x0, x1);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a1 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, a1, c2, a3);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, a1, c2, a3)]);
        sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x1, a1, c2, a3);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x1, a1, c2, a3)]);
        sort_indices<3,2,1,0,0,1,1,1>(i1data, i1data_sorted, x1.size(), a1.size(), c2.size(), a3.size());
        zgemm3m_("T", "N", x0.size(), x1.size(), a1.size()*c2.size()*a3.size(),
               1.0, i0data_sorted, a1.size()*c2.size()*a3.size(), i1data_sorted, a1.size()*c2.size()*a3.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->put_block(odata, x0, x1);
}

#endif
