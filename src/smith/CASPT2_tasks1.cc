//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks1.cc
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


#include <src/smith/CASPT2_tasks1.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task0::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: Gamma0
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task1::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: Gamma4
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(1)->get_block(x2, x1);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
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

void Task2::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  // tensor label: Gamma26
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}

void Task3::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: Gamma36
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(1)->get_block(x1, x0);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
          odata[ici0]
            += (1.0) * i0data[ici0+ci0.size()*(ix1+x1.size()*(ix0))] * fdata[ix1+x1.size()*(ix0)];
        }
      }
    }
  }
  out()->put_block(odata, ci0);
}

void Task4::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: Gamma38
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task5::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  // tensor label: Gamma42
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(1)->get_block(x2, x1);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
              odata[ici0+ci0.size()*(ix3+x3.size()*(ix0))]
                += (1.0) * i0data[ici0+ci0.size()*(ix3+x3.size()*(ix0+x0.size()*(ix2+x2.size()*(ix1))))] * fdata[ix2+x2.size()*(ix1)];
            }
          }
        }
      }
    }
  }
  out()->put_block(odata, ci0, x3, x0);
}

void Task7::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c3, a4, c1, a2);
  {
    // tensor label: I0
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, c1, a2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c3.size(), a4.size(), c1.size(), a2.size());
  }
  {
    // tensor label: I0
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<2,3,0,1,1,1,1,1>(i0data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->put_block(odata, c3, a4, c1, a2);
}

void Task8::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I0
  std::unique_ptr<double[]> odata = out()->move_block(c3, a4, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a4, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a4, c1, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c1, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c1.size(), a2.size());
    // tensor label: I1
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), x1.size());
    dgemm_("T", "N", a4.size()*c1.size()*a2.size(), c3.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a4.size()*c1.size()*a2.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), a2.size(), c3.size());
  out()->put_block(odata, c3, a4, c1, a2);
}

void Task9::Task_local::compute() {
  const Index c3 = b(0);
  const Index x1 = b(1);
  // tensor label: I1
  std::unique_ptr<double[]> odata = out()->move_block(c3, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma0
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x0)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), x0.size());
    dgemm_("T", "N", x1.size(), c3.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), c3.size());
  out()->put_block(odata, c3, x1);
}

void Task10::Task_local::compute() {
  const Index c3 = b(0);
  const Index x0 = b(1);
  // tensor label: I2
  std::unique_ptr<double[]> odata = out()->move_block(c3, x0);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x0);
    sort_indices<0,1,1,1,-2,1>(i0data, odata, c3.size(), x0.size());
  }
  out()->put_block(odata, c3, x0);
}

void Task11::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I0
  std::unique_ptr<double[]> odata = out()->move_block(c3, a4, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a4, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a4, c1, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, c1, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2, c1, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c1.size(), a4.size());
    // tensor label: I4
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), x1.size());
    dgemm_("T", "N", a2.size()*c1.size()*a4.size(), c3.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a2.size()*c1.size()*a4.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), a4.size(), c3.size());
  out()->put_block(odata, c3, a4, c1, a2);
}

void Task12::Task_local::compute() {
  const Index c3 = b(0);
  const Index x1 = b(1);
  // tensor label: I4
  std::unique_ptr<double[]> odata = out()->move_block(c3, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma0
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I5
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x0)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), x0.size());
    dgemm_("T", "N", x1.size(), c3.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), c3.size());
  out()->put_block(odata, c3, x1);
}

void Task13::Task_local::compute() {
  const Index c3 = b(0);
  const Index x0 = b(1);
  // tensor label: I5
  std::unique_ptr<double[]> odata = out()->move_block(c3, x0);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x0);
    sort_indices<0,1,1,1,1,1>(i0data, odata, c3.size(), x0.size());
  }
  out()->put_block(odata, c3, x0);
}

void Task14::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c2, a3, x0, a1);
  {
    // tensor label: I6
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, c2, a3, a1);
    sort_indices<1,2,0,3,1,1,1,1>(i0data, odata, x0.size(), c2.size(), a3.size(), a1.size());
  }
  out()->put_block(odata, c2, a3, x0, a1);
}

void Task15::Task_local::compute() {
  const Index x0 = b(0);
  const Index c2 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I6
  std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c4, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c4, a1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c4.size(), a1.size());
    // tensor label: I7
    std::unique_ptr<double[]> i1data = in(1)->get_block(c4, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c4, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c4.size(), x0.size());
    dgemm_("T", "N", c2.size()*a3.size()*a1.size(), x0.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c2.size()*a3.size()*a1.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), a1.size(), x0.size());
  out()->put_block(odata, x0, c2, a3, a1);
}

void Task16::Task_local::compute() {
  const Index c4 = b(0);
  const Index x0 = b(1);
  // tensor label: I7
  std::unique_ptr<double[]> odata = out()->move_block(c4, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c4, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma0
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I8
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c4)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), c4.size());
    dgemm_("T", "N", x0.size(), c4.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c4.size());
  out()->put_block(odata, c4, x0);
}

void Task17::Task_local::compute() {
  const Index x1 = b(0);
  const Index c4 = b(1);
  // tensor label: I8
  std::unique_ptr<double[]> odata = out()->move_block(x1, c4);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c4);
    sort_indices<0,1,1,1,-4,1>(i0data, odata, x1.size(), c4.size());
  }
  out()->put_block(odata, x1, c4);
}

void Task18::Task_local::compute() {
  const Index x0 = b(0);
  const Index c2 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I6
  std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c4, a3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c4.size(), a3.size());
    // tensor label: I10
    std::unique_ptr<double[]> i1data = in(1)->get_block(c4, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c4, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c4.size(), x0.size());
    dgemm_("T", "N", c2.size()*a1.size()*a3.size(), x0.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c2.size()*a1.size()*a3.size());
  }
  sort_indices<3,0,2,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), a3.size(), x0.size());
  out()->put_block(odata, x0, c2, a3, a1);
}

void Task19::Task_local::compute() {
  const Index c4 = b(0);
  const Index x0 = b(1);
  // tensor label: I10
  std::unique_ptr<double[]> odata = out()->move_block(c4, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c4, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma0
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I11
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c4)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), c4.size());
    dgemm_("T", "N", x0.size(), c4.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c4.size());
  out()->put_block(odata, c4, x0);
}

void Task20::Task_local::compute() {
  const Index x1 = b(0);
  const Index c4 = b(1);
  // tensor label: I11
  std::unique_ptr<double[]> odata = out()->move_block(x1, c4);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c4);
    sort_indices<0,1,1,1,2,1>(i0data, odata, x1.size(), c4.size());
  }
  out()->put_block(odata, x1, c4);
}

void Task21::Task_local::compute() {
  const Index x0 = b(0);
  const Index c2 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I6
  std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: Gamma4
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size());
    // tensor label: I13
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a3, c2, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a3.size(), c2.size(), a1.size());
    dgemm_("T", "N", x0.size(), a3.size()*c2.size()*a1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<0,2,1,3,1,1,1,1>(odata_sorted, odata, x0.size(), a3.size(), c2.size(), a1.size());
  out()->put_block(odata, x0, c2, a3, a1);
}

void Task22::Task_local::compute() {
  const Index x3 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I13
  std::unique_ptr<double[]> odata = out()->move_block(x3, a3, c2, a1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), a3.size(), c2.size(), a1.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x3, a1, c2, a3);
    sort_indices<0,3,2,1,1,1,2,1>(i1data, odata, x3.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, x3, a3, c2, a1);
}

void Task23::Task_local::compute() {
  const Index x0 = b(0);
  const Index c2 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I6
  std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma0
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I17
    std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x1, a3, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x1, a3, a1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), x1.size(), a3.size(), a1.size());
    dgemm_("T", "N", x0.size(), c2.size()*a3.size()*a1.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x0.size(), c2.size(), a3.size(), a1.size());
  out()->put_block(odata, x0, c2, a3, a1);
}

void Task24::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I17
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a1);
    dscal_(x1.size()*a3.size()*c2.size()*a1.size(), e0_, i0data.get(), 1);
    sort_indices<2,0,1,3,1,1,1,1>(i0data, odata, x1.size(), a3.size(), c2.size(), a1.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x1, a1, c2, a3);
    dscal_(x1.size()*a1.size()*c2.size()*a3.size(), e0_, i1data.get(), 1);
    sort_indices<2,0,3,1,1,1,-2,1>(i1data, odata, x1.size(), a1.size(), c2.size(), a3.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i2data = in(1)->get_block(x1, a3, c2, a1);
    sort_indices<2,0,1,3,1,1,-2,1>(i2data, odata, x1.size(), a3.size(), c2.size(), a1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i3data = in(1)->get_block(x1, a1, c2, a3);
    sort_indices<2,0,3,1,1,1,4,1>(i3data, odata, x1.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task25::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I17
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c4, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c4, a1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c4.size(), a1.size());
    // tensor label: I18
    std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c4)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c2.size(), c4.size());
    dgemm_("T", "N", x1.size()*a3.size()*a1.size(), c2.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, x1.size()*a3.size()*a1.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), a3.size(), a1.size(), c2.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task26::Task_local::compute() {
  const Index c2 = b(0);
  const Index c4 = b(1);
  // tensor label: I18
  std::unique_ptr<double[]> odata = out()->move_block(c2, c4);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c4);
    sort_indices<0,1,1,1,1,1>(i0data, odata, c2.size(), c4.size());
  }
  out()->put_block(odata, c2, c4);
}

void Task27::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I17
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c4, a3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c4.size(), a3.size());
    // tensor label: I21
    std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c4)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c2.size(), c4.size());
    dgemm_("T", "N", x1.size()*a1.size()*a3.size(), c2.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, x1.size()*a1.size()*a3.size());
  }
  sort_indices<3,0,2,1,1,1,1,1>(odata_sorted, odata, x1.size(), a1.size(), a3.size(), c2.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task28::Task_local::compute() {
  const Index c2 = b(0);
  const Index c4 = b(1);
  // tensor label: I21
  std::unique_ptr<double[]> odata = out()->move_block(c2, c4);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c4);
    sort_indices<0,1,1,1,-2,1>(i0data, odata, c2.size(), c4.size());
  }
  out()->put_block(odata, c2, c4);
}

void Task29::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I17
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a3.size());
    // tensor label: I24
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a1)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a1.size());
    dgemm_("T", "N", x1.size()*c2.size()*a3.size(), a1.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, x1.size()*c2.size()*a3.size());
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, x1.size(), c2.size(), a3.size(), a1.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task30::Task_local::compute() {
  const Index a4 = b(0);
  const Index a1 = b(1);
  // tensor label: I24
  std::unique_ptr<double[]> odata = out()->move_block(a4, a1);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1);
    sort_indices<0,1,1,1,2,1>(i0data, odata, a4.size(), a1.size());
  }
  out()->put_block(odata, a4, a1);
}

void Task31::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I17
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a4)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a4.size());
    // tensor label: I27
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a1)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a1.size());
    dgemm_("T", "N", x1.size()*a3.size()*c2.size(), a1.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, x1.size()*a3.size()*c2.size());
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x1.size(), a3.size(), c2.size(), a1.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task32::Task_local::compute() {
  const Index a4 = b(0);
  const Index a1 = b(1);
  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(a4, a1);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, a4.size(), a1.size());
  }
  out()->put_block(odata, a4, a1);
}

void Task33::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I17
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a1.size());
    // tensor label: I30
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a3)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a3.size());
    dgemm_("T", "N", x1.size()*c2.size()*a1.size(), a3.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, x1.size()*c2.size()*a1.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, x1.size(), c2.size(), a1.size(), a3.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task34::Task_local::compute() {
  const Index a4 = b(0);
  const Index a3 = b(1);
  // tensor label: I30
  std::unique_ptr<double[]> odata = out()->move_block(a4, a3);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a3);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, a4.size(), a3.size());
  }
  out()->put_block(odata, a4, a3);
}

void Task35::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I17
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a4)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a4.size());
    // tensor label: I33
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a3)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a3.size());
    dgemm_("T", "N", x1.size()*a1.size()*c2.size(), a3.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, x1.size()*a1.size()*c2.size());
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, x1.size(), a1.size(), c2.size(), a3.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task36::Task_local::compute() {
  const Index a4 = b(0);
  const Index a3 = b(1);
  // tensor label: I33
  std::unique_ptr<double[]> odata = out()->move_block(a4, a3);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a3);
    sort_indices<0,1,1,1,2,1>(i0data, odata, a4.size(), a3.size());
  }
  out()->put_block(odata, a4, a3);
}

void Task37::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c3, a4, c1, a2);
  {
    // tensor label: I38
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
    sort_indices<2,1,0,3,1,1,1,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
  }
  out()->put_block(odata, c3, a4, c1, a2);
}

void Task38::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I38
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
    sort_indices<0,1,2,3,1,1,-4,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<0,3,2,1,1,1,8,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->put_block(odata, c1, a4, c3, a2);
}

void Task39::Task_local::compute() {
  target_ = 0.0;
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index a4 = b(3);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
  sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
  // tensor label: I45
  std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a4, c3, a2)]);
  sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
  target_ += ddot_(c1.size()*a4.size()*c3.size()*a2.size(), i0data_sorted, 1, i1data_sorted, 1);
}

void Task40::Task_local::compute() {
  target_ = 0.0;
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I45
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
  }
  out()->put_block(odata, c1, a4, c3, a2);
}

void Task41::Task_local::compute() {
  target_ = 0.0;
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: v2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
  // tensor label: I47
  std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3, a2, c1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3, a2, c1)]);
  sort_indices<3,2,1,0,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size(), a2.size(), c1.size());
  target_ += ddot_(a4.size()*c3.size()*a2.size()*c1.size(), i0data_sorted, 1, i1data_sorted, 1);
}

void Task42::Task_local::compute() {
  target_ = 0.0;
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I47
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<3,2,1,0,1,1,2,1>(i0data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->put_block(odata, a4, c3, a2, c1);
}

void Task43::Task_local::compute() {
  target_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: Gamma0
  std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
  // tensor label: I49
  std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1)]);
  sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size());
  target_ += ddot_(x0.size()*x1.size(), i0data_sorted, 1, i1data_sorted, 1);
}

void Task44::Task_local::compute() {
  target_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I49
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a1 : *range_[2]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a1.size());
        // tensor label: I50
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x0)]);
        sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x0.size());
        dgemm_("T", "N", x1.size(), x0.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x0, x1);
}

void Task45::Task_local::compute() {
  target_ = 0.0;
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I50
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, x0.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task46::Task_local::compute() {
  target_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I49
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a1 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
        sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
        // tensor label: I53
        std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a1, c2, a3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a1, c2, a3)]);
        sort_indices<3,2,1,0,0,1,1,1>(i1data, i1data_sorted, x1.size(), a1.size(), c2.size(), a3.size());
        dgemm_("T", "N", x0.size(), x1.size(), a1.size()*c2.size()*a3.size(),
               1.0, i0data_sorted, a1.size()*c2.size()*a3.size(), i1data_sorted, a1.size()*c2.size()*a3.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->put_block(odata, x0, x1);
}

void Task47::Task_local::compute() {
  target_ = 0.0;
  const Index x1 = b(0);
  const Index a1 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  // tensor label: I53
  std::unique_ptr<double[]> odata = out()->move_block(x1, a1, c2, a3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x1.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, x1, a1, c2, a3);
}

void Task48::Task_local::compute() {
  target_ = 0.0;
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index a4 = b(2);
  const Index c1 = b(3);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
  // tensor label: I55
  std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3, a2, c1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3, a2, c1)]);
  sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size(), a2.size(), c1.size());
  target_ += ddot_(a4.size()*c3.size()*a2.size()*c1.size(), i0data_sorted, 1, i1data_sorted, 1);
}

void Task49::Task_local::compute() {
  target_ = 0.0;
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I55
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<3,2,1,0,1,1,-1,1>(i0data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->put_block(odata, a4, c3, a2, c1);
}

