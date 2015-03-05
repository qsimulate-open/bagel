//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks1.cc
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

#include <src/smith/MRCI_tasks1.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task0::Task_local::compute() {
  const Index x1 = b(0);
  const Index x4 = b(1);
  const Index x0 = b(2);
  const Index x5 = b(3);
  const Index x2 = b(4);
  const Index x3 = b(5);
  // tensor label: Gamma0
  std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x1);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(1)->get_block(x3, x2);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x1, x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i1)))]
                  += (1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i1+x1.size()*(i3+x3.size()*(i2)))))] * fdata[i3+x3.size()*(i2)];
              }
            }
          }
        }
      }
    }
  }
  out()->put_block(odata, x5, x0, x4, x1);
}

void Task1::Task_local::compute() {
  const Index x1 = b(0);
  const Index x6 = b(1);
  const Index x0 = b(2);
  const Index x7 = b(3);
  const Index x2 = b(4);
  const Index x3 = b(5);
  const Index x4 = b(6);
  const Index x5 = b(7);
  // tensor label: Gamma2
  std::unique_ptr<double[]> odata = out()->move_block(x7, x0, x6, x1);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(2)->get_block(x5, x4, x3, x2);
  if (x3 == x4) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x1, x5, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i4 = 0; i4 != x4.size(); ++i4) {
                  odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i1)))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i1+x1.size()*(i5+x5.size()*(i2)))))] * fdata[i5+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x6, x1, x5, x4, x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i1)))]
                      += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i1+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))))] * fdata[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  out()->put_block(odata, x7, x0, x6, x1);
}

void Task2::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x3 = b(3);
  // tensor label: Gamma7
  std::unique_ptr<double[]> odata = out()->move_block(x3, x1, x2, x0);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x2, x0);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), x1.size(), x2.size(), x0.size());
  }
  out()->put_block(odata, x3, x1, x2, x0);
}

void Task3::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  // tensor label: Gamma1
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}

void Task4::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x4 = b(3);
  const Index x0 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma3
  std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x1, x3, x2);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x1, x3, x2);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x0.size(), x4.size(), x1.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, x5, x0, x4, x1, x3, x2);
}

void Task5::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index x0 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma4
  std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x3, x2, x1);
  {
    if (x2 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i3+x3.size()*(i3+x2.size()*(i1)))))]
                  += (1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i1)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x0, x4, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x5, x0, x4, x3, x2, x1);
}

void Task6::Task_local::compute() {
  const Index x1 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x4 = b(3);
  const Index x0 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma5
  std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x2, x3, x1);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x2, x3, x1);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x5.size(), x0.size(), x4.size(), x2.size(), x3.size(), x1.size());
  }
  out()->put_block(odata, x5, x0, x4, x2, x3, x1);
}

void Task8::Task_local::compute() {
  const Index x1 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(x1, a2, x0, a1);
  {
    // tensor label: I0
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, a1, a2);
    sort_indices<1,3,0,2,1,1,1,1>(i0data, odata, x0.size(), x1.size(), a1.size(), a2.size());
  }
  out()->put_block(odata, x1, a2, x0, a1);
}

void Task9::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);
  // tensor label: I0
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a2.size());
      // tensor label: I1
      std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x0, x4, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x0, x4, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), x0.size(), x4.size(), x1.size());
      dgemm_("T", "N", a1.size()*a2.size(), x0.size()*x1.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x1, a1, a2);
}

void Task10::Task_local::compute() {
  const Index x5 = b(0);
  const Index x0 = b(1);
  const Index x4 = b(2);
  const Index x1 = b(3);
  // tensor label: I1
  std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x1);
  {
    // tensor label: Gamma0
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x1);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x5.size(), x0.size(), x4.size(), x1.size());
  }
  out()->put_block(odata, x5, x0, x4, x1);
}

void Task11::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);
  // tensor label: I0
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      // tensor label: Gamma2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x0, x6, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x6.size(), x1.size());
      // tensor label: I6
      std::unique_ptr<double[]> i1data = in(1)->get_block(x7, a1, x6, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x7, a1, x6, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x7.size(), a1.size(), x6.size(), a2.size());
      dgemm_("T", "N", x0.size()*x1.size(), a1.size()*a2.size(), x7.size()*x6.size(),
             1.0, i0data_sorted, x7.size()*x6.size(), i1data_sorted, x7.size()*x6.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), a2.size());
  out()->put_block(odata, x0, x1, a1, a2);
}

void Task12::Task_local::compute() {
  const Index x7 = b(0);
  const Index a1 = b(1);
  const Index x6 = b(2);
  const Index a2 = b(3);
  // tensor label: I6
  std::unique_ptr<double[]> odata = out()->move_block(x7, a1, x6, a2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a1, x6, a2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x7.size(), a1.size(), x6.size(), a2.size());
  }
  out()->put_block(odata, x7, a1, x6, a2);
}

void Task13::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);
  // tensor label: I0
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma7
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x2, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x2, x0)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x2.size(), x0.size());
      // tensor label: I20
      std::unique_ptr<double[]> i1data = in(1)->get_block(a1, a2, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, a2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a1.size(), a2.size(), x3.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), a1.size()*a2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a1.size(), a2.size());
  out()->put_block(odata, x0, x1, a1, a2);
}

void Task14::Task_local::compute() {
  const Index a1 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I20
  std::unique_ptr<double[]> odata = out()->move_block(a1, a2, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, a2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, a2, x3, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, x2, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, x2, a4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), x2.size(), a4.size());
      // tensor label: I21
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a1, a3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a1, a3, a2)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), a1.size(), a3.size(), a2.size());
      dgemm_("T", "N", x3.size()*x2.size(), a1.size()*a2.size(), a4.size()*a3.size(),
             1.0, i0data_sorted, a4.size()*a3.size(), i1data_sorted, a4.size()*a3.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), a1.size(), a2.size());
  out()->put_block(odata, a1, a2, x3, x2);
}

void Task15::Task_local::compute() {
  const Index a4 = b(0);
  const Index a1 = b(1);
  const Index a3 = b(2);
  const Index a2 = b(3);
  // tensor label: I21
  std::unique_ptr<double[]> odata = out()->move_block(a4, a1, a3, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1, a3, a2);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, a4.size(), a1.size(), a3.size(), a2.size());
  }
  out()->put_block(odata, a4, a1, a3, a2);
}

void Task16::Task_local::compute() {
  const Index x1 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(x1, a2, x0, a1);
  {
    // tensor label: I2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, a1, x0, x1);
    sort_indices<3,0,2,1,1,1,1,1>(i0data, odata, a2.size(), a1.size(), x0.size(), x1.size());
  }
  {
    // tensor label: I2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a1, a2, x1, x0);
    sort_indices<2,1,3,0,1,1,1,1>(i0data, odata, a1.size(), a2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x1, a2, x0, a1);
}

void Task17::Task_local::compute() {
  const Index a2 = b(0);
  const Index a1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I2
  std::unique_ptr<double[]> odata = out()->move_block(a2, a1, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, a1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, a1, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma1
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: I3
      std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x3, a1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x3, a1, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, a2.size(), x3.size(), a1.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), a2.size()*a1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a2.size(), a1.size());
  out()->put_block(odata, a2, a1, x0, x1);
}

void Task18::Task_local::compute() {
  const Index a2 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index x2 = b(3);
  // tensor label: I3
  std::unique_ptr<double[]> odata = out()->move_block(a2, x3, a1, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, x3, a1, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x3, a1, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a3.size());
    // tensor label: I4
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, a2)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), a2.size());
    dgemm_("T", "N", x3.size()*a1.size()*x2.size(), a2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x3.size()*a1.size()*x2.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size(), x2.size(), a2.size());
  out()->put_block(odata, a2, x3, a1, x2);
}

void Task19::Task_local::compute() {
  const Index a3 = b(0);
  const Index a2 = b(1);
  // tensor label: I4
  std::unique_ptr<double[]> odata = out()->move_block(a3, a2);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
    sort_indices<0,1,1,1,2,1>(i0data, odata, a3.size(), a2.size());
  }
  out()->put_block(odata, a3, a2);
}

void Task20::Task_local::compute() {
  const Index a2 = b(0);
  const Index a1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I2
  std::unique_ptr<double[]> odata = out()->move_block(a2, a1, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, a1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, a1, x0, x1), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2, x3, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, a2, x3, x2)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), a2.size(), x3.size(), x2.size());
        // tensor label: I8
        std::unique_ptr<double[]> i1data = in(1)->get_block(a1, a3, x0, x1, x3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, a3, x0, x1, x3, x2)]);
        sort_indices<1,4,5,0,2,3,0,1,1,1>(i1data, i1data_sorted, a1.size(), a3.size(), x0.size(), x1.size(), x3.size(), x2.size());
        dgemm_("T", "N", a2.size(), a1.size()*x0.size()*x1.size(), a3.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, a3.size()*x3.size()*x2.size(), i1data_sorted, a3.size()*x3.size()*x2.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, a2.size(), a1.size(), x0.size(), x1.size());
  out()->put_block(odata, a2, a1, x0, x1);
}

void Task21::Task_local::compute() {
  const Index a1 = b(0);
  const Index a3 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  const Index x3 = b(4);
  const Index x2 = b(5);
  // tensor label: I8
  std::unique_ptr<double[]> odata = out()->move_block(a1, a3, x0, x1, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, a3, x0, x1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, a3, x0, x1, x3, x2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: Gamma3
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x1, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x1, x3, x2)]);
      sort_indices<0,2,1,3,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x1.size(), x3.size(), x2.size());
      // tensor label: I9
      std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a1, x4, a3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a1, x4, a3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), a1.size(), x4.size(), a3.size());
      dgemm_("T", "N", x0.size()*x1.size()*x3.size()*x2.size(), a1.size()*a3.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, x0.size()*x1.size()*x3.size()*x2.size());
    }
  }
  sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x3.size(), x2.size(), a1.size(), a3.size());
  out()->put_block(odata, a1, a3, x0, x1, x3, x2);
}

void Task22::Task_local::compute() {
  const Index x5 = b(0);
  const Index a1 = b(1);
  const Index x4 = b(2);
  const Index a3 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata = out()->move_block(x5, a1, x4, a3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x5.size(), a1.size(), x4.size(), a3.size());
  }
  out()->put_block(odata, x5, a1, x4, a3);
}

void Task23::Task_local::compute() {
  const Index a2 = b(0);
  const Index a1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I2
  std::unique_ptr<double[]> odata = out()->move_block(a2, a1, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, a1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, a1, x0, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: Gamma4
          std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x2, x1)]);
          sort_indices<0,2,3,4,1,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
          // tensor label: I11
          std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, a2, x5, a1, x4);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, a2, x5, a1, x4)]);
          sort_indices<3,5,0,1,2,4,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), a2.size(), x5.size(), a1.size(), x4.size());
          dgemm_("T", "N", x0.size()*x1.size(), a2.size()*a1.size(), x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, x0.size()*x1.size());
        }
      }
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a2.size(), a1.size());
  out()->put_block(odata, a2, a1, x0, x1);
}

void Task24::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index a2 = b(2);
  const Index x5 = b(3);
  const Index a1 = b(4);
  const Index x4 = b(5);
  // tensor label: I11
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, a2, x5, a1, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, a2, x5, a1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, a2, x5, a1, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a3.size());
    // tensor label: I12
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x3, x2, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x3, x2, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), x3.size(), x2.size(), a2.size());
    dgemm_("T", "N", x5.size()*a1.size()*x4.size(), x3.size()*x2.size()*a2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x5.size()*a1.size()*x4.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), x3.size(), x2.size(), a2.size());
  out()->put_block(odata, x3, x2, a2, x5, a1, x4);
}

void Task25::Task_local::compute() {
  const Index a3 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a2 = b(3);
  // tensor label: I12
  std::unique_ptr<double[]> odata = out()->move_block(a3, x3, x2, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x3, x2, a2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, a3.size(), x3.size(), x2.size(), a2.size());
  }
  out()->put_block(odata, a3, x3, x2, a2);
}

void Task26::Task_local::compute() {
  const Index a2 = b(0);
  const Index a1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I2
  std::unique_ptr<double[]> odata = out()->move_block(a2, a1, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, a1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, a1, x0, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        for (auto& x3 : *range_[1]) {
          // tensor label: Gamma5
          std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x2, x3, x1);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x2, x3, x1)]);
          sort_indices<0,2,3,4,1,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x2.size(), x3.size(), x1.size());
          // tensor label: I14
          std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a2, x2, x5, a1, x4);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a2, x2, x5, a1, x4)]);
          sort_indices<3,5,2,0,1,4,0,1,1,1>(i1data, i1data_sorted, x3.size(), a2.size(), x2.size(), x5.size(), a1.size(), x4.size());
          dgemm_("T", "N", x0.size()*x1.size(), a2.size()*a1.size(), x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, x0.size()*x1.size());
        }
      }
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a2.size(), a1.size());
  out()->put_block(odata, a2, a1, x0, x1);
}

void Task27::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  const Index x2 = b(2);
  const Index x5 = b(3);
  const Index a1 = b(4);
  const Index x4 = b(5);
  // tensor label: I14
  std::unique_ptr<double[]> odata = out()->move_block(x3, a2, x2, x5, a1, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a2, x2, x5, a1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a2, x2, x5, a1, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a3.size());
    // tensor label: I15
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a2, a3, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a2, a3, x2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a2.size(), a3.size(), x2.size());
    dgemm_("T", "N", x5.size()*a1.size()*x4.size(), x3.size()*a2.size()*x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x5.size()*a1.size()*x4.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), x3.size(), a2.size(), x2.size());
  out()->put_block(odata, x3, a2, x2, x5, a1, x4);
}

void Task28::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  const Index a3 = b(2);
  const Index x2 = b(3);
  // tensor label: I15
  std::unique_ptr<double[]> odata = out()->move_block(x3, a2, a3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, a3, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), a2.size(), a3.size(), x2.size());
  }
  out()->put_block(odata, x3, a2, a3, x2);
}

void Task29::Task_local::compute() {
  const Index a2 = b(0);
  const Index a1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I2
  std::unique_ptr<double[]> odata = out()->move_block(a2, a1, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, a1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, a1, x0, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: Gamma3
          std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x1, x3, x2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x1, x3, x2)]);
          sort_indices<0,2,4,5,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x1.size(), x3.size(), x2.size());
          // tensor label: I17
          std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, a2, x5, a1, x4);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, a2, x5, a1, x4)]);
          sort_indices<3,5,0,1,2,4,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), a2.size(), x5.size(), a1.size(), x4.size());
          dgemm_("T", "N", x0.size()*x1.size(), a2.size()*a1.size(), x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, x0.size()*x1.size());
        }
      }
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a2.size(), a1.size());
  out()->put_block(odata, a2, a1, x0, x1);
}

void Task30::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index a2 = b(2);
  const Index x5 = b(3);
  const Index a1 = b(4);
  const Index x4 = b(5);
  // tensor label: I17
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, a2, x5, a1, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, a2, x5, a1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, a2, x5, a1, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a3.size());
    // tensor label: I18
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, a3, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, a3, a2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), a3.size(), a2.size());
    dgemm_("T", "N", x5.size()*a1.size()*x4.size(), x3.size()*x2.size()*a2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x5.size()*a1.size()*x4.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), x3.size(), x2.size(), a2.size());
  out()->put_block(odata, x3, x2, a2, x5, a1, x4);
}

void Task31::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index a3 = b(2);
  const Index a2 = b(3);
  // tensor label: I18
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, a3, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, a3, a2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x2.size(), a3.size(), a2.size());
  }
  out()->put_block(odata, x3, x2, a3, a2);
}

void Task33::Task_local::compute() {
  const Index x1 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(x1, a2, x0, a1);
  {
    // tensor label: I22
    std::unique_ptr<double[]> i0data = in(0)->get_block(a1, a2, x0, x1);
    sort_indices<3,1,2,0,1,1,1,1>(i0data, odata, a1.size(), a2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x1, a2, x0, a1);
}

void Task34::Task_local::compute() {
  const Index a1 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I22
  std::unique_ptr<double[]> odata = out()->move_block(a1, a2, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, a2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, a2, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma1
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: I23
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, x2, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, x2, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), x2.size(), a2.size());
      dgemm_("T", "N", x0.size()*x1.size(), a1.size()*a2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), a2.size());
  out()->put_block(odata, a1, a2, x0, x1);
}

void Task35::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index x2 = b(2);
  const Index a2 = b(3);
  // tensor label: I23
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1, x2, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), a1.size(), x2.size(), a2.size());
  }
  out()->put_block(odata, x3, a1, x2, a2);
}

void Task37::Task_local::compute() {
  const Index x1 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(x1, a2, x0, a1);
  {
    // tensor label: I24
    std::unique_ptr<double[]> i0data = in(0)->get_block(a1, a2, x0, x1);
    sort_indices<3,1,2,0,1,1,1,1>(i0data, odata, a1.size(), a2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x1, a2, x0, a1);
}

void Task38::Task_local::compute() {
  const Index a1 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I24
  std::unique_ptr<double[]> odata = out()->move_block(a1, a2, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, a2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, a2, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma1
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: I25
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, x2, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, x2, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), x2.size(), a2.size());
      dgemm_("T", "N", x0.size()*x1.size(), a1.size()*a2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), a2.size());
  out()->put_block(odata, a1, a2, x0, x1);
}

void Task39::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index x2 = b(2);
  const Index a2 = b(3);
  // tensor label: I25
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1, x2, a2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), a1.size(), x2.size(), a2.size());
  }
  out()->put_block(odata, x3, a1, x2, a2);
}

#endif
