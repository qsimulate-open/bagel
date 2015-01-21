//
// BAGEL - Parallel electron correlation program.
// Filename: CAS_test_tasks.cc
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


#include <src/smith/CAS_test_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CAS_test;

void Task1::Task_local::compute() {
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


void Task2::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  // tensor label: Gamma2
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task3::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x4 = b(3);
  const Index x0 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma9
  std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x1, x3, x2);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x1, x3, x2);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x0.size(), x4.size(), x1.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, x5, x0, x4, x1, x3, x2);
}


void Task4::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x1 = b(4);
  const Index x3 = b(5);
  const Index x2 = b(6);
  // tensor label: Gamma12
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x4, x1);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(1)->get_block(x3, x2);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x1, x3, x2);
    for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
      for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
        for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix4+x4.size()*(ix1))))]
                    += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix4+x4.size()*(ix1+x1.size()*(ix3+x3.size()*(ix2))))))] * fdata[ix3+x3.size()*(ix2)];
                }
              }
            }
          }
        }
      }
    }
  }
  out()->put_block(odata, ci0, x5, x0, x4, x1);
}


void Task5::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  // tensor label: Gamma13
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}


void Task6::Task_local::compute() {
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


void Task7::Task_local::compute() {
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


void Task8::Task_local::compute() {
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


void Task9::Task_local::compute() {
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
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

      // tensor label: I6
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

      dgemm_("T", "N", a1.size()*a2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }

  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x1, a1, a2);
}


void Task10::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);

  // tensor label: I6
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    dscal_(x3.size()*x0.size()*x2.size()*x1.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
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

  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

      // tensor label: I8
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

      dgemm_("T", "N", a1.size()*a2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }

  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x1, a1, a2);
}


void Task12::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);

  // tensor label: I8
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task13::Task_local::compute() {
  const Index x1 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);

  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(x1, a2, x0, a1);
  {
    // tensor label: I2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, a1, a2);
    sort_indices<1,3,0,2,1,1,1,1>(i0data, odata, x0.size(), x1.size(), a1.size(), a2.size());
  }
  {
    // tensor label: I2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0, a2, a1);
    sort_indices<0,2,1,3,1,1,1,1>(i0data, odata, x1.size(), x0.size(), a2.size(), a1.size());
  }
  out()->put_block(odata, x1, a2, x0, a1);
}


void Task14::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);

  // tensor label: I2
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

  for (auto& a3 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), a2.size());

    // tensor label: I3
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, a1, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, a1, a3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), a1.size(), a3.size());

    dgemm_("T", "N", a2.size(), x0.size()*x1.size()*a1.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, a2.size());
  }

  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), x0.size(), x1.size(), a1.size());
  out()->put_block(odata, x0, x1, a1, a2);
}


void Task15::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a3 = b(3);

  // tensor label: I3
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a3), 0.0);

  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a3.size());

      // tensor label: I4
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

      dgemm_("T", "N", a1.size()*a3.size(), x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*a3.size());
    }
  }

  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x1, a1, a3);
}


void Task16::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);

  // tensor label: I4
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task17::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a1 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);

  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
  sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());

  // tensor label: I10
  std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, a1, a2);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, a1, a2)]);
  sort_indices<3,1,2,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), a1.size(), a2.size());

  energy_ += ddot_(x0.size()*x1.size()*a1.size()*a2.size(), i0data_sorted, 1, i1data_sorted, 1);
}


void Task18::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);

  // tensor label: I10
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a2.size());

      // tensor label: I11
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


void Task19::Task_local::compute() {
  energy_ = 0.0;
  const Index x5 = b(0);
  const Index x0 = b(1);
  const Index x4 = b(2);
  const Index x1 = b(3);

  // tensor label: I11
  std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x1);
  {
    // tensor label: Gamma0
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x1);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x5.size(), x0.size(), x4.size(), x1.size());
  }
  out()->put_block(odata, x5, x0, x4, x1);
}


void Task20::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);

  // tensor label: I10
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

  for (auto& a3 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), a2.size());

    // tensor label: I14
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, a1, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, a1, a3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), a1.size(), a3.size());

    dgemm_("T", "N", a2.size(), x0.size()*x1.size()*a1.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, a2.size());
  }

  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), x0.size(), x1.size(), a1.size());
  out()->put_block(odata, x0, x1, a1, a2);
}


void Task21::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a3 = b(3);

  // tensor label: I14
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a3), 0.0);

  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a3.size());

      // tensor label: I15
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

      dgemm_("T", "N", a1.size()*a3.size(), x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*a3.size());
    }
  }

  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x1, a1, a3);
}


void Task22::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);

  // tensor label: I15
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task23::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);

  // tensor label: I10
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

      // tensor label: I18
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

      dgemm_("T", "N", a1.size()*a2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }

  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x1, a1, a2);
}


void Task24::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);

  // tensor label: I18
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    dscal_(x3.size()*x0.size()*x2.size()*x1.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task25::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);

  // tensor label: I10
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

      // tensor label: I21
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

      dgemm_("T", "N", a1.size()*a2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }

  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x1, a1, a2);
}


void Task26::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);

  // tensor label: I21
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task27::Task_local::compute(){
  correction_ = 0.0;
  const Index x0 = b(0);
  const Index a1 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);

  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
  sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());

  // tensor label: I23
  std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, a1, a2);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, a1, a2)]);
  sort_indices<3,1,2,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), a1.size(), a2.size());

  correction_ += ddot_(x0.size()*x1.size()*a1.size()*a2.size(), i0data_sorted, 1, i1data_sorted, 1);
}


void Task28::Task_local::compute(){
  correction_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);

  // tensor label: I23
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

      // tensor label: I24
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

      dgemm_("T", "N", a1.size()*a2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }

  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x1, a1, a2);
}


void Task29::Task_local::compute(){
  correction_ = 0.0;
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);

  // tensor label: I24
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task31::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);

  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3);
  {
    // tensor label: I25
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2);
    sort_indices<1,0,1,1,1,1>(i0data, odata, x3.size(), x2.size());
  }
  out()->put_block(odata, x2, x3);
}


void Task32::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);

  // tensor label: I25
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2), 0.0);

  for (auto& x0 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& x1 : *range_[1]) {
        for (auto& a2 : *range_[2]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());

          // tensor label: I26
          std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, x3, x2, a1, a2);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, x3, x2, a1, a2)]);
          sort_indices<0,4,1,5,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), x3.size(), x2.size(), a1.size(), a2.size());

          dgemm_("T", "N", 1, x3.size()*x2.size(), x0.size()*x1.size()*a1.size()*a2.size(),
                 1.0, i0data_sorted, x0.size()*x1.size()*a1.size()*a2.size(), i1data_sorted, x0.size()*x1.size()*a1.size()*a2.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }

  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size());
  out()->put_block(odata, x3, x2);
}


void Task33::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index a1 = b(4);
  const Index a2 = b(5);

  // tensor label: I26
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, x3, x2, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, x3, x2, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, x3, x2, a1, a2), 0.0);

  for (auto& x4 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a2)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a2.size());

      // tensor label: I27
      std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x0, x4, x1, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x0, x4, x1, x3, x2)]);
      sort_indices<2,0,1,3,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x0.size(), x4.size(), x1.size(), x3.size(), x2.size());

      dgemm_("T", "N", a1.size()*a2.size(), x0.size()*x1.size()*x3.size()*x2.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }

  sort_indices<2,3,4,5,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size(), x3.size(), x2.size());
  out()->put_block(odata, x0, x1, x3, x2, a1, a2);
}


void Task34::Task_local::compute() {
  const Index x5 = b(0);
  const Index x0 = b(1);
  const Index x4 = b(2);
  const Index x1 = b(3);
  const Index x3 = b(4);
  const Index x2 = b(5);

  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x1, x3, x2);
  {
    // tensor label: Gamma9
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x1, x3, x2);
    sort_indices<0,1,2,3,4,5,1,1,1,2>(i0data, odata, x5.size(), x0.size(), x4.size(), x1.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, x5, x0, x4, x1, x3, x2);
}


void Task35::Task_local::compute() {
  const Index a2 = b(0);
  const Index a3 = b(1);

  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(a2, a3);
  {
    // tensor label: I28
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
    sort_indices<1,0,1,1,1,1>(i0data, odata, a3.size(), a2.size());
  }
  out()->put_block(odata, a2, a3);
}


void Task36::Task_local::compute() {
  const Index a3 = b(0);
  const Index a2 = b(1);

  // tensor label: I28
  std::unique_ptr<double[]> odata = out()->move_block(a3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, a2), 0.0);

  for (auto& x0 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());

        // tensor label: I29
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, a1, a3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, a1, a3)]);
        sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), a1.size(), a3.size());

        dgemm_("T", "N", a2.size(), a3.size(), x0.size()*x1.size()*a1.size(),
               1.0, i0data_sorted, x0.size()*x1.size()*a1.size(), i1data_sorted, x0.size()*x1.size()*a1.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }

  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a2.size(), a3.size());
  out()->put_block(odata, a3, a2);
}


void Task37::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a3 = b(3);

  // tensor label: I29
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a3), 0.0);

  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a3)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a3.size());

      // tensor label: I30
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

      dgemm_("T", "N", a1.size()*a3.size(), x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*a3.size());
    }
  }

  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x1, a1, a3);
}


void Task38::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);

  // tensor label: I30
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task40::Task_local::compute() {
  const Index x1 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);

  // tensor label: Den1
  std::unique_ptr<double[]> odata = out()->move_block(x1, a2, x0, a1);
  {
    // tensor label: I31
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, a1, a2);
    sort_indices<1,3,0,2,1,1,1,1>(i0data, odata, x0.size(), x1.size(), a1.size(), a2.size());
  }
  out()->put_block(odata, x1, a2, x0, a1);
}


void Task41::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);

  // tensor label: I31
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

      // tensor label: I32
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

      dgemm_("T", "N", a1.size()*a2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }

  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x1, a1, a2);
}


void Task42::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);

  // tensor label: I32
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task44::Task_local::compute() {
  const Index ci0 = b(0);

  // tensor label: deci
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  {
    // tensor label: I33
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0);
    sort_indices<0,1,1,1,1>(i0data, odata, ci0.size());
  }
  out()->put_block(odata, ci0);
}


void Task45::Task_local::compute() {
  const Index ci0 = b(0);

  // tensor label: I33
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);

  for (auto& x0 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& x1 : *range_[1]) {
        for (auto& a2 : *range_[2]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());

          // tensor label: I34
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1, a1, a2);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1, a1, a2)]);
          sort_indices<1,3,2,4,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size(), a1.size(), a2.size());

          dgemm_("T", "N", 1, ci0.size(), x0.size()*x1.size()*a1.size()*a2.size(),
                 1.0, i0data_sorted, x0.size()*x1.size()*a1.size()*a2.size(), i1data_sorted, x0.size()*x1.size()*a1.size()*a2.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }

  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}


void Task46::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  const Index a2 = b(4);

  // tensor label: I34
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, a2), 0.0);

  for (auto& x4 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a2)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a2.size());

      // tensor label: I35
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x0, x4, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x0, x4, x1)]);
      sort_indices<3,1,0,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x0.size(), x4.size(), x1.size());

      dgemm_("T", "N", a1.size()*a2.size(), ci0.size()*x0.size()*x1.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }

  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), ci0.size(), x0.size(), x1.size());
  out()->put_block(odata, ci0, x0, x1, a1, a2);
}


void Task47::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x1 = b(4);

  // tensor label: I35
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x4, x1);
  {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x1);
    sort_indices<0,1,2,3,4,1,1,1,2>(i0data, odata, ci0.size(), x5.size(), x0.size(), x4.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x0, x4, x1);
}


void Task48::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  const Index a2 = b(4);

  // tensor label: I34
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, a2), 0.0);

  for (auto& a3 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), a2.size());

    // tensor label: I38
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1, a1, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1, a1, a3)]);
    sort_indices<4,0,1,2,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size(), a1.size(), a3.size());

    dgemm_("T", "N", a2.size(), ci0.size()*x0.size()*x1.size()*a1.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, a2.size());
  }

  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a2.size(), ci0.size(), x0.size(), x1.size(), a1.size());
  out()->put_block(odata, ci0, x0, x1, a1, a2);
}


void Task49::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  const Index a3 = b(4);

  // tensor label: I38
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, a3), 0.0);

  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a3)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a3.size());

      // tensor label: I39
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
      sort_indices<3,1,0,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());

      dgemm_("T", "N", a1.size()*a3.size(), ci0.size()*x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*a3.size());
    }
  }

  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), ci0.size(), x0.size(), x1.size());
  out()->put_block(odata, ci0, x0, x1, a1, a3);
}


void Task50::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);

  // tensor label: I39
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma13
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}


void Task51::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  const Index a2 = b(4);

  // tensor label: I34
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, a2), 0.0);

  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

      // tensor label: I49
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
      sort_indices<3,1,0,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());

      dgemm_("T", "N", a1.size()*a2.size(), ci0.size()*x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }

  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), ci0.size(), x0.size(), x1.size());
  out()->put_block(odata, ci0, x0, x1, a1, a2);
}


void Task52::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);

  // tensor label: I49
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma13
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    dscal_(ci0.size()*x3.size()*x0.size()*x2.size()*x1.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}


void Task53::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  const Index a2 = b(4);

  // tensor label: I34
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, a2), 0.0);

  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

      // tensor label: I55
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
      sort_indices<3,1,0,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());

      dgemm_("T", "N", a1.size()*a2.size(), ci0.size()*x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }

  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), ci0.size(), x0.size(), x1.size());
  out()->put_block(odata, ci0, x0, x1, a1, a2);
}


void Task54::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);

  // tensor label: I55
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma13
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}


void Task55::Task_local::compute() {
  const Index ci0 = b(0);

  // tensor label: I33
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);

  for (auto& x5 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& x4 : *range_[1]) {
        for (auto& a2 : *range_[2]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a2)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a2.size());

          // tensor label: I41
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x4, a1, a2);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x4, a1, a2)]);
          sort_indices<1,3,2,4,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x4.size(), a1.size(), a2.size());

          dgemm_("T", "N", 1, ci0.size(), x5.size()*x4.size()*a1.size()*a2.size(),
                 1.0, i0data_sorted, x5.size()*x4.size()*a1.size()*a2.size(), i1data_sorted, x5.size()*x4.size()*a1.size()*a2.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }

  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}


void Task56::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index a1 = b(3);
  const Index a2 = b(4);

  // tensor label: I41
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x5, x4, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x5, x4, a1, a2), 0.0);

  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());

      // tensor label: I42
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x0, x4, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x0, x4, x1)]);
      sort_indices<4,2,0,1,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x0.size(), x4.size(), x1.size());

      dgemm_("T", "N", a1.size()*a2.size(), ci0.size()*x5.size()*x4.size(), x0.size()*x1.size(),
             1.0, i0data_sorted, x0.size()*x1.size(), i1data_sorted, x0.size()*x1.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }

  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), ci0.size(), x5.size(), x4.size());
  out()->put_block(odata, ci0, x5, x4, a1, a2);
}


void Task57::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x1 = b(4);

  // tensor label: I42
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x4, x1);
  {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x1);
    sort_indices<0,1,2,3,4,1,1,1,2>(i0data, odata, ci0.size(), x5.size(), x0.size(), x4.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x0, x4, x1);
}


void Task58::Task_local::compute() {
  const Index ci0 = b(0);

  // tensor label: I33
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);

  for (auto& x3 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& x2 : *range_[1]) {
        for (auto& a2 : *range_[2]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

          // tensor label: I44
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, a1, a2);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, a1, a2)]);
          sort_indices<1,3,2,4,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), a1.size(), a2.size());

          dgemm_("T", "N", 1, ci0.size(), x3.size()*x2.size()*a1.size()*a2.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*a1.size()*a2.size(), i1data_sorted, x3.size()*x2.size()*a1.size()*a2.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }

  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}


void Task59::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a1 = b(3);
  const Index a2 = b(4);

  // tensor label: I44
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, a2), 0.0);

  for (auto& a3 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), a2.size());

    // tensor label: I45
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, a1, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, a1, a3)]);
    sort_indices<4,0,1,2,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), a1.size(), a3.size());

    dgemm_("T", "N", a2.size(), ci0.size()*x3.size()*x2.size()*a1.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, a2.size());
  }

  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a2.size(), ci0.size(), x3.size(), x2.size(), a1.size());
  out()->put_block(odata, ci0, x3, x2, a1, a2);
}


void Task60::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a1 = b(3);
  const Index a3 = b(4);

  // tensor label: I45
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, a3), 0.0);

  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a3)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a3.size());

      // tensor label: I46
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
      sort_indices<4,2,0,1,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());

      dgemm_("T", "N", a1.size()*a3.size(), ci0.size()*x3.size()*x2.size(), x0.size()*x1.size(),
             1.0, i0data_sorted, x0.size()*x1.size(), i1data_sorted, x0.size()*x1.size(),
             1.0, odata_sorted, a1.size()*a3.size());
    }
  }

  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), ci0.size(), x3.size(), x2.size());
  out()->put_block(odata, ci0, x3, x2, a1, a3);
}


void Task61::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);

  // tensor label: I46
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma13
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}


void Task62::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a1 = b(3);
  const Index a2 = b(4);

  // tensor label: I44
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, a2), 0.0);

  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());

      // tensor label: I52
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
      sort_indices<4,2,0,1,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());

      dgemm_("T", "N", a1.size()*a2.size(), ci0.size()*x3.size()*x2.size(), x0.size()*x1.size(),
             1.0, i0data_sorted, x0.size()*x1.size(), i1data_sorted, x0.size()*x1.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }

  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), ci0.size(), x3.size(), x2.size());
  out()->put_block(odata, ci0, x3, x2, a1, a2);
}


void Task63::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);

  // tensor label: I52
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma13
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    dscal_(ci0.size()*x3.size()*x0.size()*x2.size()*x1.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}


void Task64::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a1 = b(3);
  const Index a2 = b(4);

  // tensor label: I44
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, a2), 0.0);

  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());

      // tensor label: I58
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
      sort_indices<4,2,0,1,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());

      dgemm_("T", "N", a1.size()*a2.size(), ci0.size()*x3.size()*x2.size(), x0.size()*x1.size(),
             1.0, i0data_sorted, x0.size()*x1.size(), i1data_sorted, x0.size()*x1.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }

  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), ci0.size(), x3.size(), x2.size());
  out()->put_block(odata, ci0, x3, x2, a1, a2);
}


void Task65::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);

  // tensor label: I58
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma13
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}


