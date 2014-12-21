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
  const Index x0 = b(0);
  const Index x1 = b(1);
  // scalar
  // tensor label: Gamma0
  std::unique_ptr<double[]> odata = out()->move_block();
  // associated with merged
  std::unique_ptr<double[]> fdata = in(1)->get_block(x1, x0);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        odata[0]
          += (1.0) * i0data[i1+x1.size()*(i0)] * fdata[i1+x1.size()*(i0)];
      }
    }
  }
  out()->put_block(odata);
}


void Task2::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: Gamma2
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task3::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: Gamma6
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


void Task4::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  // tensor label: Gamma14
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task5::Task_local::compute() {
  const Index x1 = b(0);
  const Index x4 = b(1);
  const Index x0 = b(2);
  const Index x5 = b(3);
  const Index x2 = b(4);
  const Index x3 = b(5);
  // tensor label: Gamma16
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


void Task6::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x4 = b(3);
  const Index x0 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma67
  std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x1, x3, x2);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x1, x3, x2);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x0.size(), x4.size(), x1.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, x5, x0, x4, x1, x3, x2);
}


void Task7::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: Gamma72
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


void Task8::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: Gamma74
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}


void Task9::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  // tensor label: Gamma78
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


void Task10::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  // tensor label: Gamma86
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}


void Task11::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x1 = b(4);
  const Index x3 = b(5);
  const Index x2 = b(6);
  // tensor label: Gamma88
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


void Task12::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);

  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c3, a4, c1, a2);
  {
    // tensor label: I0
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
    sort_indices<2,1,0,3,1,1,1,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
  }
  out()->put_block(odata, c3, a4, c1, a2);
}


void Task13::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);

  // tensor label: I0
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
    dscal_(c1.size()*a4.size()*c3.size()*a2.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,1,1,4,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a4);
    dscal_(c1.size()*a2.size()*c3.size()*a4.size(), e0_, i1data.get(), 1);
    sort_indices<0,3,2,1,1,1,-8,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i2data = in(1)->get_block(c1, a4, c3, a2);
    sort_indices<0,1,2,3,1,1,-4,1>(i2data, odata, c1.size(), a4.size(), c3.size(), a2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i3data = in(1)->get_block(c1, a2, c3, a4);
    sort_indices<0,3,2,1,1,1,8,1>(i3data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->put_block(odata, c1, a4, c3, a2);
}


void Task14::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);

  // tensor label: I0
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, c3, a2), 0.0);

  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());

  // tensor label: I1
  std::unique_ptr<double[]> i1data = in(1)->get_block();
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size()]);
  sort_indices<0,1,1,1>(i1data, i1data_sorted);

  dgemm_("T", "N", c1.size()*a4.size()*c3.size()*a2.size(), 1, 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a4.size()*c3.size()*a2.size());

  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->put_block(odata, c1, a4, c3, a2);
}


void Task15::Task_local::compute() {

  // tensor label: I1
  std::unique_ptr<double[]> odata = out()->move_block();
  {
    // scalar
    // tensor label: Gamma0
    std::unique_ptr<double[]> i0data = in(0)->get_block();
    sort_indices<1,1,-4,1>(i0data, odata);
  }
  out()->put_block(odata);
}


void Task16::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);

  // tensor label: I0
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, c3, a2), 0.0);

  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());

  // tensor label: I3
  std::unique_ptr<double[]> i1data = in(1)->get_block();
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size()]);
  sort_indices<0,1,1,1>(i1data, i1data_sorted);

  dgemm_("T", "N", c1.size()*a2.size()*c3.size()*a4.size(), 1, 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a2.size()*c3.size()*a4.size());

  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->put_block(odata, c1, a4, c3, a2);
}


void Task17::Task_local::compute() {

  // tensor label: I3
  std::unique_ptr<double[]> odata = out()->move_block();
  {
    // scalar
    // tensor label: Gamma0
    std::unique_ptr<double[]> i0data = in(0)->get_block();
    sort_indices<1,1,8,1>(i0data, odata);
  }
  out()->put_block(odata);
}


void Task18::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);

  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c3, a4, c1, a2);
  {
    // tensor label: I4
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, a2, c3);
    sort_indices<3,1,0,2,1,1,1,1>(i0data, odata, c1.size(), a4.size(), a2.size(), c3.size());
  }
  {
    // tensor label: I4
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, a4, c1);
    sort_indices<0,2,3,1,1,1,1,1>(i0data, odata, c3.size(), a2.size(), a4.size(), c1.size());
  }
  out()->put_block(odata, c3, a4, c1, a2);
}


void Task19::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index a2 = b(2);
  const Index c3 = b(3);

  // tensor label: I4
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, a2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, a2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, a2, c3), 0.0);

  for (auto& c5 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c5);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, c5)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c3.size(), c5.size());

    // tensor label: I5
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a4, c5, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a4, c5, a2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a4.size(), c5.size(), a2.size());

    dgemm_("T", "N", c3.size(), c1.size()*a4.size()*a2.size(), c5.size(),
           1.0, i0data_sorted, c5.size(), i1data_sorted, c5.size(),
           1.0, odata_sorted, c3.size());
  }

  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c3.size(), c1.size(), a4.size(), a2.size());
  out()->put_block(odata, c1, a4, a2, c3);
}


void Task20::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c5 = b(2);
  const Index a2 = b(3);

  // tensor label: I5
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c5, a2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c5, a2);
    sort_indices<0,1,2,3,1,1,4,1>(i0data, odata, c1.size(), a4.size(), c5.size(), a2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c5, a4);
    sort_indices<0,3,2,1,1,1,-8,1>(i1data, odata, c1.size(), a2.size(), c5.size(), a4.size());
  }
  out()->put_block(odata, c1, a4, c5, a2);
}


void Task21::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index a2 = b(2);
  const Index c3 = b(3);

  // tensor label: I4
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, a2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, a2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, a2, c3), 0.0);

  for (auto& a5 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a5, a4)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a5.size(), a4.size());

    // tensor label: I9
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a5, c3, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a5, c3, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a5.size(), c3.size(), a2.size());

    dgemm_("T", "N", a4.size(), c1.size()*c3.size()*a2.size(), a5.size(),
           1.0, i0data_sorted, a5.size(), i1data_sorted, a5.size(),
           1.0, odata_sorted, a4.size());
  }

  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), c3.size(), a2.size());
  out()->put_block(odata, c1, a4, a2, c3);
}


void Task22::Task_local::compute() {
  const Index c1 = b(0);
  const Index a5 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);

  // tensor label: I9
  std::unique_ptr<double[]> odata = out()->move_block(c1, a5, c3, a2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a5, c3, a2);
    sort_indices<0,1,2,3,1,1,-4,1>(i0data, odata, c1.size(), a5.size(), c3.size(), a2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a5);
    sort_indices<0,3,2,1,1,1,8,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a5.size());
  }
  out()->put_block(odata, c1, a5, c3, a2);
}


void Task23::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index a2 = b(2);
  const Index c3 = b(3);

  // tensor label: I4
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, a2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, a2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, a2, c3), 0.0);

  for (auto& x0 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c3.size(), x0.size());

    // tensor label: I13
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4, c1, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4, c1, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size(), c1.size(), a2.size());

    dgemm_("T", "N", c3.size(), a4.size()*c1.size()*a2.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, c3.size());
  }

  sort_indices<2,1,3,0,1,1,1,1>(odata_sorted, odata, c3.size(), a4.size(), c1.size(), a2.size());
  out()->put_block(odata, c1, a4, a2, c3);
}


void Task24::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);

  // tensor label: I13
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c1, a2), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c1, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c1.size(), a2.size());

    // tensor label: I14
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a4.size()*c1.size()*a2.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a4.size()*c1.size()*a2.size());
  }

  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), a2.size(), x0.size());
  out()->put_block(odata, x0, a4, c1, a2);
}


void Task25::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I14
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-2,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task26::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);

  // tensor label: I13
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c1, a2), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, c1, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2, c1, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c1.size(), a4.size());

    // tensor label: I17
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a2.size()*c1.size()*a4.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a2.size()*c1.size()*a4.size());
  }

  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), a4.size(), x0.size());
  out()->put_block(odata, x0, a4, c1, a2);
}


void Task27::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I17
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task28::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);

  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c2, a3, x0, a1);
  {
    // tensor label: I18
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, c2, a3, a1);
    sort_indices<1,2,0,3,1,1,1,1>(i0data, odata, x0.size(), c2.size(), a3.size(), a1.size());
  }
  out()->put_block(odata, c2, a3, x0, a1);
}


void Task29::Task_local::compute() {
  const Index x0 = b(0);
  const Index c2 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);

  // tensor label: I18
  std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, c4)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), c4.size());

      // tensor label: I19
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, c2, a3, c4, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, c2, a3, c4, a1)]);
      sort_indices<0,4,1,2,3,5,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), c2.size(), a3.size(), c4.size(), a1.size());

      dgemm_("T", "N", 1, x0.size()*c2.size()*a3.size()*a1.size(), x1.size()*c4.size(),
             1.0, i0data_sorted, x1.size()*c4.size(), i1data_sorted, x1.size()*c4.size(),
             1.0, odata_sorted, 1);
    }
  }

  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x0.size(), c2.size(), a3.size(), a1.size());
  out()->put_block(odata, x0, c2, a3, a1);
}


void Task30::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index c4 = b(4);
  const Index a1 = b(5);

  // tensor label: I19
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a3, c4, a1), 0.0);

  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c4, a1);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c4, a1)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c4.size(), a1.size());

  // tensor label: I20
  std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

  dgemm_("T", "N", c2.size()*a3.size()*c4.size()*a1.size(), x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a3.size()*c4.size()*a1.size());

  sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), c4.size(), a1.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c2, a3, c4, a1);
}


void Task31::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I20
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-4,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task32::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index c4 = b(4);
  const Index a1 = b(5);

  // tensor label: I19
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a3, c4, a1), 0.0);

  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c4, a3);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c4, a3)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c4.size(), a3.size());

  // tensor label: I23
  std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

  dgemm_("T", "N", c2.size()*a1.size()*c4.size()*a3.size(), x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a1.size()*c4.size()*a3.size());

  sort_indices<4,5,0,3,2,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), c4.size(), a3.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c2, a3, c4, a1);
}


void Task33::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I23
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,2,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task34::Task_local::compute() {
  const Index x0 = b(0);
  const Index c2 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);

  // tensor label: I18
  std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a1.size());

    // tensor label: I25
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size());

    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }

  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size());
  out()->put_block(odata, x0, c2, a3, a1);
}


void Task35::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);

  // tensor label: I25
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0);
  {
    // tensor label: Gamma6
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, x3.size(), x0.size());
  }
  out()->put_block(odata, x3, x0);
}


void Task36::Task_local::compute() {
  const Index x0 = b(0);
  const Index c2 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);

  // tensor label: I18
  std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a3.size());

    // tensor label: I27
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size());

    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }

  sort_indices<3,1,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size());
  out()->put_block(odata, x0, c2, a3, a1);
}


void Task37::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);

  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0);
  {
    // tensor label: Gamma6
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0);
    sort_indices<0,1,1,1,2,1>(i0data, odata, x3.size(), x0.size());
  }
  out()->put_block(odata, x3, x0);
}


void Task38::Task_local::compute() {
  const Index x0 = b(0);
  const Index c2 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);

  // tensor label: I18
  std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

  for (auto& c4 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, c4)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), c4.size());

    // tensor label: I29
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a3, c4, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a3, c4, a1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a3.size(), c4.size(), a1.size());

    dgemm_("T", "N", c2.size(), x0.size()*a3.size()*a1.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c2.size());
  }

  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, c2.size(), x0.size(), a3.size(), a1.size());
  out()->put_block(odata, x0, c2, a3, a1);
}


void Task39::Task_local::compute() {
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c4 = b(2);
  const Index a1 = b(3);

  // tensor label: I29
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c4, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c4, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c4, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c4.size(), a1.size());

    // tensor label: I30
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a3.size()*c4.size()*a1.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a3.size()*c4.size()*a1.size());
  }

  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c4.size(), a1.size(), x0.size());
  out()->put_block(odata, x0, a3, c4, a1);
}


void Task40::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I30
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task41::Task_local::compute() {
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c4 = b(2);
  const Index a1 = b(3);

  // tensor label: I29
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c4, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c4, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c4.size(), a3.size());

    // tensor label: I33
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a1.size()*c4.size()*a3.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a1.size()*c4.size()*a3.size());
  }

  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c4.size(), a3.size(), x0.size());
  out()->put_block(odata, x0, a3, c4, a1);
}


void Task42::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I33
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-2,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task43::Task_local::compute() {
  const Index x0 = b(0);
  const Index c2 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);

  // tensor label: I18
  std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a1.size());

    // tensor label: I35
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4, c2, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4, c2, a3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size(), c2.size(), a3.size());

    dgemm_("T", "N", a1.size(), x0.size()*c2.size()*a3.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, a1.size());
  }

  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), c2.size(), a3.size());
  out()->put_block(odata, x0, c2, a3, a1);
}


void Task44::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);

  // tensor label: I35
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c2, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c2, a3), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a3.size());

    // tensor label: I36
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a4.size()*c2.size()*a3.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a4.size()*c2.size()*a3.size());
  }

  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c2.size(), a3.size(), x0.size());
  out()->put_block(odata, x0, a4, c2, a3);
}


void Task45::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I36
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,2,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task46::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);

  // tensor label: I35
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c2, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c2, a3), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a4.size());

    // tensor label: I39
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a3.size()*c2.size()*a4.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a4.size());
  }

  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a4.size(), x0.size());
  out()->put_block(odata, x0, a4, c2, a3);
}


void Task47::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I39
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task48::Task_local::compute() {
  const Index x0 = b(0);
  const Index c2 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);

  // tensor label: I18
  std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a3)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a3.size());

    // tensor label: I41
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4, c2, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4, c2, a1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size(), c2.size(), a1.size());

    dgemm_("T", "N", a3.size(), x0.size()*c2.size()*a1.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, a3.size());
  }

  sort_indices<1,2,0,3,1,1,1,1>(odata_sorted, odata, a3.size(), x0.size(), c2.size(), a1.size());
  out()->put_block(odata, x0, c2, a3, a1);
}


void Task49::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);

  // tensor label: I41
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c2, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a1.size());

    // tensor label: I42
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a4.size()*c2.size()*a1.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a4.size()*c2.size()*a1.size());
  }

  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c2.size(), a1.size(), x0.size());
  out()->put_block(odata, x0, a4, c2, a1);
}


void Task50::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I42
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task51::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);

  // tensor label: I41
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c2, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a4.size());

    // tensor label: I45
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a1.size()*c2.size()*a4.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a4.size());
  }

  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a4.size(), x0.size());
  out()->put_block(odata, x0, a4, c2, a1);
}


void Task52::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I45
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,2,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task53::Task_local::compute() {
  const Index x0 = b(0);
  const Index c2 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);

  // tensor label: I18
  std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), x1.size());

    // tensor label: I47
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, a1, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, a1, a3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), a1.size(), a3.size());

    dgemm_("T", "N", c2.size(), x0.size()*a1.size()*a3.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c2.size());
  }

  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), x0.size(), a1.size(), a3.size());
  out()->put_block(odata, x0, c2, a3, a1);
}


void Task54::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a3 = b(3);

  // tensor label: I47
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a3), 0.0);

  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a3.size());

      // tensor label: I48
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


void Task55::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);

  // tensor label: I48
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma14
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task56::Task_local::compute() {
  const Index x0 = b(0);
  const Index c2 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);

  // tensor label: I18
  std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a1.size());

    // tensor label: I60
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }

  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size());
  out()->put_block(odata, x0, c2, a3, a1);
}


void Task57::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I60
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    dscal_(x1.size()*x0.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task58::Task_local::compute() {
  const Index x0 = b(0);
  const Index c2 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);

  // tensor label: I18
  std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a3.size());

    // tensor label: I62
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }

  sort_indices<3,1,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size());
  out()->put_block(odata, x0, c2, a3, a1);
}


void Task59::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I62
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    dscal_(x1.size()*x0.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,1,1,-2,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task60::Task_local::compute() {
  const Index x0 = b(0);
  const Index c2 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);

  // tensor label: I18
  std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a1.size());

    // tensor label: I68
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }

  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size());
  out()->put_block(odata, x0, c2, a3, a1);
}


void Task61::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I68
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-2,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task62::Task_local::compute() {
  const Index x0 = b(0);
  const Index c2 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);

  // tensor label: I18
  std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a3.size());

    // tensor label: I70
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }

  sort_indices<3,1,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size());
  out()->put_block(odata, x0, c2, a3, a1);
}


void Task63::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I70
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,4,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task64::Task_local::compute() {
  const Index x1 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);

  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(x1, a2, x0, a1);
  {
    // tensor label: I49
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, a1, a2);
    sort_indices<1,3,0,2,1,1,1,1>(i0data, odata, x0.size(), x1.size(), a1.size(), a2.size());
  }
  {
    // tensor label: I49
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0, a2, a1);
    sort_indices<0,2,1,3,1,1,1,1>(i0data, odata, x1.size(), x0.size(), a2.size(), a1.size());
  }
  out()->put_block(odata, x1, a2, x0, a1);
}


void Task65::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);

  // tensor label: I49
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

  for (auto& x2 : *range_[1]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, c3)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), c3.size());

      // tensor label: I50
      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x2, x1, a1, c3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x2, x1, a1, c3, a2)]);
      sort_indices<1,4,0,2,3,5,0,1,1,1>(i1data, i1data_sorted, x0.size(), x2.size(), x1.size(), a1.size(), c3.size(), a2.size());

      dgemm_("T", "N", 1, x0.size()*x1.size()*a1.size()*a2.size(), x2.size()*c3.size(),
             1.0, i0data_sorted, x2.size()*c3.size(), i1data_sorted, x2.size()*c3.size(),
             1.0, odata_sorted, 1);
    }
  }

  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), a2.size());
  out()->put_block(odata, x0, x1, a1, a2);
}


void Task66::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  const Index c3 = b(4);
  const Index a2 = b(5);

  // tensor label: I50
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a1, c3, a2), 0.0);

  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c3, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c3.size(), a2.size());

    // tensor label: I51
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

    dgemm_("T", "N", a1.size()*c3.size()*a2.size(), x0.size()*x2.size()*x1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a1.size()*c3.size()*a2.size());
  }

  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, a1.size(), c3.size(), a2.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a1, c3, a2);
}


void Task67::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);

  // tensor label: I51
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma14
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task68::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);

  // tensor label: I49
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

  for (auto& a3 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), a2.size());

    // tensor label: I55
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


void Task69::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a3 = b(3);

  // tensor label: I55
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a3), 0.0);

  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a3.size());

      // tensor label: I56
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


void Task70::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);

  // tensor label: I56
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma14
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task71::Task_local::compute() {
  const Index x1 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);

  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(x1, a2, x0, a1);
  {
    // tensor label: I52
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, a1, a2);
    sort_indices<1,3,0,2,1,1,1,1>(i0data, odata, x0.size(), x1.size(), a1.size(), a2.size());
  }
  out()->put_block(odata, x1, a2, x0, a1);
}


void Task72::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);

  // tensor label: I52
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a2.size());

      // tensor label: I53
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


void Task73::Task_local::compute() {
  const Index x5 = b(0);
  const Index x0 = b(1);
  const Index x4 = b(2);
  const Index x1 = b(3);

  // tensor label: I53
  std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x1);
  {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x1);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x5.size(), x0.size(), x4.size(), x1.size());
  }
  out()->put_block(odata, x5, x0, x4, x1);
}


void Task74::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);

  // tensor label: I52
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

      // tensor label: I64
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


void Task75::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);

  // tensor label: I64
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma14
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    dscal_(x3.size()*x0.size()*x2.size()*x1.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task76::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);

  // tensor label: I52
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

      // tensor label: I72
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


void Task77::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);

  // tensor label: I72
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma14
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task78::Task_local::compute() {
  energy_ = 0.0;
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index a4 = b(3);

  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
  sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());

  // tensor label: I74
  std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a4, c3, a2)]);
  sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), a4.size(), c3.size(), a2.size());

  energy_ += ddot_(c1.size()*a4.size()*c3.size()*a2.size(), i0data_sorted, 1, i1data_sorted, 1);
}


void Task79::Task_local::compute() {
  energy_ = 0.0;
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);

  // tensor label: I74
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
    dscal_(c1.size()*a4.size()*c3.size()*a2.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a4);
    dscal_(c1.size()*a2.size()*c3.size()*a4.size(), e0_, i1data.get(), 1);
    sort_indices<0,3,2,1,1,1,-2,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i2data = in(1)->get_block(c1, a4, c3, a2);
    sort_indices<0,1,2,3,1,1,-2,1>(i2data, odata, c1.size(), a4.size(), c3.size(), a2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i3data = in(1)->get_block(c1, a2, c3, a4);
    sort_indices<0,3,2,1,1,1,4,1>(i3data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->put_block(odata, c1, a4, c3, a2);
}


void Task80::Task_local::compute() {
  energy_ = 0.0;
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);

  // tensor label: I74
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, c3, a2), 0.0);

  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());

  // tensor label: I75
  std::unique_ptr<double[]> i1data = in(1)->get_block();
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size()]);
  sort_indices<0,1,1,1>(i1data, i1data_sorted);

  dgemm_("T", "N", c1.size()*a4.size()*c3.size()*a2.size(), 1, 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a4.size()*c3.size()*a2.size());

  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->put_block(odata, c1, a4, c3, a2);
}


void Task81::Task_local::compute() {
  energy_ = 0.0;

  // tensor label: I75
  std::unique_ptr<double[]> odata = out()->move_block();
  {
    // scalar
    // tensor label: Gamma0
    std::unique_ptr<double[]> i0data = in(0)->get_block();
    sort_indices<1,1,-1,1>(i0data, odata);
  }
  out()->put_block(odata);
}


void Task82::Task_local::compute() {
  energy_ = 0.0;
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);

  // tensor label: I74
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, c3, a2), 0.0);

  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());

  // tensor label: I78
  std::unique_ptr<double[]> i1data = in(1)->get_block();
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size()]);
  sort_indices<0,1,1,1>(i1data, i1data_sorted);

  dgemm_("T", "N", c1.size()*a2.size()*c3.size()*a4.size(), 1, 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a2.size()*c3.size()*a4.size());

  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->put_block(odata, c1, a4, c3, a2);
}


void Task83::Task_local::compute() {
  energy_ = 0.0;

  // tensor label: I78
  std::unique_ptr<double[]> odata = out()->move_block();
  {
    // scalar
    // tensor label: Gamma0
    std::unique_ptr<double[]> i0data = in(0)->get_block();
    sort_indices<1,1,2,1>(i0data, odata);
  }
  out()->put_block(odata);
}


void Task84::Task_local::compute() {
  energy_ = 0.0;
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);

  // tensor label: I74
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, c3, a2), 0.0);

  for (auto& c5 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c5);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, c5)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c3.size(), c5.size());

    // tensor label: I81
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a4, c5, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a4, c5, a2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a4.size(), c5.size(), a2.size());

    dgemm_("T", "N", c3.size(), c1.size()*a4.size()*a2.size(), c5.size(),
           1.0, i0data_sorted, c5.size(), i1data_sorted, c5.size(),
           1.0, odata_sorted, c3.size());
  }

  sort_indices<1,2,0,3,1,1,1,1>(odata_sorted, odata, c3.size(), c1.size(), a4.size(), a2.size());
  out()->put_block(odata, c1, a4, c3, a2);
}


void Task85::Task_local::compute() {
  energy_ = 0.0;
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c5 = b(2);
  const Index a2 = b(3);

  // tensor label: I81
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c5, a2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c5, a2);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, c1.size(), a4.size(), c5.size(), a2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c5, a4);
    sort_indices<0,3,2,1,1,1,-4,1>(i1data, odata, c1.size(), a2.size(), c5.size(), a4.size());
  }
  out()->put_block(odata, c1, a4, c5, a2);
}


void Task86::Task_local::compute() {
  energy_ = 0.0;
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);

  // tensor label: I74
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, c3, a2), 0.0);

  for (auto& a5 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a5, a4)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a5.size(), a4.size());

    // tensor label: I87
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a5, c3, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a5, c3, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a5.size(), c3.size(), a2.size());

    dgemm_("T", "N", a4.size(), c1.size()*c3.size()*a2.size(), a5.size(),
           1.0, i0data_sorted, a5.size(), i1data_sorted, a5.size(),
           1.0, odata_sorted, a4.size());
  }

  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), c3.size(), a2.size());
  out()->put_block(odata, c1, a4, c3, a2);
}


void Task87::Task_local::compute() {
  energy_ = 0.0;
  const Index c1 = b(0);
  const Index a5 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);

  // tensor label: I87
  std::unique_ptr<double[]> odata = out()->move_block(c1, a5, c3, a2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a5, c3, a2);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c1.size(), a5.size(), c3.size(), a2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a5);
    sort_indices<0,3,2,1,1,1,4,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a5.size());
  }
  out()->put_block(odata, c1, a5, c3, a2);
}


void Task88::Task_local::compute() {
  energy_ = 0.0;
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);

  // tensor label: I74
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, c3, a2), 0.0);

  for (auto& x0 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c3.size(), x0.size());

    // tensor label: I93
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4, c1, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4, c1, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size(), c1.size(), a2.size());

    dgemm_("T", "N", c3.size(), a4.size()*c1.size()*a2.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, c3.size());
  }

  sort_indices<2,1,0,3,1,1,1,1>(odata_sorted, odata, c3.size(), a4.size(), c1.size(), a2.size());
  out()->put_block(odata, c1, a4, c3, a2);
}


void Task89::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);

  // tensor label: I93
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c1, a2), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c1, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c1.size(), a2.size());

    // tensor label: I94
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a4.size()*c1.size()*a2.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a4.size()*c1.size()*a2.size());
  }

  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), a2.size(), x0.size());
  out()->put_block(odata, x0, a4, c1, a2);
}


void Task90::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I94
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task91::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);

  // tensor label: I93
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c1, a2), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, c1, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2, c1, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c1.size(), a4.size());

    // tensor label: I98
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a2.size()*c1.size()*a4.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a2.size()*c1.size()*a4.size());
  }

  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), a4.size(), x0.size());
  out()->put_block(odata, x0, a4, c1, a2);
}


void Task92::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I98
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,1,2>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task93::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a1 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);

  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
  sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());

  // tensor label: I100
  std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c2, a3, a1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c2, a3, a1)]);
  sort_indices<2,1,3,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), c2.size(), a3.size(), a1.size());

  energy_ += ddot_(x0.size()*c2.size()*a3.size()*a1.size(), i0data_sorted, 1, i1data_sorted, 1);
}


void Task94::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c2 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);

  // tensor label: I100
  std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, c4)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), c4.size());

      // tensor label: I101
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, c2, a3, c4, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, c2, a3, c4, a1)]);
      sort_indices<0,4,1,2,3,5,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), c2.size(), a3.size(), c4.size(), a1.size());

      dgemm_("T", "N", 1, x0.size()*c2.size()*a3.size()*a1.size(), x1.size()*c4.size(),
             1.0, i0data_sorted, x1.size()*c4.size(), i1data_sorted, x1.size()*c4.size(),
             1.0, odata_sorted, 1);
    }
  }

  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x0.size(), c2.size(), a3.size(), a1.size());
  out()->put_block(odata, x0, c2, a3, a1);
}


void Task95::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index c4 = b(4);
  const Index a1 = b(5);

  // tensor label: I101
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a3, c4, a1), 0.0);

  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c4, a1);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c4, a1)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c4.size(), a1.size());

  // tensor label: I102
  std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

  dgemm_("T", "N", c2.size()*a3.size()*c4.size()*a1.size(), x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a3.size()*c4.size()*a1.size());

  sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), c4.size(), a1.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c2, a3, c4, a1);
}


void Task96::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I102
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task97::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index c4 = b(4);
  const Index a1 = b(5);

  // tensor label: I101
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a3, c4, a1), 0.0);

  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c4, a3);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c4, a3)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c4.size(), a3.size());

  // tensor label: I106
  std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

  dgemm_("T", "N", c2.size()*a1.size()*c4.size()*a3.size(), x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a1.size()*c4.size()*a3.size());

  sort_indices<4,5,0,3,2,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), c4.size(), a3.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c2, a3, c4, a1);
}


void Task98::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I106
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,1,2>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task99::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c2 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);

  // tensor label: I100
  std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a1.size());

    // tensor label: I109
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size());

    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }

  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size());
  out()->put_block(odata, x0, c2, a3, a1);
}


void Task100::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x0 = b(1);

  // tensor label: I109
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0);
  {
    // tensor label: Gamma6
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0);
    sort_indices<0,1,1,1,-1,4>(i0data, odata, x3.size(), x0.size());
  }
  out()->put_block(odata, x3, x0);
}


void Task101::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c2 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);

  // tensor label: I100
  std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a3.size());

    // tensor label: I112
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size());

    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }

  sort_indices<3,1,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size());
  out()->put_block(odata, x0, c2, a3, a1);
}


void Task102::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x0 = b(1);

  // tensor label: I112
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0);
  {
    // tensor label: Gamma6
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0);
    sort_indices<0,1,1,1,1,2>(i0data, odata, x3.size(), x0.size());
  }
  out()->put_block(odata, x3, x0);
}


void Task103::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c2 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);

  // tensor label: I100
  std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

  for (auto& c4 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, c4)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), c4.size());

    // tensor label: I115
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a3, c4, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a3, c4, a1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a3.size(), c4.size(), a1.size());

    dgemm_("T", "N", c2.size(), x0.size()*a3.size()*a1.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c2.size());
  }

  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, c2.size(), x0.size(), a3.size(), a1.size());
  out()->put_block(odata, x0, c2, a3, a1);
}


void Task104::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c4 = b(2);
  const Index a1 = b(3);

  // tensor label: I115
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c4, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c4, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c4, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c4.size(), a1.size());

    // tensor label: I116
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a3.size()*c4.size()*a1.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a3.size()*c4.size()*a1.size());
  }

  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c4.size(), a1.size(), x0.size());
  out()->put_block(odata, x0, a3, c4, a1);
}


void Task105::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I116
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,1,4>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task106::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c4 = b(2);
  const Index a1 = b(3);

  // tensor label: I115
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c4, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c4, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c4.size(), a3.size());

    // tensor label: I120
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a1.size()*c4.size()*a3.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a1.size()*c4.size()*a3.size());
  }

  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c4.size(), a3.size(), x0.size());
  out()->put_block(odata, x0, a3, c4, a1);
}


void Task107::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I120
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-1,2>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task108::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c2 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);

  // tensor label: I100
  std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a1.size());

    // tensor label: I123
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4, c2, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4, c2, a3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size(), c2.size(), a3.size());

    dgemm_("T", "N", a1.size(), x0.size()*c2.size()*a3.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, a1.size());
  }

  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), c2.size(), a3.size());
  out()->put_block(odata, x0, c2, a3, a1);
}


void Task109::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a4 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);

  // tensor label: I123
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c2, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c2, a3), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a3.size());

    // tensor label: I124
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a4.size()*c2.size()*a3.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a4.size()*c2.size()*a3.size());
  }

  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c2.size(), a3.size(), x0.size());
  out()->put_block(odata, x0, a4, c2, a3);
}


void Task110::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I124
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,1,2>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task111::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a4 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);

  // tensor label: I123
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c2, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c2, a3), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a4.size());

    // tensor label: I128
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a3.size()*c2.size()*a4.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a4.size());
  }

  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a4.size(), x0.size());
  out()->put_block(odata, x0, a4, c2, a3);
}


void Task112::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I128
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-1,4>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task113::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c2 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);

  // tensor label: I100
  std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a3)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a3.size());

    // tensor label: I131
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4, c2, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4, c2, a1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size(), c2.size(), a1.size());

    dgemm_("T", "N", a3.size(), x0.size()*c2.size()*a1.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, a3.size());
  }

  sort_indices<1,2,0,3,1,1,1,1>(odata_sorted, odata, a3.size(), x0.size(), c2.size(), a1.size());
  out()->put_block(odata, x0, c2, a3, a1);
}


void Task114::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a4 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);

  // tensor label: I131
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c2, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a1.size());

    // tensor label: I132
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a4.size()*c2.size()*a1.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a4.size()*c2.size()*a1.size());
  }

  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c2.size(), a1.size(), x0.size());
  out()->put_block(odata, x0, a4, c2, a1);
}


void Task115::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I132
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-1,4>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task116::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a4 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);

  // tensor label: I131
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c2, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a4.size());

    // tensor label: I136
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a1.size()*c2.size()*a4.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a4.size());
  }

  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a4.size(), x0.size());
  out()->put_block(odata, x0, a4, c2, a1);
}


void Task117::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I136
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,1,2>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task118::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c2 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);

  // tensor label: I100
  std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), x1.size());

    // tensor label: I139
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, a1, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, a1, a3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), a1.size(), a3.size());

    dgemm_("T", "N", c2.size(), x0.size()*a1.size()*a3.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c2.size());
  }

  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), x0.size(), a1.size(), a3.size());
  out()->put_block(odata, x0, c2, a3, a1);
}


void Task119::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a3 = b(3);

  // tensor label: I139
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a3), 0.0);

  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a3.size());

      // tensor label: I140
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


void Task120::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);

  // tensor label: I140
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma14
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task121::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c2 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);

  // tensor label: I100
  std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a1.size());

    // tensor label: I158
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }

  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size());
  out()->put_block(odata, x0, c2, a3, a1);
}


void Task122::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I158
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    dscal_(x1.size()*x0.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,1,1,1,4>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task123::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c2 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);

  // tensor label: I100
  std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a3.size());

    // tensor label: I161
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }

  sort_indices<3,1,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size());
  out()->put_block(odata, x0, c2, a3, a1);
}


void Task124::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I161
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    dscal_(x1.size()*x0.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,1,1,-1,2>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task125::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c2 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);

  // tensor label: I100
  std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a1.size());

    // tensor label: I171
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }

  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size());
  out()->put_block(odata, x0, c2, a3, a1);
}


void Task126::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I171
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task127::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c2 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);

  // tensor label: I100
  std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a3.size());

    // tensor label: I174
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }

  sort_indices<3,1,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size());
  out()->put_block(odata, x0, c2, a3, a1);
}


void Task128::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I174
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,2,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task129::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a1 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);

  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
  sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());

  // tensor label: I142
  std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, a1, a2);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, a1, a2)]);
  sort_indices<3,1,2,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), a1.size(), a2.size());

  energy_ += ddot_(x0.size()*x1.size()*a1.size()*a2.size(), i0data_sorted, 1, i1data_sorted, 1);
}


void Task130::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);

  // tensor label: I142
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

  for (auto& x2 : *range_[1]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, c3)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), c3.size());

      // tensor label: I143
      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x2, x1, a1, c3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x2, x1, a1, c3, a2)]);
      sort_indices<1,4,0,2,3,5,0,1,1,1>(i1data, i1data_sorted, x0.size(), x2.size(), x1.size(), a1.size(), c3.size(), a2.size());

      dgemm_("T", "N", 1, x0.size()*x1.size()*a1.size()*a2.size(), x2.size()*c3.size(),
             1.0, i0data_sorted, x2.size()*c3.size(), i1data_sorted, x2.size()*c3.size(),
             1.0, odata_sorted, 1);
    }
  }

  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), a2.size());
  out()->put_block(odata, x0, x1, a1, a2);
}


void Task131::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  const Index c3 = b(4);
  const Index a2 = b(5);

  // tensor label: I143
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a1, c3, a2), 0.0);

  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c3, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c3.size(), a2.size());

    // tensor label: I144
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

    dgemm_("T", "N", a1.size()*c3.size()*a2.size(), x0.size()*x2.size()*x1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a1.size()*c3.size()*a2.size());
  }

  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, a1.size(), c3.size(), a2.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a1, c3, a2);
}


void Task132::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);

  // tensor label: I144
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma14
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task133::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);

  // tensor label: I142
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a2.size());

      // tensor label: I147
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


void Task134::Task_local::compute() {
  energy_ = 0.0;
  const Index x5 = b(0);
  const Index x0 = b(1);
  const Index x4 = b(2);
  const Index x1 = b(3);

  // tensor label: I147
  std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x1);
  {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x1);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x5.size(), x0.size(), x4.size(), x1.size());
  }
  out()->put_block(odata, x5, x0, x4, x1);
}


void Task135::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);

  // tensor label: I142
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

  for (auto& a3 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), a2.size());

    // tensor label: I150
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


void Task136::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a3 = b(3);

  // tensor label: I150
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a3), 0.0);

  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a3.size());

      // tensor label: I151
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


void Task137::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);

  // tensor label: I151
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma14
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task138::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);

  // tensor label: I142
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

      // tensor label: I164
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


void Task139::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);

  // tensor label: I164
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma14
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    dscal_(x3.size()*x0.size()*x2.size()*x1.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task140::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);

  // tensor label: I142
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

      // tensor label: I177
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


void Task141::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);

  // tensor label: I177
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma14
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task142::Task_local::compute(){
  correction_ = 0.0;
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index a4 = b(3);

  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
  sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());

  // tensor label: I179
  std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a4, c3, a2)]);
  sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), a4.size(), c3.size(), a2.size());

  correction_ += ddot_(c1.size()*a4.size()*c3.size()*a2.size(), i0data_sorted, 1, i1data_sorted, 1);
}


void Task143::Task_local::compute(){
  correction_ = 0.0;
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);

  // tensor label: I179
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<0,3,2,1,1,1,2,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->put_block(odata, c1, a4, c3, a2);
}


void Task144::Task_local::compute(){
  correction_ = 0.0;
  const Index x0 = b(0);
  const Index a1 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);

  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
  sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());

  // tensor label: I183
  std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a3, c2, a1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a3, c2, a1)]);
  sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), a3.size(), c2.size(), a1.size());

  correction_ += ddot_(x0.size()*a3.size()*c2.size()*a1.size(), i0data_sorted, 1, i1data_sorted, 1);
}


void Task145::Task_local::compute(){
  correction_ = 0.0;
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);

  // tensor label: I183
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c2, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a1.size());

    // tensor label: I184
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }

  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size());
  out()->put_block(odata, x0, a3, c2, a1);
}


void Task146::Task_local::compute(){
  correction_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I184
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-1,4>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task147::Task_local::compute(){
  correction_ = 0.0;
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);

  // tensor label: I183
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c2, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a3.size());

    // tensor label: I187
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }

  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size());
  out()->put_block(odata, x0, a3, c2, a1);
}


void Task148::Task_local::compute(){
  correction_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I187
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,1,2>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task149::Task_local::compute(){
  correction_ = 0.0;
  const Index x0 = b(0);
  const Index a1 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);

  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
  sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());

  // tensor label: I189
  std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, a1, a2);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, a1, a2)]);
  sort_indices<3,1,2,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), a1.size(), a2.size());

  correction_ += ddot_(x0.size()*x1.size()*a1.size()*a2.size(), i0data_sorted, 1, i1data_sorted, 1);
}


void Task150::Task_local::compute(){
  correction_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);

  // tensor label: I189
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

      // tensor label: I190
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


void Task151::Task_local::compute(){
  correction_ = 0.0;
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);

  // tensor label: I190
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma14
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task153::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);

  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: I191
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<1,0,1,1,1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x0, x1);
}


void Task154::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I191
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);

  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        for (auto& a4 : *range_[2]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());

          // tensor label: I192
          std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, c1, a4, c3, a2);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, c1, a4, c3, a2)]);
          sort_indices<2,5,4,3,0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), c1.size(), a4.size(), c3.size(), a2.size());

          dgemm_("T", "N", 1, x1.size()*x0.size(), c1.size()*a4.size()*c3.size()*a2.size(),
                 1.0, i0data_sorted, c1.size()*a4.size()*c3.size()*a2.size(), i1data_sorted, c1.size()*a4.size()*c3.size()*a2.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }

  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x1, x0);
}


void Task155::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a4 = b(3);
  const Index c3 = b(4);
  const Index a2 = b(5);

  // tensor label: I192
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a4, c3, a2), 0.0);

  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());

  // tensor label: I193
  std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

  dgemm_("T", "N", c1.size()*a4.size()*c3.size()*a2.size(), x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a4.size()*c3.size()*a2.size());

  sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c1, a4, c3, a2);
}


void Task156::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I193
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task157::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a4 = b(3);
  const Index c3 = b(4);
  const Index a2 = b(5);

  // tensor label: I192
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a4, c3, a2), 0.0);

  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());

  // tensor label: I196
  std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

  dgemm_("T", "N", c1.size()*a2.size()*c3.size()*a4.size(), x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a2.size()*c3.size()*a4.size());

  sort_indices<4,5,0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c1, a4, c3, a2);
}


void Task158::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I196
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,2,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task159::Task_local::compute() {
  const Index c5 = b(0);
  const Index c3 = b(1);

  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(c5, c3);
  {
    // tensor label: I197
    std::unique_ptr<double[]> i0data = in(0)->get_block(c5, c3);
    sort_indices<0,1,1,1,1,1>(i0data, odata, c5.size(), c3.size());
  }
  out()->put_block(odata, c5, c3);
}


void Task160::Task_local::compute() {
  const Index c5 = b(0);
  const Index c3 = b(1);

  // tensor label: I197
  std::unique_ptr<double[]> odata = out()->move_block(c5, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c5, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c5, c3), 0.0);

  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& a4 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());

        // tensor label: I198
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a4, c5, a2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a4, c5, a2)]);
        sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), a4.size(), c5.size(), a2.size());

        dgemm_("T", "N", c3.size(), c5.size(), c1.size()*a4.size()*a2.size(),
               1.0, i0data_sorted, c1.size()*a4.size()*a2.size(), i1data_sorted, c1.size()*a4.size()*a2.size(),
               1.0, odata_sorted, c3.size());
      }
    }
  }

  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c3.size(), c5.size());
  out()->put_block(odata, c5, c3);
}


void Task161::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c5 = b(2);
  const Index a2 = b(3);

  // tensor label: I198
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c5, a2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c5, a2);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, c1.size(), a4.size(), c5.size(), a2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c5, a4);
    sort_indices<0,3,2,1,1,1,-4,1>(i1data, odata, c1.size(), a2.size(), c5.size(), a4.size());
  }
  out()->put_block(odata, c1, a4, c5, a2);
}


void Task162::Task_local::compute() {
  const Index a4 = b(0);
  const Index a5 = b(1);

  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(a4, a5);
  {
    // tensor label: I201
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, a4);
    sort_indices<1,0,1,1,1,1>(i0data, odata, a5.size(), a4.size());
  }
  out()->put_block(odata, a4, a5);
}


void Task163::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);

  // tensor label: I201
  std::unique_ptr<double[]> odata = out()->move_block(a5, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a5, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a5, a4), 0.0);

  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());

        // tensor label: I202
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a5, c3, a2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a5, c3, a2)]);
        sort_indices<0,3,2,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), a5.size(), c3.size(), a2.size());

        dgemm_("T", "N", a4.size(), a5.size(), c1.size()*c3.size()*a2.size(),
               1.0, i0data_sorted, c1.size()*c3.size()*a2.size(), i1data_sorted, c1.size()*c3.size()*a2.size(),
               1.0, odata_sorted, a4.size());
      }
    }
  }

  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a4.size(), a5.size());
  out()->put_block(odata, a5, a4);
}


void Task164::Task_local::compute() {
  const Index c1 = b(0);
  const Index a5 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);

  // tensor label: I202
  std::unique_ptr<double[]> odata = out()->move_block(c1, a5, c3, a2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a5, c3, a2);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c1.size(), a5.size(), c3.size(), a2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a5);
    sort_indices<0,3,2,1,1,1,4,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a5.size());
  }
  out()->put_block(odata, c1, a5, c3, a2);
}


void Task165::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);

  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3);
  {
    // tensor label: I205
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, c3);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x0.size(), c3.size());
  }
  out()->put_block(odata, x0, c3);
}


void Task166::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);

  // tensor label: I205
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c3), 0.0);

  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& a4 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());

        // tensor label: I206
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4, c1, a2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4, c1, a2)]);
        sort_indices<2,3,1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size(), c1.size(), a2.size());

        dgemm_("T", "N", c3.size(), x0.size(), a4.size()*c1.size()*a2.size(),
               1.0, i0data_sorted, a4.size()*c1.size()*a2.size(), i1data_sorted, a4.size()*c1.size()*a2.size(),
               1.0, odata_sorted, c3.size());
      }
    }
  }

  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c3.size(), x0.size());
  out()->put_block(odata, x0, c3);
}


void Task167::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);

  // tensor label: I206
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c1, a2), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c1, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c1.size(), a2.size());

    // tensor label: I207
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a4.size()*c1.size()*a2.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a4.size()*c1.size()*a2.size());
  }

  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), a2.size(), x0.size());
  out()->put_block(odata, x0, a4, c1, a2);
}


void Task168::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I207
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task169::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);

  // tensor label: I206
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c1, a2), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, c1, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2, c1, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c1.size(), a4.size());

    // tensor label: I210
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a2.size()*c1.size()*a4.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a2.size()*c1.size()*a4.size());
  }

  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), a4.size(), x0.size());
  out()->put_block(odata, x0, a4, c1, a2);
}


void Task170::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I210
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,1,2>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task171::Task_local::compute() {
  const Index c4 = b(0);
  const Index x1 = b(1);

  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(c4, x1);
  {
    // tensor label: I211
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c4);
    sort_indices<1,0,1,1,1,1>(i0data, odata, x1.size(), c4.size());
  }
  out()->put_block(odata, c4, x1);
}


void Task172::Task_local::compute() {
  const Index x1 = b(0);
  const Index c4 = b(1);

  // tensor label: I211
  std::unique_ptr<double[]> odata = out()->move_block(x1, c4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, c4), 0.0);

  for (auto& x0 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& c2 : *range_[0]) {
        for (auto& a3 : *range_[2]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());

          // tensor label: I212
          std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, c2, a3, c4, a1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, c2, a3, c4, a1)]);
          sort_indices<1,5,2,3,0,4,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), c2.size(), a3.size(), c4.size(), a1.size());

          dgemm_("T", "N", 1, x1.size()*c4.size(), x0.size()*c2.size()*a3.size()*a1.size(),
                 1.0, i0data_sorted, x0.size()*c2.size()*a3.size()*a1.size(), i1data_sorted, x0.size()*c2.size()*a3.size()*a1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }

  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), c4.size());
  out()->put_block(odata, x1, c4);
}


void Task173::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index c4 = b(4);
  const Index a1 = b(5);

  // tensor label: I212
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a3, c4, a1), 0.0);

  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c4, a1);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c4, a1)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c4.size(), a1.size());

  // tensor label: I213
  std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

  dgemm_("T", "N", c2.size()*a3.size()*c4.size()*a1.size(), x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a3.size()*c4.size()*a1.size());

  sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), c4.size(), a1.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c2, a3, c4, a1);
}


void Task174::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I213
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task175::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index c4 = b(4);
  const Index a1 = b(5);

  // tensor label: I212
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a3, c4, a1), 0.0);

  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c4, a3);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c4, a3)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c4.size(), a3.size());

  // tensor label: I216
  std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

  dgemm_("T", "N", c2.size()*a1.size()*c4.size()*a3.size(), x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a1.size()*c4.size()*a3.size());

  sort_indices<4,5,0,3,2,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), c4.size(), a3.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c2, a3, c4, a1);
}


void Task176::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I216
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,1,2>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task177::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);

  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(x1, x2);
  {
    // tensor label: I217
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x1);
    sort_indices<1,0,1,1,1,1>(i0data, odata, x2.size(), x1.size());
  }
  out()->put_block(odata, x1, x2);
}


void Task178::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);

  // tensor label: I217
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1), 0.0);

  for (auto& x0 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& c2 : *range_[0]) {
        for (auto& a3 : *range_[2]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());

          // tensor label: I218
          std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x2, x1, a3, c2, a1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x2, x1, a3, c2, a1)]);
          sort_indices<0,5,4,3,1,2,0,1,1,1>(i1data, i1data_sorted, x0.size(), x2.size(), x1.size(), a3.size(), c2.size(), a1.size());

          dgemm_("T", "N", 1, x2.size()*x1.size(), x0.size()*a3.size()*c2.size()*a1.size(),
                 1.0, i0data_sorted, x0.size()*a3.size()*c2.size()*a1.size(), i1data_sorted, x0.size()*a3.size()*c2.size()*a1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }

  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size());
  out()->put_block(odata, x2, x1);
}


void Task179::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a3 = b(3);
  const Index c2 = b(4);
  const Index a1 = b(5);

  // tensor label: I218
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a3, c2, a1), 0.0);

  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a1.size());

    // tensor label: I219
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), x0.size()*x2.size()*x1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }

  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a3, c2, a1);
}


void Task180::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);

  // tensor label: I219
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma14
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,4>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task181::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a3 = b(3);
  const Index c2 = b(4);
  const Index a1 = b(5);

  // tensor label: I218
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a3, c2, a1), 0.0);

  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a3.size());

    // tensor label: I222
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), x0.size()*x2.size()*x1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }

  sort_indices<3,4,5,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a3, c2, a1);
}


void Task182::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);

  // tensor label: I222
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma14
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task183::Task_local::compute() {
  const Index c4 = b(0);
  const Index c2 = b(1);

  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(c4, c2);
  {
    // tensor label: I223
    std::unique_ptr<double[]> i0data = in(0)->get_block(c4, c2);
    sort_indices<0,1,1,1,1,1>(i0data, odata, c4.size(), c2.size());
  }
  out()->put_block(odata, c4, c2);
}


void Task184::Task_local::compute() {
  const Index c4 = b(0);
  const Index c2 = b(1);

  // tensor label: I223
  std::unique_ptr<double[]> odata = out()->move_block(c4, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c4, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c4, c2), 0.0);

  for (auto& x0 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());

        // tensor label: I224
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a3, c4, a1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a3, c4, a1)]);
        sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x0.size(), a3.size(), c4.size(), a1.size());

        dgemm_("T", "N", c2.size(), c4.size(), x0.size()*a3.size()*a1.size(),
               1.0, i0data_sorted, x0.size()*a3.size()*a1.size(), i1data_sorted, x0.size()*a3.size()*a1.size(),
               1.0, odata_sorted, c2.size());
      }
    }
  }

  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c2.size(), c4.size());
  out()->put_block(odata, c4, c2);
}


void Task185::Task_local::compute() {
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c4 = b(2);
  const Index a1 = b(3);

  // tensor label: I224
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c4, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c4, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c4, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c4.size(), a1.size());

    // tensor label: I225
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a3.size()*c4.size()*a1.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a3.size()*c4.size()*a1.size());
  }

  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c4.size(), a1.size(), x0.size());
  out()->put_block(odata, x0, a3, c4, a1);
}


void Task186::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I225
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,1,4>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task187::Task_local::compute() {
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c4 = b(2);
  const Index a1 = b(3);

  // tensor label: I224
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c4, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c4, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c4.size(), a3.size());

    // tensor label: I228
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a1.size()*c4.size()*a3.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a1.size()*c4.size()*a3.size());
  }

  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c4.size(), a3.size(), x0.size());
  out()->put_block(odata, x0, a3, c4, a1);
}


void Task188::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I228
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-1,2>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task189::Task_local::compute() {
  const Index a1 = b(0);
  const Index a4 = b(1);

  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(a1, a4);
  {
    // tensor label: I229
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1);
    sort_indices<1,0,1,1,1,1>(i0data, odata, a4.size(), a1.size());
  }
  out()->put_block(odata, a1, a4);
}


void Task190::Task_local::compute() {
  const Index a4 = b(0);
  const Index a1 = b(1);

  // tensor label: I229
  std::unique_ptr<double[]> odata = out()->move_block(a4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a1), 0.0);

  for (auto& x0 : *range_[1]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());

        // tensor label: I230
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4, c2, a3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4, c2, a3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size(), c2.size(), a3.size());

        dgemm_("T", "N", a1.size(), a4.size(), x0.size()*c2.size()*a3.size(),
               1.0, i0data_sorted, x0.size()*c2.size()*a3.size(), i1data_sorted, x0.size()*c2.size()*a3.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }

  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a1.size(), a4.size());
  out()->put_block(odata, a4, a1);
}


void Task191::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);

  // tensor label: I230
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c2, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c2, a3), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a3.size());

    // tensor label: I231
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a4.size()*c2.size()*a3.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a4.size()*c2.size()*a3.size());
  }

  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c2.size(), a3.size(), x0.size());
  out()->put_block(odata, x0, a4, c2, a3);
}


void Task192::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I231
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,1,2>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task193::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);

  // tensor label: I230
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c2, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c2, a3), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a4.size());

    // tensor label: I234
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a3.size()*c2.size()*a4.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a4.size());
  }

  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a4.size(), x0.size());
  out()->put_block(odata, x0, a4, c2, a3);
}


void Task194::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I234
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-1,4>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task195::Task_local::compute() {
  const Index a3 = b(0);
  const Index a4 = b(1);

  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(a3, a4);
  {
    // tensor label: I235
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a3);
    sort_indices<1,0,1,1,1,1>(i0data, odata, a4.size(), a3.size());
  }
  out()->put_block(odata, a3, a4);
}


void Task196::Task_local::compute() {
  const Index a4 = b(0);
  const Index a3 = b(1);

  // tensor label: I235
  std::unique_ptr<double[]> odata = out()->move_block(a4, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a3), 0.0);

  for (auto& x0 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& c2 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());

        // tensor label: I236
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4, c2, a1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4, c2, a1)]);
        sort_indices<0,3,2,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size(), c2.size(), a1.size());

        dgemm_("T", "N", a3.size(), a4.size(), x0.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, x0.size()*c2.size()*a1.size(), i1data_sorted, x0.size()*c2.size()*a1.size(),
               1.0, odata_sorted, a3.size());
      }
    }
  }

  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a3.size(), a4.size());
  out()->put_block(odata, a4, a3);
}


void Task197::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);

  // tensor label: I236
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c2, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a1.size());

    // tensor label: I237
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a4.size()*c2.size()*a1.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a4.size()*c2.size()*a1.size());
  }

  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c2.size(), a1.size(), x0.size());
  out()->put_block(odata, x0, a4, c2, a1);
}


void Task198::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I237
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-1,4>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task199::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);

  // tensor label: I236
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c2, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a4.size());

    // tensor label: I240
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a1.size()*c2.size()*a4.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a4.size());
  }

  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a4.size(), x0.size());
  out()->put_block(odata, x0, a4, c2, a1);
}


void Task200::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I240
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,1,2>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task201::Task_local::compute() {
  const Index x1 = b(0);
  const Index c2 = b(1);

  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(x1, c2);
  {
    // tensor label: I241
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c2);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x1.size(), c2.size());
  }
  out()->put_block(odata, x1, c2);
}


void Task202::Task_local::compute() {
  const Index x1 = b(0);
  const Index c2 = b(1);

  // tensor label: I241
  std::unique_ptr<double[]> odata = out()->move_block(x1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, c2), 0.0);

  for (auto& x0 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());

        // tensor label: I242
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, a1, a3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, a1, a3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), a1.size(), a3.size());

        dgemm_("T", "N", c2.size(), x1.size(), x0.size()*a1.size()*a3.size(),
               1.0, i0data_sorted, x0.size()*a1.size()*a3.size(), i1data_sorted, x0.size()*a1.size()*a3.size(),
               1.0, odata_sorted, c2.size());
      }
    }
  }

  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c2.size(), x1.size());
  out()->put_block(odata, x1, c2);
}


void Task203::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a3 = b(3);

  // tensor label: I242
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a3), 0.0);

  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a3)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a3.size());

      // tensor label: I243
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


void Task204::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);

  // tensor label: I243
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma14
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task205::Task_local::compute() {
  const Index c3 = b(0);
  const Index x2 = b(1);

  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(c3, x2);
  {
    // tensor label: I244
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c3);
    sort_indices<1,0,1,1,1,1>(i0data, odata, x2.size(), c3.size());
  }
  out()->put_block(odata, c3, x2);
}


void Task206::Task_local::compute() {
  const Index x2 = b(0);
  const Index c3 = b(1);

  // tensor label: I244
  std::unique_ptr<double[]> odata = out()->move_block(x2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, c3), 0.0);

  for (auto& x0 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& x1 : *range_[1]) {
        for (auto& a2 : *range_[2]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());

          // tensor label: I245
          std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x2, x1, a1, c3, a2);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x2, x1, a1, c3, a2)]);
          sort_indices<0,3,2,5,1,4,0,1,1,1>(i1data, i1data_sorted, x0.size(), x2.size(), x1.size(), a1.size(), c3.size(), a2.size());

          dgemm_("T", "N", 1, x2.size()*c3.size(), x0.size()*x1.size()*a1.size()*a2.size(),
                 1.0, i0data_sorted, x0.size()*x1.size()*a1.size()*a2.size(), i1data_sorted, x0.size()*x1.size()*a1.size()*a2.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }

  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x2.size(), c3.size());
  out()->put_block(odata, x2, c3);
}


void Task207::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  const Index c3 = b(4);
  const Index a2 = b(5);

  // tensor label: I245
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a1, c3, a2), 0.0);

  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c3, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c3.size(), a2.size());

    // tensor label: I246
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

    dgemm_("T", "N", a1.size()*c3.size()*a2.size(), x0.size()*x2.size()*x1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a1.size()*c3.size()*a2.size());
  }

  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, a1.size(), c3.size(), a2.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a1, c3, a2);
}


void Task208::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);

  // tensor label: I246
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma14
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task209::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);

  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3);
  {
    // tensor label: I247
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2);
    sort_indices<1,0,1,1,1,1>(i0data, odata, x3.size(), x2.size());
  }
  out()->put_block(odata, x2, x3);
}


void Task210::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);

  // tensor label: I247
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

          // tensor label: I248
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


void Task211::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index a1 = b(4);
  const Index a2 = b(5);

  // tensor label: I248
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, x3, x2, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, x3, x2, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, x3, x2, a1, a2), 0.0);

  for (auto& x4 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a2)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a2.size());

      // tensor label: I249
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


void Task212::Task_local::compute() {
  const Index x5 = b(0);
  const Index x0 = b(1);
  const Index x4 = b(2);
  const Index x1 = b(3);
  const Index x3 = b(4);
  const Index x2 = b(5);

  // tensor label: I249
  std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x1, x3, x2);
  {
    // tensor label: Gamma67
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x1, x3, x2);
    sort_indices<0,1,2,3,4,5,1,1,1,2>(i0data, odata, x5.size(), x0.size(), x4.size(), x1.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, x5, x0, x4, x1, x3, x2);
}


void Task213::Task_local::compute() {
  const Index a2 = b(0);
  const Index a3 = b(1);

  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(a2, a3);
  {
    // tensor label: I250
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
    sort_indices<1,0,1,1,1,1>(i0data, odata, a3.size(), a2.size());
  }
  out()->put_block(odata, a2, a3);
}


void Task214::Task_local::compute() {
  const Index a3 = b(0);
  const Index a2 = b(1);

  // tensor label: I250
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

        // tensor label: I251
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


void Task215::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a3 = b(3);

  // tensor label: I251
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a3), 0.0);

  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a3)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a3.size());

      // tensor label: I252
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


void Task216::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);

  // tensor label: I252
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma14
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task218::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);

  // tensor label: Den1
  std::unique_ptr<double[]> odata = out()->move_block(c3, a4, c1, a2);
  {
    // tensor label: I253
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
    sort_indices<2,1,0,3,1,1,1,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
  }
  out()->put_block(odata, c3, a4, c1, a2);
}


void Task219::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);

  // tensor label: I253
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<0,3,2,1,1,1,4,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->put_block(odata, c1, a4, c3, a2);
}


void Task220::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);

  // tensor label: Den1
  std::unique_ptr<double[]> odata = out()->move_block(c2, a3, x0, a1);
  {
    // tensor label: I255
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c2, a1);
    sort_indices<2,1,0,3,1,1,1,1>(i0data, odata, x0.size(), a3.size(), c2.size(), a1.size());
  }
  out()->put_block(odata, c2, a3, x0, a1);
}


void Task221::Task_local::compute() {
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);

  // tensor label: I255
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c2, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a1.size());

    // tensor label: I256
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }

  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size());
  out()->put_block(odata, x0, a3, c2, a1);
}


void Task222::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I256
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-1,2>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task223::Task_local::compute() {
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);

  // tensor label: I255
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c2, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a3.size());

    // tensor label: I258
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }

  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size());
  out()->put_block(odata, x0, a3, c2, a1);
}


void Task224::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);

  // tensor label: I258
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}


void Task225::Task_local::compute() {
  const Index x1 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);

  // tensor label: Den1
  std::unique_ptr<double[]> odata = out()->move_block(x1, a2, x0, a1);
  {
    // tensor label: I259
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, a1, a2);
    sort_indices<1,3,0,2,1,1,1,1>(i0data, odata, x0.size(), x1.size(), a1.size(), a2.size());
  }
  out()->put_block(odata, x1, a2, x0, a1);
}


void Task226::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);

  // tensor label: I259
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

      // tensor label: I260
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


void Task227::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);

  // tensor label: I260
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma14
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task229::Task_local::compute() {
  const Index ci0 = b(0);

  // tensor label: deci
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  {
    // tensor label: I261
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0);
    sort_indices<0,1,1,1,1>(i0data, odata, ci0.size());
  }
  out()->put_block(odata, ci0);
}


void Task230::Task_local::compute() {
  const Index ci0 = b(0);

  // tensor label: I261
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);

  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        for (auto& a4 : *range_[2]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());

          // tensor label: I262
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, c1, a4, c3, a2);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, c1, a4, c3, a2)]);
          sort_indices<1,4,3,2,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), c1.size(), a4.size(), c3.size(), a2.size());

          dgemm_("T", "N", 1, ci0.size(), c1.size()*a4.size()*c3.size()*a2.size(),
                 1.0, i0data_sorted, c1.size()*a4.size()*c3.size()*a2.size(), i1data_sorted, c1.size()*a4.size()*c3.size()*a2.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }

  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}


void Task231::Task_local::compute() {
  const Index ci0 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  const Index a2 = b(4);

  // tensor label: I262
  std::unique_ptr<double[]> odata = out()->move_block(ci0, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, c1, a4, c3, a2), 0.0);

  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());

  // tensor label: I263
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0)]);
  sort_indices<0,0,1,1,1>(i1data, i1data_sorted, ci0.size());

  dgemm_("T", "N", c1.size()*a4.size()*c3.size()*a2.size(), ci0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a4.size()*c3.size()*a2.size());

  sort_indices<4,0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size(), ci0.size());
  out()->put_block(odata, ci0, c1, a4, c3, a2);
}


void Task232::Task_local::compute() {
  const Index ci0 = b(0);

  // tensor label: I263
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  {
    // tensor label: Gamma72
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0);
    sort_indices<0,1,1,-1,1>(i0data, odata, ci0.size());
  }
  {
    // tensor label: Gamma72
    std::unique_ptr<double[]> i1data = in(0)->get_block(ci0);
    sort_indices<0,1,1,-1,1>(i1data, odata, ci0.size());
  }
  out()->put_block(odata, ci0);
}


void Task233::Task_local::compute() {
  const Index ci0 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  const Index a2 = b(4);

  // tensor label: I262
  std::unique_ptr<double[]> odata = out()->move_block(ci0, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, c1, a4, c3, a2), 0.0);

  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());

  // tensor label: I266
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0)]);
  sort_indices<0,0,1,1,1>(i1data, i1data_sorted, ci0.size());

  dgemm_("T", "N", c1.size()*a2.size()*c3.size()*a4.size(), ci0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a2.size()*c3.size()*a4.size());

  sort_indices<4,0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size(), ci0.size());
  out()->put_block(odata, ci0, c1, a4, c3, a2);
}


void Task234::Task_local::compute() {
  const Index ci0 = b(0);

  // tensor label: I266
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  {
    // tensor label: Gamma72
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0);
    sort_indices<0,1,1,2,1>(i0data, odata, ci0.size());
  }
  {
    // tensor label: Gamma72
    std::unique_ptr<double[]> i1data = in(0)->get_block(ci0);
    sort_indices<0,1,1,2,1>(i1data, odata, ci0.size());
  }
  out()->put_block(odata, ci0);
}


void Task235::Task_local::compute() {
  const Index ci0 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  const Index a2 = b(4);

  // tensor label: I262
  std::unique_ptr<double[]> odata = out()->move_block(ci0, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, c1, a4, c3, a2), 0.0);

  for (auto& x0 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c3.size(), x0.size());

    // tensor label: I269
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, a4, c1, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, a4, c1, a2)]);
    sort_indices<1,0,2,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), a4.size(), c1.size(), a2.size());

    dgemm_("T", "N", c3.size(), ci0.size()*a4.size()*c1.size()*a2.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, c3.size());
  }

  sort_indices<1,3,2,0,4,1,1,1,1>(odata_sorted, odata, c3.size(), ci0.size(), a4.size(), c1.size(), a2.size());
  out()->put_block(odata, ci0, c1, a4, c3, a2);
}


void Task236::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a4 = b(2);
  const Index c1 = b(3);
  const Index a2 = b(4);

  // tensor label: I269
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a4, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a4, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a4, c1, a2), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c1, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c1.size(), a2.size());

    // tensor label: I270
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

    dgemm_("T", "N", a4.size()*c1.size()*a2.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a4.size()*c1.size()*a2.size());
  }

  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), a2.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, a4, c1, a2);
}


void Task237::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);

  // tensor label: I270
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma74
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}


void Task238::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a4 = b(2);
  const Index c1 = b(3);
  const Index a2 = b(4);

  // tensor label: I269
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a4, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a4, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a4, c1, a2), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, c1, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2, c1, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c1.size(), a4.size());

    // tensor label: I274
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

    dgemm_("T", "N", a2.size()*c1.size()*a4.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a2.size()*c1.size()*a4.size());
  }

  sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), a4.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, a4, c1, a2);
}


void Task239::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);

  // tensor label: I274
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma74
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}


void Task240::Task_local::compute() {
  const Index ci0 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  const Index a2 = b(4);

  // tensor label: I262
  std::unique_ptr<double[]> odata = out()->move_block(ci0, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, c1, a4, c3, a2), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c3.size(), x1.size());

    // tensor label: I336
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, a4, c1, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, a4, c1, a2)]);
    sort_indices<1,0,2,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), a4.size(), c1.size(), a2.size());

    dgemm_("T", "N", c3.size(), ci0.size()*a4.size()*c1.size()*a2.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c3.size());
  }

  sort_indices<1,3,2,0,4,1,1,1,1>(odata_sorted, odata, c3.size(), ci0.size(), a4.size(), c1.size(), a2.size());
  out()->put_block(odata, ci0, c1, a4, c3, a2);
}


void Task241::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index a4 = b(2);
  const Index c1 = b(3);
  const Index a2 = b(4);

  // tensor label: I336
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a4, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a4, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a4, c1, a2), 0.0);

  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a4, c1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a4, c1, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a4.size(), c1.size(), a2.size());

    // tensor label: I337
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

    dgemm_("T", "N", a4.size()*c1.size()*a2.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a4.size()*c1.size()*a2.size());
  }

  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), a2.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, a4, c1, a2);
}


void Task242::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);

  // tensor label: I337
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma74
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}


void Task243::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index a4 = b(2);
  const Index c1 = b(3);
  const Index a2 = b(4);

  // tensor label: I336
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a4, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a4, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a4, c1, a2), 0.0);

  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a2, c1, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a2, c1, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a2.size(), c1.size(), a4.size());

    // tensor label: I341
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

    dgemm_("T", "N", a2.size()*c1.size()*a4.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a2.size()*c1.size()*a4.size());
  }

  sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), a4.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, a4, c1, a2);
}


void Task244::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);

  // tensor label: I341
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma74
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}


void Task245::Task_local::compute() {
  const Index ci0 = b(0);

  // tensor label: I261
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);

  for (auto& x0 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& c2 : *range_[0]) {
        for (auto& a3 : *range_[2]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());

          // tensor label: I276
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, c2, a3, a1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, c2, a3, a1)]);
          sort_indices<1,4,2,3,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), c2.size(), a3.size(), a1.size());

          dgemm_("T", "N", 1, ci0.size(), x0.size()*c2.size()*a3.size()*a1.size(),
                 1.0, i0data_sorted, x0.size()*c2.size()*a3.size()*a1.size(), i1data_sorted, x0.size()*c2.size()*a3.size()*a1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }

  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}


void Task246::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index a1 = b(4);

  // tensor label: I276
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c2, a3, a1), 0.0);

  for (auto& c4 : *range_[0]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, c4)]);
      sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), c4.size());

      // tensor label: I277
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0, c2, a3, c4, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0, c2, a3, c4, a1)]);
      sort_indices<5,1,0,2,3,4,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size(), c2.size(), a3.size(), c4.size(), a1.size());

      dgemm_("T", "N", 1, ci0.size()*x0.size()*c2.size()*a3.size()*a1.size(), x1.size()*c4.size(),
             1.0, i0data_sorted, x1.size()*c4.size(), i1data_sorted, x1.size()*c4.size(),
             1.0, odata_sorted, 1);
    }
  }

  sort_indices<0,1,2,3,4,1,1,1,1>(odata_sorted, odata, ci0.size(), x0.size(), c2.size(), a3.size(), a1.size());
  out()->put_block(odata, ci0, x0, c2, a3, a1);
}


void Task247::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a3 = b(4);
  const Index c4 = b(5);
  const Index a1 = b(6);

  // tensor label: I277
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a3, c4, a1), 0.0);

  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c4, a1);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c4, a1)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c4.size(), a1.size());

  // tensor label: I278
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
  sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

  dgemm_("T", "N", c2.size()*a3.size()*c4.size()*a1.size(), ci0.size()*x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a3.size()*c4.size()*a1.size());

  sort_indices<4,5,6,0,1,2,3,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), c4.size(), a1.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, c2, a3, c4, a1);
}


void Task248::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);

  // tensor label: I278
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma74
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}


void Task249::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a3 = b(4);
  const Index c4 = b(5);
  const Index a1 = b(6);

  // tensor label: I277
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a3, c4, a1), 0.0);

  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c4, a3);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c4, a3)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c4.size(), a3.size());

  // tensor label: I282
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
  sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

  dgemm_("T", "N", c2.size()*a1.size()*c4.size()*a3.size(), ci0.size()*x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a1.size()*c4.size()*a3.size());

  sort_indices<4,5,6,0,3,2,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), c4.size(), a3.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, c2, a3, c4, a1);
}


void Task250::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);

  // tensor label: I282
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma74
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}


void Task251::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index a1 = b(4);

  // tensor label: I276
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c2, a3, a1), 0.0);

  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a1.size());

    // tensor label: I285
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size());

    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), ci0.size()*x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }

  sort_indices<3,4,1,0,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, c2, a3, a1);
}


void Task252::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);

  // tensor label: I285
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0);
  {
    // tensor label: Gamma78
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0);
    sort_indices<0,1,2,1,1,-1,4>(i0data, odata, ci0.size(), x3.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x0);
}


void Task253::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index a1 = b(4);

  // tensor label: I276
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c2, a3, a1), 0.0);

  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a3.size());

    // tensor label: I288
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size());

    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), ci0.size()*x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }

  sort_indices<3,4,1,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, c2, a3, a1);
}


void Task254::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);

  // tensor label: I288
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0);
  {
    // tensor label: Gamma78
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0);
    sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x3.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x0);
}


void Task255::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index a1 = b(4);

  // tensor label: I276
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c2, a3, a1), 0.0);

  for (auto& c4 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, c4)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), c4.size());

    // tensor label: I291
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, a3, c4, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, a3, c4, a1)]);
    sort_indices<3,0,1,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), a3.size(), c4.size(), a1.size());

    dgemm_("T", "N", c2.size(), ci0.size()*x0.size()*a3.size()*a1.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c2.size());
  }

  sort_indices<1,2,0,3,4,1,1,1,1>(odata_sorted, odata, c2.size(), ci0.size(), x0.size(), a3.size(), a1.size());
  out()->put_block(odata, ci0, x0, c2, a3, a1);
}


void Task256::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a3 = b(2);
  const Index c4 = b(3);
  const Index a1 = b(4);

  // tensor label: I291
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a3, c4, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c4, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c4, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c4.size(), a1.size());

    // tensor label: I292
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

    dgemm_("T", "N", a3.size()*c4.size()*a1.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a3.size()*c4.size()*a1.size());
  }

  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c4.size(), a1.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, a3, c4, a1);
}


void Task257::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);

  // tensor label: I292
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma74
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,1,4>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}


void Task258::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a3 = b(2);
  const Index c4 = b(3);
  const Index a1 = b(4);

  // tensor label: I291
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a3, c4, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c4, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c4.size(), a3.size());

    // tensor label: I296
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

    dgemm_("T", "N", a1.size()*c4.size()*a3.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a1.size()*c4.size()*a3.size());
  }

  sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c4.size(), a3.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, a3, c4, a1);
}


void Task259::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);

  // tensor label: I296
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma74
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}


void Task260::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index a1 = b(4);

  // tensor label: I276
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c2, a3, a1), 0.0);

  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a1.size());

    // tensor label: I299
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, a4, c2, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, a4, c2, a3)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), a4.size(), c2.size(), a3.size());

    dgemm_("T", "N", a1.size(), ci0.size()*x0.size()*c2.size()*a3.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, a1.size());
  }

  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x0.size(), c2.size(), a3.size());
  out()->put_block(odata, ci0, x0, c2, a3, a1);
}


void Task261::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a4 = b(2);
  const Index c2 = b(3);
  const Index a3 = b(4);

  // tensor label: I299
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a4, c2, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a4, c2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a4, c2, a3), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a3.size());

    // tensor label: I300
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

    dgemm_("T", "N", a4.size()*c2.size()*a3.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a4.size()*c2.size()*a3.size());
  }

  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c2.size(), a3.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, a4, c2, a3);
}


void Task262::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);

  // tensor label: I300
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma74
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}


void Task263::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a4 = b(2);
  const Index c2 = b(3);
  const Index a3 = b(4);

  // tensor label: I299
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a4, c2, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a4, c2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a4, c2, a3), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a4.size());

    // tensor label: I304
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

    dgemm_("T", "N", a3.size()*c2.size()*a4.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a4.size());
  }

  sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a4.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, a4, c2, a3);
}


void Task264::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);

  // tensor label: I304
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma74
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,4>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}


void Task265::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index a1 = b(4);

  // tensor label: I276
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c2, a3, a1), 0.0);

  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a3)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a3.size());

    // tensor label: I307
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, a4, c2, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, a4, c2, a1)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), a4.size(), c2.size(), a1.size());

    dgemm_("T", "N", a3.size(), ci0.size()*x0.size()*c2.size()*a1.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, a3.size());
  }

  sort_indices<1,2,3,0,4,1,1,1,1>(odata_sorted, odata, a3.size(), ci0.size(), x0.size(), c2.size(), a1.size());
  out()->put_block(odata, ci0, x0, c2, a3, a1);
}


void Task266::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a4 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);

  // tensor label: I307
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a4, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a4, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a4, c2, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a1.size());

    // tensor label: I308
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

    dgemm_("T", "N", a4.size()*c2.size()*a1.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a4.size()*c2.size()*a1.size());
  }

  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c2.size(), a1.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, a4, c2, a1);
}


void Task267::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);

  // tensor label: I308
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma74
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,4>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}


void Task268::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a4 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);

  // tensor label: I307
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a4, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a4, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a4, c2, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a4.size());

    // tensor label: I312
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

    dgemm_("T", "N", a1.size()*c2.size()*a4.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a4.size());
  }

  sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a4.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, a4, c2, a1);
}


void Task269::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);

  // tensor label: I312
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma74
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}


void Task270::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index a1 = b(4);

  // tensor label: I276
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c2, a3, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), x1.size());

    // tensor label: I315
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1, a1, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1, a1, a3)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size(), a1.size(), a3.size());

    dgemm_("T", "N", c2.size(), ci0.size()*x0.size()*a1.size()*a3.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c2.size());
  }

  sort_indices<1,2,0,4,3,1,1,1,1>(odata_sorted, odata, c2.size(), ci0.size(), x0.size(), a1.size(), a3.size());
  out()->put_block(odata, ci0, x0, c2, a3, a1);
}


void Task271::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  const Index a3 = b(4);

  // tensor label: I315
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, a3), 0.0);

  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a3)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a3.size());

      // tensor label: I316
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


void Task272::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);

  // tensor label: I316
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma86
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}


void Task273::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index a1 = b(4);

  // tensor label: I276
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c2, a3, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a1.size());

    // tensor label: I397
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }

  sort_indices<3,4,1,0,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, c2, a3, a1);
}


void Task274::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);

  // tensor label: I397
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma74
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    dscal_(ci0.size()*x1.size()*x0.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,1,1,1,4>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}


void Task275::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index a1 = b(4);

  // tensor label: I276
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c2, a3, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a3.size());

    // tensor label: I400
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }

  sort_indices<3,4,1,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, c2, a3, a1);
}


void Task276::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);

  // tensor label: I400
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma74
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    dscal_(ci0.size()*x1.size()*x0.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,1,1,-1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}


void Task277::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index a1 = b(4);

  // tensor label: I276
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c2, a3, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a1.size());

    // tensor label: I415
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }

  sort_indices<3,4,1,0,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, c2, a3, a1);
}


void Task278::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);

  // tensor label: I415
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma74
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}


void Task279::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index a1 = b(4);

  // tensor label: I276
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c2, a3, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a3.size());

    // tensor label: I418
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }

  sort_indices<3,4,1,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, c2, a3, a1);
}


void Task280::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);

  // tensor label: I418
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma74
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,2,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}


void Task281::Task_local::compute() {
  const Index ci0 = b(0);

  // tensor label: I261
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

          // tensor label: I318
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


void Task282::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  const Index a2 = b(4);

  // tensor label: I318
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, a2), 0.0);

  for (auto& c3 : *range_[0]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, c3)]);
      sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x2.size(), c3.size());

      // tensor label: I319
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x2, x1, a1, c3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x2, x1, a1, c3, a2)]);
      sort_indices<5,2,0,1,3,4,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x2.size(), x1.size(), a1.size(), c3.size(), a2.size());

      dgemm_("T", "N", 1, ci0.size()*x0.size()*x1.size()*a1.size()*a2.size(), x2.size()*c3.size(),
             1.0, i0data_sorted, x2.size()*c3.size(), i1data_sorted, x2.size()*c3.size(),
             1.0, odata_sorted, 1);
    }
  }

  sort_indices<0,1,2,3,4,1,1,1,1>(odata_sorted, odata, ci0.size(), x0.size(), x1.size(), a1.size(), a2.size());
  out()->put_block(odata, ci0, x0, x1, a1, a2);
}


void Task283::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a1 = b(4);
  const Index c3 = b(5);
  const Index a2 = b(6);

  // tensor label: I319
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x2, x1, a1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x2, x1, a1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x2, x1, a1, c3, a2), 0.0);

  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c3, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c3.size(), a2.size());

    // tensor label: I320
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
    sort_indices<1,0,2,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());

    dgemm_("T", "N", a1.size()*c3.size()*a2.size(), ci0.size()*x0.size()*x2.size()*x1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a1.size()*c3.size()*a2.size());
  }

  sort_indices<3,4,5,6,0,1,2,1,1,1,1>(odata_sorted, odata, a1.size(), c3.size(), a2.size(), ci0.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x0, x2, x1, a1, c3, a2);
}


void Task284::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);

  // tensor label: I320
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma86
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}


void Task285::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  const Index a2 = b(4);

  // tensor label: I318
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, a2), 0.0);

  for (auto& x4 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a2)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a2.size());

      // tensor label: I323
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


void Task286::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x1 = b(4);

  // tensor label: I323
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x4, x1);
  {
    // tensor label: Gamma88
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x1);
    sort_indices<0,1,2,3,4,1,1,1,2>(i0data, odata, ci0.size(), x5.size(), x0.size(), x4.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x0, x4, x1);
}


void Task287::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  const Index a2 = b(4);

  // tensor label: I318
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, a2), 0.0);

  for (auto& a3 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), a2.size());

    // tensor label: I326
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


void Task288::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  const Index a3 = b(4);

  // tensor label: I326
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, a3), 0.0);

  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a3)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a3.size());

      // tensor label: I327
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


void Task289::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);

  // tensor label: I327
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma86
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}


void Task290::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  const Index a2 = b(4);

  // tensor label: I318
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, a2), 0.0);

  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

      // tensor label: I403
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


void Task291::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);

  // tensor label: I403
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma86
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    dscal_(ci0.size()*x3.size()*x0.size()*x2.size()*x1.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}


void Task292::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  const Index a2 = b(4);

  // tensor label: I318
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, a2), 0.0);

  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

      // tensor label: I421
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


void Task293::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);

  // tensor label: I421
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma86
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}


void Task294::Task_local::compute() {
  const Index ci0 = b(0);

  // tensor label: I261
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);

  for (auto& x1 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& c2 : *range_[0]) {
        for (auto& a3 : *range_[2]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a3)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a3.size());

          // tensor label: I343
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, c2, a3, a1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, c2, a3, a1)]);
          sort_indices<1,4,2,3,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), c2.size(), a3.size(), a1.size());

          dgemm_("T", "N", 1, ci0.size(), x1.size()*c2.size()*a3.size()*a1.size(),
                 1.0, i0data_sorted, x1.size()*c2.size()*a3.size()*a1.size(), i1data_sorted, x1.size()*c2.size()*a3.size()*a1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }

  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}


void Task295::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index a1 = b(4);

  // tensor label: I343
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c2, a3, a1), 0.0);

  for (auto& c4 : *range_[0]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, c4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, c4)]);
      sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), c4.size());

      // tensor label: I344
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0, c2, a3, c4, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0, c2, a3, c4, a1)]);
      sort_indices<5,2,0,1,3,4,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size(), c2.size(), a3.size(), c4.size(), a1.size());

      dgemm_("T", "N", 1, ci0.size()*x1.size()*c2.size()*a3.size()*a1.size(), x0.size()*c4.size(),
             1.0, i0data_sorted, x0.size()*c4.size(), i1data_sorted, x0.size()*c4.size(),
             1.0, odata_sorted, 1);
    }
  }

  sort_indices<0,1,2,3,4,1,1,1,1>(odata_sorted, odata, ci0.size(), x1.size(), c2.size(), a3.size(), a1.size());
  out()->put_block(odata, ci0, x1, c2, a3, a1);
}


void Task296::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a3 = b(4);
  const Index c4 = b(5);
  const Index a1 = b(6);

  // tensor label: I344
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a3, c4, a1), 0.0);

  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c4, a1);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c4, a1)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c4.size(), a1.size());

  // tensor label: I345
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
  sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

  dgemm_("T", "N", c2.size()*a3.size()*c4.size()*a1.size(), ci0.size()*x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a3.size()*c4.size()*a1.size());

  sort_indices<4,5,6,0,1,2,3,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), c4.size(), a1.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, c2, a3, c4, a1);
}


void Task297::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);

  // tensor label: I345
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma74
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}


void Task298::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a3 = b(4);
  const Index c4 = b(5);
  const Index a1 = b(6);

  // tensor label: I344
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a3, c4, a1), 0.0);

  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c4, a3);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c4, a3)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c4.size(), a3.size());

  // tensor label: I349
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
  sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

  dgemm_("T", "N", c2.size()*a1.size()*c4.size()*a3.size(), ci0.size()*x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a1.size()*c4.size()*a3.size());

  sort_indices<4,5,6,0,3,2,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), c4.size(), a3.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, c2, a3, c4, a1);
}


void Task299::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);

  // tensor label: I349
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma74
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}


void Task300::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index a1 = b(4);

  // tensor label: I343
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c2, a3, a1), 0.0);

  for (auto& c4 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, c4)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), c4.size());

    // tensor label: I358
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, a3, c4, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, a3, c4, a1)]);
    sort_indices<3,0,1,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), a3.size(), c4.size(), a1.size());

    dgemm_("T", "N", c2.size(), ci0.size()*x1.size()*a3.size()*a1.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c2.size());
  }

  sort_indices<1,2,0,3,4,1,1,1,1>(odata_sorted, odata, c2.size(), ci0.size(), x1.size(), a3.size(), a1.size());
  out()->put_block(odata, ci0, x1, c2, a3, a1);
}


void Task301::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index c4 = b(3);
  const Index a1 = b(4);

  // tensor label: I358
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a3, c4, a1), 0.0);

  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c4, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a3, c4, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), c4.size(), a1.size());

    // tensor label: I359
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

    dgemm_("T", "N", a3.size()*c4.size()*a1.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a3.size()*c4.size()*a1.size());
  }

  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c4.size(), a1.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, a3, c4, a1);
}


void Task302::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);

  // tensor label: I359
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma74
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,1,4>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}


void Task303::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index c4 = b(3);
  const Index a1 = b(4);

  // tensor label: I358
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a3, c4, a1), 0.0);

  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c4, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c4.size(), a3.size());

    // tensor label: I363
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

    dgemm_("T", "N", a1.size()*c4.size()*a3.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a1.size()*c4.size()*a3.size());
  }

  sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c4.size(), a3.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, a3, c4, a1);
}


void Task304::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);

  // tensor label: I363
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma74
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}


void Task305::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index a1 = b(4);

  // tensor label: I343
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c2, a3, a1), 0.0);

  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a1.size());

    // tensor label: I366
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, a4, c2, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, a4, c2, a3)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), a4.size(), c2.size(), a3.size());

    dgemm_("T", "N", a1.size(), ci0.size()*x1.size()*c2.size()*a3.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, a1.size());
  }

  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x1.size(), c2.size(), a3.size());
  out()->put_block(odata, ci0, x1, c2, a3, a1);
}


void Task306::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index a4 = b(2);
  const Index c2 = b(3);
  const Index a3 = b(4);

  // tensor label: I366
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a4, c2, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a4, c2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a4, c2, a3), 0.0);

  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a4, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a4, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a4.size(), c2.size(), a3.size());

    // tensor label: I367
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

    dgemm_("T", "N", a4.size()*c2.size()*a3.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a4.size()*c2.size()*a3.size());
  }

  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c2.size(), a3.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, a4, c2, a3);
}


void Task307::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);

  // tensor label: I367
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma74
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}


void Task308::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index a4 = b(2);
  const Index c2 = b(3);
  const Index a3 = b(4);

  // tensor label: I366
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a4, c2, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a4, c2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a4, c2, a3), 0.0);

  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c2, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a3, c2, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), c2.size(), a4.size());

    // tensor label: I371
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

    dgemm_("T", "N", a3.size()*c2.size()*a4.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a4.size());
  }

  sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a4.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, a4, c2, a3);
}


void Task309::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);

  // tensor label: I371
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma74
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,4>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}


void Task310::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index a1 = b(4);

  // tensor label: I343
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c2, a3, a1), 0.0);

  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a3)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a3.size());

    // tensor label: I374
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, a4, c2, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, a4, c2, a1)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), a4.size(), c2.size(), a1.size());

    dgemm_("T", "N", a3.size(), ci0.size()*x1.size()*c2.size()*a1.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, a3.size());
  }

  sort_indices<1,2,3,0,4,1,1,1,1>(odata_sorted, odata, a3.size(), ci0.size(), x1.size(), c2.size(), a1.size());
  out()->put_block(odata, ci0, x1, c2, a3, a1);
}


void Task311::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index a4 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);

  // tensor label: I374
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a4, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a4, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a4, c2, a1), 0.0);

  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a4, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a4, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a4.size(), c2.size(), a1.size());

    // tensor label: I375
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

    dgemm_("T", "N", a4.size()*c2.size()*a1.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a4.size()*c2.size()*a1.size());
  }

  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c2.size(), a1.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, a4, c2, a1);
}


void Task312::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);

  // tensor label: I375
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma74
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,4>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}


void Task313::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index a4 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);

  // tensor label: I374
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a4, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a4, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a4, c2, a1), 0.0);

  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a4.size());

    // tensor label: I379
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

    dgemm_("T", "N", a1.size()*c2.size()*a4.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a4.size());
  }

  sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a4.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, a4, c2, a1);
}


void Task314::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);

  // tensor label: I379
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma74
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}


void Task315::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index a1 = b(4);

  // tensor label: I343
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c2, a3, a1), 0.0);

  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), c2.size(), a1.size());

    // tensor label: I406
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }

  sort_indices<3,4,1,0,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, c2, a3, a1);
}


void Task316::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);

  // tensor label: I406
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma74
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    dscal_(ci0.size()*x1.size()*x0.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,1,1,1,4>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}


void Task317::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index a1 = b(4);

  // tensor label: I343
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c2, a3, a1), 0.0);

  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());

    // tensor label: I409
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }

  sort_indices<3,4,1,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, c2, a3, a1);
}


void Task318::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);

  // tensor label: I409
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma74
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    dscal_(ci0.size()*x1.size()*x0.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,1,1,-1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}


void Task319::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index a1 = b(4);

  // tensor label: I343
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c2, a3, a1), 0.0);

  for (auto& x0 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), c2.size(), a1.size());

    // tensor label: I424
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }

  sort_indices<3,4,1,0,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, c2, a3, a1);
}


void Task320::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);

  // tensor label: I424
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma74
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}


void Task321::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index a1 = b(4);

  // tensor label: I343
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c2, a3, a1), 0.0);

  for (auto& x0 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());

    // tensor label: I427
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }

  sort_indices<3,4,1,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, c2, a3, a1);
}


void Task322::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);

  // tensor label: I427
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma74
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,2,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}


void Task323::Task_local::compute() {
  const Index ci0 = b(0);

  // tensor label: I261
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);

  for (auto& x3 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& c2 : *range_[0]) {
        for (auto& a3 : *range_[2]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a3)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a3.size());

          // tensor label: I351
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, a3, c2, a1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, a3, c2, a1)]);
          sort_indices<1,4,3,2,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), a3.size(), c2.size(), a1.size());

          dgemm_("T", "N", 1, ci0.size(), x3.size()*a3.size()*c2.size()*a1.size(),
                 1.0, i0data_sorted, x3.size()*a3.size()*c2.size()*a1.size(), i1data_sorted, x3.size()*a3.size()*c2.size()*a1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }

  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}


void Task324::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);

  // tensor label: I351
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, a3, c2, a1), 0.0);

  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), c2.size(), a1.size());

    // tensor label: I352
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size());

    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), ci0.size()*x3.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }

  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), ci0.size(), x3.size());
  out()->put_block(odata, ci0, x3, a3, c2, a1);
}


void Task325::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);

  // tensor label: I352
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0);
  {
    // tensor label: Gamma78
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0);
    sort_indices<0,1,2,1,1,-1,4>(i0data, odata, ci0.size(), x3.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x0);
}


void Task326::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);

  // tensor label: I351
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, a3, c2, a1), 0.0);

  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());

    // tensor label: I355
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size());

    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), ci0.size()*x3.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }

  sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), ci0.size(), x3.size());
  out()->put_block(odata, ci0, x3, a3, c2, a1);
}


void Task327::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);

  // tensor label: I355
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0);
  {
    // tensor label: Gamma78
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0);
    sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x3.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x0);
}


void Task328::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);

  // tensor label: I351
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, a3, c2, a1), 0.0);

  for (auto& x2 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), x2.size());

    // tensor label: I382
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, a1, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, a1, a3)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), a1.size(), a3.size());

    dgemm_("T", "N", c2.size(), ci0.size()*x3.size()*a1.size()*a3.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, c2.size());
  }

  sort_indices<1,2,4,0,3,1,1,1,1>(odata_sorted, odata, c2.size(), ci0.size(), x3.size(), a1.size(), a3.size());
  out()->put_block(odata, ci0, x3, a3, c2, a1);
}


void Task329::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a1 = b(3);
  const Index a3 = b(4);

  // tensor label: I382
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, a3), 0.0);

  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a3)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a3.size());

      // tensor label: I383
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


void Task330::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);

  // tensor label: I383
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma86
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}


void Task331::Task_local::compute() {
  const Index ci0 = b(0);

  // tensor label: I261
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

          // tensor label: I385
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


void Task332::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a1 = b(3);
  const Index a2 = b(4);

  // tensor label: I385
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, a2), 0.0);

  for (auto& c3 : *range_[0]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, c3)]);
      sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), c3.size());

      // tensor label: I386
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, a1, c3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, a1, c3, a2)]);
      sort_indices<5,3,0,1,2,4,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), a1.size(), c3.size(), a2.size());

      dgemm_("T", "N", 1, ci0.size()*x3.size()*x2.size()*a1.size()*a2.size(), x1.size()*c3.size(),
             1.0, i0data_sorted, x1.size()*c3.size(), i1data_sorted, x1.size()*c3.size(),
             1.0, odata_sorted, 1);
    }
  }

  sort_indices<0,1,2,3,4,1,1,1,1>(odata_sorted, odata, ci0.size(), x3.size(), x2.size(), a1.size(), a2.size());
  out()->put_block(odata, ci0, x3, x2, a1, a2);
}


void Task333::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a1 = b(4);
  const Index c3 = b(5);
  const Index a2 = b(6);

  // tensor label: I386
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, a1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, x1, a1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, x1, a1, c3, a2), 0.0);

  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c3, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c3.size(), a2.size());

    // tensor label: I387
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());

    dgemm_("T", "N", a1.size()*c3.size()*a2.size(), ci0.size()*x3.size()*x2.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a1.size()*c3.size()*a2.size());
  }

  sort_indices<3,4,5,6,0,1,2,1,1,1,1>(odata_sorted, odata, a1.size(), c3.size(), a2.size(), ci0.size(), x3.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x3, x2, x1, a1, c3, a2);
}


void Task334::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);

  // tensor label: I387
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma86
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}


void Task335::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a1 = b(3);
  const Index a2 = b(4);

  // tensor label: I385
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, a2), 0.0);

  for (auto& a3 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), a2.size());

    // tensor label: I393
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


void Task336::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a1 = b(3);
  const Index a3 = b(4);

  // tensor label: I393
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, a3), 0.0);

  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a3)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a3.size());

      // tensor label: I394
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


void Task337::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);

  // tensor label: I394
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma86
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}


void Task338::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a1 = b(3);
  const Index a2 = b(4);

  // tensor label: I385
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, a2), 0.0);

  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());

      // tensor label: I412
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


void Task339::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);

  // tensor label: I412
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma86
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    dscal_(ci0.size()*x3.size()*x0.size()*x2.size()*x1.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}


void Task340::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a1 = b(3);
  const Index a2 = b(4);

  // tensor label: I385
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, a2), 0.0);

  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());

      // tensor label: I430
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


void Task341::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);

  // tensor label: I430
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma86
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}


void Task342::Task_local::compute() {
  const Index ci0 = b(0);

  // tensor label: I261
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

          // tensor label: I389
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


void Task343::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index a1 = b(3);
  const Index a2 = b(4);

  // tensor label: I389
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x5, x4, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x5, x4, a1, a2), 0.0);

  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());

      // tensor label: I390
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


void Task344::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x1 = b(4);

  // tensor label: I390
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x4, x1);
  {
    // tensor label: Gamma88
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x1);
    sort_indices<0,1,2,3,4,1,1,1,2>(i0data, odata, ci0.size(), x5.size(), x0.size(), x4.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x0, x4, x1);
}


