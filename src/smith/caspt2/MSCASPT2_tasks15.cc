//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_tasks15.cc
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

#include <src/smith/caspt2/MSCASPT2_tasks15.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MSCASPT2;

void Task700::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I1166
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a1 : *range_[2]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c2, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a3, c2, a1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), c2.size(), a1.size());
        // tensor label: l2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a1, c2, a3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a1, c2, a3)]);
        sort_indices<3,2,1,0,0,1,-2,1>(i1data, i1data_sorted, x1.size(), a1.size(), c2.size(), a3.size());
        dgemm_("T", "N", x0.size(), x1.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->add_block(odata, x0, x1);
}

void Task701::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I1166
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
        // tensor label: l2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a1, c2, a3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a1, c2, a3)]);
        sort_indices<1,2,3,0,0,1,4,1>(i1data, i1data_sorted, x1.size(), a1.size(), c2.size(), a3.size());
        dgemm_("T", "N", x0.size(), x1.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->add_block(odata, x0, x1);
}

void Task702::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I1166
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: h1
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size());
      // tensor label: I1233
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c2, a1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c2, a1, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), c2.size(), a1.size(), x0.size());
      dgemm_("T", "N", 1, x1.size()*x0.size(), c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->add_block(odata, x0, x1);
}

void Task703::Task_local::compute() {
  const Index x1 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1233
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, c2, a1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, c2, a1, x0), 0.0);
  {
    // tensor label: l2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
    sort_indices<3,2,1,0,1,1,-2,1>(i0data, odata, x0.size(), a1.size(), c2.size(), x1.size());
  }
  {
    // tensor label: l2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x1, a1, c2, x0);
    sort_indices<0,2,1,3,1,1,-2,1>(i1data, odata, x1.size(), a1.size(), c2.size(), x0.size());
  }
  out()->add_block(odata, x1, c2, a1, x0);
}

void Task704::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I1166
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: h1
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size());
      // tensor label: I1236
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, a2, c1)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), a2.size(), c1.size());
      dgemm_("T", "N", 1, x1.size()*x0.size(), a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->add_block(odata, x0, x1);
}

void Task705::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1236
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0, a2, c1)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, a2, c1), 0.0);
  {
    // tensor label: l2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
    sort_indices<3,2,1,0,1,1,4,1>(i0data, odata, c1.size(), a2.size(), x0.size(), x1.size());
  }
  {
    // tensor label: l2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, x1, x0);
    sort_indices<2,3,1,0,1,1,4,1>(i1data, odata, c1.size(), a2.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, x1, x0, a2, c1);
}

void Task706::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1120
  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        for (auto& x1 : *range_[1]) {
          // tensor label: Gamma170
          std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x3, x0, x2, x1)]);
          sort_indices<1,2,3,4,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
          // tensor label: I1172
          std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, x3, x2);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, x3, x2)]);
          sort_indices<2,1,3,0,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), x3.size(), x2.size());
          dgemm_("T", "N", ci0.size(), 1, x1.size()*x0.size()*x3.size()*x2.size(),
                 1.0, i0data_sorted, x1.size()*x0.size()*x3.size()*x2.size(), i1data_sorted, x1.size()*x0.size()*x3.size()*x2.size(),
                 1.0, odata_sorted, ci0.size());
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->add_block(odata, ci0);
}

void Task707::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1172
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());
      // tensor label: l2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a1, x1, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a1, x1, a2)]);
      sort_indices<1,3,0,2,0,1,2,1>(i1data, i1data_sorted, x0.size(), a1.size(), x1.size(), a2.size());
      dgemm_("T", "N", x3.size()*x2.size(), x1.size()*x0.size(), a2.size()*a1.size(),
             1.0, i0data_sorted, a2.size()*a1.size(), i1data_sorted, a2.size()*a1.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<3,2,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x0.size(), x1.size());
  out()->add_block(odata, x1, x0, x3, x2);
}

void Task708::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1172
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());
      // tensor label: l2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, x2, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, x2, a2)]);
      sort_indices<1,3,0,2,0,1,2,1>(i1data, i1data_sorted, x3.size(), a1.size(), x2.size(), a2.size());
      dgemm_("T", "N", x0.size()*x1.size(), x2.size()*x3.size(), a2.size()*a1.size(),
             1.0, i0data_sorted, a2.size()*a1.size(), i1data_sorted, a2.size()*a1.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x3.size(), x2.size());
  out()->add_block(odata, x1, x0, x3, x2);
}

void Task709::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1172
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size());
    // tensor label: l2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a1, x1, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a1, x1, x2)]);
    sort_indices<1,0,2,3,0,1,2,1>(i1data, i1data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
    dgemm_("T", "N", x3.size(), x2.size()*x1.size()*x0.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x3.size());
  }
  sort_indices<2,1,0,3,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), x1.size(), x2.size());
  out()->add_block(odata, x1, x0, x3, x2);
}

void Task710::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1172
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size());
    // tensor label: l2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, x2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, x2, x1)]);
    sort_indices<1,0,2,3,0,1,2,1>(i1data, i1data_sorted, x3.size(), a1.size(), x2.size(), x1.size());
    dgemm_("T", "N", x0.size(), x1.size()*x2.size()*x3.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x3.size(), x2.size(), x1.size());
  out()->add_block(odata, x1, x0, x3, x2);
}

void Task711::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1120
  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        for (auto& x3 : *range_[1]) {
          for (auto& x2 : *range_[1]) {
            for (auto& x1 : *range_[1]) {
              // tensor label: Gamma341
              std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x4, x0, x3, x2, x1);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x5, x4, x0, x3, x2, x1)]);
              sort_indices<1,2,3,4,5,6,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x5.size(), x4.size(), x0.size(), x3.size(), x2.size(), x1.size());
              // tensor label: I1178
              std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x4, x5, x0, x1, x2);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x4, x5, x0, x1, x2)]);
              sort_indices<2,1,3,0,5,4,0,1,1,1>(i1data, i1data_sorted, x3.size(), x4.size(), x5.size(), x0.size(), x1.size(), x2.size());
              dgemm_("T", "N", ci0.size(), 1, x3.size()*x4.size()*x5.size()*x0.size()*x1.size()*x2.size(),
                     1.0, i0data_sorted, x3.size()*x4.size()*x5.size()*x0.size()*x1.size()*x2.size(), i1data_sorted, x3.size()*x4.size()*x5.size()*x0.size()*x1.size()*x2.size(),
                     1.0, odata_sorted, ci0.size());
            }
          }
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->add_block(odata, ci0);
}

void Task712::Task_local::compute() {
  const Index x3 = b(0);
  const Index x4 = b(1);
  const Index x5 = b(2);
  const Index x0 = b(3);
  const Index x1 = b(4);
  const Index x2 = b(5);
  // tensor label: I1178
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x4, x5, x0, x1, x2)]);
  std::fill_n(odata.get(), out()->get_size(x3, x4, x5, x0, x1, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x4, x5, x0, x1, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x4, x5, x0, x1, x2), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x0, x1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x0, x1, x2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x0.size(), x1.size(), x2.size());
    // tensor label: l2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, c1, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, c1, x3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), c1.size(), x3.size());
    dgemm_("T", "N", x0.size()*x1.size()*x2.size(), x3.size()*x4.size()*x5.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, x0.size()*x1.size()*x2.size());
  }
  sort_indices<5,4,3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, x3, x4, x5, x0, x1, x2);
}

void Task713::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1120
  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: Gamma132
          std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x1, x0, x2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x3, x1, x0, x2)]);
          sort_indices<1,2,3,4,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x3.size(), x1.size(), x0.size(), x2.size());
          // tensor label: I1193
          std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x3, x0, x1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x3, x0, x1)]);
          sort_indices<1,3,2,0,0,1,1,1>(i1data, i1data_sorted, x2.size(), x3.size(), x0.size(), x1.size());
          dgemm_("T", "N", ci0.size(), 1, x2.size()*x3.size()*x0.size()*x1.size(),
                 1.0, i0data_sorted, x2.size()*x3.size()*x0.size()*x1.size(), i1data_sorted, x2.size()*x3.size()*x0.size()*x1.size(),
                 1.0, odata_sorted, ci0.size());
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->add_block(odata, ci0);
}

void Task714::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1193
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, x3, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x2, x3, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x3, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, x0, x1), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x0, x1, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x0, x1, a1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, c2.size(), x0.size(), x1.size(), a1.size());
      // tensor label: l2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, c2, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, c2, x2)]);
      sort_indices<2,1,0,3,0,1,-1,1>(i1data, i1data_sorted, x3.size(), a1.size(), c2.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), x2.size()*x3.size(), c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<3,2,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x3.size(), x2.size());
  out()->add_block(odata, x2, x3, x0, x1);
}

void Task715::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1120
  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        for (auto& x1 : *range_[1]) {
          // tensor label: Gamma122
          std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x0, x1);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x3, x2, x0, x1)]);
          sort_indices<1,2,3,4,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x3.size(), x2.size(), x0.size(), x1.size());
          // tensor label: I1205
          std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x3, x0, x1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x3, x0, x1)]);
          sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), x3.size(), x0.size(), x1.size());
          dgemm_("T", "N", ci0.size(), 1, x2.size()*x3.size()*x0.size()*x1.size(),
                 1.0, i0data_sorted, x2.size()*x3.size()*x0.size()*x1.size(), i1data_sorted, x2.size()*x3.size()*x0.size()*x1.size(),
                 1.0, odata_sorted, ci0.size());
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->add_block(odata, ci0);
}

void Task716::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1205
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, x3, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x2, x3, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x3, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, x0, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x0, x1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x0, x1, a2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), x0.size(), x1.size(), a2.size());
      // tensor label: l2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, x3, x2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size(), x3.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), x2.size()*x3.size(), a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<3,2,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x3.size(), x2.size());
  out()->add_block(odata, x2, x3, x0, x1);
}

void Task717::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1205
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, x3, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x2, x3, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x3, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, x0, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), x0.size());
    // tensor label: l2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, c1, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, c1, x1)]);
    sort_indices<2,0,1,3,0,1,2,1>(i1data, i1data_sorted, x3.size(), x2.size(), c1.size(), x1.size());
    dgemm_("T", "N", x0.size(), x1.size()*x2.size()*x3.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<2,1,0,3,1,1,1,1>(odata_sorted, odata, x0.size(), x3.size(), x2.size(), x1.size());
  out()->add_block(odata, x2, x3, x0, x1);
}

void Task718::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1120
  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        for (auto& x3 : *range_[1]) {
          for (auto& x1 : *range_[1]) {
            for (auto& x0 : *range_[1]) {
              // tensor label: Gamma161
              std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x2, x4, x3, x1, x0);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x5, x2, x4, x3, x1, x0)]);
              sort_indices<1,2,3,4,5,6,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
              // tensor label: I1217
              std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x4, x5, x0, x1, x2);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x4, x5, x0, x1, x2)]);
              sort_indices<2,5,1,0,4,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x4.size(), x5.size(), x0.size(), x1.size(), x2.size());
              dgemm_("T", "N", ci0.size(), 1, x3.size()*x4.size()*x5.size()*x0.size()*x1.size()*x2.size(),
                     1.0, i0data_sorted, x3.size()*x4.size()*x5.size()*x0.size()*x1.size()*x2.size(), i1data_sorted, x3.size()*x4.size()*x5.size()*x0.size()*x1.size()*x2.size(),
                     1.0, odata_sorted, ci0.size());
            }
          }
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->add_block(odata, ci0);
}

void Task719::Task_local::compute() {
  const Index x3 = b(0);
  const Index x4 = b(1);
  const Index x5 = b(2);
  const Index x0 = b(3);
  const Index x1 = b(4);
  const Index x2 = b(5);
  // tensor label: I1217
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x4, x5, x0, x1, x2)]);
  std::fill_n(odata.get(), out()->get_size(x3, x4, x5, x0, x1, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x4, x5, x0, x1, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x4, x5, x0, x1, x2), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, x2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1, x2, a1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size(), x2.size(), a1.size());
    // tensor label: l2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a1, x4, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a1, x4, x3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), a1.size(), x4.size(), x3.size());
    dgemm_("T", "N", x0.size()*x1.size()*x2.size(), x3.size()*x4.size()*x5.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x0.size()*x1.size()*x2.size());
  }
  sort_indices<5,4,3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, x3, x4, x5, x0, x1, x2);
}

#endif
