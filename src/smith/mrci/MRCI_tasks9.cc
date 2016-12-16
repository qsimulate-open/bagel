//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_tasks9.cc
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

#include <src/smith/mrci/MRCI_tasks9.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task400::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a1, c2, c3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, a1, c2, c3)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x2.size(), a1.size(), c2.size(), c3.size());
      // tensor label: I666
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x1, x2, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x1, x2, x0)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), x1.size(), x2.size(), x0.size());
      dgemm_("T", "N", a1.size()*c2.size(), x1.size()*x0.size(), c3.size()*x2.size(),
             1.0, i0data_sorted, c3.size()*x2.size(), i1data_sorted, c3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*c2.size());
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), x1.size(), x0.size());
  out()->add_block(odata, c2, x1, x0, a1);
}

void Task401::Task_local::compute() {
  const Index c3 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index x0 = b(3);
  // tensor label: I666
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, x1, x2, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, x1, x2, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x1, x2, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x1, x2, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma24
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x1, x3, x2, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x1, x3, x2, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x1.size(), x3.size(), x2.size(), x0.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, c3, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, c3, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), c3.size(), x3.size());
        dgemm_("T", "N", x1.size()*x2.size()*x0.size(), c3.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x1.size()*x2.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x2.size(), x0.size(), c3.size());
  out()->add_block(odata, c3, x1, x2, x0);
}

void Task402::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& x5 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a1, c2, x5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a1, c2, x5)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, c3.size(), a1.size(), c2.size(), x5.size());
      // tensor label: I669
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x1, x5, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x1, x5, x0)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), x1.size(), x5.size(), x0.size());
      dgemm_("T", "N", a1.size()*c2.size(), x1.size()*x0.size(), c3.size()*x5.size(),
             1.0, i0data_sorted, c3.size()*x5.size(), i1data_sorted, c3.size()*x5.size(),
             1.0, odata_sorted, a1.size()*c2.size());
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), x1.size(), x0.size());
  out()->add_block(odata, c2, x1, x0, a1);
}

void Task403::Task_local::compute() {
  const Index c3 = b(0);
  const Index x1 = b(1);
  const Index x5 = b(2);
  const Index x0 = b(3);
  // tensor label: I669
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, x1, x5, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, x1, x5, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x1, x5, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x1, x5, x0), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma220
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x5, x4, x0, x3, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x5, x4, x0, x3, x2)]);
        sort_indices<2,4,5,0,1,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), x5.size(), x4.size(), x0.size(), x3.size(), x2.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, c3, x3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, c3, x3, x2)]);
        sort_indices<0,2,3,1,0,1,1,2>(i1data, i1data_sorted, x4.size(), c3.size(), x3.size(), x2.size());
        dgemm_("T", "N", x1.size()*x5.size()*x0.size(), c3.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x1.size()*x5.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x5.size(), x0.size(), c3.size());
  out()->add_block(odata, c3, x1, x5, x0);
}

void Task404::Task_local::compute() {
  const Index c3 = b(0);
  const Index x1 = b(1);
  const Index x5 = b(2);
  const Index x0 = b(3);
  // tensor label: I669
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, x1, x5, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, x1, x5, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x1, x5, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x1, x5, x0), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma222
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x5, x4, x3, x2, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x5, x4, x3, x2, x0)]);
        sort_indices<2,3,4,0,1,5,0,1,1,1>(i0data, i0data_sorted, x1.size(), x5.size(), x4.size(), x3.size(), x2.size(), x0.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x3, x2, c3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x3, x2, c3)]);
        sort_indices<0,1,2,3,0,1,1,2>(i1data, i1data_sorted, x4.size(), x3.size(), x2.size(), c3.size());
        dgemm_("T", "N", x1.size()*x5.size()*x0.size(), c3.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x1.size()*x5.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x5.size(), x0.size(), c3.size());
  out()->add_block(odata, c3, x1, x5, x0);
}

void Task405::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& x5 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c3, x5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c3, x5)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c3.size(), x5.size());
      // tensor label: I672
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x5, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x5, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), x5.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), x1.size()*x0.size(), c3.size()*x5.size(),
             1.0, i0data_sorted, c3.size()*x5.size(), i1data_sorted, c3.size()*x5.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x1.size(), x0.size());
  out()->add_block(odata, c2, x1, x0, a1);
}

void Task406::Task_local::compute() {
  const Index c3 = b(0);
  const Index x5 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I672
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, x5, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, x5, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x5, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x5, x1, x0), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma221
        std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x5, x3, x2, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x4, x5, x3, x2, x1, x0)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x4.size(), x5.size(), x3.size(), x2.size(), x1.size(), x0.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, c3, x3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, c3, x3, x2)]);
        sort_indices<0,2,3,1,0,1,-1,2>(i1data, i1data_sorted, x4.size(), c3.size(), x3.size(), x2.size());
        dgemm_("T", "N", x5.size()*x1.size()*x0.size(), c3.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x5.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x1.size(), x0.size(), c3.size());
  out()->add_block(odata, c3, x5, x1, x0);
}

void Task407::Task_local::compute() {
  const Index c3 = b(0);
  const Index x5 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I672
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, x5, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, x5, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x5, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x5, x1, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma104
        std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x5, x4, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x5, x4, x3, x1, x0)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x2.size(), x5.size(), x4.size(), x3.size(), x1.size(), x0.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x3, x2, c3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x3, x2, c3)]);
        sort_indices<2,0,1,3,0,1,-1,2>(i1data, i1data_sorted, x4.size(), x3.size(), x2.size(), c3.size());
        dgemm_("T", "N", x5.size()*x1.size()*x0.size(), c3.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x5.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x1.size(), x0.size(), c3.size());
  out()->add_block(odata, c3, x5, x1, x0);
}

void Task408::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& c3 : *range_[0]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c3, x4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, c3, x4)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c3.size(), x4.size());
        // tensor label: I699
        std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c3, x5, x0, x1, x4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c3, x5, x0, x1, x4)]);
        sort_indices<2,1,5,0,3,4,0,1,1,1>(i1data, i1data_sorted, c2.size(), c3.size(), x5.size(), x0.size(), x1.size(), x4.size());
        dgemm_("T", "N", a1.size(), c2.size()*x0.size()*x1.size(), c3.size()*x5.size()*x4.size(),
               1.0, i0data_sorted, c3.size()*x5.size()*x4.size(), i1data_sorted, c3.size()*x5.size()*x4.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,3,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), x0.size(), x1.size());
  out()->add_block(odata, c2, x1, x0, a1);
}

void Task409::Task_local::compute() {
  const Index c2 = b(0);
  const Index c3 = b(1);
  const Index x5 = b(2);
  const Index x0 = b(3);
  const Index x1 = b(4);
  const Index x4 = b(5);
  // tensor label: I699
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c3, x5, x0, x1, x4)]);
  std::fill_n(odata.get(), out()->get_size(c2, c3, x5, x0, x1, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c3, x5, x0, x1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c3, x5, x0, x1, x4), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma230
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x1, x4, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x1, x4, x3, x2)]);
      sort_indices<4,5,0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x1.size(), x4.size(), x3.size(), x2.size());
      // tensor label: I700
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c3, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c3, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), c3.size(), x3.size(), x2.size());
      dgemm_("T", "N", x5.size()*x0.size()*x1.size()*x4.size(), c2.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x5.size()*x0.size()*x1.size()*x4.size());
    }
  }
  sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, x5.size(), x0.size(), x1.size(), x4.size(), c2.size(), c3.size());
  out()->add_block(odata, c2, c3, x5, x0, x1, x4);
}

void Task410::Task_local::compute() {
  const Index c2 = b(0);
  const Index c3 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I700
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c3, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c2, c3, x3, x2), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c3, x3, x2);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, c2.size(), c3.size(), x3.size(), x2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x3, x2, c2, c3);
    sort_indices<2,3,0,1,1,1,-1,2>(i1data, odata, x3.size(), x2.size(), c2.size(), c3.size());
  }
  out()->add_block(odata, c2, c3, x3, x2);
}

void Task411::Task_local::compute() {
  const Index c2 = b(0);
  const Index c3 = b(1);
  const Index x5 = b(2);
  const Index x0 = b(3);
  const Index x1 = b(4);
  const Index x4 = b(5);
  // tensor label: I699
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c3, x5, x0, x1, x4)]);
  std::fill_n(odata.get(), out()->get_size(c2, c3, x5, x0, x1, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c3, x5, x0, x1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c3, x5, x0, x1, x4), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: Gamma232
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x2, x4, x1, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x2, x4, x1, x3)]);
      sort_indices<2,5,0,1,3,4,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x2.size(), x4.size(), x1.size(), x3.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x3, x2, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x3, x2, c3)]);
      sort_indices<2,1,0,3,0,1,-1,2>(i1data, i1data_sorted, c2.size(), x3.size(), x2.size(), c3.size());
      dgemm_("T", "N", x5.size()*x0.size()*x4.size()*x1.size(), c2.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x5.size()*x0.size()*x4.size()*x1.size());
    }
  }
  sort_indices<4,5,0,1,3,2,1,1,1,1>(odata_sorted, odata, x5.size(), x0.size(), x4.size(), x1.size(), c2.size(), c3.size());
  out()->add_block(odata, c2, c3, x5, x0, x1, x4);
}

void Task412::Task_local::compute() {
  const Index c2 = b(0);
  const Index c3 = b(1);
  const Index x5 = b(2);
  const Index x0 = b(3);
  const Index x1 = b(4);
  const Index x4 = b(5);
  // tensor label: I699
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c3, x5, x0, x1, x4)]);
  std::fill_n(odata.get(), out()->get_size(c2, c3, x5, x0, x1, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c3, x5, x0, x1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c3, x5, x0, x1, x4), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma234
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x3, x4, x1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x3, x4, x1, x2)]);
      sort_indices<2,5,0,1,3,4,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x3.size(), x4.size(), x1.size(), x2.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c3, c2, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c3, c2, x2)]);
      sort_indices<0,3,1,2,0,1,1,2>(i1data, i1data_sorted, x3.size(), c3.size(), c2.size(), x2.size());
      dgemm_("T", "N", x5.size()*x0.size()*x4.size()*x1.size(), c3.size()*c2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x5.size()*x0.size()*x4.size()*x1.size());
    }
  }
  sort_indices<5,4,0,1,3,2,1,1,1,1>(odata_sorted, odata, x5.size(), x0.size(), x4.size(), x1.size(), c3.size(), c2.size());
  out()->add_block(odata, c2, c3, x5, x0, x1, x4);
}

void Task413::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: Gamma230
          std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x1, x4, x3, x2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x1, x4, x3, x2)]);
          sort_indices<0,3,4,5,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x1.size(), x4.size(), x3.size(), x2.size());
          // tensor label: I702
          std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x3, x2, x5, c2, x4);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x3, x2, x5, c2, x4)]);
          sort_indices<3,5,1,2,0,4,0,1,1,1>(i1data, i1data_sorted, a1.size(), x3.size(), x2.size(), x5.size(), c2.size(), x4.size());
          dgemm_("T", "N", x0.size()*x1.size(), a1.size()*c2.size(), x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, x0.size()*x1.size());
        }
      }
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), c2.size());
  out()->add_block(odata, c2, x1, x0, a1);
}

void Task414::Task_local::compute() {
  const Index a1 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x5 = b(3);
  const Index c2 = b(4);
  const Index x4 = b(5);
  // tensor label: I702
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x3, x2, x5, c2, x4)]);
  std::fill_n(odata.get(), out()->get_size(a1, x3, x2, x5, c2, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x3, x2, x5, c2, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x3, x2, x5, c2, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a3, c2, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a3, c2, x4)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a3.size(), c2.size(), x4.size());
    // tensor label: I703
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, a1, x3, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, a1, x3, x2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), a1.size(), x3.size(), x2.size());
    dgemm_("T", "N", x5.size()*c2.size()*x4.size(), a1.size()*x3.size()*x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x5.size()*c2.size()*x4.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), c2.size(), x4.size(), a1.size(), x3.size(), x2.size());
  out()->add_block(odata, a1, x3, x2, x5, c2, x4);
}

void Task415::Task_local::compute() {
  const Index a3 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I703
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, a1, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(a3, a1, x3, x2), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a1, x3, x2);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, a3.size(), a1.size(), x3.size(), x2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x3, x2, a3, a1);
    sort_indices<2,3,0,1,1,1,1,2>(i1data, odata, x3.size(), x2.size(), a3.size(), a1.size());
  }
  out()->add_block(odata, a3, a1, x3, x2);
}

void Task416::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: Gamma233
          std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x3, x1, x4, x2, x0);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x3, x1, x4, x2, x0)]);
          sort_indices<0,1,3,4,2,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x3.size(), x1.size(), x4.size(), x2.size(), x0.size());
          // tensor label: I708
          std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, a1, x5, c2, x4);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, a1, x5, c2, x4)]);
          sort_indices<3,0,5,1,2,4,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), a1.size(), x5.size(), c2.size(), x4.size());
          dgemm_("T", "N", x1.size()*x0.size(), a1.size()*c2.size(), x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, x1.size()*x0.size());
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a1.size(), c2.size());
  out()->add_block(odata, c2, x1, x0, a1);
}

void Task417::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index a1 = b(2);
  const Index x5 = b(3);
  const Index c2 = b(4);
  const Index x4 = b(5);
  // tensor label: I708
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x2, a1, x5, c2, x4)]);
  std::fill_n(odata.get(), out()->get_size(x3, x2, a1, x5, c2, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, a1, x5, c2, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, a1, x5, c2, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a3, c2, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a3, c2, x4)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a3.size(), c2.size(), x4.size());
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x3, x2, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x3, x2, a1)]);
    sort_indices<0,1,2,3,0,1,1,2>(i1data, i1data_sorted, a3.size(), x3.size(), x2.size(), a1.size());
    dgemm_("T", "N", x5.size()*c2.size()*x4.size(), x3.size()*x2.size()*a1.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x5.size()*c2.size()*x4.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), c2.size(), x4.size(), x3.size(), x2.size(), a1.size());
  out()->add_block(odata, x3, x2, a1, x5, c2, x4);
}

void Task418::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        for (auto& x3 : *range_[1]) {
          // tensor label: Gamma235
          std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2, x1, x4, x3, x0);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x2, x1, x4, x3, x0)]);
          sort_indices<0,1,3,4,2,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x2.size(), x1.size(), x4.size(), x3.size(), x0.size());
          // tensor label: I714
          std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, x2, x5, c2, x4);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, x2, x5, c2, x4)]);
          sort_indices<3,2,5,0,1,4,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), x2.size(), x5.size(), c2.size(), x4.size());
          dgemm_("T", "N", x1.size()*x0.size(), a1.size()*c2.size(), x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, x1.size()*x0.size());
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a1.size(), c2.size());
  out()->add_block(odata, c2, x1, x0, a1);
}

void Task419::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index x2 = b(2);
  const Index x5 = b(3);
  const Index c2 = b(4);
  const Index x4 = b(5);
  // tensor label: I714
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a1, x2, x5, c2, x4)]);
  std::fill_n(odata.get(), out()->get_size(x3, a1, x2, x5, c2, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1, x2, x5, c2, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1, x2, x5, c2, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a3, c2, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a3, c2, x4)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a3.size(), c2.size(), x4.size());
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, a3, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, a3, x2)]);
    sort_indices<2,0,1,3,0,1,-1,2>(i1data, i1data_sorted, x3.size(), a1.size(), a3.size(), x2.size());
    dgemm_("T", "N", x5.size()*c2.size()*x4.size(), x3.size()*a1.size()*x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x5.size()*c2.size()*x4.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), c2.size(), x4.size(), x3.size(), a1.size(), x2.size());
  out()->add_block(odata, x3, a1, x2, x5, c2, x4);
}

void Task420::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& x5 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a1, x5, x4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a1, x5, x4)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a1.size(), x5.size(), x4.size());
        // tensor label: I729
        std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c3, x5, x4, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c3, x5, x4, x1, x0)]);
        sort_indices<1,2,3,0,4,5,0,1,1,1>(i1data, i1data_sorted, c2.size(), c3.size(), x5.size(), x4.size(), x1.size(), x0.size());
        dgemm_("T", "N", a1.size(), c2.size()*x1.size()*x0.size(), c3.size()*x5.size()*x4.size(),
               1.0, i0data_sorted, c3.size()*x5.size()*x4.size(), i1data_sorted, c3.size()*x5.size()*x4.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), x1.size(), x0.size());
  out()->add_block(odata, c2, x1, x0, a1);
}

void Task421::Task_local::compute() {
  const Index c2 = b(0);
  const Index c3 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I729
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c3, x5, x4, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c2, c3, x5, x4, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c3, x5, x4, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c3, x5, x4, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma240
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x2, x1, x0)]);
      sort_indices<2,3,0,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I730
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c3, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c3, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), c3.size(), x3.size(), x2.size());
      dgemm_("T", "N", x5.size()*x4.size()*x1.size()*x0.size(), c2.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x5.size()*x4.size()*x1.size()*x0.size());
    }
  }
  sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x1.size(), x0.size(), c2.size(), c3.size());
  out()->add_block(odata, c2, c3, x5, x4, x1, x0);
}

void Task422::Task_local::compute() {
  const Index c2 = b(0);
  const Index c3 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I730
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c3, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c2, c3, x3, x2), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c3, x3, x2);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, c2.size(), c3.size(), x3.size(), x2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x3, x2, c2, c3);
    sort_indices<2,3,0,1,1,1,1,2>(i1data, odata, x3.size(), x2.size(), c2.size(), c3.size());
  }
  out()->add_block(odata, c2, c3, x3, x2);
}

void Task423::Task_local::compute() {
  const Index c2 = b(0);
  const Index c3 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I729
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c3, x5, x4, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c2, c3, x5, x4, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c3, x5, x4, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c3, x5, x4, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma24
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x1, x3, x2, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x1, x3, x2, x0)]);
      sort_indices<3,4,0,1,2,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x1.size(), x3.size(), x2.size(), x0.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x3, x2, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x3, x2, c3)]);
      sort_indices<1,2,0,3,0,1,1,2>(i1data, i1data_sorted, c2.size(), x3.size(), x2.size(), c3.size());
      dgemm_("T", "N", x5.size()*x4.size()*x1.size()*x0.size(), c2.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x5.size()*x4.size()*x1.size()*x0.size());
    }
  }
  sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x1.size(), x0.size(), c2.size(), c3.size());
  out()->add_block(odata, c2, c3, x5, x4, x1, x0);
}

void Task424::Task_local::compute() {
  const Index c2 = b(0);
  const Index c3 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I729
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c3, x5, x4, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c2, c3, x5, x4, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c3, x5, x4, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c3, x5, x4, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma244
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x0, x1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x0, x1, x2)]);
      sort_indices<2,5,0,1,3,4,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x0.size(), x1.size(), x2.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c3, c2, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c3, c2, x2)]);
      sort_indices<0,3,1,2,0,1,-1,2>(i1data, i1data_sorted, x3.size(), c3.size(), c2.size(), x2.size());
      dgemm_("T", "N", x5.size()*x4.size()*x0.size()*x1.size(), c3.size()*c2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x5.size()*x4.size()*x0.size()*x1.size());
    }
  }
  sort_indices<5,4,0,1,3,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x0.size(), x1.size(), c3.size(), c2.size());
  out()->add_block(odata, c2, c3, x5, x4, x1, x0);
}

void Task425::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: Gamma240
          std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x2, x1, x0);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x2, x1, x0)]);
          sort_indices<0,1,2,3,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x2.size(), x1.size(), x0.size());
          // tensor label: I732
          std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x3, x2, c2, x5, x4);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x3, x2, c2, x5, x4)]);
          sort_indices<4,5,1,2,0,3,0,1,1,1>(i1data, i1data_sorted, a1.size(), x3.size(), x2.size(), c2.size(), x5.size(), x4.size());
          dgemm_("T", "N", x1.size()*x0.size(), a1.size()*c2.size(), x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, x1.size()*x0.size());
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a1.size(), c2.size());
  out()->add_block(odata, c2, x1, x0, a1);
}

void Task426::Task_local::compute() {
  const Index a1 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index c2 = b(3);
  const Index x5 = b(4);
  const Index x4 = b(5);
  // tensor label: I732
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x3, x2, c2, x5, x4)]);
  std::fill_n(odata.get(), out()->get_size(a1, x3, x2, c2, x5, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x3, x2, c2, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x3, x2, c2, x5, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, x5, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, x5, x4)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), x5.size(), x4.size());
    // tensor label: I733
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, a1, x3, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, a1, x3, x2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), a1.size(), x3.size(), x2.size());
    dgemm_("T", "N", c2.size()*x5.size()*x4.size(), a1.size()*x3.size()*x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, c2.size()*x5.size()*x4.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, c2.size(), x5.size(), x4.size(), a1.size(), x3.size(), x2.size());
  out()->add_block(odata, a1, x3, x2, c2, x5, x4);
}

void Task427::Task_local::compute() {
  const Index a3 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I733
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, a1, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(a3, a1, x3, x2), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a1, x3, x2);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, a3.size(), a1.size(), x3.size(), x2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x3, x2, a3, a1);
    sort_indices<2,3,0,1,1,1,-1,2>(i1data, odata, x3.size(), x2.size(), a3.size(), a1.size());
  }
  out()->add_block(odata, a3, a1, x3, x2);
}

void Task428::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: Gamma24
          std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x1, x3, x2, x0);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x1, x3, x2, x0)]);
          sort_indices<0,1,3,4,2,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x1.size(), x3.size(), x2.size(), x0.size());
          // tensor label: I738
          std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, a1, c2, x5, x4);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, a1, c2, x5, x4)]);
          sort_indices<4,5,0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), a1.size(), c2.size(), x5.size(), x4.size());
          dgemm_("T", "N", x1.size()*x0.size(), a1.size()*c2.size(), x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, x1.size()*x0.size());
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a1.size(), c2.size());
  out()->add_block(odata, c2, x1, x0, a1);
}

void Task429::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index a1 = b(2);
  const Index c2 = b(3);
  const Index x5 = b(4);
  const Index x4 = b(5);
  // tensor label: I738
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x2, a1, c2, x5, x4)]);
  std::fill_n(odata.get(), out()->get_size(x3, x2, a1, c2, x5, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, a1, c2, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, a1, c2, x5, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, x5, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, x5, x4)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), x5.size(), x4.size());
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x3, x2, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x3, x2, a1)]);
    sort_indices<0,1,2,3,0,1,-1,2>(i1data, i1data_sorted, a3.size(), x3.size(), x2.size(), a1.size());
    dgemm_("T", "N", c2.size()*x5.size()*x4.size(), x3.size()*x2.size()*a1.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, c2.size()*x5.size()*x4.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, c2.size(), x5.size(), x4.size(), x3.size(), x2.size(), a1.size());
  out()->add_block(odata, x3, x2, a1, c2, x5, x4);
}

void Task430::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: Gamma244
          std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x0, x1, x2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x0, x1, x2)]);
          sort_indices<0,1,2,5,3,4,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x0.size(), x1.size(), x2.size());
          // tensor label: I744
          std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, x2, c2, x5, x4);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, x2, c2, x5, x4)]);
          sort_indices<4,5,0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), x2.size(), c2.size(), x5.size(), x4.size());
          dgemm_("T", "N", x0.size()*x1.size(), a1.size()*c2.size(), x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, x0.size()*x1.size());
        }
      }
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), c2.size());
  out()->add_block(odata, c2, x1, x0, a1);
}

void Task431::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index x2 = b(2);
  const Index c2 = b(3);
  const Index x5 = b(4);
  const Index x4 = b(5);
  // tensor label: I744
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a1, x2, c2, x5, x4)]);
  std::fill_n(odata.get(), out()->get_size(x3, a1, x2, c2, x5, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1, x2, c2, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1, x2, c2, x5, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, x5, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, x5, x4)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), x5.size(), x4.size());
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, a3, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, a3, x2)]);
    sort_indices<2,0,1,3,0,1,1,2>(i1data, i1data_sorted, x3.size(), a1.size(), a3.size(), x2.size());
    dgemm_("T", "N", c2.size()*x5.size()*x4.size(), x3.size()*a1.size()*x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, c2.size()*x5.size()*x4.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, c2.size(), x5.size(), x4.size(), x3.size(), a1.size(), x2.size());
  out()->add_block(odata, x3, a1, x2, c2, x5, x4);
}

void Task432::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a1, x6, x5);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a1, x6, x5)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x7.size(), a1.size(), x6.size(), x5.size());
        // tensor label: I759
        std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x7, x0, x6, x5, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x7, x0, x6, x5, x1)]);
        sort_indices<1,3,4,0,2,5,0,1,1,1>(i1data, i1data_sorted, c2.size(), x7.size(), x0.size(), x6.size(), x5.size(), x1.size());
        dgemm_("T", "N", a1.size(), c2.size()*x0.size()*x1.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,3,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), x0.size(), x1.size());
  out()->add_block(odata, c2, x1, x0, a1);
}

void Task433::Task_local::compute() {
  const Index c2 = b(0);
  const Index x7 = b(1);
  const Index x0 = b(2);
  const Index x6 = b(3);
  const Index x5 = b(4);
  const Index x1 = b(5);
  // tensor label: I759
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x7, x0, x6, x5, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x7, x0, x6, x5, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x7, x0, x6, x5, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x7, x0, x6, x5, x1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma250
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x5, x1, x4, x3, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x0, x6, x5, x1, x4, x3, x2)]);
        sort_indices<5,6,7,0,1,2,3,4,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x6.size(), x5.size(), x1.size(), x4.size(), x3.size(), x2.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x4, x3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x4, x3, x2)]);
        sort_indices<1,2,3,0,0,1,1,2>(i1data, i1data_sorted, c2.size(), x4.size(), x3.size(), x2.size());
        dgemm_("T", "N", x7.size()*x0.size()*x6.size()*x5.size()*x1.size(), c2.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x7.size()*x0.size()*x6.size()*x5.size()*x1.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x7.size(), x0.size(), x6.size(), x5.size(), x1.size(), c2.size());
  out()->add_block(odata, c2, x7, x0, x6, x5, x1);
}

void Task434::Task_local::compute() {
  const Index c2 = b(0);
  const Index x7 = b(1);
  const Index x0 = b(2);
  const Index x6 = b(3);
  const Index x5 = b(4);
  const Index x1 = b(5);
  // tensor label: I759
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x7, x0, x6, x5, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x7, x0, x6, x5, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x7, x0, x6, x5, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x7, x0, x6, x5, x1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma251
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x5, x4, x3, x1, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x0, x6, x5, x4, x3, x1, x2)]);
        sort_indices<4,5,7,0,1,2,3,6,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x6.size(), x5.size(), x4.size(), x3.size(), x1.size(), x2.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x3, c2, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x3, c2, x2)]);
        sort_indices<0,1,3,2,0,1,1,2>(i1data, i1data_sorted, x4.size(), x3.size(), c2.size(), x2.size());
        dgemm_("T", "N", x7.size()*x0.size()*x6.size()*x5.size()*x1.size(), c2.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x7.size()*x0.size()*x6.size()*x5.size()*x1.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x7.size(), x0.size(), x6.size(), x5.size(), x1.size(), c2.size());
  out()->add_block(odata, c2, x7, x0, x6, x5, x1);
}

void Task435::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x2, c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, x2, c2, a1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, a3.size(), x2.size(), c2.size(), a1.size());
      // tensor label: I765
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), x1.size()*x0.size(), a3.size()*x2.size(),
             1.0, i0data_sorted, a3.size()*x2.size(), i1data_sorted, a3.size()*x2.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x1.size(), x0.size());
  out()->add_block(odata, c2, x1, x0, a1);
}

void Task436::Task_local::compute() {
  const Index a3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I765
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a3, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma252
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2, x4, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x2, x4, x3, x1, x0)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a3, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a3, x4, x3)]);
        sort_indices<0,2,3,1,0,1,-1,1>(i1data, i1data_sorted, x5.size(), a3.size(), x4.size(), x3.size());
        dgemm_("T", "N", x2.size()*x1.size()*x0.size(), a3.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), a3.size());
  out()->add_block(odata, a3, x2, x1, x0);
}

void Task437::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x2, a3, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x2, a3, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), x2.size(), a3.size(), a1.size());
      // tensor label: I768
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x0, x1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x0, x1, x2)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a3.size(), x0.size(), x1.size(), x2.size());
      dgemm_("T", "N", c2.size()*a1.size(), x0.size()*x1.size(), a3.size()*x2.size(),
             1.0, i0data_sorted, a3.size()*x2.size(), i1data_sorted, a3.size()*x2.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x0.size(), x1.size());
  out()->add_block(odata, c2, x1, x0, a1);
}

void Task438::Task_local::compute() {
  const Index a3 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: I768
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, x0, x1, x2)]);
  std::fill_n(odata.get(), out()->get_size(a3, x0, x1, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x0, x1, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x0, x1, x2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma31
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x1, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x1, x2)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x1.size(), x2.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a3, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a3, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x5.size(), a3.size(), x4.size(), x3.size());
        dgemm_("T", "N", x0.size()*x1.size()*x2.size(), a3.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x0.size()*x1.size()*x2.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), a3.size());
  out()->add_block(odata, a3, x0, x1, x2);
}

void Task439::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a3, c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a3, c2, a1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a3.size(), c2.size(), a1.size());
      // tensor label: I807
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x5, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x5, x1, x0)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), x5.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), x1.size()*x0.size(), a3.size()*x5.size(),
             1.0, i0data_sorted, a3.size()*x5.size(), i1data_sorted, a3.size()*x5.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x1.size(), x0.size());
  out()->add_block(odata, c2, x1, x0, a1);
}

void Task440::Task_local::compute() {
  const Index a3 = b(0);
  const Index x5 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I807
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, x5, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a3, x5, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x5, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x5, x1, x0), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma240
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x2, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x2, x1, x0)]);
        sort_indices<1,2,3,0,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x2.size(), x1.size(), x0.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x4, x3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x4, x3, x2)]);
        sort_indices<1,2,3,0,0,1,-1,2>(i1data, i1data_sorted, a3.size(), x4.size(), x3.size(), x2.size());
        dgemm_("T", "N", x5.size()*x1.size()*x0.size(), a3.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x5.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x1.size(), x0.size(), a3.size());
  out()->add_block(odata, a3, x5, x1, x0);
}

void Task441::Task_local::compute() {
  const Index a3 = b(0);
  const Index x5 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I807
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, x5, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a3, x5, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x5, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x5, x1, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma252
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2, x4, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x2, x4, x3, x1, x0)]);
        sort_indices<1,2,3,0,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x3, a3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x3, a3, x2)]);
        sort_indices<3,0,1,2,0,1,-1,2>(i1data, i1data_sorted, x4.size(), x3.size(), a3.size(), x2.size());
        dgemm_("T", "N", x5.size()*x1.size()*x0.size(), a3.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x5.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x1.size(), x0.size(), a3.size());
  out()->add_block(odata, a3, x5, x1, x0);
}

void Task442::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, c2, a3)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c2.size(), a3.size());
      // tensor label: I810
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x5, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x5, x0, x1)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), x5.size(), x0.size(), x1.size());
      dgemm_("T", "N", a1.size()*c2.size(), x0.size()*x1.size(), a3.size()*x5.size(),
             1.0, i0data_sorted, a3.size()*x5.size(), i1data_sorted, a3.size()*x5.size(),
             1.0, odata_sorted, a1.size()*c2.size());
    }
  }
  sort_indices<1,3,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), x0.size(), x1.size());
  out()->add_block(odata, c2, x1, x0, a1);
}

void Task443::Task_local::compute() {
  const Index a3 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I810
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, x5, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a3, x5, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x5, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x5, x0, x1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma230
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x1, x4, x3, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x1, x4, x3, x2)]);
        sort_indices<3,4,5,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x1.size(), x4.size(), x3.size(), x2.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x4, x3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x4, x3, x2)]);
        sort_indices<1,2,3,0,0,1,1,2>(i1data, i1data_sorted, a3.size(), x4.size(), x3.size(), x2.size());
        dgemm_("T", "N", x5.size()*x0.size()*x1.size(), a3.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x5.size()*x0.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x0.size(), x1.size(), a3.size());
  out()->add_block(odata, a3, x5, x0, x1);
}

void Task444::Task_local::compute() {
  const Index a3 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I810
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, x5, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a3, x5, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x5, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x5, x0, x1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma31
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x1, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x1, x2)]);
        sort_indices<2,3,5,0,1,4,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x1.size(), x2.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x3, a3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x3, a3, x2)]);
        sort_indices<0,1,3,2,0,1,1,2>(i1data, i1data_sorted, x4.size(), x3.size(), a3.size(), x2.size());
        dgemm_("T", "N", x5.size()*x0.size()*x1.size(), a3.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x5.size()*x0.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x0.size(), x1.size(), a3.size());
  out()->add_block(odata, a3, x5, x0, x1);
}

void Task445::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        for (auto& x3 : *range_[1]) {
          // tensor label: Gamma276
          std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x2, x1, x3);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x2, x1, x3)]);
          sort_indices<0,2,3,5,1,4,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x2.size(), x1.size(), x3.size());
          // tensor label: I837
          std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x3, x2, x5, a1, x4);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x3, x2, x5, a1, x4)]);
          sort_indices<3,5,2,1,0,4,0,1,1,1>(i1data, i1data_sorted, c2.size(), x3.size(), x2.size(), x5.size(), a1.size(), x4.size());
          dgemm_("T", "N", x0.size()*x1.size(), c2.size()*a1.size(), x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, x0.size()*x1.size());
        }
      }
    }
  }
  sort_indices<2,1,0,3,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c2.size(), a1.size());
  out()->add_block(odata, c2, x1, x0, a1);
}

void Task446::Task_local::compute() {
  const Index c2 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x5 = b(3);
  const Index a1 = b(4);
  const Index x4 = b(5);
  // tensor label: I837
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x3, x2, x5, a1, x4)]);
  std::fill_n(odata.get(), out()->get_size(c2, x3, x2, x5, a1, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x3, x2, x5, a1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x3, x2, x5, a1, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a3.size());
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x3, a3, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x3, a3, x2)]);
    sort_indices<2,0,1,3,0,1,-2,1>(i1data, i1data_sorted, c2.size(), x3.size(), a3.size(), x2.size());
    dgemm_("T", "N", x5.size()*a1.size()*x4.size(), c2.size()*x3.size()*x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x5.size()*a1.size()*x4.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), c2.size(), x3.size(), x2.size());
  out()->add_block(odata, c2, x3, x2, x5, a1, x4);
}

void Task447::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: Gamma568
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x1, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x1, x4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x1.size(), x4.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a1, c2, x4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a1, c2, x4)]);
      sort_indices<0,3,1,2,0,1,-1,1>(i1data, i1data_sorted, x5.size(), a1.size(), c2.size(), x4.size());
      dgemm_("T", "N", x0.size()*x1.size(), a1.size()*c2.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), c2.size());
  out()->add_block(odata, c2, x1, x0, a1);
}

void Task448::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: Gamma569
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x1.size(), x0.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, x5, x4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, x5, x4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size(), x5.size(), x4.size());
      dgemm_("T", "N", x1.size()*x0.size(), c2.size()*a1.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a1.size());
  out()->add_block(odata, c2, x1, x0, a1);
}

void Task449::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      // tensor label: Gamma572
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x1, x6);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x0, x1, x6)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x1.size(), x6.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x7, a1, c2, x6);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x7, a1, c2, x6)]);
      sort_indices<0,3,1,2,0,1,-1,2>(i1data, i1data_sorted, x7.size(), a1.size(), c2.size(), x6.size());
      dgemm_("T", "N", x0.size()*x1.size(), a1.size()*c2.size(), x7.size()*x6.size(),
             1.0, i0data_sorted, x7.size()*x6.size(), i1data_sorted, x7.size()*x6.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), c2.size());
  out()->add_block(odata, c2, x1, x0, a1);
}

#endif
