//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_tasks4.cc
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

#include <src/smith/mrci/MRCI_tasks4.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task150::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: Gamma5
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x3, x1, x0)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
    // tensor label: I16
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), x3.size());
    dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, x2.size()*x1.size()*x0.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task151::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  // tensor label: I16
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(c1, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c1, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c1, x3)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c1.size(), x3.size());
      // tensor label: h1
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2)]);
      sort_indices<1,0,0,1,2,1>(i1data, i1data_sorted, a3.size(), c2.size());
      dgemm_("T", "N", c1.size()*x3.size(), 1, a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, c1.size()*x3.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c1.size(), x3.size());
  out()->add_block(odata, c1, x3);
}

void Task152::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  // tensor label: I16
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(c1, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, c2, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, c2, x3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), c2.size(), x3.size());
      // tensor label: h1
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2)]);
      sort_indices<0,1,0,1,-1,1>(i1data, i1data_sorted, a3.size(), c2.size());
      dgemm_("T", "N", c1.size()*x3.size(), 1, a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, c1.size()*x3.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c1.size(), x3.size());
  out()->add_block(odata, c1, x3);
}

void Task153::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  // tensor label: I16
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(c1, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      for (auto& c2 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a3, c2, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a3, c2, x3)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a3.size(), c2.size(), x3.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c4, a3, c2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c4, a3, c2)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), c4.size(), a3.size(), c2.size());
        dgemm_("T", "N", x3.size(), c1.size(), c4.size()*a3.size()*c2.size(),
               1.0, i0data_sorted, c4.size()*a3.size()*c2.size(), i1data_sorted, c4.size()*a3.size()*c2.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), c1.size());
  out()->add_block(odata, c1, x3);
}

void Task154::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  // tensor label: I16
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(c1, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      for (auto& c4 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c4, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c4, x3)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c4.size(), x3.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c4, a3, c2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c4, a3, c2)]);
        sort_indices<3,2,1,0,0,1,-2,1>(i1data, i1data_sorted, c1.size(), c4.size(), a3.size(), c2.size());
        dgemm_("T", "N", x3.size(), c1.size(), c4.size()*a3.size()*c2.size(),
               1.0, i0data_sorted, c4.size()*a3.size()*c2.size(), i1data_sorted, c4.size()*a3.size()*c2.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), c1.size());
  out()->add_block(odata, c1, x3);
}

void Task155::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  // tensor label: I16
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(c1, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c2, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c2, a3)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c2.size(), a3.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x3, a3, c2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x3, a3, c2)]);
        sort_indices<0,3,2,1,0,1,4,1>(i1data, i1data_sorted, a4.size(), x3.size(), a3.size(), c2.size());
        dgemm_("T", "N", c1.size(), x3.size(), a4.size()*a3.size()*c2.size(),
               1.0, i0data_sorted, a4.size()*a3.size()*c2.size(), i1data_sorted, a4.size()*a3.size()*c2.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c1.size(), x3.size());
  out()->add_block(odata, c1, x3);
}

void Task156::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  // tensor label: I16
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(c1, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a4 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, c2, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, c2, a4)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), c2.size(), a4.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x3, a3, c2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x3, a3, c2)]);
        sort_indices<2,3,0,1,0,1,-2,1>(i1data, i1data_sorted, a4.size(), x3.size(), a3.size(), c2.size());
        dgemm_("T", "N", c1.size(), x3.size(), a4.size()*a3.size()*c2.size(),
               1.0, i0data_sorted, a4.size()*a3.size()*c2.size(), i1data_sorted, a4.size()*a3.size()*c2.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c1.size(), x3.size());
  out()->add_block(odata, c1, x3);
}

void Task157::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: Gamma7
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x3, x2, x4, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x3, x2, x4, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x3.size(), x2.size(), x4.size(), x1.size(), x0.size());
        // tensor label: I22
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x5, c1, x4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x5, c1, x4)]);
        sort_indices<1,0,3,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x5.size(), c1.size(), x4.size());
        dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x3.size()*x5.size()*x4.size(),
               1.0, i0data_sorted, x3.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x5.size()*x4.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task158::Task_local::compute() {
  const Index x3 = b(0);
  const Index x5 = b(1);
  const Index c1 = b(2);
  const Index x4 = b(3);
  // tensor label: I22
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x5, c1, x4)]);
  std::fill_n(odata.get(), out()->get_size(x3, x5, c1, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x5, c1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x5, c1, x4), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, c1, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, c1, x4)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), c1.size(), x4.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x3)]);
    sort_indices<0,1,0,1,-1,1>(i1data, i1data_sorted, a2.size(), x3.size());
    dgemm_("T", "N", x5.size()*c1.size()*x4.size(), x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x5.size()*c1.size()*x4.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), c1.size(), x4.size(), x3.size());
  out()->add_block(odata, x3, x5, c1, x4);
}

void Task159::Task_local::compute() {
  const Index x3 = b(0);
  const Index x5 = b(1);
  const Index c1 = b(2);
  const Index x4 = b(3);
  // tensor label: I22
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x5, c1, x4)]);
  std::fill_n(odata.get(), out()->get_size(x3, x5, c1, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x5, c1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x5, c1, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a3, c2, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a3, c2, x4)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a3.size(), c2.size(), x4.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x3, c1, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x3, c1, c2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, a3.size(), x3.size(), c1.size(), c2.size());
      dgemm_("T", "N", x5.size()*x4.size(), x3.size()*c1.size(), a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, x5.size()*x4.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), c1.size());
  out()->add_block(odata, x3, x5, c1, x4);
}

void Task160::Task_local::compute() {
  const Index x3 = b(0);
  const Index x5 = b(1);
  const Index c1 = b(2);
  const Index x4 = b(3);
  // tensor label: I22
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x5, c1, x4)]);
  std::fill_n(odata.get(), out()->get_size(x3, x5, c1, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x5, c1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x5, c1, x4), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, c1, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, c1, a3)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), c1.size(), a3.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x4, a2, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x4, a2, x3)]);
      sort_indices<2,0,1,3,0,1,-1,1>(i1data, i1data_sorted, a3.size(), x4.size(), a2.size(), x3.size());
      dgemm_("T", "N", x5.size()*c1.size(), x4.size()*x3.size(), a3.size()*a2.size(),
             1.0, i0data_sorted, a3.size()*a2.size(), i1data_sorted, a3.size()*a2.size(),
             1.0, odata_sorted, x5.size()*c1.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), c1.size(), x4.size(), x3.size());
  out()->add_block(odata, x3, x5, c1, x4);
}

void Task161::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      for (auto& x6 : *range_[1]) {
        for (auto& x4 : *range_[1]) {
          for (auto& x3 : *range_[1]) {
            // tensor label: Gamma97
            std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x7, x5, x6, x4, x3, x1, x0);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x7, x5, x6, x4, x3, x1, x0)]);
            sort_indices<1,2,3,4,5,0,6,7,0,1,1,1>(i0data, i0data_sorted, x2.size(), x7.size(), x5.size(), x6.size(), x4.size(), x3.size(), x1.size(), x0.size());
            // tensor label: I300
            std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x3, c1, x7, x6);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x3, c1, x7, x6)]);
            sort_indices<4,0,5,1,2,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x3.size(), c1.size(), x7.size(), x6.size());
            dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x5.size()*x4.size()*x3.size()*x7.size()*x6.size(),
                   1.0, i0data_sorted, x5.size()*x4.size()*x3.size()*x7.size()*x6.size(), i1data_sorted, x5.size()*x4.size()*x3.size()*x7.size()*x6.size(),
                   1.0, odata_sorted, x2.size()*x1.size()*x0.size());
          }
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task162::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index c1 = b(3);
  const Index x7 = b(4);
  const Index x6 = b(5);
  // tensor label: I300
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x4, x3, c1, x7, x6)]);
  std::fill_n(odata.get(), out()->get_size(x5, x4, x3, c1, x7, x6), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x5, x4, x3, c1, x7, x6)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x5, x4, x3, c1, x7, x6), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x7, c2, x6);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x7, c2, x6)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x7.size(), c2.size(), x6.size());
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x5, c2, x4, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, c2, x4, x3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), c2.size(), x4.size(), x3.size());
    dgemm_("T", "N", c1.size()*x7.size()*x6.size(), x5.size()*x4.size()*x3.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, c1.size()*x7.size()*x6.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), x7.size(), x6.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, x5, x4, x3, c1, x7, x6);
}

void Task163::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x6 : *range_[1]) {
        for (auto& x5 : *range_[1]) {
          for (auto& x4 : *range_[1]) {
            // tensor label: Gamma98
            std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x7, x3, x6, x5, x4, x1, x0);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x7, x3, x6, x5, x4, x1, x0)]);
            sort_indices<1,2,3,4,5,0,6,7,0,1,1,1>(i0data, i0data_sorted, x2.size(), x7.size(), x3.size(), x6.size(), x5.size(), x4.size(), x1.size(), x0.size());
            // tensor label: I303
            std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x3, c1, x7, x6);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x3, c1, x7, x6)]);
            sort_indices<4,2,5,0,1,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x3.size(), c1.size(), x7.size(), x6.size());
            dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x5.size()*x4.size()*x3.size()*x7.size()*x6.size(),
                   1.0, i0data_sorted, x5.size()*x4.size()*x3.size()*x7.size()*x6.size(), i1data_sorted, x5.size()*x4.size()*x3.size()*x7.size()*x6.size(),
                   1.0, odata_sorted, x2.size()*x1.size()*x0.size());
          }
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task164::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index c1 = b(3);
  const Index x7 = b(4);
  const Index x6 = b(5);
  // tensor label: I303
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x4, x3, c1, x7, x6)]);
  std::fill_n(odata.get(), out()->get_size(x5, x4, x3, c1, x7, x6), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x5, x4, x3, c1, x7, x6)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x5, x4, x3, c1, x7, x6), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x7, c2, x6);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x7, c2, x6)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x7.size(), c2.size(), x6.size());
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x3, c2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x3, c2)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x3.size(), c2.size());
    dgemm_("T", "N", c1.size()*x7.size()*x6.size(), x5.size()*x4.size()*x3.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, c1.size()*x7.size()*x6.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), x7.size(), x6.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, x5, x4, x3, c1, x7, x6);
}

void Task165::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        for (auto& x4 : *range_[1]) {
          for (auto& x3 : *range_[1]) {
            // tensor label: Gamma100
            std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, x2, x5, x4, x3, x1, x0);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, x2, x5, x4, x3, x1, x0)]);
            sort_indices<0,1,3,4,5,2,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x2.size(), x5.size(), x4.size(), x3.size(), x1.size(), x0.size());
            // tensor label: I309
            std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x4, x3, x7, x6, x5);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x4, x3, x7, x6, x5)]);
            sort_indices<3,4,5,1,2,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), x4.size(), x3.size(), x7.size(), x6.size(), x5.size());
            dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, i0data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, odata_sorted, x2.size()*x1.size()*x0.size());
          }
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task166::Task_local::compute() {
  const Index c1 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x7 = b(3);
  const Index x6 = b(4);
  const Index x5 = b(5);
  // tensor label: I309
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x4, x3, x7, x6, x5)]);
  std::fill_n(odata.get(), out()->get_size(c1, x4, x3, x7, x6, x5), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x4, x3, x7, x6, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x4, x3, x7, x6, x5), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, c2, x5);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, c2, x5)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), c2.size(), x5.size());
    // tensor label: I310
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c2, x4, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c2, x4, x3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), c2.size(), x4.size(), x3.size());
    dgemm_("T", "N", x7.size()*x6.size()*x5.size(), c1.size()*x4.size()*x3.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, x7.size()*x6.size()*x5.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x7.size(), x6.size(), x5.size(), c1.size(), x4.size(), x3.size());
  out()->add_block(odata, c1, x4, x3, x7, x6, x5);
}

void Task167::Task_local::compute() {
  const Index c1 = b(0);
  const Index c2 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I310
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c2, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(c1, c2, x4, x3), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c2, x4, x3);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, c1.size(), c2.size(), x4.size(), x3.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x4, x3, c1, c2);
    sort_indices<2,3,0,1,1,1,-1,2>(i1data, odata, x4.size(), x3.size(), c1.size(), c2.size());
  }
  out()->add_block(odata, c1, c2, x4, x3);
}

void Task168::Task_local::compute() {
  const Index c1 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x7 = b(3);
  const Index x6 = b(4);
  const Index x5 = b(5);
  // tensor label: I309
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x4, x3, x7, x6, x5)]);
  std::fill_n(odata.get(), out()->get_size(c1, x4, x3, x7, x6, x5), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x4, x3, x7, x6, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x4, x3, x7, x6, x5), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x7, x6);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x7, x6)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x7.size(), x6.size());
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x5, x4, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x5, x4, x3)]);
    sort_indices<0,1,2,3,0,1,1,2>(i1data, i1data_sorted, a2.size(), x5.size(), x4.size(), x3.size());
    dgemm_("T", "N", c1.size()*x7.size()*x6.size(), x5.size()*x4.size()*x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, c1.size()*x7.size()*x6.size());
  }
  sort_indices<0,4,5,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), x7.size(), x6.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, c1, x4, x3, x7, x6, x5);
}

void Task169::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        for (auto& x5 : *range_[1]) {
          for (auto& x4 : *range_[1]) {
            // tensor label: Gamma101
            std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, x3, x5, x2, x4, x1, x0);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, x3, x5, x2, x4, x1, x0)]);
            sort_indices<0,1,2,3,5,4,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x3.size(), x5.size(), x2.size(), x4.size(), x1.size(), x0.size());
            // tensor label: I312
            std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x4, x3, x7, x6, x5);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x4, x3, x7, x6, x5)]);
            sort_indices<3,4,2,5,1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), x4.size(), x3.size(), x7.size(), x6.size(), x5.size());
            dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, i0data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, odata_sorted, x2.size()*x1.size()*x0.size());
          }
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task170::Task_local::compute() {
  const Index c1 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x7 = b(3);
  const Index x6 = b(4);
  const Index x5 = b(5);
  // tensor label: I312
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x4, x3, x7, x6, x5)]);
  std::fill_n(odata.get(), out()->get_size(c1, x4, x3, x7, x6, x5), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x4, x3, x7, x6, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x4, x3, x7, x6, x5), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, c2, x5);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, c2, x5)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), c2.size(), x5.size());
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x4, x3, c2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x4, x3, c2)]);
    sort_indices<3,0,1,2,0,1,-1,2>(i1data, i1data_sorted, c1.size(), x4.size(), x3.size(), c2.size());
    dgemm_("T", "N", x7.size()*x6.size()*x5.size(), c1.size()*x4.size()*x3.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, x7.size()*x6.size()*x5.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x7.size(), x6.size(), x5.size(), c1.size(), x4.size(), x3.size());
  out()->add_block(odata, c1, x4, x3, x7, x6, x5);
}

void Task171::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        for (auto& x5 : *range_[1]) {
          for (auto& x3 : *range_[1]) {
            // tensor label: Gamma102
            std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, x4, x5, x2, x3, x1, x0);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, x4, x5, x2, x3, x1, x0)]);
            sort_indices<0,1,2,3,5,4,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x4.size(), x5.size(), x2.size(), x3.size(), x1.size(), x0.size());
            // tensor label: I315
            std::unique_ptr<double[]> i1data = in(1)->get_block(x4, c1, x3, x7, x6, x5);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, c1, x3, x7, x6, x5)]);
            sort_indices<3,4,0,5,2,1,0,1,1,1>(i1data, i1data_sorted, x4.size(), c1.size(), x3.size(), x7.size(), x6.size(), x5.size());
            dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, i0data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, odata_sorted, x2.size()*x1.size()*x0.size());
          }
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task172::Task_local::compute() {
  const Index x4 = b(0);
  const Index c1 = b(1);
  const Index x3 = b(2);
  const Index x7 = b(3);
  const Index x6 = b(4);
  const Index x5 = b(5);
  // tensor label: I315
  std::unique_ptr<double[]> odata(new double[out()->get_size(x4, c1, x3, x7, x6, x5)]);
  std::fill_n(odata.get(), out()->get_size(x4, c1, x3, x7, x6, x5), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, c1, x3, x7, x6, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, c1, x3, x7, x6, x5), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, c2, x5);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, c2, x5)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), c2.size(), x5.size());
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x4, c2, c1, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, c2, c1, x3)]);
    sort_indices<1,0,2,3,0,1,1,2>(i1data, i1data_sorted, x4.size(), c2.size(), c1.size(), x3.size());
    dgemm_("T", "N", x7.size()*x6.size()*x5.size(), x4.size()*c1.size()*x3.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, x7.size()*x6.size()*x5.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x7.size(), x6.size(), x5.size(), x4.size(), c1.size(), x3.size());
  out()->add_block(odata, x4, c1, x3, x7, x6, x5);
}

void Task173::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma104
        std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x5, x4, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x5, x4, x3, x1, x0)]);
        sort_indices<1,2,3,0,4,5,0,1,1,1>(i0data, i0data_sorted, x2.size(), x5.size(), x4.size(), x3.size(), x1.size(), x0.size());
        // tensor label: I321
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x3, c1, x5);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x3, c1, x5)]);
        sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x4.size(), x3.size(), c1.size(), x5.size());
        dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x4.size()*x3.size()*x5.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x5.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task174::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index c1 = b(2);
  const Index x5 = b(3);
  // tensor label: I321
  std::unique_ptr<double[]> odata(new double[out()->get_size(x4, x3, c1, x5)]);
  std::fill_n(odata.get(), out()->get_size(x4, x3, c1, x5), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, x3, c1, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x3, c1, x5), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c1, x5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c1, x5)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c1.size(), x5.size());
      // tensor label: I322
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, x4, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, x4, x3)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), x4.size(), x3.size());
      dgemm_("T", "N", c1.size()*x5.size(), x4.size()*x3.size(), a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, c1.size()*x5.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, x4, x3, c1, x5);
}

void Task175::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I322
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x4, x3), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, c2, x4, x3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, a3.size(), c2.size(), x4.size(), x3.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x4, x3, a3, c2);
    sort_indices<2,3,0,1,1,1,1,1>(i1data, odata, x4.size(), x3.size(), a3.size(), c2.size());
  }
  out()->add_block(odata, a3, c2, x4, x3);
}

void Task176::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index c1 = b(2);
  const Index x5 = b(3);
  // tensor label: I321
  std::unique_ptr<double[]> odata(new double[out()->get_size(x4, x3, c1, x5)]);
  std::fill_n(odata.get(), out()->get_size(x4, x3, c1, x5), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, x3, c1, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x3, c1, x5), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, c2, x5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, c2, x5)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), c2.size(), x5.size());
      // tensor label: I325
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, x4, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, x4, x3)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), x4.size(), x3.size());
      dgemm_("T", "N", c1.size()*x5.size(), x4.size()*x3.size(), a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, c1.size()*x5.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, x4, x3, c1, x5);
}

void Task177::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I325
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x4, x3), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, c2, x4, x3);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, a3.size(), c2.size(), x4.size(), x3.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x4, x3, a3, c2);
    sort_indices<2,3,0,1,1,1,-1,2>(i1data, odata, x4.size(), x3.size(), a3.size(), c2.size());
  }
  out()->add_block(odata, a3, c2, x4, x3);
}

void Task178::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index c1 = b(2);
  const Index x5 = b(3);
  // tensor label: I321
  std::unique_ptr<double[]> odata(new double[out()->get_size(x4, x3, c1, x5)]);
  std::fill_n(odata.get(), out()->get_size(x4, x3, c1, x5), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, x3, c1, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x3, c1, x5), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c1, x5)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c1.size(), x5.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x4, c3, a2, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, c3, a2, x3)]);
      sort_indices<1,2,0,3,0,1,-1,2>(i1data, i1data_sorted, x4.size(), c3.size(), a2.size(), x3.size());
      dgemm_("T", "N", c1.size()*x5.size(), x4.size()*x3.size(), c3.size()*a2.size(),
             1.0, i0data_sorted, c3.size()*a2.size(), i1data_sorted, c3.size()*a2.size(),
             1.0, odata_sorted, c1.size()*x5.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, x4, x3, c1, x5);
}

void Task179::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: Gamma107
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x5, x2, x4, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x5, x2, x4, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x3.size(), x5.size(), x2.size(), x4.size(), x1.size(), x0.size());
        // tensor label: I330
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x3, c1, x5);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x3, c1, x5)]);
        sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, x4.size(), x3.size(), c1.size(), x5.size());
        dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x4.size()*x3.size()*x5.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x5.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task180::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index c1 = b(2);
  const Index x5 = b(3);
  // tensor label: I330
  std::unique_ptr<double[]> odata(new double[out()->get_size(x4, x3, c1, x5)]);
  std::fill_n(odata.get(), out()->get_size(x4, x3, c1, x5), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, x3, c1, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x3, c1, x5), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, c2, x5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, c2, x5)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), c2.size(), x5.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x4, x3, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x4, x3, c2)]);
      sort_indices<0,3,1,2,0,1,-1,2>(i1data, i1data_sorted, a3.size(), x4.size(), x3.size(), c2.size());
      dgemm_("T", "N", c1.size()*x5.size(), x4.size()*x3.size(), a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, c1.size()*x5.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, x4, x3, c1, x5);
}

void Task181::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma109
        std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x5, x2, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x4, x5, x2, x3, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x4.size(), x5.size(), x2.size(), x3.size(), x1.size(), x0.size());
        // tensor label: I336
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x3, c1, x5);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x3, c1, x5)]);
        sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x4.size(), x3.size(), c1.size(), x5.size());
        dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x4.size()*x3.size()*x5.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x5.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task182::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index c1 = b(2);
  const Index x5 = b(3);
  // tensor label: I336
  std::unique_ptr<double[]> odata(new double[out()->get_size(x4, x3, c1, x5)]);
  std::fill_n(odata.get(), out()->get_size(x4, x3, c1, x5), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, x3, c1, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x3, c1, x5), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x5)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x5.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x4, c3, a2, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, c3, a2, x3)]);
      sort_indices<2,1,0,3,0,1,1,2>(i1data, i1data_sorted, x4.size(), c3.size(), a2.size(), x3.size());
      dgemm_("T", "N", c1.size()*x5.size(), x4.size()*x3.size(), c3.size()*a2.size(),
             1.0, i0data_sorted, c3.size()*a2.size(), i1data_sorted, c3.size()*a2.size(),
             1.0, odata_sorted, c1.size()*x5.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, x4, x3, c1, x5);
}

void Task183::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      for (auto& x6 : *range_[1]) {
        for (auto& x4 : *range_[1]) {
          for (auto& x3 : *range_[1]) {
            // tensor label: Gamma114
            std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x5, x2, x6, x4, x3, x1, x0);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x5, x2, x6, x4, x3, x1, x0)]);
            sort_indices<0,1,3,4,5,2,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x5.size(), x2.size(), x6.size(), x4.size(), x3.size(), x1.size(), x0.size());
            // tensor label: I351
            std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x3, x7, c1, x6);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x3, x7, c1, x6)]);
            sort_indices<3,0,5,1,2,4,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x3.size(), x7.size(), c1.size(), x6.size());
            dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x5.size()*x4.size()*x3.size()*x7.size()*x6.size(),
                   1.0, i0data_sorted, x5.size()*x4.size()*x3.size()*x7.size()*x6.size(), i1data_sorted, x5.size()*x4.size()*x3.size()*x7.size()*x6.size(),
                   1.0, odata_sorted, x2.size()*x1.size()*x0.size());
          }
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task184::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x7 = b(3);
  const Index c1 = b(4);
  const Index x6 = b(5);
  // tensor label: I351
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x4, x3, x7, c1, x6)]);
  std::fill_n(odata.get(), out()->get_size(x5, x4, x3, x7, c1, x6), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x5, x4, x3, x7, c1, x6)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x5, x4, x3, x7, c1, x6), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a2, c1, x6);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a2, c1, x6)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), a2.size(), c1.size(), x6.size());
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x5, x4, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x5, x4, x3)]);
    sort_indices<0,1,2,3,0,1,-1,2>(i1data, i1data_sorted, a2.size(), x5.size(), x4.size(), x3.size());
    dgemm_("T", "N", x7.size()*c1.size()*x6.size(), x5.size()*x4.size()*x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x7.size()*c1.size()*x6.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x7.size(), c1.size(), x6.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, x5, x4, x3, x7, c1, x6);
}

void Task185::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x6 : *range_[1]) {
        for (auto& x5 : *range_[1]) {
          for (auto& x4 : *range_[1]) {
            // tensor label: Gamma115
            std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x3, x2, x6, x5, x4, x1, x0);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x3, x2, x6, x5, x4, x1, x0)]);
            sort_indices<0,1,3,4,5,2,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x3.size(), x2.size(), x6.size(), x5.size(), x4.size(), x1.size(), x0.size());
            // tensor label: I354
            std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x3, x7, c1, x6);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x3, x7, c1, x6)]);
            sort_indices<3,2,5,0,1,4,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x3.size(), x7.size(), c1.size(), x6.size());
            dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x5.size()*x4.size()*x3.size()*x7.size()*x6.size(),
                   1.0, i0data_sorted, x5.size()*x4.size()*x3.size()*x7.size()*x6.size(), i1data_sorted, x5.size()*x4.size()*x3.size()*x7.size()*x6.size(),
                   1.0, odata_sorted, x2.size()*x1.size()*x0.size());
          }
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task186::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x7 = b(3);
  const Index c1 = b(4);
  const Index x6 = b(5);
  // tensor label: I354
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x4, x3, x7, c1, x6)]);
  std::fill_n(odata.get(), out()->get_size(x5, x4, x3, x7, c1, x6), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x5, x4, x3, x7, c1, x6)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x5, x4, x3, x7, c1, x6), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a2, c1, x6);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a2, c1, x6)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), a2.size(), c1.size(), x6.size());
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, a2, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, a2, x3)]);
    sort_indices<2,0,1,3,0,1,-1,2>(i1data, i1data_sorted, x5.size(), x4.size(), a2.size(), x3.size());
    dgemm_("T", "N", x7.size()*c1.size()*x6.size(), x5.size()*x4.size()*x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x7.size()*c1.size()*x6.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x7.size(), c1.size(), x6.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, x5, x4, x3, x7, c1, x6);
}

void Task187::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        for (auto& x4 : *range_[1]) {
          for (auto& x3 : *range_[1]) {
            // tensor label: Gamma119
            std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, x5, x4, x2, x3, x1, x0);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, x5, x4, x2, x3, x1, x0)]);
            sort_indices<0,1,2,3,5,4,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
            // tensor label: I366
            std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x3, c1, x7, x6);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x3, c1, x7, x6)]);
            sort_indices<4,5,0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x3.size(), c1.size(), x7.size(), x6.size());
            dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x5.size()*x4.size()*x3.size()*x7.size()*x6.size(),
                   1.0, i0data_sorted, x5.size()*x4.size()*x3.size()*x7.size()*x6.size(), i1data_sorted, x5.size()*x4.size()*x3.size()*x7.size()*x6.size(),
                   1.0, odata_sorted, x2.size()*x1.size()*x0.size());
          }
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task188::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index c1 = b(3);
  const Index x7 = b(4);
  const Index x6 = b(5);
  // tensor label: I366
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x4, x3, c1, x7, x6)]);
  std::fill_n(odata.get(), out()->get_size(x5, x4, x3, c1, x7, x6), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x5, x4, x3, c1, x7, x6)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x5, x4, x3, c1, x7, x6), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x7, x6);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x7, x6)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x7.size(), x6.size());
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, a2, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, a2, x3)]);
    sort_indices<2,0,1,3,0,1,1,2>(i1data, i1data_sorted, x5.size(), x4.size(), a2.size(), x3.size());
    dgemm_("T", "N", c1.size()*x7.size()*x6.size(), x5.size()*x4.size()*x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, c1.size()*x7.size()*x6.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), x7.size(), x6.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, x5, x4, x3, c1, x7, x6);
}

void Task189::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x6 : *range_[1]) {
        for (auto& x5 : *range_[1]) {
          for (auto& x4 : *range_[1]) {
            // tensor label: Gamma122
            std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x3, x6, x5, x2, x4, x1, x0);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x3, x6, x5, x2, x4, x1, x0)]);
            sort_indices<0,1,2,3,5,4,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x3.size(), x6.size(), x5.size(), x2.size(), x4.size(), x1.size(), x0.size());
            // tensor label: I375
            std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x4, x3, x7, x6, x5);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x4, x3, x7, x6, x5)]);
            sort_indices<3,2,4,5,1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), x4.size(), x3.size(), x7.size(), x6.size(), x5.size());
            dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, i0data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, odata_sorted, x2.size()*x1.size()*x0.size());
          }
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task190::Task_local::compute() {
  const Index c1 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x7 = b(3);
  const Index x6 = b(4);
  const Index x5 = b(5);
  // tensor label: I375
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x4, x3, x7, x6, x5)]);
  std::fill_n(odata.get(), out()->get_size(c1, x4, x3, x7, x6, x5), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x4, x3, x7, x6, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x4, x3, x7, x6, x5), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a2, x6, x5);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a2, x6, x5)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), a2.size(), x6.size(), x5.size());
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x4, a2, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x4, a2, x3)]);
    sort_indices<2,0,1,3,0,1,-1,1>(i1data, i1data_sorted, c1.size(), x4.size(), a2.size(), x3.size());
    dgemm_("T", "N", x7.size()*x6.size()*x5.size(), c1.size()*x4.size()*x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x7.size()*x6.size()*x5.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x7.size(), x6.size(), x5.size(), c1.size(), x4.size(), x3.size());
  out()->add_block(odata, c1, x4, x3, x7, x6, x5);
}

void Task191::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: Gamma550
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, x2, x5, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, x2, x5, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x2.size(), x5.size(), x1.size(), x0.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x7, x6, c1, x5);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x7, x6, c1, x5)]);
        sort_indices<0,1,3,2,0,1,-1,1>(i1data, i1data_sorted, x7.size(), x6.size(), c1.size(), x5.size());
        dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task192::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x9 : *range_[1]) {
    for (auto& x8 : *range_[1]) {
      for (auto& x7 : *range_[1]) {
        // tensor label: Gamma551
        std::unique_ptr<double[]> i0data = in(0)->get_block(x9, x8, x2, x7, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x9, x8, x2, x7, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x9.size(), x8.size(), x2.size(), x7.size(), x1.size(), x0.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x9, x8, c1, x7);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x9, x8, c1, x7)]);
        sort_indices<0,1,3,2,0,1,-1,2>(i1data, i1data_sorted, x9.size(), x8.size(), c1.size(), x7.size());
        dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x9.size()*x8.size()*x7.size(),
               1.0, i0data_sorted, x9.size()*x8.size()*x7.size(), i1data_sorted, x9.size()*x8.size()*x7.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task193::Task_local::compute() {
  const Index c3 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, x0, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(c3, x0, c1, a2), 0.0);
  {
    // tensor label: I27
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c3, x0, a2);
    sort_indices<1,2,0,3,1,1,1,1>(i0data, odata, c1.size(), c3.size(), x0.size(), a2.size());
  }
  out()->add_block(odata, c3, x0, c1, a2);
}

void Task194::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size());
    // tensor label: I28
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c3, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c3, x1, x0)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), c3.size(), x1.size(), x0.size());
    dgemm_("T", "N", a2.size(), c1.size()*c3.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), c3.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task195::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I28
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x3, x0, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x3, x0, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), x3.size(), x0.size(), x2.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x3, c3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x3, c3, x2)]);
      sort_indices<1,3,0,2,0,1,-2,1>(i1data, i1data_sorted, c1.size(), x3.size(), c3.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), c1.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), c3.size());
  out()->add_block(odata, c1, c3, x1, x0);
}

void Task196::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  // tensor label: h1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size());
  // tensor label: I31
  std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), x0.size());
  dgemm_("T", "N", c1.size()*a2.size(), c3.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a2.size());
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task197::Task_local::compute() {
  const Index c3 = b(0);
  const Index x0 = b(1);
  // tensor label: I31
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma10
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, c3, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, c3, x1)]);
        sort_indices<0,1,3,2,0,1,2,1>(i1data, i1data_sorted, x3.size(), x2.size(), c3.size(), x1.size());
        dgemm_("T", "N", x0.size(), c3.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c3.size());
  out()->add_block(odata, c3, x0);
}

void Task198::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  // tensor label: h1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size());
  // tensor label: I34
  std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), x0.size());
  dgemm_("T", "N", c3.size()*a2.size(), c1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c3.size()*a2.size());
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task199::Task_local::compute() {
  const Index c1 = b(0);
  const Index x0 = b(1);
  // tensor label: I34
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma10
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, c1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, c1, x1)]);
        sort_indices<0,1,3,2,0,1,-1,1>(i1data, i1data_sorted, x3.size(), x2.size(), c1.size(), x1.size());
        dgemm_("T", "N", x0.size(), c1.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c1.size());
  out()->add_block(odata, c1, x0);
}

#endif
