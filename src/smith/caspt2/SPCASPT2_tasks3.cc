//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: SPCASPT2_tasks3.cc
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

#include <src/smith/caspt2/SPCASPT2_tasks3.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::SPCASPT2;

void Task100::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a2)]);
  std::fill_n(odata.get(), out()->get_size(x3, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& c1 : *range_[0]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, c1, x4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, c1, x4)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), c1.size(), x4.size());
        // tensor label: I28
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x4, x5, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x4, x5, x3)]);
        sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), x4.size(), x5.size(), x3.size());
        dgemm_("T", "N", a2.size(), x3.size(), c1.size()*x4.size()*x5.size(),
               1.0, i0data_sorted, c1.size()*x4.size()*x5.size(), i1data_sorted, c1.size()*x4.size()*x5.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a2.size(), x3.size());
  out()->add_block(odata, x3, a2);
}

void Task101::Task_local::compute() {
  const Index c1 = b(0);
  const Index x4 = b(1);
  const Index x5 = b(2);
  const Index x3 = b(3);
  // tensor label: I28
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x4, x5, x3)]);
  std::fill_n(odata.get(), out()->get_size(c1, x4, x5, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x4, x5, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x4, x5, x3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: Gamma9
        std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x4, x1, x0, x5, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x4, x1, x0, x5, x3)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x2.size(), x4.size(), x1.size(), x0.size(), x5.size(), x3.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, c1, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, c1, x2)]);
        sort_indices<3,1,0,2,0,1,-1,1>(i1data, i1data_sorted, x0.size(), x1.size(), c1.size(), x2.size());
        dgemm_("T", "N", x4.size()*x5.size()*x3.size(), c1.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x4.size()*x5.size()*x3.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x4.size(), x5.size(), x3.size(), c1.size());
  out()->add_block(odata, c1, x4, x5, x3);
}

void Task102::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a2)]);
  std::fill_n(odata.get(), out()->get_size(x3, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a2), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& x5 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x5, x4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x5, x4)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x5.size(), x4.size());
        // tensor label: I31
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x5, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x5, x4, x3)]);
        sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), x5.size(), x4.size(), x3.size());
        dgemm_("T", "N", a2.size(), x3.size(), c1.size()*x5.size()*x4.size(),
               1.0, i0data_sorted, c1.size()*x5.size()*x4.size(), i1data_sorted, c1.size()*x5.size()*x4.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a2.size(), x3.size());
  out()->add_block(odata, x3, a2);
}

void Task103::Task_local::compute() {
  const Index c1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I31
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x5, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(c1, x5, x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x5, x4, x3), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma6
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x1, x0, x2, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x1, x0, x2, x3)]);
        sort_indices<2,3,4,0,1,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x1.size(), x0.size(), x2.size(), x3.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, c1, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, c1, x2)]);
        sort_indices<1,0,3,2,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), c1.size(), x2.size());
        dgemm_("T", "N", x5.size()*x4.size()*x3.size(), c1.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x5.size()*x4.size()*x3.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), c1.size());
  out()->add_block(odata, c1, x5, x4, x3);
}

void Task104::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a2)]);
  std::fill_n(odata.get(), out()->get_size(x3, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a2)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a2.size());
        // tensor label: I187
        std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x5, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x5, x4, x3)]);
        sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a1.size(), x5.size(), x4.size(), x3.size());
        dgemm_("T", "N", a2.size(), x3.size(), a1.size()*x5.size()*x4.size(),
               1.0, i0data_sorted, a1.size()*x5.size()*x4.size(), i1data_sorted, a1.size()*x5.size()*x4.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a2.size(), x3.size());
  out()->add_block(odata, x3, a2);
}

void Task105::Task_local::compute() {
  const Index a1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I187
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x5, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(a1, x5, x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x5, x4, x3), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma62
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x2, x1, x4, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x2, x1, x4, x3)]);
        sort_indices<1,2,3,0,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x2.size(), x1.size(), x4.size(), x3.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a1, x1, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a1, x1, x2)]);
        sort_indices<0,3,2,1,0,1,2,1>(i1data, i1data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
        dgemm_("T", "N", x5.size()*x4.size()*x3.size(), a1.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x5.size()*x4.size()*x3.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), a1.size());
  out()->add_block(odata, a1, x5, x4, x3);
}

void Task106::Task_local::compute() {
  const Index a2 = b(0);
  const Index x1 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a2, x1), 0.0);
  {
    // tensor label: I33
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2);
    sort_indices<1,0,1,1,1,1>(i0data, odata, x1.size(), a2.size());
  }
  out()->add_block(odata, a2, x1);
}

void Task107::Task_local::compute() {
  const Index x1 = b(0);
  const Index a2 = b(1);
  // tensor label: I33
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x1, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, a2), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& c3 : *range_[0]) {
      for (auto& c1 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        // tensor label: I34
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c3, x0, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c3, x0, x1)]);
        sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), c3.size(), x0.size(), x1.size());
        dgemm_("T", "N", a2.size(), x1.size(), c1.size()*c3.size()*x0.size(),
               1.0, i0data_sorted, c1.size()*c3.size()*x0.size(), i1data_sorted, c1.size()*c3.size()*x0.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a2.size(), x1.size());
  out()->add_block(odata, x1, a2);
}

void Task108::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I34
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, x1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: Gamma11
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x2, x1, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x2, x1, x3)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), x2.size(), x1.size(), x3.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x3, c3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x3, c3, x2)]);
      sort_indices<3,1,0,2,0,1,-2,1>(i1data, i1data_sorted, c1.size(), x3.size(), c3.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), c1.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c1.size(), c3.size());
  out()->add_block(odata, c1, c3, x0, x1);
}

void Task109::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1), 0.0);
  {
    // tensor label: I36
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, c1);
    sort_indices<0,1,1,1,1,1>(i0data, odata, a2.size(), c1.size());
  }
  out()->add_block(odata, a2, c1);
}

void Task110::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  // tensor label: I36
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
      // tensor label: I37
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x0)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), x0.size());
      dgemm_("T", "N", a2.size()*c1.size(), 1, c3.size()*x0.size(),
             1.0, i0data_sorted, c3.size()*x0.size(), i1data_sorted, c3.size()*x0.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size());
  out()->add_block(odata, a2, c1);
}

void Task111::Task_local::compute() {
  const Index c3 = b(0);
  const Index x0 = b(1);
  // tensor label: I37
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma12
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, c3, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, c3, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), c3.size(), x1.size());
        dgemm_("T", "N", x0.size(), c3.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c3.size());
  out()->add_block(odata, c3, x0);
}

void Task112::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c3), 0.0);
  {
    // tensor label: I39
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2);
    sort_indices<1,0,1,1,1,1>(i0data, odata, c3.size(), a2.size());
  }
  out()->add_block(odata, a2, c3);
}

void Task113::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  // tensor label: I39
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a2)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a2), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& c1 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
      // tensor label: I40
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x0)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), x0.size());
      dgemm_("T", "N", c3.size()*a2.size(), 1, c1.size()*x0.size(),
             1.0, i0data_sorted, c1.size()*x0.size(), i1data_sorted, c1.size()*x0.size(),
             1.0, odata_sorted, c3.size()*a2.size());
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c3.size());
  out()->add_block(odata, c3, a2);
}

void Task114::Task_local::compute() {
  const Index c1 = b(0);
  const Index x0 = b(1);
  // tensor label: I40
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma13
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

void Task115::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  // tensor label: I39
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a2)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
      // tensor label: I196
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c1)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), c1.size());
      dgemm_("T", "N", c3.size()*a2.size(), 1, a4.size()*c1.size(),
             1.0, i0data_sorted, a4.size()*c1.size(), i1data_sorted, a4.size()*c1.size(),
             1.0, odata_sorted, c3.size()*a2.size());
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c3.size());
  out()->add_block(odata, c3, a2);
}

void Task116::Task_local::compute() {
  const Index a4 = b(0);
  const Index c1 = b(1);
  // tensor label: I196
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, c1)]);
  std::fill_n(odata.get(), out()->get_size(a4, c1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma65
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a4, c1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a4, c1, x0)]);
      sort_indices<0,3,1,2,0,1,2,1>(i1data, i1data_sorted, x1.size(), a4.size(), c1.size(), x0.size());
      dgemm_("T", "N", 1, a4.size()*c1.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size());
  out()->add_block(odata, a4, c1);
}

void Task117::Task_local::compute() {
  const Index a4 = b(0);
  const Index c1 = b(1);
  // tensor label: I196
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, c1)]);
  std::fill_n(odata.get(), out()->get_size(a4, c1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma67
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a4, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a4, x1, x0)]);
      sort_indices<2,3,0,1,0,1,-2,1>(i1data, i1data_sorted, c1.size(), a4.size(), x1.size(), x0.size());
      dgemm_("T", "N", 1, c1.size()*a4.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size());
  out()->add_block(odata, a4, c1);
}

void Task118::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x2), 0.0);
  {
    // tensor label: I42
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x1);
    sort_indices<1,0,1,1,1,1>(i0data, odata, x2.size(), x1.size());
  }
  out()->add_block(odata, x1, x2);
}

void Task119::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  // tensor label: I42
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: Gamma14
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x3, x2, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size(), x2.size(), x1.size());
      // tensor label: I43
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size());
      dgemm_("T", "N", x2.size()*x1.size(), 1, x3.size()*x0.size(),
             1.0, i0data_sorted, x3.size()*x0.size(), i1data_sorted, x3.size()*x0.size(),
             1.0, odata_sorted, x2.size()*x1.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size());
  out()->add_block(odata, x2, x1);
}

void Task120::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  // tensor label: I43
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(x3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x0), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c1 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
        sort_indices<2,1,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        // tensor label: I44
        std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, c1, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, c1, x3)]);
        sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), c1.size(), x3.size());
        dgemm_("T", "N", x0.size(), x3.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x3.size());
  out()->add_block(odata, x3, x0);
}

void Task121::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I44
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a2, c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, c1, x3), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x3);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c3.size(), a2.size(), c1.size(), x3.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, x3);
    sort_indices<2,1,0,3,1,1,2,1>(i1data, odata, c1.size(), a2.size(), c3.size(), x3.size());
  }
  out()->add_block(odata, c3, a2, c1, x3);
}

void Task122::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  // tensor label: I42
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma81
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: I252
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size());
      dgemm_("T", "N", x2.size()*x1.size(), 1, x3.size()*x0.size(),
             1.0, i0data_sorted, x3.size()*x0.size(), i1data_sorted, x3.size()*x0.size(),
             1.0, odata_sorted, x2.size()*x1.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size());
  out()->add_block(odata, x2, x1);
}

void Task123::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  // tensor label: I252
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(x3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x0), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a1 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
        sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
        // tensor label: I253
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a3, c2, a1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a3, c2, a1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, x3.size(), a3.size(), c2.size(), a1.size());
        dgemm_("T", "N", x0.size(), x3.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x3.size());
  out()->add_block(odata, x3, x0);
}

void Task124::Task_local::compute() {
  const Index x3 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I253
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a3, c2, a1)]);
  std::fill_n(odata.get(), out()->get_size(x3, a3, c2, a1), 0.0);
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
  out()->add_block(odata, x3, a3, c2, a1);
}

void Task125::Task_local::compute() {
  const Index c4 = b(0);
  const Index c1 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(c4, c1)]);
  std::fill_n(odata.get(), out()->get_size(c4, c1), 0.0);
  {
    // tensor label: I48
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c4);
    sort_indices<1,0,1,1,1,1>(i0data, odata, c1.size(), c4.size());
  }
  out()->add_block(odata, c4, c1);
}

void Task126::Task_local::compute() {
  const Index c1 = b(0);
  const Index c4 = b(1);
  // tensor label: I48
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c4)]);
  std::fill_n(odata.get(), out()->get_size(c1, c4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c4), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c3, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a2, c3, x1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c3.size(), x1.size());
        // tensor label: I49
        std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, c1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, c1, x1)]);
        sort_indices<1,0,3,2,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), c1.size(), x1.size());
        dgemm_("T", "N", c4.size(), c1.size(), c3.size()*a2.size()*x1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*x1.size(), i1data_sorted, c3.size()*a2.size()*x1.size(),
               1.0, odata_sorted, c4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c4.size(), c1.size());
  out()->add_block(odata, c1, c4);
}

void Task127::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x1 = b(3);
  // tensor label: I49
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a2, c1, x1)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, c1, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a2, c1, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a2, c1, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x0)]);
    sort_indices<3,0,1,2,0,1,-1,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
    dgemm_("T", "N", x1.size(), c3.size()*a2.size()*c1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), c1.size(), a2.size(), c3.size());
  out()->add_block(odata, c3, a2, c1, x1);
}

void Task128::Task_local::compute() {
  const Index c1 = b(0);
  const Index c4 = b(1);
  // tensor label: I48
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c4)]);
  std::fill_n(odata.get(), out()->get_size(c1, c4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c4), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c4, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c4, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c4.size(), x1.size());
        // tensor label: I52
        std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, c1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, c1, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), c1.size(), x1.size());
        dgemm_("T", "N", c4.size(), c1.size(), c3.size()*a2.size()*x1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*x1.size(), i1data_sorted, c3.size()*a2.size()*x1.size(),
               1.0, odata_sorted, c4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c4.size(), c1.size());
  out()->add_block(odata, c1, c4);
}

void Task129::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x1 = b(3);
  // tensor label: I52
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a2, c1, x1)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, c1, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a2, c1, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a2, c1, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma17
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x0)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
    dgemm_("T", "N", x1.size(), c3.size()*a2.size()*c1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), c1.size(), a2.size(), c3.size());
  out()->add_block(odata, c3, a2, c1, x1);
}

void Task130::Task_local::compute() {
  const Index c4 = b(0);
  const Index c3 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(c4, c3)]);
  std::fill_n(odata.get(), out()->get_size(c4, c3), 0.0);
  {
    // tensor label: I54
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c4);
    sort_indices<1,0,1,1,1,1>(i0data, odata, c3.size(), c4.size());
  }
  out()->add_block(odata, c4, c3);
}

void Task131::Task_local::compute() {
  const Index c3 = b(0);
  const Index c4 = b(1);
  // tensor label: I54
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, c4)]);
  std::fill_n(odata.get(), out()->get_size(c3, c4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c4), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c1, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a2, c1, x1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c1.size(), x1.size());
        // tensor label: I55
        std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, c1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, c1, x1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), c1.size(), x1.size());
        dgemm_("T", "N", c4.size(), c3.size(), a2.size()*c1.size()*x1.size(),
               1.0, i0data_sorted, a2.size()*c1.size()*x1.size(), i1data_sorted, a2.size()*c1.size()*x1.size(),
               1.0, odata_sorted, c4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c4.size(), c3.size());
  out()->add_block(odata, c3, c4);
}

void Task132::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x1 = b(3);
  // tensor label: I55
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a2, c1, x1)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, c1, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a2, c1, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a2, c1, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma17
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x0)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
    dgemm_("T", "N", x1.size(), c3.size()*a2.size()*c1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), c1.size(), a2.size(), c3.size());
  out()->add_block(odata, c3, a2, c1, x1);
}

void Task133::Task_local::compute() {
  const Index c3 = b(0);
  const Index c4 = b(1);
  // tensor label: I54
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, c4)]);
  std::fill_n(odata.get(), out()->get_size(c3, c4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c4), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c4, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c4, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c4.size(), x1.size());
        // tensor label: I61
        std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, c1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, c1, x1)]);
        sort_indices<2,1,3,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), c1.size(), x1.size());
        dgemm_("T", "N", c4.size(), c3.size(), a2.size()*c1.size()*x1.size(),
               1.0, i0data_sorted, a2.size()*c1.size()*x1.size(), i1data_sorted, a2.size()*c1.size()*x1.size(),
               1.0, odata_sorted, c4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c4.size(), c3.size());
  out()->add_block(odata, c3, c4);
}

void Task134::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x1 = b(3);
  // tensor label: I61
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a2, c1, x1)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, c1, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a2, c1, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a2, c1, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma17
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x0)]);
    sort_indices<3,0,1,2,0,1,-2,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
    dgemm_("T", "N", x1.size(), c3.size()*a2.size()*c1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), c1.size(), a2.size(), c3.size());
  out()->add_block(odata, c3, a2, c1, x1);
}

void Task135::Task_local::compute() {
  const Index a2 = b(0);
  const Index a4 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, a4)]);
  std::fill_n(odata.get(), out()->get_size(a2, a4), 0.0);
  {
    // tensor label: I57
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, a4);
    sort_indices<0,1,1,1,1,1>(i0data, odata, a2.size(), a4.size());
  }
  out()->add_block(odata, a2, a4);
}

void Task136::Task_local::compute() {
  const Index a2 = b(0);
  const Index a4 = b(1);
  // tensor label: I57
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, a4)]);
  std::fill_n(odata.get(), out()->get_size(a2, a4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, a4), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& c1 : *range_[0]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, c1, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4, c1, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), c1.size(), x1.size());
        // tensor label: I58
        std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, c1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, c1, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), c1.size(), x1.size());
        dgemm_("T", "N", a4.size(), a2.size(), c3.size()*c1.size()*x1.size(),
               1.0, i0data_sorted, c3.size()*c1.size()*x1.size(), i1data_sorted, c3.size()*c1.size()*x1.size(),
               1.0, odata_sorted, a4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a4.size(), a2.size());
  out()->add_block(odata, a2, a4);
}

void Task137::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x1 = b(3);
  // tensor label: I58
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a2, c1, x1)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, c1, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a2, c1, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a2, c1, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma17
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x0)]);
    sort_indices<3,0,1,2,0,1,-1,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
    dgemm_("T", "N", x1.size(), c3.size()*a2.size()*c1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), c1.size(), a2.size(), c3.size());
  out()->add_block(odata, c3, a2, c1, x1);
}

void Task138::Task_local::compute() {
  const Index a2 = b(0);
  const Index a4 = b(1);
  // tensor label: I57
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, a4)]);
  std::fill_n(odata.get(), out()->get_size(a2, a4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, a4), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& c3 : *range_[0]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x1.size());
        // tensor label: I64
        std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, c1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, c1, x1)]);
        sort_indices<2,0,3,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), c1.size(), x1.size());
        dgemm_("T", "N", a4.size(), a2.size(), c3.size()*c1.size()*x1.size(),
               1.0, i0data_sorted, c3.size()*c1.size()*x1.size(), i1data_sorted, c3.size()*c1.size()*x1.size(),
               1.0, odata_sorted, a4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a4.size(), a2.size());
  out()->add_block(odata, a2, a4);
}

void Task139::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x1 = b(3);
  // tensor label: I64
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a2, c1, x1)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, c1, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a2, c1, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a2, c1, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x0)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
    dgemm_("T", "N", x1.size(), c3.size()*a2.size()*c1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), c1.size(), a2.size(), c3.size());
  out()->add_block(odata, c3, a2, c1, x1);
}

void Task140::Task_local::compute() {
  const Index x1 = b(0);
  const Index c1 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, c1)]);
  std::fill_n(odata.get(), out()->get_size(x1, c1), 0.0);
  {
    // tensor label: I66
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c1);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x1.size(), c1.size());
  }
  out()->add_block(odata, x1, c1);
}

void Task141::Task_local::compute() {
  const Index x1 = b(0);
  const Index c1 = b(1);
  // tensor label: I66
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, c1)]);
  std::fill_n(odata.get(), out()->get_size(x1, c1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, c1), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& c3 : *range_[0]) {
      for (auto& a2 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
        sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        // tensor label: I67
        std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c3, x0, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c3, x0, x1)]);
        sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), c3.size(), x0.size(), x1.size());
        dgemm_("T", "N", c1.size(), x1.size(), a2.size()*c3.size()*x0.size(),
               1.0, i0data_sorted, a2.size()*c3.size()*x0.size(), i1data_sorted, a2.size()*c3.size()*x0.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c1.size(), x1.size());
  out()->add_block(odata, x1, c1);
}

void Task142::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I67
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c3, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c3, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c3, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x0, x1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: Gamma22
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x2, x3, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x2, x3, x1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), x2.size(), x3.size(), x1.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a2, c3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a2, c3, x2)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a2.size(), c3.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), a2.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a2.size(), c3.size());
  out()->add_block(odata, a2, c3, x0, x1);
}

void Task143::Task_local::compute() {
  const Index x1 = b(0);
  const Index c1 = b(1);
  // tensor label: I66
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, c1)]);
  std::fill_n(odata.get(), out()->get_size(x1, c1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, c1), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& c3 : *range_[0]) {
      for (auto& a2 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
        sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        // tensor label: I73
        std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, x0, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, x0, x1)]);
        sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), x0.size(), x1.size());
        dgemm_("T", "N", c1.size(), x1.size(), c3.size()*a2.size()*x0.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*x0.size(), i1data_sorted, c3.size()*a2.size()*x0.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c1.size(), x1.size());
  out()->add_block(odata, x1, c1);
}

void Task144::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I73
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a2, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a2, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma13
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,-1,1>(i1data, i1data_sorted, c3.size(), a2.size(), x3.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), c3.size()*a2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c3.size(), a2.size());
  out()->add_block(odata, c3, a2, x0, x1);
}

void Task145::Task_local::compute() {
  const Index x1 = b(0);
  const Index c3 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, c3)]);
  std::fill_n(odata.get(), out()->get_size(x1, c3), 0.0);
  {
    // tensor label: I69
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c3);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x1.size(), c3.size());
  }
  out()->add_block(odata, x1, c3);
}

void Task146::Task_local::compute() {
  const Index x1 = b(0);
  const Index c3 = b(1);
  // tensor label: I69
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, c3)]);
  std::fill_n(odata.get(), out()->get_size(x1, c3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, c3), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c1 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
        sort_indices<3,1,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        // tensor label: I70
        std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c1, x0, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c1, x0, x1)]);
        sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), c1.size(), x0.size(), x1.size());
        dgemm_("T", "N", c3.size(), x1.size(), a2.size()*c1.size()*x0.size(),
               1.0, i0data_sorted, a2.size()*c1.size()*x0.size(), i1data_sorted, a2.size()*c1.size()*x0.size(),
               1.0, odata_sorted, c3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c3.size(), x1.size());
  out()->add_block(odata, x1, c3);
}

void Task147::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I70
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c1, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma13
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: I71
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a2, c1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a2, c1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a2.size(), c1.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), a2.size()*c1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a2.size(), c1.size());
  out()->add_block(odata, a2, c1, x0, x1);
}

void Task148::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x2 = b(3);
  // tensor label: I71
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a2, c1, x2)]);
  std::fill_n(odata.get(), out()->get_size(x3, a2, c1, x2), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), a2.size(), c1.size(), x2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, x3, x2);
    sort_indices<2,1,0,3,1,1,2,1>(i1data, odata, c1.size(), a2.size(), x3.size(), x2.size());
  }
  out()->add_block(odata, x3, a2, c1, x2);
}

void Task149::Task_local::compute() {
  const Index x1 = b(0);
  const Index a4 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, a4)]);
  std::fill_n(odata.get(), out()->get_size(x1, a4), 0.0);
  {
    // tensor label: I78
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x1);
    sort_indices<1,0,1,1,1,1>(i0data, odata, a4.size(), x1.size());
  }
  out()->add_block(odata, x1, a4);
}

#endif
