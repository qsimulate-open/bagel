//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_tasks7.cc
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

#include <src/smith/relmrci/RelMRCI_tasks7.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelMRCI;

void Task300::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I63
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: Gamma5
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x2, x3, x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x2, x3, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
      // tensor label: I67
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x2, a2, c1, x3);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x2, a2, c1, x3)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), a2.size(), c1.size(), x3.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), a2.size()*c1.size(), x2.size()*x3.size(),
             1.0, i0data_sorted, x2.size()*x3.size(), i1data_sorted, x2.size()*x3.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a2.size(), c1.size());
  out()->add_block(odata, c1, x1, x0, a2);
}

void Task301::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I67
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a2, c1, x3);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c3, a2, c1, x3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c1.size(), x3.size());
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x2, c3);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x2, c3)]);
    sort_indices<1,0,0,1,-1,1>(i1data, i1data_sorted, x2.size(), c3.size());
    zgemm3m_("T", "N", a2.size()*c1.size()*x3.size(), x2.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, a2.size()*c1.size()*x3.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x3.size(), x2.size());
  out()->add_block(odata, x2, a2, c1, x3);
}

void Task302::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I67
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a2, c3, x3);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a2, c3, x3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x3.size());
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x2, c3);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x2, c3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x2.size(), c3.size());
    zgemm3m_("T", "N", c1.size()*a2.size()*x3.size(), x2.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, c1.size()*a2.size()*x3.size());
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x3.size(), x2.size());
  out()->add_block(odata, x2, a2, c1, x3);
}

void Task303::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I67
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c4, a2, c3, x3);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c4, a2, c3, x3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c3.size(), x3.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x2, c4, c1, c3);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x2, c4, c1, c3)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), c4.size(), c1.size(), c3.size());
      zgemm3m_("T", "N", a2.size()*x3.size(), x2.size()*c1.size(), c4.size()*c3.size(),
             1.0, i0data_sorted, c4.size()*c3.size(), i1data_sorted, c4.size()*c3.size(),
             1.0, odata_sorted, a2.size()*x3.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, a2.size(), x3.size(), x2.size(), c1.size());
  out()->add_block(odata, x2, a2, c1, x3);
}

void Task304::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I67
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a2, c4, x3);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c3, a2, c4, x3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c4.size(), x3.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x2, c4, c1, c3);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x2, c4, c1, c3)]);
      sort_indices<3,1,0,2,0,1,-1,1>(i1data, i1data_sorted, x2.size(), c4.size(), c1.size(), c3.size());
      zgemm3m_("T", "N", a2.size()*x3.size(), x2.size()*c1.size(), c4.size()*c3.size(),
             1.0, i0data_sorted, c4.size()*c3.size(), i1data_sorted, c4.size()*c3.size(),
             1.0, odata_sorted, a2.size()*x3.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, a2.size(), x3.size(), x2.size(), c1.size());
  out()->add_block(odata, x2, a2, c1, x3);
}

void Task305::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I67
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c4, a3, c1, x3);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c4, a3, c1, x3)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a3.size(), c1.size(), x3.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x2, c4, a3, a2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x2, c4, a3, a2)]);
      sort_indices<1,2,0,3,0,1,-1,1>(i1data, i1data_sorted, x2.size(), c4.size(), a3.size(), a2.size());
      zgemm3m_("T", "N", c1.size()*x3.size(), x2.size()*a2.size(), c4.size()*a3.size(),
             1.0, i0data_sorted, c4.size()*a3.size(), i1data_sorted, c4.size()*a3.size(),
             1.0, odata_sorted, c1.size()*x3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), x3.size(), x2.size(), a2.size());
  out()->add_block(odata, x2, a2, c1, x3);
}

void Task306::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I67
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a4, c1, x3);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c3, a4, c1, x3)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), c1.size(), x3.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x2, a2, a4, c3);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x2, a2, a4, c3)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, x2.size(), a2.size(), a4.size(), c3.size());
      zgemm3m_("T", "N", c1.size()*x3.size(), x2.size()*a2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*x3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), x3.size(), x2.size(), a2.size());
  out()->add_block(odata, x2, a2, c1, x3);
}

void Task307::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I67
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a3, c4, x3);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a3, c4, x3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), c4.size(), x3.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x2, c4, a3, a2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x2, c4, a3, a2)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), c4.size(), a3.size(), a2.size());
      zgemm3m_("T", "N", c1.size()*x3.size(), x2.size()*a2.size(), c4.size()*a3.size(),
             1.0, i0data_sorted, c4.size()*a3.size(), i1data_sorted, c4.size()*a3.size(),
             1.0, odata_sorted, c1.size()*x3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), x3.size(), x2.size(), a2.size());
  out()->add_block(odata, x2, a2, c1, x3);
}

void Task308::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I67
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a4, c3, x3);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a4, c3, x3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x3.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x2, a2, a4, c3);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x2, a2, a4, c3)]);
      sort_indices<2,3,0,1,0,1,-1,1>(i1data, i1data_sorted, x2.size(), a2.size(), a4.size(), c3.size());
      zgemm3m_("T", "N", c1.size()*x3.size(), x2.size()*a2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*x3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), x3.size(), x2.size(), a2.size());
  out()->add_block(odata, x2, a2, c1, x3);
}

void Task309::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I67
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a4, c3, a2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a4, c3, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, x3, x2, c3);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, x3, x2, c3)]);
      sort_indices<0,3,1,2,0,1,-1,1>(i1data, i1data_sorted, a4.size(), x3.size(), x2.size(), c3.size());
      zgemm3m_("T", "N", c1.size()*a2.size(), x3.size()*x2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x3.size(), x2.size());
  out()->add_block(odata, x2, a2, c1, x3);
}

void Task310::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I67
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a2, c3, a4);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a2, c3, a4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, x3, x2, c3);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, x3, x2, c3)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a4.size(), x3.size(), x2.size(), c3.size());
      zgemm3m_("T", "N", c1.size()*a2.size(), x3.size()*x2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x3.size(), x2.size());
  out()->add_block(odata, x2, a2, c1, x3);
}

void Task311::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I63
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma24
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I73
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a2, x3, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size(), x3.size(), x2.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), c1.size()*a2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), a2.size());
  out()->add_block(odata, c1, x1, x0, a2);
}

void Task312::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I73
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, x3, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, x3, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a2, x3, x2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c3, a2, x3, x2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), x3.size(), x2.size());
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, c3);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, c3)]);
    sort_indices<1,0,0,1,-1,1>(i1data, i1data_sorted, c1.size(), c3.size());
    zgemm3m_("T", "N", a2.size()*x3.size()*x2.size(), c1.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, a2.size()*x3.size()*x2.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a2.size(), x3.size(), x2.size(), c1.size());
  out()->add_block(odata, c1, a2, x3, x2);
}

void Task313::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I73
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, x3, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, x3, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a3, x3, x2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a3, x3, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), x3.size(), x2.size());
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a3, a2);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a3, a2)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), a2.size());
    zgemm3m_("T", "N", c1.size()*x3.size()*x2.size(), a2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, c1.size()*x3.size()*x2.size());
  }
  sort_indices<0,3,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), x3.size(), x2.size(), a2.size());
  out()->add_block(odata, c1, a2, x3, x2);
}

void Task314::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I73
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, x3, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, x3, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a3, x2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(a3, x2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), x2.size());
    // tensor label: I89
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, a3, c1, a2);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, a3, c1, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a3.size(), c1.size(), a2.size());
    zgemm3m_("T", "N", x2.size(), x3.size()*c1.size()*a2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x2.size());
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, x2.size(), x3.size(), c1.size(), a2.size());
  out()->add_block(odata, c1, a2, x3, x2);
}

void Task315::Task_local::compute() {
  const Index x3 = b(0);
  const Index a3 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I89
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x3, a3, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x3, a3, c1, a2), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, a3, c1, a2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), a3.size(), c1.size(), a2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(x3, a2, c1, a3);
    sort_indices<0,3,2,1,1,1,-1,1>(i1data, odata, x3.size(), a2.size(), c1.size(), a3.size());
  }
  out()->add_block(odata, x3, a3, c1, a2);
}

void Task316::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I73
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, x3, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, x3, x2), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c4, a3, x3, x2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c4, a3, x3, x2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a3.size(), x3.size(), x2.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, c4, a3, a2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, c4, a3, a2)]);
      sort_indices<1,2,0,3,0,1,-1,1>(i1data, i1data_sorted, c1.size(), c4.size(), a3.size(), a2.size());
      zgemm3m_("T", "N", x3.size()*x2.size(), c1.size()*a2.size(), c4.size()*a3.size(),
             1.0, i0data_sorted, c4.size()*a3.size(), i1data_sorted, c4.size()*a3.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c1.size(), a2.size());
  out()->add_block(odata, c1, a2, x3, x2);
}

void Task317::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I73
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, x3, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, x3, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a4, x3, x2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c3, a4, x3, x2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), x3.size(), x2.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a2, a4, c3);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a2, a4, c3)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size(), a4.size(), c3.size());
      zgemm3m_("T", "N", x3.size()*x2.size(), c1.size()*a2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c1.size(), a2.size());
  out()->add_block(odata, c1, a2, x3, x2);
}

void Task318::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I73
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, x3, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, x3, x2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a4, c3, a2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a4, c3, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
      // tensor label: I631
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, c3, x3, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, c3, x3, x2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size(), x3.size(), x2.size());
      zgemm3m_("T", "N", c1.size()*a2.size(), x3.size()*x2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x3.size(), x2.size());
  out()->add_block(odata, c1, a2, x3, x2);
}

void Task319::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I631
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a4, c3, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(a4, c3, x3, x2), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a4, c3, x3, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, a4.size(), c3.size(), x3.size(), x2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(x3, x2, a4, c3);
    sort_indices<2,3,0,1,1,1,-1,1>(i1data, odata, x3.size(), x2.size(), a4.size(), c3.size());
  }
  out()->add_block(odata, a4, c3, x3, x2);
}

void Task320::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I73
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, x3, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, x3, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a2, c3, a4);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a2, c3, a4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
      // tensor label: I634
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, c3, x3, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, c3, x3, x2)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size(), x3.size(), x2.size());
      zgemm3m_("T", "N", c1.size()*a2.size(), x3.size()*x2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x3.size(), x2.size());
  out()->add_block(odata, c1, a2, x3, x2);
}

void Task321::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I634
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a4, c3, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(a4, c3, x3, x2), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a4, c3, x3, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, a4.size(), c3.size(), x3.size(), x2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(x3, x2, a4, c3);
    sort_indices<2,3,0,1,1,1,1,1>(i1data, odata, x3.size(), x2.size(), a4.size(), c3.size());
  }
  out()->add_block(odata, a4, c3, x3, x2);
}

void Task322::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I73
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, x3, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, x3, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a3, c4, a2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a3, c4, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), c4.size(), a2.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, c4, a3, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, c4, a3, x2)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), c4.size(), a3.size(), x2.size());
      zgemm3m_("T", "N", c1.size()*a2.size(), x3.size()*x2.size(), c4.size()*a3.size(),
             1.0, i0data_sorted, c4.size()*a3.size(), i1data_sorted, c4.size()*a3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x3.size(), x2.size());
  out()->add_block(odata, c1, a2, x3, x2);
}

void Task323::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I73
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, x3, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, x3, x2), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a2, c4, a3);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a2, c4, a3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c4.size(), a3.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, c4, a3, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, c4, a3, x2)]);
      sort_indices<1,2,0,3,0,1,-1,1>(i1data, i1data_sorted, x3.size(), c4.size(), a3.size(), x2.size());
      zgemm3m_("T", "N", c1.size()*a2.size(), x3.size()*x2.size(), c4.size()*a3.size(),
             1.0, i0data_sorted, c4.size()*a3.size(), i1data_sorted, c4.size()*a3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x3.size(), x2.size());
  out()->add_block(odata, c1, a2, x3, x2);
}

void Task324::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I73
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, x3, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, x3, x2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, a4, c3, a2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, a4, c3, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a4.size(), c3.size(), a2.size());
      // tensor label: I679
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, x2, c1, c3);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, x2, c1, c3)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, a4.size(), x2.size(), c1.size(), c3.size());
      zgemm3m_("T", "N", x3.size()*a2.size(), x2.size()*c1.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, x3.size()*a2.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, x3.size(), a2.size(), x2.size(), c1.size());
  out()->add_block(odata, c1, a2, x3, x2);
}

void Task325::Task_local::compute() {
  const Index a4 = b(0);
  const Index x2 = b(1);
  const Index c1 = b(2);
  const Index c3 = b(3);
  // tensor label: I679
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a4, x2, c1, c3)]);
  std::fill_n(odata.get(), out()->get_size(a4, x2, c1, c3), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a4, x2, c1, c3);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, a4.size(), x2.size(), c1.size(), c3.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(c1, x2, a4, c3);
    sort_indices<2,1,0,3,1,1,1,1>(i1data, odata, c1.size(), x2.size(), a4.size(), c3.size());
  }
  out()->add_block(odata, a4, x2, c1, c3);
}

void Task326::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I73
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, x3, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, x3, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, a2, c3, a4);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, a2, c3, a4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c3.size(), a4.size());
      // tensor label: I682
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, x2, c1, c3);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, x2, c1, c3)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a4.size(), x2.size(), c1.size(), c3.size());
      zgemm3m_("T", "N", x3.size()*a2.size(), x2.size()*c1.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, x3.size()*a2.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, x3.size(), a2.size(), x2.size(), c1.size());
  out()->add_block(odata, c1, a2, x3, x2);
}

void Task327::Task_local::compute() {
  const Index a4 = b(0);
  const Index x2 = b(1);
  const Index c1 = b(2);
  const Index c3 = b(3);
  // tensor label: I682
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a4, x2, c1, c3)]);
  std::fill_n(odata.get(), out()->get_size(a4, x2, c1, c3), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a4, x2, c1, c3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, a4.size(), x2.size(), c1.size(), c3.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(c1, x2, a4, c3);
    sort_indices<2,1,0,3,1,1,-1,1>(i1data, odata, c1.size(), x2.size(), a4.size(), c3.size());
  }
  out()->add_block(odata, a4, x2, c1, c3);
}

void Task328::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I73
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, x3, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, x3, x2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, a4, c1, a3);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, a4, c1, a3)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a4.size(), c1.size(), a3.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, x2, a3, a2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, x2, a3, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), x2.size(), a3.size(), a2.size());
      zgemm3m_("T", "N", x3.size()*c1.size(), x2.size()*a2.size(), a4.size()*a3.size(),
             1.0, i0data_sorted, a4.size()*a3.size(), i1data_sorted, a4.size()*a3.size(),
             1.0, odata_sorted, x3.size()*c1.size());
    }
  }
  sort_indices<1,3,0,2,1,1,1,1>(odata_sorted, odata, x3.size(), c1.size(), x2.size(), a2.size());
  out()->add_block(odata, c1, a2, x3, x2);
}

void Task329::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I73
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, x3, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, x3, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, a3, c1, a4);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, a3, c1, a4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c1.size(), a4.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, x2, a3, a2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, x2, a3, a2)]);
      sort_indices<2,0,1,3,0,1,-1,1>(i1data, i1data_sorted, a4.size(), x2.size(), a3.size(), a2.size());
      zgemm3m_("T", "N", x3.size()*c1.size(), x2.size()*a2.size(), a4.size()*a3.size(),
             1.0, i0data_sorted, a4.size()*a3.size(), i1data_sorted, a4.size()*a3.size(),
             1.0, odata_sorted, x3.size()*c1.size());
    }
  }
  sort_indices<1,3,0,2,1,1,1,1>(odata_sorted, odata, x3.size(), c1.size(), x2.size(), a2.size());
  out()->add_block(odata, c1, a2, x3, x2);
}

void Task330::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I63
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, x2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, x2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), x2.size());
    // tensor label: I79
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a2, x2, x1, x0);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a2, x2, x1, x0)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), x2.size(), x1.size(), x0.size());
    zgemm3m_("T", "N", c1.size(), a2.size()*x1.size()*x0.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, c1.size());
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), x0.size());
  out()->add_block(odata, c1, x1, x0, a2);
}

void Task331::Task_local::compute() {
  const Index a2 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I79
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a2, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a2, x2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a2, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma26
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x2, x4, x3, x1, x0);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x2, x4, x3, x1, x0)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x5, a2, x4, x3);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x5, a2, x4, x3)]);
        sort_indices<0,2,3,1,0,1,-1,1>(i1data, i1data_sorted, x5.size(), a2.size(), x4.size(), x3.size());
        zgemm3m_("T", "N", x2.size()*x1.size()*x0.size(), a2.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), a2.size());
  out()->add_block(odata, a2, x2, x1, x0);
}

void Task332::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I63
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  // tensor label: Gamma27
  std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
  std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
  // tensor label: I82
  std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a2);
  std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a2)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size());
  zgemm3m_("T", "N", x1.size()*x0.size(), c1.size()*a2.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, x1.size()*x0.size());
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), a2.size());
  out()->add_block(odata, c1, x1, x0, a2);
}

void Task333::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  // tensor label: I82
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a4, c3, a2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a4, c3, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
      // tensor label: h1
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, c3);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, c3)]);
      sort_indices<0,1,0,1,-2,1>(i1data, i1data_sorted, a4.size(), c3.size());
      zgemm3m_("T", "N", c1.size()*a2.size(), 1, a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size());
  out()->add_block(odata, c1, a2);
}

void Task334::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  // tensor label: I82
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a2, c3, a4);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a2, c3, a4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
      // tensor label: h1
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, c3);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, c3)]);
      sort_indices<1,0,0,1,2,1>(i1data, i1data_sorted, a4.size(), c3.size());
      zgemm3m_("T", "N", c1.size()*a2.size(), 1, a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size());
  out()->add_block(odata, c1, a2);
}

void Task335::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  // tensor label: I82
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      for (auto& c5 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a4, c5, a2);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c3, a4, c5, a2)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), c5.size(), a2.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, c5, a4, c3);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, c5, a4, c3)]);
        sort_indices<3,2,1,0,0,1,-2,1>(i1data, i1data_sorted, c1.size(), c5.size(), a4.size(), c3.size());
        zgemm3m_("T", "N", a2.size(), c1.size(), c5.size()*a4.size()*c3.size(),
               1.0, i0data_sorted, c5.size()*a4.size()*c3.size(), i1data_sorted, c5.size()*a4.size()*c3.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size());
  out()->add_block(odata, c1, a2);
}

void Task336::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  // tensor label: I82
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& c5 : *range_[0]) {
      for (auto& a4 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a2, c5, a4);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c3, a2, c5, a4)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c5.size(), a4.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, c5, a4, c3);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, c5, a4, c3)]);
        sort_indices<3,1,2,0,0,1,2,1>(i1data, i1data_sorted, c1.size(), c5.size(), a4.size(), c3.size());
        zgemm3m_("T", "N", a2.size(), c1.size(), c5.size()*a4.size()*c3.size(),
               1.0, i0data_sorted, c5.size()*a4.size()*c3.size(), i1data_sorted, c5.size()*a4.size()*c3.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size());
  out()->add_block(odata, c1, a2);
}

void Task337::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  // tensor label: I82
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2), 0.0);
  for (auto& a5 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      for (auto& a4 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a5, c3, a4);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a5, c3, a4)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), a5.size(), c3.size(), a4.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a5, a2, a4, c3);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a5, a2, a4, c3)]);
        sort_indices<0,3,2,1,0,1,2,1>(i1data, i1data_sorted, a5.size(), a2.size(), a4.size(), c3.size());
        zgemm3m_("T", "N", c1.size(), a2.size(), a5.size()*a4.size()*c3.size(),
               1.0, i0data_sorted, a5.size()*a4.size()*c3.size(), i1data_sorted, a5.size()*a4.size()*c3.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size());
  out()->add_block(odata, c1, a2);
}

void Task338::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  // tensor label: I82
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      for (auto& a5 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a4, c3, a5);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a4, c3, a5)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a5.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a5, a2, a4, c3);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a5, a2, a4, c3)]);
        sort_indices<2,3,0,1,0,1,-2,1>(i1data, i1data_sorted, a5.size(), a2.size(), a4.size(), c3.size());
        zgemm3m_("T", "N", c1.size(), a2.size(), a5.size()*a4.size()*c3.size(),
               1.0, i0data_sorted, a5.size()*a4.size()*c3.size(), i1data_sorted, a5.size()*a4.size()*c3.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size());
  out()->add_block(odata, c1, a2);
}

void Task339::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I63
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, a2, x2, c3);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, a2, x2, c3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), x2.size(), c3.size());
        // tensor label: I543
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, c3, x3, x2, x1, x0);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, c3, x3, x2, x1, x0)]);
        sort_indices<2,3,1,0,4,5,0,1,1,1>(i1data, i1data_sorted, c1.size(), c3.size(), x3.size(), x2.size(), x1.size(), x0.size());
        zgemm3m_("T", "N", a2.size(), c1.size()*x1.size()*x0.size(), c3.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, c3.size()*x3.size()*x2.size(), i1data_sorted, c3.size()*x3.size()*x2.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x1.size(), x0.size());
  out()->add_block(odata, c1, x1, x0, a2);
}

void Task340::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I543
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, c3, x3, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x3, x2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c3, x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x3, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: Gamma84
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x5, x2, x4, x1, x0);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x5, x2, x4, x1, x0)]);
      sort_indices<1,3,0,2,4,5,0,1,1,1>(i0data, i0data_sorted, x3.size(), x5.size(), x2.size(), x4.size(), x1.size(), x0.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, x5, c3, x4);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, x5, c3, x4)]);
      sort_indices<1,3,0,2,0,1,-2,1>(i1data, i1data_sorted, c1.size(), x5.size(), c3.size(), x4.size());
      zgemm3m_("T", "N", x3.size()*x2.size()*x1.size()*x0.size(), c1.size()*c3.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, x3.size()*x2.size()*x1.size()*x0.size());
    }
  }
  sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size(), c1.size(), c3.size());
  out()->add_block(odata, c1, c3, x3, x2, x1, x0);
}

void Task341::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I63
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x4, a2, x3, x2);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x4, a2, x3, x2)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x4.size(), a2.size(), x3.size(), x2.size());
        // tensor label: I546
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, x4, x3, x2, x1, x0);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, x4, x3, x2, x1, x0)]);
        sort_indices<1,2,3,0,4,5,0,1,1,1>(i1data, i1data_sorted, c1.size(), x4.size(), x3.size(), x2.size(), x1.size(), x0.size());
        zgemm3m_("T", "N", a2.size(), c1.size()*x1.size()*x0.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x1.size(), x0.size());
  out()->add_block(odata, c1, x1, x0, a2);
}

void Task342::Task_local::compute() {
  const Index c1 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I546
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x4, x3, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x4, x3, x2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x4, x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x4, x3, x2, x1, x0), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: Gamma179
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x6, x4, x5, x3, x2, x1, x0);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x7, x6, x4, x5, x3, x2, x1, x0)]);
        sort_indices<0,1,3,2,4,5,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x4.size(), x5.size(), x3.size(), x2.size(), x1.size(), x0.size());
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x7, x6, c1, x5);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x7, x6, c1, x5)]);
        sort_indices<0,1,3,2,0,1,1,2>(i1data, i1data_sorted, x7.size(), x6.size(), c1.size(), x5.size());
        zgemm3m_("T", "N", x4.size()*x3.size()*x2.size()*x1.size()*x0.size(), c1.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, x4.size()*x3.size()*x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x4.size(), x3.size(), x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x4, x3, x2, x1, x0);
}

void Task343::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I63
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x4, x3, x2, a2);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x4, x3, x2, a2)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x4.size(), x3.size(), x2.size(), a2.size());
        // tensor label: I549
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, x2, x4, x3, x1, x0);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, x2, x4, x3, x1, x0)]);
        sort_indices<2,3,1,0,4,5,0,1,1,1>(i1data, i1data_sorted, c1.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
        zgemm3m_("T", "N", a2.size(), c1.size()*x1.size()*x0.size(), x2.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x2.size()*x4.size()*x3.size(), i1data_sorted, x2.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x1.size(), x0.size());
  out()->add_block(odata, c1, x1, x0, a2);
}

void Task344::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I549
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x2, x4, x3, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x4, x3, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x2, x4, x3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x4, x3, x1, x0), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: Gamma77
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x6, x2, x5, x4, x3, x1, x0);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x7, x6, x2, x5, x4, x3, x1, x0)]);
        sort_indices<0,1,3,2,4,5,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x2.size(), x5.size(), x4.size(), x3.size(), x1.size(), x0.size());
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x7, x6, c1, x5);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x7, x6, c1, x5)]);
        sort_indices<0,1,3,2,0,1,1,2>(i1data, i1data_sorted, x7.size(), x6.size(), c1.size(), x5.size());
        zgemm3m_("T", "N", x2.size()*x4.size()*x3.size()*x1.size()*x0.size(), c1.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, x2.size()*x4.size()*x3.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x2.size(), x4.size(), x3.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x4, x3, x1, x0);
}

void Task345::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I63
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x2, c3, c1, a2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x2, c3, c1, a2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x2.size(), c3.size(), c1.size(), a2.size());
      // tensor label: I552
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, x2, x1, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, x2, x1, x0)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), x2.size(), x1.size(), x0.size());
      zgemm3m_("T", "N", c1.size()*a2.size(), x1.size()*x0.size(), c3.size()*x2.size(),
             1.0, i0data_sorted, c3.size()*x2.size(), i1data_sorted, c3.size()*x2.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), x0.size());
  out()->add_block(odata, c1, x1, x0, a2);
}

void Task346::Task_local::compute() {
  const Index c3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I552
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, x2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma4
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x4, x2, x3, x1, x0);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x4, x2, x3, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x5, x4, c3, x3);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x5, x4, c3, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), c3.size(), x3.size());
        zgemm3m_("T", "N", x2.size()*x1.size()*x0.size(), c3.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c3.size());
  out()->add_block(odata, c3, x2, x1, x0);
}

void Task347::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I63
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x2, a2, c1, c3);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x2, a2, c1, c3)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x2.size(), a2.size(), c1.size(), c3.size());
      // tensor label: I555
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, x2, x1, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, x2, x1, x0)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), x2.size(), x1.size(), x0.size());
      zgemm3m_("T", "N", a2.size()*c1.size(), x1.size()*x0.size(), c3.size()*x2.size(),
             1.0, i0data_sorted, c3.size()*x2.size(), i1data_sorted, c3.size()*x2.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x1.size(), x0.size());
  out()->add_block(odata, c1, x1, x0, a2);
}

void Task348::Task_local::compute() {
  const Index c3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I555
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, x2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma4
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x4, x2, x3, x1, x0);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x4, x2, x3, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x5, x4, c3, x3);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x5, x4, c3, x3)]);
        sort_indices<0,1,3,2,0,1,-1,1>(i1data, i1data_sorted, x5.size(), x4.size(), c3.size(), x3.size());
        zgemm3m_("T", "N", x2.size()*x1.size()*x0.size(), c3.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c3.size());
  out()->add_block(odata, c3, x2, x1, x0);
}

void Task349::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I63
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& x5 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a2, c1, x5);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c3, a2, c1, x5)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c1.size(), x5.size());
      // tensor label: I558
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, x5, x1, x0);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, x5, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), x5.size(), x1.size(), x0.size());
      zgemm3m_("T", "N", a2.size()*c1.size(), x1.size()*x0.size(), c3.size()*x5.size(),
             1.0, i0data_sorted, c3.size()*x5.size(), i1data_sorted, c3.size()*x5.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x1.size(), x0.size());
  out()->add_block(odata, c1, x1, x0, a2);
}

#endif
