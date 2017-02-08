//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_contract_tasks.cc
// Copyright (C) 2017 Toru Shiozaki
//
// Author: Jae Woo Park <jwpk1201@northwestern.edu>
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


#include <src/smith/caspt2/MSCASPT2.h>
#include <src/smith/caspt2/MSCASPT2_contract_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MSCASPT2;

void Task901::Task_local::compute() {
  const Index ci0 = b(0);

  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);

  std::unique_ptr<double[]> i0data = in(0)->get_block(ci0);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0)]);
  sort_indices<0,0,1,1,1>(i0data, i0data_sorted, ci0.size());

  std::unique_ptr<double[]> i1data = in(1)->get_block();
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size()]);
  sort_indices<0,1,1,1>(i1data, i1data_sorted);
  dgemm_("T", "N", ci0.size(), 1, 1, 1.0, i0data_sorted, 1,
         i1data_sorted, 1, 1.0, odata_sorted, ci0.size());

  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->add_block(odata, ci0);
}

void Task902::Task_local::compute() {
  const Index ci0 = b(0);

  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);

  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x0, x1)]);
      sort_indices<1,2,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x0.size(), x1.size());

      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size());
      dgemm_("T", "N", ci0.size(), 1, x0.size()*x1.size(),
             1.0, i0data_sorted, x0.size()*x1.size(),
             i1data_sorted, x0.size()*x1.size(), 1.0, odata_sorted, ci0.size());
    }
  }

  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->add_block(odata, ci0);
}

void Task903::Task_local::compute() {
  const Index ci0 = b(0);

  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);

  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        for (auto& x3 : *range_[1]) {
          std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1, x2, x3);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x0, x1, x2, x3)]);
          sort_indices<1,2,3,4,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x0.size(), x1.size(), x2.size(), x3.size());

          std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, x2, x3);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, x2, x3)]);
          sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), x2.size(), x3.size());
          dgemm_("T", "N", ci0.size(), 1, x0.size()*x1.size()*x2.size()*x3.size(),
                 1.0, i0data_sorted, x0.size()*x1.size()*x2.size()*x3.size(),
                 i1data_sorted, x0.size()*x1.size()*x2.size()*x3.size(), 1.0, odata_sorted, ci0.size());
        }
      }
    }
  }

  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->add_block(odata, ci0);
}


void Task904::Task_local::compute() {
  const Index ci0 = b(0);

  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);

  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        for (auto& x3 : *range_[1]) {
          for (auto& x4 : *range_[1]) {
            for (auto& x5 : *range_[1]) {
              std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1, x2, x3, x4, x5);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x0, x1, x2, x3, x4, x5)]);
              sort_indices<1,2,3,4,5,6,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x0.size(), x1.size(), x2.size(), x3.size(), x4.size(), x5.size());

              std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, x2, x3, x4, x5);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, x2, x3, x4, x5)]);
              sort_indices<0,1,2,3,4,5,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), x2.size(), x3.size(), x4.size(), x5.size());
              dgemm_("T", "N", ci0.size(), 1, x0.size()*x1.size()*x2.size()*x3.size()*x4.size()*x5.size(),
                     1.0, i0data_sorted, x0.size()*x1.size()*x2.size()*x3.size()*x4.size()*x5.size(),
                     i1data_sorted, x0.size()*x1.size()*x2.size()*x3.size()*x4.size()*x5.size(), 1.0, odata_sorted, ci0.size());
            }
          }
        }
      }
    }
  }

  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->add_block(odata, ci0);
}

void Task905::Task_local::compute() {
  const Index ci0 = b(0);

  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);

  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        for (auto& x3 : *range_[1]) {
          for (auto& x4 : *range_[1]) {
            for (auto& x5 : *range_[1]) {
              std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1, x2, x3, x4, x5);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x0, x1, x2, x3, x4, x5)]);
              sort_indices<1,2,3,4,5,6,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x0.size(), x1.size(), x2.size(), x3.size(), x4.size(), x5.size());

              std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, x2, x3, x4, x5);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, x2, x3, x4, x5)]);
              sort_indices<0,1,2,3,4,5,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), x2.size(), x3.size(), x4.size(), x5.size());
              dgemm_("T", "N", ci0.size(), 1, x0.size()*x1.size()*x2.size()*x3.size()*x4.size()*x5.size(),
                     1.0, i0data_sorted, x0.size()*x1.size()*x2.size()*x3.size()*x4.size()*x5.size(),
                     i1data_sorted, x0.size()*x1.size()*x2.size()*x3.size()*x4.size()*x5.size(), 1.0, odata_sorted, ci0.size());
            }
          }
        }
      }
    }
  }

  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->add_block(odata, ci0);
}

void Task914::Task_local::compute() {
  const Index ci0 = b(0);

  int jx0 = 0;
  for (auto& x0 : *range_[1]) {
    int jx1 = 0;
    const int x0offset = x0.size() * jx0;
    for (auto& x1 : *range_[1]) {
      const int x1offset = x1.size() * jx1;
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
      const int lena = ciwfn_->det()->lena();
      const int lenb = ciwfn_->det()->lenb();

      int jci1 = 0;
      for (auto& ci1 : *range_[3]) {
        const int ci1offset = ci1.size() * jci1;
        std::unique_ptr<double[]> odata(new double[out()->get_size(ci1)]);
        std::fill_n(odata.get(), out()->get_size(ci1), 0.0);

        for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
          int x0p = ix0 + x0offset;
          for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
            int x1p = ix1 + x1offset;
            for (auto& iter : ciwfn_->det()->phia(x0p, x1p)) {
              size_t iaJ = iter.source;
              size_t iaK = iter.target;
              double sign = static_cast<double>(iter.sign);
              for (size_t ib = 0; ib != lenb; ++ib) {
                size_t iK = ib+iaK*lenb;
                size_t iJ = ib+iaJ*lenb;
                if (iK < ci0offset_ || iK >= (ci0offset_ + ci0.size())) continue;
                if (iJ < ci1offset || iJ >= (ci1offset + ci1.size())) continue;

                odata[iJ-ci1offset] += sign * i0data[(iK-ci0offset_)+ci0.size()*(ix0+x0.size()*ix1)];
              }
            }

            for (auto& iter : ciwfn_->det()->phib(x0p, x1p)) {
              size_t ibJ = iter.source;
              size_t ibK = iter.target;
              double sign = static_cast<double>(iter.sign);
              for (size_t ia = 0; ia != lena; ++ia) {
                size_t iK = ibK+ia*lenb;
                size_t iJ = ibJ+ia*lenb;
                if (iK < ci0offset_ || iK >= (ci0offset_ + ci0.size())) continue;
                if (iJ < ci1offset || iJ >= (ci1offset + ci1.size())) continue;

                odata[iJ-ci1offset] += sign * i0data[(iK-ci0offset_)+ci0.size()*(ix0+x0.size()*ix1)];
              }
            }
          }
        }
        out()->add_block(odata, ci1);
        ++jci1;
      }
      ++jx1;
    }
    ++jx0;
  }
}


void Task915::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);

  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(ci0, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1), 0.0);

  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        for (auto& x5 : *range_[1]) {
          std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x2, x3, x4, x5);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x2, x3, x4, x5)]);
          sort_indices<1,2,3,4,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x2.size(), x3.size(), x4.size(), x5.size());

          std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, x2, x3, x4, x5);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, x2, x3, x4, x5)]);
          sort_indices<2,3,4,5,0,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), x2.size(), x3.size(), x4.size(), x5.size());

          dgemm_("T", "N", ci0.size(), x0.size()*x1.size(), x2.size()*x3.size()*x4.size()*x5.size(),
                 1.0, i0data_sorted, x2.size()*x3.size()*x4.size()*x5.size(),
                 i1data_sorted, x2.size()*x3.size()*x4.size()*x5.size(), 1.0, odata_sorted, ci0.size());
        }
      }
    }
  }

  sort_indices<0,1,2,1,1,1,1>(odata_sorted, odata, ci0.size(), x0.size(), x1.size());
  out()->add_block(odata, ci0, x0, x1);
}


void Task916::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  const Index x3 = b(4);

  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1, x2, x3);

  for (auto& x4 : *range_[1]) {
    {
      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x4, x4, x1, x2, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x4, x4, x1, x2, x3)]);
      sort_indices<0,3,4,5,1,2,0,1,1,1>(i1data, i1data_sorted, x0.size(), x4.size(), x4.size(), x1.size(), x2.size(), x3.size());
      for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
            for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
              for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0] -= i0data[ici0+ci0.size()*(ix0+x0.size()*(ix1+x1.size()*(ix2+x2.size()*ix3)))]
                               * i1data_sorted[ix0+x0.size()*(ix1+x1.size()*(ix2+x2.size()*(ix3+x3.size()*(ix4+x4.size()*ix4))))];
                }
              }
            }
          }
        }
      }
    }
    {
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x4, x0, x1, x4, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x4, x0, x1, x4, x3)]);
      sort_indices<2,3,0,5,1,4,0,1,1,1>(i1data, i1data_sorted, x2.size(), x4.size(), x0.size(), x1.size(), x4.size(), x3.size());
      for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
            for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
              for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0] -= i0data[ici0+ci0.size()*(ix0+x0.size()*(ix1+x1.size()*(ix2+x2.size()*ix3)))]
                               * i1data_sorted[ix0+x0.size()*(ix1+x1.size()*(ix2+x2.size()*(ix3+x3.size()*(ix4+x4.size()*ix4))))];
                }
              }
            }
          }
        }
      }
    }
  }
  out()->add_block(odata, ci0);
}

#endif
