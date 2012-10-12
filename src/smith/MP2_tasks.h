//
// Newint - Parallel electron correlation program.
// Filename: MP2_tasks.h
// Copyright (C) 2012 Shiozaki group
//
// Author: Shiozaki group <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The Newint package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#ifndef __SRC_SMITH_MP2_TASKS_H 
#define __SRC_SMITH_MP2_TASKS_H 

#include <memory>
#include <algorithm>
#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <vector>

namespace bagel {
namespace SMITH {
namespace MP2{

template <typename T>
class Task0 : public Task<T> {
  protected:
    std::shared_ptr<Tensor<T> > r_;
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;

    void compute_() {
      r_->zero();
    };  

  public:
    Task0(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      r_ =  t[0];
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
    };  
    ~Task0() {}; 
};

template <typename T>
class Task1 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > r;
    std::shared_ptr<Tensor<T> > I0;

    void compute_() {
      for (auto& a2 : virt_) {
        for (auto& c1 : closed_) {
          for (auto& a4 : virt_) {
            for (auto& c3 : closed_) {
              std::vector<size_t> ohash = {c3.key(), a4.key(), c1.key(), a2.key()};
              std::unique_ptr<double[]> odata = r->move_block(ohash);
              {
                std::vector<size_t> i0hash = {c1.key(), a4.key(), c3.key(), a2.key()};
                std::unique_ptr<double[]> i0data = I0->get_block(i0hash);
                sort_indices<2,1,0,3,1,1,1,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
              }
              r->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task1(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      r = t[0];
      I0 = t[1];
    };
    ~Task1() {};
};

template <typename T>
class Task2 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I0;
    std::shared_ptr<Tensor<T> > v2;

    void compute_() {
      for (auto& a2 : virt_) {
        for (auto& c3 : closed_) {
          for (auto& a4 : virt_) {
            for (auto& c1 : closed_) {
              std::vector<size_t> ohash = {c1.key(), a4.key(), c3.key(), a2.key()};
              std::unique_ptr<double[]> odata = I0->move_block(ohash);
              {
                std::vector<size_t> i0hash = {c1.key(), a4.key(), c3.key(), a2.key()};
                std::unique_ptr<double[]> i0data = v2->get_block(i0hash);
                sort_indices<0,1,2,3,1,1,-4,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
              }
              {
                std::vector<size_t> i1hash = {c1.key(), a2.key(), c3.key(), a4.key()};
                std::unique_ptr<double[]> i1data = v2->get_block(i1hash);
                sort_indices<0,3,2,1,1,1,8,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a4.size());
              }
              I0->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task2(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I0 = t[0];
      v2 = t[1];
    };
    ~Task2() {};
};

#if 0
template <typename T>
class Task3 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I0;
    std::shared_ptr<Tensor<T> > t2;
    std::shared_ptr<Tensor<T> > I1;

    void compute_() {
      for (auto& a2 : virt_) {
        for (auto& c3 : closed_) {
          for (auto& a4 : virt_) {
            for (auto& c1 : closed_) {
              std::vector<size_t> ohash = {c1.key(), a4.key(), c3.key(), a2.key()};
              std::unique_ptr<double[]> odata = I0->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I0->get_size(ohash)]);
              std::fill(odata_sorted.get(), odata_sorted.get()+I0->get_size(ohash), 0.0);

              std::vector<size_t> i0hash = {c1.key(), a4.key(), c3.key(), a2.key()};
              std::unique_ptr<double[]> i0data = t2->get_block(i0hash);
              std::unique_ptr<double[]> i0data_sorted(new double[t2->get_size(i0hash)]);
              sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());

              std::vector<size_t> i1hash = {};
              std::unique_ptr<double[]> i1data = I1->get_block(i1hash);
              std::unique_ptr<double[]> i1data_sorted(new double[I1->get_size(i1hash)]);
              sort_indices<0,1,1,1>(i1data, i1data_sorted);

              dgemm_("T", "N", c1.size()*a4.size()*c3.size()*a2.size(), 1, 1,
                     1.0, i0data_sorted, 1, i1data_sorted, 1,
                     1.0, odata_sorted, c1.size()*a4.size()*c3.size()*a2.size());

              sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
              I0->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task3(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I0 = t[0];
      t2 = t[1];
      I1 = t[2];
    };
    ~Task3() {};
};

template <typename T>
class Task4 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I1;
    std::shared_ptr<Tensor<T> > Gamma;

    void compute_() {
      std::vector<size_t> ohash = {};
      std::unique_ptr<double[]> odata = I1->move_block(ohash);
      {
        std::vector<size_t> i0hash = {};
        std::unique_ptr<double[]> i0data = Gamma->get_block(i0hash);
        sort_indices<1,1,-4,1>(i0data, odata);
      }
      I1->put_block(ohash, odata);
    };

  public:
    Task4(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I1 = t[0];
      Gamma = t[1];
    };
    ~Task4() {};
};

template <typename T>
class Task5 : public Task<T> {  // associated with gamma
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > Gamma;
    std::shared_ptr<Tensor<T> > rdm1;
    std::shared_ptr<Tensor<T> > f1;

    void compute_() {
      std::vector<size_t> ohash = {};
      std::unique_ptr<double[]> odata = Gamma->move_block(ohash);
      // associated with merged
      for (auto& x0 : active_) {
        for (auto& x1 : active_) {
          std::vector<size_t> fhash = {x1.key(), x0.key()};
          std::unique_ptr<double[]> fdata = f1->get_block(fhash);
          {
            std::vector<size_t> i0hash = {x1.key(), x0.key()};
            std::unique_ptr<double[]> data = rdm1->get_block(i0hash);
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                odata[]
                  += (1.0) * data[i1+x1.size()*(i0)] * fdata[i1+x1.size()*(i0)];
              }
            }
          }
        }
      }
      Gamma->put_block(ohash, odata);
    };  


  public:
    Task5(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      Gamma   = t[0];
      rdm1    = t[1];
      f1      = t[2];
    };
    ~Task5() {};
};

template <typename T>
class Task6 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I0;
    std::shared_ptr<Tensor<T> > t2;
    std::shared_ptr<Tensor<T> > I3;

    void compute_() {
      for (auto& a2 : virt_) {
        for (auto& c3 : closed_) {
          for (auto& a4 : virt_) {
            for (auto& c1 : closed_) {
              std::vector<size_t> ohash = {c1.key(), a4.key(), c3.key(), a2.key()};
              std::unique_ptr<double[]> odata = I0->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I0->get_size(ohash)]);
              std::fill(odata_sorted.get(), odata_sorted.get()+I0->get_size(ohash), 0.0);

              std::vector<size_t> i0hash = {c1.key(), a2.key(), c3.key(), a4.key()};
              std::unique_ptr<double[]> i0data = t2->get_block(i0hash);
              std::unique_ptr<double[]> i0data_sorted(new double[t2->get_size(i0hash)]);
              sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());

              std::vector<size_t> i1hash = {};
              std::unique_ptr<double[]> i1data = I3->get_block(i1hash);
              std::unique_ptr<double[]> i1data_sorted(new double[I3->get_size(i1hash)]);
              sort_indices<0,1,1,1>(i1data, i1data_sorted);

              dgemm_("T", "N", c1.size()*a2.size()*c3.size()*a4.size(), 1, 1,
                     1.0, i0data_sorted, 1, i1data_sorted, 1,
                     1.0, odata_sorted, c1.size()*a2.size()*c3.size()*a4.size());

              sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
              I0->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task6(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I0 = t[0];
      t2 = t[1];
      I3 = t[2];
    };
    ~Task6() {};
};

template <typename T>
class Task7 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I3;
    std::shared_ptr<Tensor<T> > Gamma;

    void compute_() {
      std::vector<size_t> ohash = {};
      std::unique_ptr<double[]> odata = I3->move_block(ohash);
      {
        std::vector<size_t> i0hash = {};
        std::unique_ptr<double[]> i0data = Gamma->get_block(i0hash);
        sort_indices<1,1,8,1>(i0data, odata);
      }
      I3->put_block(ohash, odata);
    };

  public:
    Task7(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I3 = t[0];
      Gamma = t[1];
    };
    ~Task7() {};
};

template <typename T>
class Task8 : public Task<T> {  // associated with gamma
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > Gamma;
    std::shared_ptr<Tensor<T> > rdm1;
    std::shared_ptr<Tensor<T> > f1;

    void compute_() {
      std::vector<size_t> ohash = {};
      std::unique_ptr<double[]> odata = Gamma->move_block(ohash);
      // associated with merged
      for (auto& x0 : active_) {
        for (auto& x1 : active_) {
          std::vector<size_t> fhash = {x1.key(), x0.key()};
          std::unique_ptr<double[]> fdata = f1->get_block(fhash);
          {
            std::vector<size_t> i0hash = {x1.key(), x0.key()};
            std::unique_ptr<double[]> data = rdm1->get_block(i0hash);
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                odata[]
                  += (1.0) * data[i1+x1.size()*(i0)] * fdata[i1+x1.size()*(i0)];
              }
            }
          }
        }
      }
      Gamma->put_block(ohash, odata);
    };  


  public:
    Task8(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      Gamma   = t[0];
      rdm1    = t[1];
      f1      = t[2];
    };
    ~Task8() {};
};

#endif
template <typename T>
class Task9 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > r;
    std::shared_ptr<Tensor<T> > I4;

    void compute_() {
      for (auto& a2 : virt_) {
        for (auto& c1 : closed_) {
          for (auto& a4 : virt_) {
            for (auto& c3 : closed_) {
              std::vector<size_t> ohash = {c3.key(), a4.key(), c1.key(), a2.key()};
              std::unique_ptr<double[]> odata = r->move_block(ohash);
              {
                std::vector<size_t> i0hash = {c1.key(), a4.key(), a2.key(), c3.key()};
                std::unique_ptr<double[]> i0data = I4->get_block(i0hash);
                sort_indices<3,1,0,2,1,1,1,1>(i0data, odata, c1.size(), a4.size(), a2.size(), c3.size());
              }
              {
                std::vector<size_t> i0hash = {c3.key(), a2.key(), a4.key(), c1.key()};
                std::unique_ptr<double[]> i0data = I4->get_block(i0hash);
                sort_indices<0,2,3,1,1,1,1,1>(i0data, odata, c3.size(), a2.size(), a4.size(), c1.size());
              }
              r->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task9(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      r = t[0];
      I4 = t[1];
    };
    ~Task9() {};
};

template <typename T>
class Task10 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I4;
    std::shared_ptr<Tensor<T> > f1;
    std::shared_ptr<Tensor<T> > I5;

    void compute_() {
      for (auto& c3 : closed_) {
        for (auto& a2 : virt_) {
          for (auto& a4 : virt_) {
            for (auto& c1 : closed_) {
              std::vector<size_t> ohash = {c1.key(), a4.key(), a2.key(), c3.key()};
              std::unique_ptr<double[]> odata = I4->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I4->get_size(ohash)]);
              std::fill(odata_sorted.get(), odata_sorted.get()+I4->get_size(ohash), 0.0);

              for (auto& c5 : closed_) {
                std::vector<size_t> i0hash = {c3.key(), c5.key()};
                std::unique_ptr<double[]> i0data = f1->get_block(i0hash);
                std::unique_ptr<double[]> i0data_sorted(new double[f1->get_size(i0hash)]);
                sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c3.size(), c5.size());

                std::vector<size_t> i1hash = {c1.key(), a4.key(), c5.key(), a2.key()};
                std::unique_ptr<double[]> i1data = I5->get_block(i1hash);
                std::unique_ptr<double[]> i1data_sorted(new double[I5->get_size(i1hash)]);
                sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a4.size(), c5.size(), a2.size());

                dgemm_("T", "N", c3.size(), c1.size()*a4.size()*a2.size(), c5.size(),
                       1.0, i0data_sorted, c5.size(), i1data_sorted, c5.size(),
                       1.0, odata_sorted, c3.size());
              }

              sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c3.size(), c1.size(), a4.size(), a2.size());
              I4->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task10(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I4 = t[0];
      f1 = t[1];
      I5 = t[2];
    };
    ~Task10() {};
};

template <typename T>
class Task11 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I5;
    std::shared_ptr<Tensor<T> > t2;

    void compute_() {
      for (auto& a2 : virt_) {
        for (auto& c5 : closed_) {
          for (auto& a4 : virt_) {
            for (auto& c1 : closed_) {
              std::vector<size_t> ohash = {c1.key(), a4.key(), c5.key(), a2.key()};
              std::unique_ptr<double[]> odata = I5->move_block(ohash);
              {
                std::vector<size_t> i0hash = {c1.key(), a4.key(), c5.key(), a2.key()};
                std::unique_ptr<double[]> i0data = t2->get_block(i0hash);
                sort_indices<0,1,2,3,1,1,4,1>(i0data, odata, c1.size(), a4.size(), c5.size(), a2.size());
              }
              {
                std::vector<size_t> i1hash = {c1.key(), a2.key(), c5.key(), a4.key()};
                std::unique_ptr<double[]> i1data = t2->get_block(i1hash);
                sort_indices<0,3,2,1,1,1,-8,1>(i1data, odata, c1.size(), a2.size(), c5.size(), a4.size());
              }
              I5->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task11(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I5 = t[0];
      t2 = t[1];
    };
    ~Task11() {};
};

template <typename T>
class Task12 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I4;
    std::shared_ptr<Tensor<T> > f1;
    std::shared_ptr<Tensor<T> > I9;

    void compute_() {
      for (auto& c3 : closed_) {
        for (auto& a2 : virt_) {
          for (auto& a4 : virt_) {
            for (auto& c1 : closed_) {
              std::vector<size_t> ohash = {c1.key(), a4.key(), a2.key(), c3.key()};
              std::unique_ptr<double[]> odata = I4->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I4->get_size(ohash)]);
              std::fill(odata_sorted.get(), odata_sorted.get()+I4->get_size(ohash), 0.0);

              for (auto& a5 : virt_) {
                std::vector<size_t> i0hash = {a5.key(), a4.key()};
                std::unique_ptr<double[]> i0data = f1->get_block(i0hash);
                std::unique_ptr<double[]> i0data_sorted(new double[f1->get_size(i0hash)]);
                sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a5.size(), a4.size());

                std::vector<size_t> i1hash = {c1.key(), a5.key(), c3.key(), a2.key()};
                std::unique_ptr<double[]> i1data = I9->get_block(i1hash);
                std::unique_ptr<double[]> i1data_sorted(new double[I9->get_size(i1hash)]);
                sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a5.size(), c3.size(), a2.size());

                dgemm_("T", "N", a4.size(), c1.size()*c3.size()*a2.size(), a5.size(),
                       1.0, i0data_sorted, a5.size(), i1data_sorted, a5.size(),
                       1.0, odata_sorted, a4.size());
              }

              sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), c3.size(), a2.size());
              I4->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task12(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I4 = t[0];
      f1 = t[1];
      I9 = t[2];
    };
    ~Task12() {};
};

template <typename T>
class Task13 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I9;
    std::shared_ptr<Tensor<T> > t2;

    void compute_() {
      for (auto& a2 : virt_) {
        for (auto& c3 : closed_) {
          for (auto& a5 : virt_) {
            for (auto& c1 : closed_) {
              std::vector<size_t> ohash = {c1.key(), a5.key(), c3.key(), a2.key()};
              std::unique_ptr<double[]> odata = I9->move_block(ohash);
              {
                std::vector<size_t> i0hash = {c1.key(), a5.key(), c3.key(), a2.key()};
                std::unique_ptr<double[]> i0data = t2->get_block(i0hash);
                sort_indices<0,1,2,3,1,1,-4,1>(i0data, odata, c1.size(), a5.size(), c3.size(), a2.size());
              }
              {
                std::vector<size_t> i1hash = {c1.key(), a2.key(), c3.key(), a5.key()};
                std::unique_ptr<double[]> i1data = t2->get_block(i1hash);
                sort_indices<0,3,2,1,1,1,8,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a5.size());
              }
              I9->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task13(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I9 = t[0];
      t2 = t[1];
    };
    ~Task13() {};
};

template <typename T>
class Task14 : public EnergyTask<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I14;
    std::shared_ptr<Tensor<T> > t2;
    std::shared_ptr<Tensor<T> > I15;

    void compute_() {
      this->energy_ = 0.0;

      for (auto& a4 : virt_) {
        for (auto& c3 : closed_) {
          for (auto& a2 : virt_) {
            for (auto& c1 : closed_) {
              std::vector<size_t> i0hash = {c3.key(), a4.key(), c1.key(), a2.key()};
              std::unique_ptr<double[]> i0data = t2->get_block(i0hash);
              std::unique_ptr<double[]> i0data_sorted(new double[t2->get_size(i0hash)]);
              sort_indices<1,0,3,2,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), c1.size(), a2.size());

              std::vector<size_t> i1hash = {c1.key(), a4.key(), c3.key(), a2.key()};
              std::unique_ptr<double[]> i1data = I15->get_block(i1hash);
              std::unique_ptr<double[]> i1data_sorted(new double[I15->get_size(i1hash)]);
              sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), a4.size(), c3.size(), a2.size());

              this->energy_ += ddot_(c1.size()*a4.size()*c3.size()*a2.size(), i0data_sorted, 1, i1data_sorted, 1);
            }
          }
        }
      }

    };

  public:
    Task14(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : EnergyTask<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I14 = t[0];
      t2 = t[1];
      I15 = t[2];
    };
    ~Task14() {};
};

template <typename T>
class Task15 : public EnergyTask<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I15;
    std::shared_ptr<Tensor<T> > v2;
    std::shared_ptr<Tensor<T> > r;

    void compute_() {
      this->energy_ = 0.0;
      for (auto& a2 : virt_) {
        for (auto& c3 : closed_) {
          for (auto& a4 : virt_) {
            for (auto& c1 : closed_) {
              std::vector<size_t> ohash = {c1.key(), a4.key(), c3.key(), a2.key()};
              std::unique_ptr<double[]> odata = I15->move_block(ohash);
              {
                std::vector<size_t> i0hash = {c1.key(), a4.key(), c3.key(), a2.key()};
                std::unique_ptr<double[]> i0data = v2->get_block(i0hash);
                sort_indices<0,1,2,3,1,1,-4,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
              }
              {
                std::vector<size_t> i1hash = {c1.key(), a2.key(), c3.key(), a4.key()};
                std::unique_ptr<double[]> i1data = v2->get_block(i1hash);
                sort_indices<0,3,2,1,1,1,8,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a4.size());
              }
              {
                std::vector<size_t> i2hash = {c1.key(), a4.key(), c3.key(), a2.key()};
                std::unique_ptr<double[]> i2data = r->get_block(i2hash);
                sort_indices<0,1,2,3,1,1,-4,1>(i2data, odata, c1.size(), a4.size(), c3.size(), a2.size());
              }
              {
                std::vector<size_t> i3hash = {c1.key(), a2.key(), c3.key(), a4.key()};
                std::unique_ptr<double[]> i3data = r->get_block(i3hash);
                sort_indices<0,3,2,1,1,1,8,1>(i3data, odata, c1.size(), a2.size(), c3.size(), a4.size());
              }
              I15->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task15(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : EnergyTask<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I15 = t[0];
      v2 = t[1];
      r = t[2];
    };
    ~Task15() {};
};


}
}
}
#endif

