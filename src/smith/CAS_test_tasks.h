//
// BAGEL - Parallel electron correlation program.
// Filename: CAS_test_tasks.h
// Copyright (C) 2012 Shiozaki group
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


#ifndef __SRC_SMITH_CAS_test_TASKS_H 
#define __SRC_SMITH_CAS_test_TASKS_H 

#include <memory>
#include <algorithm>
#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <vector>

namespace bagel {
namespace SMITH {
namespace CAS_test{

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
class Task1 : public Task<T> {  // associated with gamma
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > Gamma0;
    std::shared_ptr<Tensor<T> > rdm2;
    std::shared_ptr<Tensor<T> > f1;

    void compute_() {
      for (auto& x0 : active_) {
        for (auto& x3 : active_) {
          std::vector<size_t> ohash = {x3.key(), x0.key()};
          std::unique_ptr<double[]> odata = Gamma0->move_block(ohash);
          // associated with merged
          for (auto& x1 : active_) {
            for (auto& x2 : active_) {
              std::vector<size_t> fhash = {x2.key(), x1.key()};
              std::unique_ptr<double[]> fdata = f1->get_block(fhash);
              {
                std::vector<size_t> i0hash = {x3.key(), x0.key(), x2.key(), x1.key()};
                std::unique_ptr<double[]> i0data = rdm2->get_block(i0hash);
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
            }
          }
          Gamma0->put_block(ohash, odata);
        }
      }
    };  


  public:
    Task1(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      Gamma0  = t[0];
      rdm2    = t[1];
      f1      = t[2];
    };
    ~Task1() {};
};

template <typename T>
class Task2 : public Task<T> {  // associated with gamma
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > Gamma2;
    std::shared_ptr<Tensor<T> > rdm1;

    void compute_() {
      for (auto& x0 : active_) {
        for (auto& x1 : active_) {
          std::vector<size_t> ohash = {x1.key(), x0.key()};
          std::unique_ptr<double[]> odata = Gamma2->move_block(ohash);
          {
            std::vector<size_t> i0hash = {x1.key(), x0.key()};
            std::unique_ptr<double[]> i0data = rdm1->get_block(i0hash);
            sort_indices<0,1,1,1,1,1>(i0data, odata, x1.size(), x0.size());
          }
          Gamma2->put_block(ohash, odata);
        }
      }
    };  


  public:
    Task2(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      Gamma2  = t[0];
      rdm1    = t[1];
    };
    ~Task2() {};
};

template <typename T>
class Task3 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > r;
    std::shared_ptr<Tensor<T> > I0;

    void compute_() {
      for (auto& a2 : virt_) {
        for (auto& c1 : closed_) {
          for (auto& a3 : virt_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), a3.key(), c1.key(), a2.key()};
              std::unique_ptr<double[]> odata = r->move_block(ohash);
              {
                std::vector<size_t> i0hash = {x0.key(), a3.key(), c1.key(), a2.key()};
                std::unique_ptr<double[]> i0data = I0->get_block(i0hash);
                sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x0.size(), a3.size(), c1.size(), a2.size());
              }
              r->put_block(ohash, odata);
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
      r = t[0];
      I0 = t[1];
    };
    ~Task3() {};
};

template <typename T>
class Task4 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I0;
    std::shared_ptr<Tensor<T> > t2;
    std::shared_ptr<Tensor<T> > I1;

    void compute_() {
      for (auto& a2 : virt_) {
        for (auto& c1 : closed_) {
          for (auto& a3 : virt_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), a3.key(), c1.key(), a2.key()};
              std::unique_ptr<double[]> odata = I0->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I0->get_size(ohash)]);
              std::fill_n(odata_sorted.get(), I0->get_size(ohash), 0.0);

              for (auto& x3 : active_) {
                std::vector<size_t> i0hash = {x3.key(), a3.key(), c1.key(), a2.key()};
                std::unique_ptr<double[]> i0data = t2->get_block(i0hash);
                std::unique_ptr<double[]> i0data_sorted(new double[t2->get_size(i0hash)]);
                sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c1.size(), a2.size());

                std::vector<size_t> i1hash = {x3.key(), x0.key()};
                std::unique_ptr<double[]> i1data = I1->get_block(i1hash);
                std::unique_ptr<double[]> i1data_sorted(new double[I1->get_size(i1hash)]);
                sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size());

                dgemm_("T", "N", a3.size()*c1.size()*a2.size(), x0.size(), x3.size(),
                       1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
                       1.0, odata_sorted, a3.size()*c1.size()*a2.size());
              }

              sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c1.size(), a2.size(), x0.size());
              I0->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task4(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I0 = t[0];
      t2 = t[1];
      I1 = t[2];
    };
    ~Task4() {};
};

template <typename T>
class Task5 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I1;
    std::shared_ptr<Tensor<T> > Gamma0;

    void compute_() {
      for (auto& x0 : active_) {
        for (auto& x3 : active_) {
          std::vector<size_t> ohash = {x3.key(), x0.key()};
          std::unique_ptr<double[]> odata = I1->move_block(ohash);
          {
            std::vector<size_t> i0hash = {x3.key(), x0.key()};
            std::unique_ptr<double[]> i0data = Gamma0->get_block(i0hash);
            sort_indices<0,1,1,1,2,1>(i0data, odata, x3.size(), x0.size());
          }
          I1->put_block(ohash, odata);
        }
      }
    };

  public:
    Task5(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I1 = t[0];
      Gamma0 = t[1];
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
        for (auto& c1 : closed_) {
          for (auto& a3 : virt_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), a3.key(), c1.key(), a2.key()};
              std::unique_ptr<double[]> odata = I0->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I0->get_size(ohash)]);
              std::fill_n(odata_sorted.get(), I0->get_size(ohash), 0.0);

              for (auto& x3 : active_) {
                std::vector<size_t> i0hash = {x3.key(), a2.key(), c1.key(), a3.key()};
                std::unique_ptr<double[]> i0data = t2->get_block(i0hash);
                std::unique_ptr<double[]> i0data_sorted(new double[t2->get_size(i0hash)]);
                sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c1.size(), a3.size());

                std::vector<size_t> i1hash = {x3.key(), x0.key()};
                std::unique_ptr<double[]> i1data = I3->get_block(i1hash);
                std::unique_ptr<double[]> i1data_sorted(new double[I3->get_size(i1hash)]);
                sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size());

                dgemm_("T", "N", a2.size()*c1.size()*a3.size(), x0.size(), x3.size(),
                       1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
                       1.0, odata_sorted, a2.size()*c1.size()*a3.size());
              }

              sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), a3.size(), x0.size());
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
    std::shared_ptr<Tensor<T> > Gamma0;

    void compute_() {
      for (auto& x0 : active_) {
        for (auto& x3 : active_) {
          std::vector<size_t> ohash = {x3.key(), x0.key()};
          std::unique_ptr<double[]> odata = I3->move_block(ohash);
          {
            std::vector<size_t> i0hash = {x3.key(), x0.key()};
            std::unique_ptr<double[]> i0data = Gamma0->get_block(i0hash);
            sort_indices<0,1,1,1,-1,1>(i0data, odata, x3.size(), x0.size());
          }
          I3->put_block(ohash, odata);
        }
      }
    };

  public:
    Task7(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I3 = t[0];
      Gamma0 = t[1];
    };
    ~Task7() {};
};

template <typename T>
class Task8 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I0;
    std::shared_ptr<Tensor<T> > f1;
    std::shared_ptr<Tensor<T> > I5;

    void compute_() {
      for (auto& a2 : virt_) {
        for (auto& c1 : closed_) {
          for (auto& a3 : virt_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), a3.key(), c1.key(), a2.key()};
              std::unique_ptr<double[]> odata = I0->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I0->get_size(ohash)]);
              std::fill_n(odata_sorted.get(), I0->get_size(ohash), 0.0);

              for (auto& c4 : closed_) {
                std::vector<size_t> i0hash = {c1.key(), c4.key()};
                std::unique_ptr<double[]> i0data = f1->get_block(i0hash);
                std::unique_ptr<double[]> i0data_sorted(new double[f1->get_size(i0hash)]);
                sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), c4.size());

                std::vector<size_t> i1hash = {x0.key(), a3.key(), c4.key(), a2.key()};
                std::unique_ptr<double[]> i1data = I5->get_block(i1hash);
                std::unique_ptr<double[]> i1data_sorted(new double[I5->get_size(i1hash)]);
                sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a3.size(), c4.size(), a2.size());

                dgemm_("T", "N", c1.size(), x0.size()*a3.size()*a2.size(), c4.size(),
                       1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
                       1.0, odata_sorted, c1.size());
              }

              sort_indices<1,2,0,3,1,1,1,1>(odata_sorted, odata, c1.size(), x0.size(), a3.size(), a2.size());
              I0->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task8(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I0 = t[0];
      f1 = t[1];
      I5 = t[2];
    };
    ~Task8() {};
};

template <typename T>
class Task9 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I5;
    std::shared_ptr<Tensor<T> > t2;
    std::shared_ptr<Tensor<T> > I6;

    void compute_() {
      for (auto& a2 : virt_) {
        for (auto& c4 : closed_) {
          for (auto& a3 : virt_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), a3.key(), c4.key(), a2.key()};
              std::unique_ptr<double[]> odata = I5->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I5->get_size(ohash)]);
              std::fill_n(odata_sorted.get(), I5->get_size(ohash), 0.0);

              for (auto& x1 : active_) {
                std::vector<size_t> i0hash = {x1.key(), a3.key(), c4.key(), a2.key()};
                std::unique_ptr<double[]> i0data = t2->get_block(i0hash);
                std::unique_ptr<double[]> i0data_sorted(new double[t2->get_size(i0hash)]);
                sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c4.size(), a2.size());

                std::vector<size_t> i1hash = {x1.key(), x0.key()};
                std::unique_ptr<double[]> i1data = I6->get_block(i1hash);
                std::unique_ptr<double[]> i1data_sorted(new double[I6->get_size(i1hash)]);
                sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

                dgemm_("T", "N", a3.size()*c4.size()*a2.size(), x0.size(), x1.size(),
                       1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                       1.0, odata_sorted, a3.size()*c4.size()*a2.size());
              }

              sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c4.size(), a2.size(), x0.size());
              I5->put_block(ohash, odata);
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
      I5 = t[0];
      t2 = t[1];
      I6 = t[2];
    };
    ~Task9() {};
};

template <typename T>
class Task10 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I6;
    std::shared_ptr<Tensor<T> > Gamma2;

    void compute_() {
      for (auto& x0 : active_) {
        for (auto& x1 : active_) {
          std::vector<size_t> ohash = {x1.key(), x0.key()};
          std::unique_ptr<double[]> odata = I6->move_block(ohash);
          {
            std::vector<size_t> i0hash = {x1.key(), x0.key()};
            std::unique_ptr<double[]> i0data = Gamma2->get_block(i0hash);
            sort_indices<0,1,1,1,-2,1>(i0data, odata, x1.size(), x0.size());
          }
          I6->put_block(ohash, odata);
        }
      }
    };

  public:
    Task10(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I6 = t[0];
      Gamma2 = t[1];
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
    std::shared_ptr<Tensor<T> > I9;

    void compute_() {
      for (auto& a2 : virt_) {
        for (auto& c4 : closed_) {
          for (auto& a3 : virt_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), a3.key(), c4.key(), a2.key()};
              std::unique_ptr<double[]> odata = I5->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I5->get_size(ohash)]);
              std::fill_n(odata_sorted.get(), I5->get_size(ohash), 0.0);

              for (auto& x1 : active_) {
                std::vector<size_t> i0hash = {x1.key(), a2.key(), c4.key(), a3.key()};
                std::unique_ptr<double[]> i0data = t2->get_block(i0hash);
                std::unique_ptr<double[]> i0data_sorted(new double[t2->get_size(i0hash)]);
                sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c4.size(), a3.size());

                std::vector<size_t> i1hash = {x1.key(), x0.key()};
                std::unique_ptr<double[]> i1data = I9->get_block(i1hash);
                std::unique_ptr<double[]> i1data_sorted(new double[I9->get_size(i1hash)]);
                sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

                dgemm_("T", "N", a2.size()*c4.size()*a3.size(), x0.size(), x1.size(),
                       1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                       1.0, odata_sorted, a2.size()*c4.size()*a3.size());
              }

              sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c4.size(), a3.size(), x0.size());
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
      I9 = t[2];
    };
    ~Task11() {};
};

template <typename T>
class Task12 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I9;
    std::shared_ptr<Tensor<T> > Gamma2;

    void compute_() {
      for (auto& x0 : active_) {
        for (auto& x1 : active_) {
          std::vector<size_t> ohash = {x1.key(), x0.key()};
          std::unique_ptr<double[]> odata = I9->move_block(ohash);
          {
            std::vector<size_t> i0hash = {x1.key(), x0.key()};
            std::unique_ptr<double[]> i0data = Gamma2->get_block(i0hash);
            sort_indices<0,1,1,1,1,1>(i0data, odata, x1.size(), x0.size());
          }
          I9->put_block(ohash, odata);
        }
      }
    };

  public:
    Task12(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I9 = t[0];
      Gamma2 = t[1];
    };
    ~Task12() {};
};

template <typename T>
class Task13 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I0;
    std::shared_ptr<Tensor<T> > f1;
    std::shared_ptr<Tensor<T> > I11;

    void compute_() {
      for (auto& a2 : virt_) {
        for (auto& c1 : closed_) {
          for (auto& a3 : virt_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), a3.key(), c1.key(), a2.key()};
              std::unique_ptr<double[]> odata = I0->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I0->get_size(ohash)]);
              std::fill_n(odata_sorted.get(), I0->get_size(ohash), 0.0);

              for (auto& a4 : virt_) {
                std::vector<size_t> i0hash = {a4.key(), a2.key()};
                std::unique_ptr<double[]> i0data = f1->get_block(i0hash);
                std::unique_ptr<double[]> i0data_sorted(new double[f1->get_size(i0hash)]);
                sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a2.size());

                std::vector<size_t> i1hash = {x0.key(), a4.key(), c1.key(), a3.key()};
                std::unique_ptr<double[]> i1data = I11->get_block(i1hash);
                std::unique_ptr<double[]> i1data_sorted(new double[I11->get_size(i1hash)]);
                sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size(), c1.size(), a3.size());

                dgemm_("T", "N", a2.size(), x0.size()*c1.size()*a3.size(), a4.size(),
                       1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
                       1.0, odata_sorted, a2.size());
              }

              sort_indices<1,3,2,0,1,1,1,1>(odata_sorted, odata, a2.size(), x0.size(), c1.size(), a3.size());
              I0->put_block(ohash, odata);
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
      I0 = t[0];
      f1 = t[1];
      I11 = t[2];
    };
    ~Task13() {};
};

template <typename T>
class Task14 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I11;
    std::shared_ptr<Tensor<T> > t2;
    std::shared_ptr<Tensor<T> > I12;

    void compute_() {
      for (auto& a3 : virt_) {
        for (auto& c1 : closed_) {
          for (auto& a4 : virt_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), a4.key(), c1.key(), a3.key()};
              std::unique_ptr<double[]> odata = I11->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I11->get_size(ohash)]);
              std::fill_n(odata_sorted.get(), I11->get_size(ohash), 0.0);

              for (auto& x1 : active_) {
                std::vector<size_t> i0hash = {x1.key(), a4.key(), c1.key(), a3.key()};
                std::unique_ptr<double[]> i0data = t2->get_block(i0hash);
                std::unique_ptr<double[]> i0data_sorted(new double[t2->get_size(i0hash)]);
                sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c1.size(), a3.size());

                std::vector<size_t> i1hash = {x1.key(), x0.key()};
                std::unique_ptr<double[]> i1data = I12->get_block(i1hash);
                std::unique_ptr<double[]> i1data_sorted(new double[I12->get_size(i1hash)]);
                sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

                dgemm_("T", "N", a4.size()*c1.size()*a3.size(), x0.size(), x1.size(),
                       1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                       1.0, odata_sorted, a4.size()*c1.size()*a3.size());
              }

              sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), a3.size(), x0.size());
              I11->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task14(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I11 = t[0];
      t2 = t[1];
      I12 = t[2];
    };
    ~Task14() {};
};

template <typename T>
class Task15 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I12;
    std::shared_ptr<Tensor<T> > Gamma2;

    void compute_() {
      for (auto& x0 : active_) {
        for (auto& x1 : active_) {
          std::vector<size_t> ohash = {x1.key(), x0.key()};
          std::unique_ptr<double[]> odata = I12->move_block(ohash);
          {
            std::vector<size_t> i0hash = {x1.key(), x0.key()};
            std::unique_ptr<double[]> i0data = Gamma2->get_block(i0hash);
            sort_indices<0,1,1,1,-1,1>(i0data, odata, x1.size(), x0.size());
          }
          I12->put_block(ohash, odata);
        }
      }
    };

  public:
    Task15(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I12 = t[0];
      Gamma2 = t[1];
    };
    ~Task15() {};
};

template <typename T>
class Task16 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I11;
    std::shared_ptr<Tensor<T> > t2;
    std::shared_ptr<Tensor<T> > I15;

    void compute_() {
      for (auto& a3 : virt_) {
        for (auto& c1 : closed_) {
          for (auto& a4 : virt_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), a4.key(), c1.key(), a3.key()};
              std::unique_ptr<double[]> odata = I11->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I11->get_size(ohash)]);
              std::fill_n(odata_sorted.get(), I11->get_size(ohash), 0.0);

              for (auto& x1 : active_) {
                std::vector<size_t> i0hash = {x1.key(), a3.key(), c1.key(), a4.key()};
                std::unique_ptr<double[]> i0data = t2->get_block(i0hash);
                std::unique_ptr<double[]> i0data_sorted(new double[t2->get_size(i0hash)]);
                sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c1.size(), a4.size());

                std::vector<size_t> i1hash = {x1.key(), x0.key()};
                std::unique_ptr<double[]> i1data = I15->get_block(i1hash);
                std::unique_ptr<double[]> i1data_sorted(new double[I15->get_size(i1hash)]);
                sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

                dgemm_("T", "N", a3.size()*c1.size()*a4.size(), x0.size(), x1.size(),
                       1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                       1.0, odata_sorted, a3.size()*c1.size()*a4.size());
              }

              sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a3.size(), c1.size(), a4.size(), x0.size());
              I11->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task16(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I11 = t[0];
      t2 = t[1];
      I15 = t[2];
    };
    ~Task16() {};
};

template <typename T>
class Task17 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I15;
    std::shared_ptr<Tensor<T> > Gamma2;

    void compute_() {
      for (auto& x0 : active_) {
        for (auto& x1 : active_) {
          std::vector<size_t> ohash = {x1.key(), x0.key()};
          std::unique_ptr<double[]> odata = I15->move_block(ohash);
          {
            std::vector<size_t> i0hash = {x1.key(), x0.key()};
            std::unique_ptr<double[]> i0data = Gamma2->get_block(i0hash);
            sort_indices<0,1,1,1,2,1>(i0data, odata, x1.size(), x0.size());
          }
          I15->put_block(ohash, odata);
        }
      }
    };

  public:
    Task17(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I15 = t[0];
      Gamma2 = t[1];
    };
    ~Task17() {};
};

template <typename T>
class Task18 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I0;
    std::shared_ptr<Tensor<T> > f1;
    std::shared_ptr<Tensor<T> > I17;

    void compute_() {
      for (auto& a2 : virt_) {
        for (auto& c1 : closed_) {
          for (auto& a3 : virt_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), a3.key(), c1.key(), a2.key()};
              std::unique_ptr<double[]> odata = I0->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I0->get_size(ohash)]);
              std::fill_n(odata_sorted.get(), I0->get_size(ohash), 0.0);

              for (auto& a4 : virt_) {
                std::vector<size_t> i0hash = {a4.key(), a3.key()};
                std::unique_ptr<double[]> i0data = f1->get_block(i0hash);
                std::unique_ptr<double[]> i0data_sorted(new double[f1->get_size(i0hash)]);
                sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a3.size());

                std::vector<size_t> i1hash = {x0.key(), a4.key(), c1.key(), a2.key()};
                std::unique_ptr<double[]> i1data = I17->get_block(i1hash);
                std::unique_ptr<double[]> i1data_sorted(new double[I17->get_size(i1hash)]);
                sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size(), c1.size(), a2.size());

                dgemm_("T", "N", a3.size(), x0.size()*c1.size()*a2.size(), a4.size(),
                       1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
                       1.0, odata_sorted, a3.size());
              }

              sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, a3.size(), x0.size(), c1.size(), a2.size());
              I0->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task18(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I0 = t[0];
      f1 = t[1];
      I17 = t[2];
    };
    ~Task18() {};
};

template <typename T>
class Task19 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I17;
    std::shared_ptr<Tensor<T> > t2;
    std::shared_ptr<Tensor<T> > I18;

    void compute_() {
      for (auto& a2 : virt_) {
        for (auto& c1 : closed_) {
          for (auto& a4 : virt_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), a4.key(), c1.key(), a2.key()};
              std::unique_ptr<double[]> odata = I17->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I17->get_size(ohash)]);
              std::fill_n(odata_sorted.get(), I17->get_size(ohash), 0.0);

              for (auto& x1 : active_) {
                std::vector<size_t> i0hash = {x1.key(), a4.key(), c1.key(), a2.key()};
                std::unique_ptr<double[]> i0data = t2->get_block(i0hash);
                std::unique_ptr<double[]> i0data_sorted(new double[t2->get_size(i0hash)]);
                sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c1.size(), a2.size());

                std::vector<size_t> i1hash = {x1.key(), x0.key()};
                std::unique_ptr<double[]> i1data = I18->get_block(i1hash);
                std::unique_ptr<double[]> i1data_sorted(new double[I18->get_size(i1hash)]);
                sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

                dgemm_("T", "N", a4.size()*c1.size()*a2.size(), x0.size(), x1.size(),
                       1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                       1.0, odata_sorted, a4.size()*c1.size()*a2.size());
              }

              sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), a2.size(), x0.size());
              I17->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task19(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I17 = t[0];
      t2 = t[1];
      I18 = t[2];
    };
    ~Task19() {};
};

template <typename T>
class Task20 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I18;
    std::shared_ptr<Tensor<T> > Gamma2;

    void compute_() {
      for (auto& x0 : active_) {
        for (auto& x1 : active_) {
          std::vector<size_t> ohash = {x1.key(), x0.key()};
          std::unique_ptr<double[]> odata = I18->move_block(ohash);
          {
            std::vector<size_t> i0hash = {x1.key(), x0.key()};
            std::unique_ptr<double[]> i0data = Gamma2->get_block(i0hash);
            sort_indices<0,1,1,1,2,1>(i0data, odata, x1.size(), x0.size());
          }
          I18->put_block(ohash, odata);
        }
      }
    };

  public:
    Task20(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I18 = t[0];
      Gamma2 = t[1];
    };
    ~Task20() {};
};

template <typename T>
class Task21 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I17;
    std::shared_ptr<Tensor<T> > t2;
    std::shared_ptr<Tensor<T> > I21;

    void compute_() {
      for (auto& a2 : virt_) {
        for (auto& c1 : closed_) {
          for (auto& a4 : virt_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), a4.key(), c1.key(), a2.key()};
              std::unique_ptr<double[]> odata = I17->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I17->get_size(ohash)]);
              std::fill_n(odata_sorted.get(), I17->get_size(ohash), 0.0);

              for (auto& x1 : active_) {
                std::vector<size_t> i0hash = {x1.key(), a2.key(), c1.key(), a4.key()};
                std::unique_ptr<double[]> i0data = t2->get_block(i0hash);
                std::unique_ptr<double[]> i0data_sorted(new double[t2->get_size(i0hash)]);
                sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c1.size(), a4.size());

                std::vector<size_t> i1hash = {x1.key(), x0.key()};
                std::unique_ptr<double[]> i1data = I21->get_block(i1hash);
                std::unique_ptr<double[]> i1data_sorted(new double[I21->get_size(i1hash)]);
                sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

                dgemm_("T", "N", a2.size()*c1.size()*a4.size(), x0.size(), x1.size(),
                       1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                       1.0, odata_sorted, a2.size()*c1.size()*a4.size());
              }

              sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), a4.size(), x0.size());
              I17->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task21(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I17 = t[0];
      t2 = t[1];
      I21 = t[2];
    };
    ~Task21() {};
};

template <typename T>
class Task22 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I21;
    std::shared_ptr<Tensor<T> > Gamma2;

    void compute_() {
      for (auto& x0 : active_) {
        for (auto& x1 : active_) {
          std::vector<size_t> ohash = {x1.key(), x0.key()};
          std::unique_ptr<double[]> odata = I21->move_block(ohash);
          {
            std::vector<size_t> i0hash = {x1.key(), x0.key()};
            std::unique_ptr<double[]> i0data = Gamma2->get_block(i0hash);
            sort_indices<0,1,1,1,-1,1>(i0data, odata, x1.size(), x0.size());
          }
          I21->put_block(ohash, odata);
        }
      }
    };

  public:
    Task22(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I21 = t[0];
      Gamma2 = t[1];
    };
    ~Task22() {};
};

template <typename T>
class Task23 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I0;
    std::shared_ptr<Tensor<T> > t2;
    std::shared_ptr<Tensor<T> > I23;

    void compute_() {
      for (auto& a2 : virt_) {
        for (auto& c1 : closed_) {
          for (auto& a3 : virt_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), a3.key(), c1.key(), a2.key()};
              std::unique_ptr<double[]> odata = I0->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I0->get_size(ohash)]);
              std::fill_n(odata_sorted.get(), I0->get_size(ohash), 0.0);

              for (auto& x1 : active_) {
                std::vector<size_t> i0hash = {x1.key(), a3.key(), c1.key(), a2.key()};
                std::unique_ptr<double[]> i0data = t2->get_block(i0hash);
                std::unique_ptr<double[]> i0data_sorted(new double[t2->get_size(i0hash)]);
                sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c1.size(), a2.size());

                std::vector<size_t> i1hash = {x1.key(), x0.key()};
                std::unique_ptr<double[]> i1data = I23->get_block(i1hash);
                std::unique_ptr<double[]> i1data_sorted(new double[I23->get_size(i1hash)]);
                sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

                dgemm_("T", "N", a3.size()*c1.size()*a2.size(), x0.size(), x1.size(),
                       1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                       1.0, odata_sorted, a3.size()*c1.size()*a2.size());
              }

              sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c1.size(), a2.size(), x0.size());
              I0->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task23(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I0 = t[0];
      t2 = t[1];
      I23 = t[2];
    };
    ~Task23() {};
};

template <typename T>
class Task24 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I23;
    std::shared_ptr<Tensor<T> > Gamma2;
    double e0_;

    void compute_() {
      for (auto& x0 : active_) {
        for (auto& x1 : active_) {
          std::vector<size_t> ohash = {x1.key(), x0.key()};
          std::unique_ptr<double[]> odata = I23->move_block(ohash);
          {
            std::vector<size_t> i0hash = {x1.key(), x0.key()};
            std::unique_ptr<double[]> i0data = Gamma2->get_block(i0hash);
            dscal_(x1.size()*x0.size(), -e0_, i0data.get(), 1);
            sort_indices<0,1,1,1,2,1>(i0data, odata, x1.size(), x0.size());
          }
          I23->put_block(ohash, odata);
        }
      }
    };

  public:
    Task24(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i, const double e) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I23 = t[0];
      Gamma2 = t[1];
      e0_ = e; 
    };
    ~Task24() {};
};

template <typename T>
class Task25 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I0;
    std::shared_ptr<Tensor<T> > t2;
    std::shared_ptr<Tensor<T> > I25;

    void compute_() {
      for (auto& a2 : virt_) {
        for (auto& c1 : closed_) {
          for (auto& a3 : virt_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), a3.key(), c1.key(), a2.key()};
              std::unique_ptr<double[]> odata = I0->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I0->get_size(ohash)]);
              std::fill_n(odata_sorted.get(), I0->get_size(ohash), 0.0);

              for (auto& x1 : active_) {
                std::vector<size_t> i0hash = {x1.key(), a2.key(), c1.key(), a3.key()};
                std::unique_ptr<double[]> i0data = t2->get_block(i0hash);
                std::unique_ptr<double[]> i0data_sorted(new double[t2->get_size(i0hash)]);
                sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c1.size(), a3.size());

                std::vector<size_t> i1hash = {x1.key(), x0.key()};
                std::unique_ptr<double[]> i1data = I25->get_block(i1hash);
                std::unique_ptr<double[]> i1data_sorted(new double[I25->get_size(i1hash)]);
                sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

                dgemm_("T", "N", a2.size()*c1.size()*a3.size(), x0.size(), x1.size(),
                       1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                       1.0, odata_sorted, a2.size()*c1.size()*a3.size());
              }

              sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), a3.size(), x0.size());
              I0->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task25(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I0 = t[0];
      t2 = t[1];
      I25 = t[2];
    };
    ~Task25() {};
};

template <typename T>
class Task26 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I25;
    std::shared_ptr<Tensor<T> > Gamma2;
    double e0_;

    void compute_() {
      for (auto& x0 : active_) {
        for (auto& x1 : active_) {
          std::vector<size_t> ohash = {x1.key(), x0.key()};
          std::unique_ptr<double[]> odata = I25->move_block(ohash);
          {
            std::vector<size_t> i0hash = {x1.key(), x0.key()};
            std::unique_ptr<double[]> i0data = Gamma2->get_block(i0hash);
            dscal_(x1.size()*x0.size(), -e0_, i0data.get(), 1);
            sort_indices<0,1,1,1,-1,1>(i0data, odata, x1.size(), x0.size());
          }
          I25->put_block(ohash, odata);
        }
      }
    };

  public:
    Task26(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i, const double e) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I25 = t[0];
      Gamma2 = t[1];
      e0_ = e; 
    };
    ~Task26() {};
};

template <typename T>
class Task27 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I0;
    std::shared_ptr<Tensor<T> > v2;
    std::shared_ptr<Tensor<T> > I27;

    void compute_() {
      for (auto& a2 : virt_) {
        for (auto& c1 : closed_) {
          for (auto& a3 : virt_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), a3.key(), c1.key(), a2.key()};
              std::unique_ptr<double[]> odata = I0->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I0->get_size(ohash)]);
              std::fill_n(odata_sorted.get(), I0->get_size(ohash), 0.0);

              for (auto& x1 : active_) {
                std::vector<size_t> i0hash = {x1.key(), a3.key(), c1.key(), a2.key()};
                std::unique_ptr<double[]> i0data = v2->get_block(i0hash);
                std::unique_ptr<double[]> i0data_sorted(new double[v2->get_size(i0hash)]);
                sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c1.size(), a2.size());

                std::vector<size_t> i1hash = {x1.key(), x0.key()};
                std::unique_ptr<double[]> i1data = I27->get_block(i1hash);
                std::unique_ptr<double[]> i1data_sorted(new double[I27->get_size(i1hash)]);
                sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

                dgemm_("T", "N", a3.size()*c1.size()*a2.size(), x0.size(), x1.size(),
                       1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                       1.0, odata_sorted, a3.size()*c1.size()*a2.size());
              }

              sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c1.size(), a2.size(), x0.size());
              I0->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task27(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I0 = t[0];
      v2 = t[1];
      I27 = t[2];
    };
    ~Task27() {};
};

template <typename T>
class Task28 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I27;
    std::shared_ptr<Tensor<T> > Gamma2;

    void compute_() {
      for (auto& x0 : active_) {
        for (auto& x1 : active_) {
          std::vector<size_t> ohash = {x1.key(), x0.key()};
          std::unique_ptr<double[]> odata = I27->move_block(ohash);
          {
            std::vector<size_t> i0hash = {x1.key(), x0.key()};
            std::unique_ptr<double[]> i0data = Gamma2->get_block(i0hash);
            sort_indices<0,1,1,1,4,1>(i0data, odata, x1.size(), x0.size());
          }
          I27->put_block(ohash, odata);
        }
      }
    };

  public:
    Task28(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I27 = t[0];
      Gamma2 = t[1];
    };
    ~Task28() {};
};

template <typename T>
class Task29 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I0;
    std::shared_ptr<Tensor<T> > v2;
    std::shared_ptr<Tensor<T> > I29;

    void compute_() {
      for (auto& a2 : virt_) {
        for (auto& c1 : closed_) {
          for (auto& a3 : virt_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), a3.key(), c1.key(), a2.key()};
              std::unique_ptr<double[]> odata = I0->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I0->get_size(ohash)]);
              std::fill_n(odata_sorted.get(), I0->get_size(ohash), 0.0);

              for (auto& x1 : active_) {
                std::vector<size_t> i0hash = {x1.key(), a2.key(), c1.key(), a3.key()};
                std::unique_ptr<double[]> i0data = v2->get_block(i0hash);
                std::unique_ptr<double[]> i0data_sorted(new double[v2->get_size(i0hash)]);
                sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c1.size(), a3.size());

                std::vector<size_t> i1hash = {x1.key(), x0.key()};
                std::unique_ptr<double[]> i1data = I29->get_block(i1hash);
                std::unique_ptr<double[]> i1data_sorted(new double[I29->get_size(i1hash)]);
                sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

                dgemm_("T", "N", a2.size()*c1.size()*a3.size(), x0.size(), x1.size(),
                       1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                       1.0, odata_sorted, a2.size()*c1.size()*a3.size());
              }

              sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), a3.size(), x0.size());
              I0->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task29(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I0 = t[0];
      v2 = t[1];
      I29 = t[2];
    };
    ~Task29() {};
};

template <typename T>
class Task30 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I29;
    std::shared_ptr<Tensor<T> > Gamma2;

    void compute_() {
      for (auto& x0 : active_) {
        for (auto& x1 : active_) {
          std::vector<size_t> ohash = {x1.key(), x0.key()};
          std::unique_ptr<double[]> odata = I29->move_block(ohash);
          {
            std::vector<size_t> i0hash = {x1.key(), x0.key()};
            std::unique_ptr<double[]> i0data = Gamma2->get_block(i0hash);
            sort_indices<0,1,1,1,-2,1>(i0data, odata, x1.size(), x0.size());
          }
          I29->put_block(ohash, odata);
        }
      }
    };

  public:
    Task30(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I29 = t[0];
      Gamma2 = t[1];
    };
    ~Task30() {};
};

template <typename T>
class Task31 : public EnergyTask<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I30;
    std::shared_ptr<Tensor<T> > t2;
    std::shared_ptr<Tensor<T> > I31;

    void compute_() {
      this->energy_ = 0.0;

      for (auto& a3 : virt_) {
        for (auto& x0 : active_) {
          for (auto& a2 : virt_) {
            for (auto& c1 : closed_) {
              std::vector<size_t> i0hash = {x0.key(), a3.key(), c1.key(), a2.key()};
              std::unique_ptr<double[]> i0data = t2->get_block(i0hash);
              std::unique_ptr<double[]> i0data_sorted(new double[t2->get_size(i0hash)]);
              sort_indices<1,0,3,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), c1.size(), a2.size());

              std::vector<size_t> i1hash = {x0.key(), a3.key(), c1.key(), a2.key()};
              std::unique_ptr<double[]> i1data = I31->get_block(i1hash);
              std::unique_ptr<double[]> i1data_sorted(new double[I31->get_size(i1hash)]);
              sort_indices<1,0,3,2,0,1,1,1>(i1data, i1data_sorted, x0.size(), a3.size(), c1.size(), a2.size());

              this->energy_ += ddot_(x0.size()*a3.size()*c1.size()*a2.size(), i0data_sorted, 1, i1data_sorted, 1);
            }
          }
        }
      }

    };

  public:
    Task31(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : EnergyTask<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I30 = t[0];
      t2 = t[1];
      I31 = t[2];
    };
    ~Task31() {};
};

template <typename T>
class Task32 : public EnergyTask<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I31;
    std::shared_ptr<Tensor<T> > v2;
    std::shared_ptr<Tensor<T> > I32;

    void compute_() {
      this->energy_ = 0.0;
      for (auto& a2 : virt_) {
        for (auto& c1 : closed_) {
          for (auto& a3 : virt_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), a3.key(), c1.key(), a2.key()};
              std::unique_ptr<double[]> odata = I31->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I31->get_size(ohash)]);
              std::fill_n(odata_sorted.get(), I31->get_size(ohash), 0.0);

              for (auto& x1 : active_) {
                std::vector<size_t> i0hash = {x1.key(), a3.key(), c1.key(), a2.key()};
                std::unique_ptr<double[]> i0data = v2->get_block(i0hash);
                std::unique_ptr<double[]> i0data_sorted(new double[v2->get_size(i0hash)]);
                sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c1.size(), a2.size());

                std::vector<size_t> i1hash = {x1.key(), x0.key()};
                std::unique_ptr<double[]> i1data = I32->get_block(i1hash);
                std::unique_ptr<double[]> i1data_sorted(new double[I32->get_size(i1hash)]);
                sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

                dgemm_("T", "N", a3.size()*c1.size()*a2.size(), x0.size(), x1.size(),
                       1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                       1.0, odata_sorted, a3.size()*c1.size()*a2.size());
              }

              sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c1.size(), a2.size(), x0.size());
              I31->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task32(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : EnergyTask<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I31 = t[0];
      v2 = t[1];
      I32 = t[2];
    };
    ~Task32() {};
};

template <typename T>
class Task33 : public EnergyTask<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I32;
    std::shared_ptr<Tensor<T> > Gamma2;

    void compute_() {
      this->energy_ = 0.0;
      for (auto& x0 : active_) {
        for (auto& x1 : active_) {
          std::vector<size_t> ohash = {x1.key(), x0.key()};
          std::unique_ptr<double[]> odata = I32->move_block(ohash);
          {
            std::vector<size_t> i0hash = {x1.key(), x0.key()};
            std::unique_ptr<double[]> i0data = Gamma2->get_block(i0hash);
            sort_indices<0,1,1,1,4,1>(i0data, odata, x1.size(), x0.size());
          }
          I32->put_block(ohash, odata);
        }
      }
    };

  public:
    Task33(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : EnergyTask<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I32 = t[0];
      Gamma2 = t[1];
    };
    ~Task33() {};
};

template <typename T>
class Task34 : public EnergyTask<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I31;
    std::shared_ptr<Tensor<T> > v2;
    std::shared_ptr<Tensor<T> > I35;

    void compute_() {
      this->energy_ = 0.0;
      for (auto& a2 : virt_) {
        for (auto& c1 : closed_) {
          for (auto& a3 : virt_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), a3.key(), c1.key(), a2.key()};
              std::unique_ptr<double[]> odata = I31->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I31->get_size(ohash)]);
              std::fill_n(odata_sorted.get(), I31->get_size(ohash), 0.0);

              for (auto& x1 : active_) {
                std::vector<size_t> i0hash = {x1.key(), a2.key(), c1.key(), a3.key()};
                std::unique_ptr<double[]> i0data = v2->get_block(i0hash);
                std::unique_ptr<double[]> i0data_sorted(new double[v2->get_size(i0hash)]);
                sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c1.size(), a3.size());

                std::vector<size_t> i1hash = {x1.key(), x0.key()};
                std::unique_ptr<double[]> i1data = I35->get_block(i1hash);
                std::unique_ptr<double[]> i1data_sorted(new double[I35->get_size(i1hash)]);
                sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

                dgemm_("T", "N", a2.size()*c1.size()*a3.size(), x0.size(), x1.size(),
                       1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                       1.0, odata_sorted, a2.size()*c1.size()*a3.size());
              }

              sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), a3.size(), x0.size());
              I31->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task34(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : EnergyTask<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I31 = t[0];
      v2 = t[1];
      I35 = t[2];
    };
    ~Task34() {};
};

template <typename T>
class Task35 : public EnergyTask<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I35;
    std::shared_ptr<Tensor<T> > Gamma2;

    void compute_() {
      this->energy_ = 0.0;
      for (auto& x0 : active_) {
        for (auto& x1 : active_) {
          std::vector<size_t> ohash = {x1.key(), x0.key()};
          std::unique_ptr<double[]> odata = I35->move_block(ohash);
          {
            std::vector<size_t> i0hash = {x1.key(), x0.key()};
            std::unique_ptr<double[]> i0data = Gamma2->get_block(i0hash);
            sort_indices<0,1,1,1,-2,1>(i0data, odata, x1.size(), x0.size());
          }
          I35->put_block(ohash, odata);
        }
      }
    };

  public:
    Task35(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : EnergyTask<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I35 = t[0];
      Gamma2 = t[1];
    };
    ~Task35() {};
};

template <typename T>
class Task36 : public EnergyTask<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I31;
    std::shared_ptr<Tensor<T> > r;
    std::shared_ptr<Tensor<T> > I38;

    void compute_() {
      this->energy_ = 0.0;
      for (auto& a2 : virt_) {
        for (auto& c1 : closed_) {
          for (auto& a3 : virt_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), a3.key(), c1.key(), a2.key()};
              std::unique_ptr<double[]> odata = I31->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I31->get_size(ohash)]);
              std::fill_n(odata_sorted.get(), I31->get_size(ohash), 0.0);

              for (auto& x1 : active_) {
                std::vector<size_t> i0hash = {x1.key(), a3.key(), c1.key(), a2.key()};
                std::unique_ptr<double[]> i0data = r->get_block(i0hash);
                std::unique_ptr<double[]> i0data_sorted(new double[r->get_size(i0hash)]);
                sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c1.size(), a2.size());

                std::vector<size_t> i1hash = {x1.key(), x0.key()};
                std::unique_ptr<double[]> i1data = I38->get_block(i1hash);
                std::unique_ptr<double[]> i1data_sorted(new double[I38->get_size(i1hash)]);
                sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

                dgemm_("T", "N", a3.size()*c1.size()*a2.size(), x0.size(), x1.size(),
                       1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                       1.0, odata_sorted, a3.size()*c1.size()*a2.size());
              }

              sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c1.size(), a2.size(), x0.size());
              I31->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task36(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : EnergyTask<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I31 = t[0];
      r = t[1];
      I38 = t[2];
    };
    ~Task36() {};
};

template <typename T>
class Task37 : public EnergyTask<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I38;
    std::shared_ptr<Tensor<T> > Gamma2;

    void compute_() {
      this->energy_ = 0.0;
      for (auto& x0 : active_) {
        for (auto& x1 : active_) {
          std::vector<size_t> ohash = {x1.key(), x0.key()};
          std::unique_ptr<double[]> odata = I38->move_block(ohash);
          {
            std::vector<size_t> i0hash = {x1.key(), x0.key()};
            std::unique_ptr<double[]> i0data = Gamma2->get_block(i0hash);
            sort_indices<0,1,1,1,2,1>(i0data, odata, x1.size(), x0.size());
          }
          I38->put_block(ohash, odata);
        }
      }
    };

  public:
    Task37(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : EnergyTask<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I38 = t[0];
      Gamma2 = t[1];
    };
    ~Task37() {};
};

template <typename T>
class Task38 : public EnergyTask<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I31;
    std::shared_ptr<Tensor<T> > r;
    std::shared_ptr<Tensor<T> > I41;

    void compute_() {
      this->energy_ = 0.0;
      for (auto& a2 : virt_) {
        for (auto& c1 : closed_) {
          for (auto& a3 : virt_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), a3.key(), c1.key(), a2.key()};
              std::unique_ptr<double[]> odata = I31->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I31->get_size(ohash)]);
              std::fill_n(odata_sorted.get(), I31->get_size(ohash), 0.0);

              for (auto& x1 : active_) {
                std::vector<size_t> i0hash = {x1.key(), a2.key(), c1.key(), a3.key()};
                std::unique_ptr<double[]> i0data = r->get_block(i0hash);
                std::unique_ptr<double[]> i0data_sorted(new double[r->get_size(i0hash)]);
                sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c1.size(), a3.size());

                std::vector<size_t> i1hash = {x1.key(), x0.key()};
                std::unique_ptr<double[]> i1data = I41->get_block(i1hash);
                std::unique_ptr<double[]> i1data_sorted(new double[I41->get_size(i1hash)]);
                sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

                dgemm_("T", "N", a2.size()*c1.size()*a3.size(), x0.size(), x1.size(),
                       1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                       1.0, odata_sorted, a2.size()*c1.size()*a3.size());
              }

              sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), a3.size(), x0.size());
              I31->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task38(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : EnergyTask<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I31 = t[0];
      r = t[1];
      I41 = t[2];
    };
    ~Task38() {};
};

template <typename T>
class Task39 : public EnergyTask<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I41;
    std::shared_ptr<Tensor<T> > Gamma2;

    void compute_() {
      this->energy_ = 0.0;
      for (auto& x0 : active_) {
        for (auto& x1 : active_) {
          std::vector<size_t> ohash = {x1.key(), x0.key()};
          std::unique_ptr<double[]> odata = I41->move_block(ohash);
          {
            std::vector<size_t> i0hash = {x1.key(), x0.key()};
            std::unique_ptr<double[]> i0data = Gamma2->get_block(i0hash);
            sort_indices<0,1,1,1,-1,1>(i0data, odata, x1.size(), x0.size());
          }
          I41->put_block(ohash, odata);
        }
      }
    };

  public:
    Task39(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : EnergyTask<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I41 = t[0];
      Gamma2 = t[1];
    };
    ~Task39() {};
};


}
}
}
#endif

