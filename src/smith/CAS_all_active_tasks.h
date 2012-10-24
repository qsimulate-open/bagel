//
// BAGEL - Parallel electron correlation program.
// Filename: CAS_all_active_tasks.h
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


#ifndef __SRC_SMITH_CAS_all_active_TASKS_H 
#define __SRC_SMITH_CAS_all_active_TASKS_H 

#include <memory>
#include <algorithm>
#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <vector>

namespace bagel {
namespace SMITH {
namespace CAS_all_active{

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
    std::shared_ptr<Tensor<T> > rdm3;
    std::shared_ptr<Tensor<T> > f1;

    void compute_() {
      for (auto& x1 : active_) {
        for (auto& x4 : active_) {
          for (auto& x0 : active_) {
            for (auto& x5 : active_) {
              std::vector<size_t> ohash = {x5.key(), x0.key(), x4.key(), x1.key()};
              std::unique_ptr<double[]> odata = Gamma0->move_block(ohash);
              // associated with merged
              for (auto& x2 : active_) {
                for (auto& x3 : active_) {
                  std::vector<size_t> fhash = {x3.key(), x2.key()};
                  std::unique_ptr<double[]> fdata = f1->get_block(fhash);
                  {
                    std::vector<size_t> i0hash = {x5.key(), x0.key(), x4.key(), x1.key(), x3.key(), x2.key()};
                    std::unique_ptr<double[]> data = rdm3->get_block(i0hash);
                    for (int i2 = 0; i2 != x2.size(); ++i2) {
                      for (int i3 = 0; i3 != x3.size(); ++i3) {
                        for (int i1 = 0; i1 != x1.size(); ++i1) {
                          for (int i4 = 0; i4 != x4.size(); ++i4) {
                            for (int i0 = 0; i0 != x0.size(); ++i0) {
                              for (int i5 = 0; i5 != x5.size(); ++i5) {
                                odata[i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i1)))]
                                  += (1.0) * data[i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i1+x1.size()*(i3+x3.size()*(i2)))))] * fdata[i3+x3.size()*(i2)];
                              }
                            }
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
        }
      }
    };  


  public:
    Task1(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      Gamma0  = t[0];
      rdm3    = t[1];
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
    std::shared_ptr<Tensor<T> > rdm2;

    void compute_() {
      for (auto& x1 : active_) {
        for (auto& x2 : active_) {
          for (auto& x0 : active_) {
            for (auto& x3 : active_) {
              std::vector<size_t> ohash = {x3.key(), x0.key(), x2.key(), x1.key()};
              std::unique_ptr<double[]> odata = Gamma2->move_block(ohash);
              {
                std::vector<size_t> i0hash = {x3.key(), x0.key(), x2.key(), x1.key()};
                std::unique_ptr<double[]> data = rdm2->get_block(i0hash);
                sort_indices<0,1,2,3,1,1,-1,1>(data, odata, x3.size(), x0.size(), x2.size(), x1.size());
              }
              Gamma2->put_block(ohash, odata);
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
      Gamma2  = t[0];
      rdm2    = t[1];
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
      for (auto& a1 : virt_) {
        for (auto& x0 : active_) {
          for (auto& a2 : virt_) {
            for (auto& x1 : active_) {
              std::vector<size_t> ohash = {x1.key(), a2.key(), x0.key(), a1.key()};
              std::unique_ptr<double[]> odata = r->move_block(ohash);
              {
                std::vector<size_t> i0hash = {x0.key(), x1.key(), a1.key(), a2.key()};
                std::unique_ptr<double[]> i0data = I0->get_block(i0hash);
                sort_indices<1,3,0,2,1,1,1,1>(i0data, odata, x0.size(), x1.size(), a1.size(), a2.size());
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
        for (auto& a1 : virt_) {
          for (auto& x1 : active_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), x1.key(), a1.key(), a2.key()};
              std::unique_ptr<double[]> odata = I0->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I0->get_size(ohash)]);
              std::fill(odata_sorted.get(), odata_sorted.get()+I0->get_size(ohash), 0.0);

              for (auto& x5 : active_) {
                for (auto& x4 : active_) {
                  std::vector<size_t> i0hash = {x5.key(), a1.key(), x4.key(), a2.key()};
                  std::unique_ptr<double[]> i0data = t2->get_block(i0hash);
                  std::unique_ptr<double[]> i0data_sorted(new double[t2->get_size(i0hash)]);
                  sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a2.size());

                  std::vector<size_t> i1hash = {x5.key(), x0.key(), x4.key(), x1.key()};
                  std::unique_ptr<double[]> i1data = I1->get_block(i1hash);
                  std::unique_ptr<double[]> i1data_sorted(new double[I1->get_size(i1hash)]);
                  sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), x0.size(), x4.size(), x1.size());

                  dgemm_("T", "N", a1.size()*a2.size(), x0.size()*x1.size(), x5.size()*x4.size(),
                         1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
                         1.0, odata_sorted, a1.size()*a2.size());
                }
              }

              sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
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
      for (auto& x1 : active_) {
        for (auto& x4 : active_) {
          for (auto& x0 : active_) {
            for (auto& x5 : active_) {
              std::vector<size_t> ohash = {x5.key(), x0.key(), x4.key(), x1.key()};
              std::unique_ptr<double[]> odata = I1->move_block(ohash);
              {
                std::vector<size_t> i0hash = {x5.key(), x0.key(), x4.key(), x1.key()};
                std::unique_ptr<double[]> i0data = Gamma0->get_block(i0hash);
                sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x5.size(), x0.size(), x4.size(), x1.size());
              }
              I1->put_block(ohash, odata);
            }
          }
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
    std::shared_ptr<Tensor<T> > I6;

    void compute_() {
      for (auto& a2 : virt_) {
        for (auto& a1 : virt_) {
          for (auto& x1 : active_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), x1.key(), a1.key(), a2.key()};
              std::unique_ptr<double[]> odata = I0->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I0->get_size(ohash)]);
              std::fill(odata_sorted.get(), odata_sorted.get()+I0->get_size(ohash), 0.0);

              for (auto& x3 : active_) {
                for (auto& x2 : active_) {
                  std::vector<size_t> i0hash = {x3.key(), a1.key(), x2.key(), a2.key()};
                  std::unique_ptr<double[]> i0data = t2->get_block(i0hash);
                  std::unique_ptr<double[]> i0data_sorted(new double[t2->get_size(i0hash)]);
                  sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

                  std::vector<size_t> i1hash = {x3.key(), x0.key(), x2.key(), x1.key()};
                  std::unique_ptr<double[]> i1data = I6->get_block(i1hash);
                  std::unique_ptr<double[]> i1data_sorted(new double[I6->get_size(i1hash)]);
                  sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

                  dgemm_("T", "N", a1.size()*a2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
                         1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
                         1.0, odata_sorted, a1.size()*a2.size());
                }
              }

              sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
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
      I6 = t[2];
    };
    ~Task6() {};
};

template <typename T>
class Task7 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I6;
    std::shared_ptr<Tensor<T> > Gamma2;
    double e0_;

    void compute_() {
      for (auto& x1 : active_) {
        for (auto& x2 : active_) {
          for (auto& x0 : active_) {
            for (auto& x3 : active_) {
              std::vector<size_t> ohash = {x3.key(), x0.key(), x2.key(), x1.key()};
              std::unique_ptr<double[]> odata = I6->move_block(ohash);
              {
                std::vector<size_t> i0hash = {x3.key(), x0.key(), x2.key(), x1.key()};
                std::unique_ptr<double[]> i0data = Gamma2->get_block(i0hash);
                dscal_(x3.size()*x0.size()*x2.size()*x1.size(), -e0_, i0data.get(), 1);
                sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
              }
              I6->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task7(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i, const double e) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I6 = t[0];
      Gamma2 = t[1];
      e0_ = e; 
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
    std::shared_ptr<Tensor<T> > v2;
    std::shared_ptr<Tensor<T> > I8;

    void compute_() {
      for (auto& a2 : virt_) {
        for (auto& a1 : virt_) {
          for (auto& x1 : active_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), x1.key(), a1.key(), a2.key()};
              std::unique_ptr<double[]> odata = I0->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I0->get_size(ohash)]);
              std::fill(odata_sorted.get(), odata_sorted.get()+I0->get_size(ohash), 0.0);

              for (auto& x3 : active_) {
                for (auto& x2 : active_) {
                  std::vector<size_t> i0hash = {x3.key(), a1.key(), x2.key(), a2.key()};
                  std::unique_ptr<double[]> i0data = v2->get_block(i0hash);
                  std::unique_ptr<double[]> i0data_sorted(new double[v2->get_size(i0hash)]);
                  sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

                  std::vector<size_t> i1hash = {x3.key(), x0.key(), x2.key(), x1.key()};
                  std::unique_ptr<double[]> i1data = I8->get_block(i1hash);
                  std::unique_ptr<double[]> i1data_sorted(new double[I8->get_size(i1hash)]);
                  sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

                  dgemm_("T", "N", a1.size()*a2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
                         1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
                         1.0, odata_sorted, a1.size()*a2.size());
                }
              }

              sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
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
      v2 = t[1];
      I8 = t[2];
    };
    ~Task8() {};
};

template <typename T>
class Task9 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I8;
    std::shared_ptr<Tensor<T> > Gamma2;

    void compute_() {
      for (auto& x1 : active_) {
        for (auto& x2 : active_) {
          for (auto& x0 : active_) {
            for (auto& x3 : active_) {
              std::vector<size_t> ohash = {x3.key(), x0.key(), x2.key(), x1.key()};
              std::unique_ptr<double[]> odata = I8->move_block(ohash);
              {
                std::vector<size_t> i0hash = {x3.key(), x0.key(), x2.key(), x1.key()};
                std::unique_ptr<double[]> i0data = Gamma2->get_block(i0hash);
                sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
              }
              I8->put_block(ohash, odata);
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
      I8 = t[0];
      Gamma2 = t[1];
    };
    ~Task9() {};
};

template <typename T>
class Task10 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > r;
    std::shared_ptr<Tensor<T> > I2;

    void compute_() {
      for (auto& a1 : virt_) {
        for (auto& x0 : active_) {
          for (auto& a2 : virt_) {
            for (auto& x1 : active_) {
              std::vector<size_t> ohash = {x1.key(), a2.key(), x0.key(), a1.key()};
              std::unique_ptr<double[]> odata = r->move_block(ohash);
              {
                std::vector<size_t> i0hash = {x0.key(), x1.key(), a1.key(), a2.key()};
                std::unique_ptr<double[]> i0data = I2->get_block(i0hash);
                sort_indices<1,3,0,2,1,1,1,1>(i0data, odata, x0.size(), x1.size(), a1.size(), a2.size());
              }
              {
                std::vector<size_t> i0hash = {x1.key(), x0.key(), a2.key(), a1.key()};
                std::unique_ptr<double[]> i0data = I2->get_block(i0hash);
                sort_indices<0,2,1,3,1,1,1,1>(i0data, odata, x1.size(), x0.size(), a2.size(), a1.size());
              }
              r->put_block(ohash, odata);
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
      r = t[0];
      I2 = t[1];
    };
    ~Task10() {};
};

template <typename T>
class Task11 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I2;
    std::shared_ptr<Tensor<T> > f1;
    std::shared_ptr<Tensor<T> > I3;

    void compute_() {
      for (auto& a2 : virt_) {
        for (auto& a1 : virt_) {
          for (auto& x1 : active_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), x1.key(), a1.key(), a2.key()};
              std::unique_ptr<double[]> odata = I2->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I2->get_size(ohash)]);
              std::fill(odata_sorted.get(), odata_sorted.get()+I2->get_size(ohash), 0.0);

              for (auto& a3 : virt_) {
                std::vector<size_t> i0hash = {a3.key(), a2.key()};
                std::unique_ptr<double[]> i0data = f1->get_block(i0hash);
                std::unique_ptr<double[]> i0data_sorted(new double[f1->get_size(i0hash)]);
                sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), a2.size());

                std::vector<size_t> i1hash = {x0.key(), x1.key(), a1.key(), a3.key()};
                std::unique_ptr<double[]> i1data = I3->get_block(i1hash);
                std::unique_ptr<double[]> i1data_sorted(new double[I3->get_size(i1hash)]);
                sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), a1.size(), a3.size());

                dgemm_("T", "N", a2.size(), x0.size()*x1.size()*a1.size(), a3.size(),
                       1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
                       1.0, odata_sorted, a2.size());
              }

              sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), x0.size(), x1.size(), a1.size());
              I2->put_block(ohash, odata);
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
      I2 = t[0];
      f1 = t[1];
      I3 = t[2];
    };
    ~Task11() {};
};

template <typename T>
class Task12 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I3;
    std::shared_ptr<Tensor<T> > t2;
    std::shared_ptr<Tensor<T> > I4;

    void compute_() {
      for (auto& a3 : virt_) {
        for (auto& a1 : virt_) {
          for (auto& x1 : active_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), x1.key(), a1.key(), a3.key()};
              std::unique_ptr<double[]> odata = I3->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I3->get_size(ohash)]);
              std::fill(odata_sorted.get(), odata_sorted.get()+I3->get_size(ohash), 0.0);

              for (auto& x3 : active_) {
                for (auto& x2 : active_) {
                  std::vector<size_t> i0hash = {x3.key(), a1.key(), x2.key(), a3.key()};
                  std::unique_ptr<double[]> i0data = t2->get_block(i0hash);
                  std::unique_ptr<double[]> i0data_sorted(new double[t2->get_size(i0hash)]);
                  sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a3.size());

                  std::vector<size_t> i1hash = {x3.key(), x0.key(), x2.key(), x1.key()};
                  std::unique_ptr<double[]> i1data = I4->get_block(i1hash);
                  std::unique_ptr<double[]> i1data_sorted(new double[I4->get_size(i1hash)]);
                  sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

                  dgemm_("T", "N", a1.size()*a3.size(), x0.size()*x1.size(), x3.size()*x2.size(),
                         1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
                         1.0, odata_sorted, a1.size()*a3.size());
                }
              }

              sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), x0.size(), x1.size());
              I3->put_block(ohash, odata);
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
      I3 = t[0];
      t2 = t[1];
      I4 = t[2];
    };
    ~Task12() {};
};

template <typename T>
class Task13 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I4;
    std::shared_ptr<Tensor<T> > Gamma2;

    void compute_() {
      for (auto& x1 : active_) {
        for (auto& x2 : active_) {
          for (auto& x0 : active_) {
            for (auto& x3 : active_) {
              std::vector<size_t> ohash = {x3.key(), x0.key(), x2.key(), x1.key()};
              std::unique_ptr<double[]> odata = I4->move_block(ohash);
              {
                std::vector<size_t> i0hash = {x3.key(), x0.key(), x2.key(), x1.key()};
                std::unique_ptr<double[]> i0data = Gamma2->get_block(i0hash);
                sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
              }
              I4->put_block(ohash, odata);
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
      I4 = t[0];
      Gamma2 = t[1];
    };
    ~Task13() {};
};

template <typename T>
class Task14 : public EnergyTask<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I9;
    std::shared_ptr<Tensor<T> > t2;
    std::shared_ptr<Tensor<T> > I10;

    void compute_() {
      this->energy_ = 0.0;

      for (auto& a2 : virt_) {
        for (auto& x1 : active_) {
          for (auto& a1 : virt_) {
            for (auto& x0 : active_) {
              std::vector<size_t> i0hash = {x1.key(), a2.key(), x0.key(), a1.key()};
              std::unique_ptr<double[]> i0data = t2->get_block(i0hash);
              std::unique_ptr<double[]> i0data_sorted(new double[t2->get_size(i0hash)]);
              sort_indices<1,0,3,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), x0.size(), a1.size());

              std::vector<size_t> i1hash = {x0.key(), x1.key(), a1.key(), a2.key()};
              std::unique_ptr<double[]> i1data = I10->get_block(i1hash);
              std::unique_ptr<double[]> i1data_sorted(new double[I10->get_size(i1hash)]);
              sort_indices<3,1,2,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), a1.size(), a2.size());

              this->energy_ += ddot_(x0.size()*x1.size()*a1.size()*a2.size(), i0data_sorted, 1, i1data_sorted, 1);
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
      I9 = t[0];
      t2 = t[1];
      I10 = t[2];
    };
    ~Task14() {};
};

template <typename T>
class Task15 : public EnergyTask<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I10;
    std::shared_ptr<Tensor<T> > v2;
    std::shared_ptr<Tensor<T> > I11;

    void compute_() {
      this->energy_ = 0.0;
      for (auto& a2 : virt_) {
        for (auto& a1 : virt_) {
          for (auto& x1 : active_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), x1.key(), a1.key(), a2.key()};
              std::unique_ptr<double[]> odata = I10->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I10->get_size(ohash)]);
              std::fill(odata_sorted.get(), odata_sorted.get()+I10->get_size(ohash), 0.0);

              for (auto& x3 : active_) {
                for (auto& x2 : active_) {
                  std::vector<size_t> i0hash = {x3.key(), a1.key(), x2.key(), a2.key()};
                  std::unique_ptr<double[]> i0data = v2->get_block(i0hash);
                  std::unique_ptr<double[]> i0data_sorted(new double[v2->get_size(i0hash)]);
                  sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

                  std::vector<size_t> i1hash = {x3.key(), x0.key(), x2.key(), x1.key()};
                  std::unique_ptr<double[]> i1data = I11->get_block(i1hash);
                  std::unique_ptr<double[]> i1data_sorted(new double[I11->get_size(i1hash)]);
                  sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

                  dgemm_("T", "N", a1.size()*a2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
                         1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
                         1.0, odata_sorted, a1.size()*a2.size());
                }
              }

              sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
              I10->put_block(ohash, odata);
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
      I10 = t[0];
      v2 = t[1];
      I11 = t[2];
    };
    ~Task15() {};
};

template <typename T>
class Task16 : public EnergyTask<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I11;
    std::shared_ptr<Tensor<T> > Gamma2;

    void compute_() {
      this->energy_ = 0.0;
      for (auto& x1 : active_) {
        for (auto& x2 : active_) {
          for (auto& x0 : active_) {
            for (auto& x3 : active_) {
              std::vector<size_t> ohash = {x3.key(), x0.key(), x2.key(), x1.key()};
              std::unique_ptr<double[]> odata = I11->move_block(ohash);
              {
                std::vector<size_t> i0hash = {x3.key(), x0.key(), x2.key(), x1.key()};
                std::unique_ptr<double[]> i0data = Gamma2->get_block(i0hash);
                sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
              }
              I11->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task16(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : EnergyTask<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I11 = t[0];
      Gamma2 = t[1];
    };
    ~Task16() {};
};

template <typename T>
class Task17 : public EnergyTask<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I10;
    std::shared_ptr<Tensor<T> > r;
    std::shared_ptr<Tensor<T> > I14;

    void compute_() {
      this->energy_ = 0.0;
      for (auto& a2 : virt_) {
        for (auto& a1 : virt_) {
          for (auto& x1 : active_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), x1.key(), a1.key(), a2.key()};
              std::unique_ptr<double[]> odata = I10->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I10->get_size(ohash)]);
              std::fill(odata_sorted.get(), odata_sorted.get()+I10->get_size(ohash), 0.0);

              for (auto& x3 : active_) {
                for (auto& x2 : active_) {
                  std::vector<size_t> i0hash = {x3.key(), a1.key(), x2.key(), a2.key()};
                  std::unique_ptr<double[]> i0data = r->get_block(i0hash);
                  std::unique_ptr<double[]> i0data_sorted(new double[r->get_size(i0hash)]);
                  sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

                  std::vector<size_t> i1hash = {x3.key(), x0.key(), x2.key(), x1.key()};
                  std::unique_ptr<double[]> i1data = I14->get_block(i1hash);
                  std::unique_ptr<double[]> i1data_sorted(new double[I14->get_size(i1hash)]);
                  sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

                  dgemm_("T", "N", a1.size()*a2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
                         1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
                         1.0, odata_sorted, a1.size()*a2.size());
                }
              }

              sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
              I10->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task17(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : EnergyTask<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I10 = t[0];
      r = t[1];
      I14 = t[2];
    };
    ~Task17() {};
};

template <typename T>
class Task18 : public EnergyTask<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I14;
    std::shared_ptr<Tensor<T> > Gamma2;

    void compute_() {
      this->energy_ = 0.0;
      for (auto& x1 : active_) {
        for (auto& x2 : active_) {
          for (auto& x0 : active_) {
            for (auto& x3 : active_) {
              std::vector<size_t> ohash = {x3.key(), x0.key(), x2.key(), x1.key()};
              std::unique_ptr<double[]> odata = I14->move_block(ohash);
              {
                std::vector<size_t> i0hash = {x3.key(), x0.key(), x2.key(), x1.key()};
                std::unique_ptr<double[]> i0data = Gamma2->get_block(i0hash);
                sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
              }
              I14->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task18(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : EnergyTask<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I14 = t[0];
      Gamma2 = t[1];
    };
    ~Task18() {};
};


}
}
}
#endif

