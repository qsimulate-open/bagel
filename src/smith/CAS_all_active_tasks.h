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
                    for (int i1 = 0; i1 != x1.size(); ++i1) {
                      for (int i4 = 0; i4 != x4.size(); ++i4) {
                        for (int i0 = 0; i0 != x0.size(); ++i0) {
                          for (int i5 = 0; i5 != x5.size(); ++i5) {
                            for (int i2 = 0; i2 != x2.size(); ++i2) {
                              for (int i3 = 0; i3 != x3.size(); ++i3) {
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
    std::shared_ptr<Tensor<T> > Gamma6;
    std::shared_ptr<Tensor<T> > rdm2;

    void compute_() {
      for (auto& x1 : active_) {
        for (auto& x2 : active_) {
          for (auto& x0 : active_) {
            for (auto& x3 : active_) {
              std::vector<size_t> ohash = {x3.key(), x0.key(), x2.key(), x1.key()};
              std::unique_ptr<double[]> odata = Gamma6->move_block(ohash);
              {
                std::vector<size_t> i0hash = {x3.key(), x0.key(), x2.key(), x1.key()};
                std::unique_ptr<double[]> data = rdm2->get_block(i0hash);
                sort_indices<0,1,2,3,1,1,-1,1>(data, odata, x3.size(), x0.size(), x2.size(), x1.size());
              }
              Gamma6->put_block(ohash, odata);
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
      Gamma6  = t[0];
      rdm2    = t[1];
    };
    ~Task2() {};
};

template <typename T>
class Task3 : public Task<T> {  // associated with gamma
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > Gamma5;
    std::shared_ptr<Tensor<T> > rdm2;
    std::shared_ptr<Tensor<T> > rdm3;

    void compute_() {
      for (auto& x1 : active_) {
        for (auto& x2 : active_) {
          for (auto& x0 : active_) {
            for (auto& x3 : active_) {
              for (auto& x4 : active_) {
                for (auto& x5 : active_) {
                  std::vector<size_t> ohash = {x5.key(), x4.key(), x3.key(), x0.key(), x2.key(), x1.key()};
                  std::unique_ptr<double[]> odata = Gamma5->move_block(ohash);
                  {
                    if (x2 == x4) {
                      std::vector<size_t> i0hash = {x5.key(), x1.key(), x3.key(), x0.key()};
                      std::unique_ptr<double[]> data = rdm2->get_block(i0hash);
                      for (int i1 = 0; i1 != x1.size(); ++i1) {
                        for (int i0 = 0; i0 != x0.size(); ++i0) {
                          for (int i3 = 0; i3 != x3.size(); ++i3) {
                            for (int i4 = 0; i4 != x4.size(); ++i4) {
                              for (int i5 = 0; i5 != x5.size(); ++i5) {
                                odata[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i0+x0.size()*(i4+x2.size()*(i1)))))]
                                  += (1.0) * data[i5+x5.size()*(i1+x1.size()*(i3+x3.size()*(i0)))];
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                  {
                    if (x3 == x4) {
                      std::vector<size_t> i0hash = {x5.key(), x0.key(), x2.key(), x1.key()};
                      std::unique_ptr<double[]> data = rdm2->get_block(i0hash);
                      for (int i1 = 0; i1 != x1.size(); ++i1) {
                        for (int i2 = 0; i2 != x2.size(); ++i2) {
                          for (int i0 = 0; i0 != x0.size(); ++i0) {
                            for (int i4 = 0; i4 != x4.size(); ++i4) {
                              for (int i5 = 0; i5 != x5.size(); ++i5) {
                                odata[i5+x5.size()*(i4+x4.size()*(i4+x3.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))]
                                  += (-1.0) * data[i5+x5.size()*(i0+x0.size()*(i2+x2.size()*(i1)))];
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                  {
                    std::vector<size_t> i0hash = {x5.key(), x4.key(), x3.key(), x0.key(), x2.key(), x1.key()};
                    std::unique_ptr<double[]> data = rdm3->get_block(i0hash);
                    sort_indices<0,1,2,3,4,5,1,1,-1,1>(data, odata, x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
                  }
                  Gamma5->put_block(ohash, odata);
                }
              }
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
      Gamma5  = t[0];
      rdm2    = t[1];
      rdm3    = t[2];
    };
    ~Task3() {};
};

template <typename T>
class Task4 : public Task<T> {  // associated with gamma
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > Gamma2;
    std::shared_ptr<Tensor<T> > rdm2;
    std::shared_ptr<Tensor<T> > rdm3;
    std::shared_ptr<Tensor<T> > rdm4;
    std::shared_ptr<Tensor<T> > f1;

    void compute_() {
      for (auto& x0 : active_) {
        for (auto& x1 : active_) {
          for (auto& x2 : active_) {
            for (auto& x5 : active_) {
              for (auto& x6 : active_) {
                for (auto& x7 : active_) {
                  std::vector<size_t> ohash = {x7.key(), x6.key(), x5.key(), x2.key(), x1.key(), x0.key()};
                  std::unique_ptr<double[]> odata = Gamma2->move_block(ohash);
                  // associated with merged
                  for (auto& x3 : active_) {
                    for (auto& x4 : active_) {
                      std::vector<size_t> fhash = {x4.key(), x3.key()};
                      std::unique_ptr<double[]> fdata = f1->get_block(fhash);
                      if (x1 == x6) {
                        std::vector<size_t> i0hash = {x7.key(), x0.key(), x5.key(), x2.key(), x4.key(), x3.key()};
                        std::unique_ptr<double[]> data = rdm3->get_block(i0hash);
                        for (int i0 = 0; i0 != x0.size(); ++i0) {
                          for (int i2 = 0; i2 != x2.size(); ++i2) {
                            for (int i5 = 0; i5 != x5.size(); ++i5) {
                              for (int i6 = 0; i6 != x6.size(); ++i6) {
                                for (int i7 = 0; i7 != x7.size(); ++i7) {
                                  for (int i3 = 0; i3 != x3.size(); ++i3) {
                                    for (int i4 = 0; i4 != x4.size(); ++i4) {
                                      odata[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i2+x2.size()*(i6+x1.size()*(i0)))))]
                                        += (1.0) * data[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i3)))))] * fdata[i4+x4.size()*(i3)];
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                      if (x1 == x2 && x4 == x6) {
                        std::vector<size_t> i0hash = {x7.key(), x3.key(), x5.key(), x0.key()};
                        std::unique_ptr<double[]> data = rdm2->get_block(i0hash);
                        for (int i0 = 0; i0 != x0.size(); ++i0) {
                          for (int i2 = 0; i2 != x2.size(); ++i2) {
                            for (int i5 = 0; i5 != x5.size(); ++i5) {
                              for (int i6 = 0; i6 != x6.size(); ++i6) {
                                for (int i7 = 0; i7 != x7.size(); ++i7) {
                                  for (int i3 = 0; i3 != x3.size(); ++i3) {
                                    for (int i4 = 0; i4 != x4.size(); ++i4) {
                                  //    const int i4 = i6;
                                      odata[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i2+x2.size()*(i2+x1.size()*(i0)))))]
                                        += (1.0) * data[i7+x7.size()*(i3+x3.size()*(i5+x5.size()*(i0)))] * fdata[i4+x4.size()*(i3)];
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                      if (x4 == x6 && x1 == x3) {
                        std::vector<size_t> i0hash = {x7.key(), x0.key(), x5.key(), x2.key()};
                        std::unique_ptr<double[]> data = rdm2->get_block(i0hash);
                        for (int i0 = 0; i0 != x0.size(); ++i0) {
                          for (int i2 = 0; i2 != x2.size(); ++i2) {
                            for (int i5 = 0; i5 != x5.size(); ++i5) {
                              for (int i6 = 0; i6 != x6.size(); ++i6) {
                                for (int i7 = 0; i7 != x7.size(); ++i7) {
                                  for (int i3 = 0; i3 != x3.size(); ++i3) {
                                    const int i4 = i6;
                                    odata[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i2+x2.size()*(i3+x1.size()*(i0)))))]
                                      += (-1.0) * data[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i2)))] * fdata[i4+x4.size()*(i3)];
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                      if (x4 == x6) {
                        std::vector<size_t> i0hash = {x7.key(), x3.key(), x5.key(), x2.key(), x1.key(), x0.key()};
                        std::unique_ptr<double[]> data = rdm3->get_block(i0hash);
                        for (int i0 = 0; i0 != x0.size(); ++i0) {
                          for (int i1 = 0; i1 != x1.size(); ++i1) {
                            for (int i2 = 0; i2 != x2.size(); ++i2) {
                              for (int i5 = 0; i5 != x5.size(); ++i5) {
                                for (int i6 = 0; i6 != x6.size(); ++i6) {
                                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                                    for (int i3 = 0; i3 != x3.size(); ++i3) {
                                      const int i4 = i6;
                                      odata[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i2+x2.size()*(i1+x1.size()*(i0)))))]
                                        += (-1.0) * data[i7+x7.size()*(i3+x3.size()*(i5+x5.size()*(i2+x2.size()*(i1+x1.size()*(i0)))))] * fdata[i4+x4.size()*(i3)];
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                      if (x5 == x6 && x1 == x2) {
                        std::vector<size_t> i0hash = {x7.key(), x0.key(), x4.key(), x3.key()};
                        std::unique_ptr<double[]> data = rdm2->get_block(i0hash);
                        for (int i0 = 0; i0 != x0.size(); ++i0) {
                          for (int i2 = 0; i2 != x2.size(); ++i2) {
                            for (int i6 = 0; i6 != x6.size(); ++i6) {
                              for (int i7 = 0; i7 != x7.size(); ++i7) {
                                for (int i3 = 0; i3 != x3.size(); ++i3) {
                                  for (int i4 = 0; i4 != x4.size(); ++i4) {
                                    odata[i7+x7.size()*(i6+x6.size()*(i6+x5.size()*(i2+x2.size()*(i2+x1.size()*(i0)))))]
                                      += (-1.0) * data[i7+x7.size()*(i0+x0.size()*(i4+x4.size()*(i3)))] * fdata[i4+x4.size()*(i3)];
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                      if (x1 == x2) {
                        std::vector<size_t> i0hash = {x7.key(), x6.key(), x5.key(), x0.key(), x4.key(), x3.key()};
                        std::unique_ptr<double[]> data = rdm3->get_block(i0hash);
                        for (int i0 = 0; i0 != x0.size(); ++i0) {
                          for (int i2 = 0; i2 != x2.size(); ++i2) {
                            for (int i5 = 0; i5 != x5.size(); ++i5) {
                              for (int i6 = 0; i6 != x6.size(); ++i6) {
                                for (int i7 = 0; i7 != x7.size(); ++i7) {
                                  for (int i3 = 0; i3 != x3.size(); ++i3) {
                                    for (int i4 = 0; i4 != x4.size(); ++i4) {
                                      odata[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i2+x2.size()*(i2+x1.size()*(i0)))))]
                                        += (-1.0) * data[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i3)))))] * fdata[i4+x4.size()*(i3)];
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                      if (x1 == x3 && x5 == x6) {
                        std::vector<size_t> i0hash = {x7.key(), x2.key(), x4.key(), x0.key()};
                        std::unique_ptr<double[]> data = rdm2->get_block(i0hash);
                        for (int i0 = 0; i0 != x0.size(); ++i0) {
                          for (int i2 = 0; i2 != x2.size(); ++i2) {
                            for (int i6 = 0; i6 != x6.size(); ++i6) {
                              for (int i7 = 0; i7 != x7.size(); ++i7) {
                                for (int i3 = 0; i3 != x3.size(); ++i3) {
                                  for (int i4 = 0; i4 != x4.size(); ++i4) {
                                    odata[i7+x7.size()*(i6+x6.size()*(i6+x5.size()*(i2+x2.size()*(i3+x1.size()*(i0)))))]
                                      += (1.0) * data[i7+x7.size()*(i2+x2.size()*(i4+x4.size()*(i0)))] * fdata[i4+x4.size()*(i3)];
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                      if (x1 == x3) {
                        std::vector<size_t> i0hash = {x7.key(), x6.key(), x5.key(), x2.key(), x4.key(), x0.key()};
                        std::unique_ptr<double[]> data = rdm3->get_block(i0hash);
                        for (int i0 = 0; i0 != x0.size(); ++i0) {
                          for (int i2 = 0; i2 != x2.size(); ++i2) {
                            for (int i5 = 0; i5 != x5.size(); ++i5) {
                              for (int i6 = 0; i6 != x6.size(); ++i6) {
                                for (int i7 = 0; i7 != x7.size(); ++i7) {
                                  for (int i3 = 0; i3 != x3.size(); ++i3) {
                                    for (int i4 = 0; i4 != x4.size(); ++i4) {
                                      odata[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i2+x2.size()*(i3+x1.size()*(i0)))))]
                                        += (1.0) * data[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i0)))))] * fdata[i4+x4.size()*(i3)];
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                      if (x5 == x6) {
                        std::vector<size_t> i0hash = {x7.key(), x2.key(), x4.key(), x3.key(), x1.key(), x0.key()};
                        std::unique_ptr<double[]> data = rdm3->get_block(i0hash);
                        for (int i0 = 0; i0 != x0.size(); ++i0) {
                          for (int i1 = 0; i1 != x1.size(); ++i1) {
                            for (int i2 = 0; i2 != x2.size(); ++i2) {
                              for (int i6 = 0; i6 != x6.size(); ++i6) {
                                for (int i7 = 0; i7 != x7.size(); ++i7) {
                                  for (int i3 = 0; i3 != x3.size(); ++i3) {
                                    for (int i4 = 0; i4 != x4.size(); ++i4) {
                                      odata[i7+x7.size()*(i6+x6.size()*(i6+x5.size()*(i2+x2.size()*(i1+x1.size()*(i0)))))]
                                        += (1.0) * data[i7+x7.size()*(i2+x2.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))] * fdata[i4+x4.size()*(i3)];
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                      {
                        std::vector<size_t> i0hash = {x7.key(), x6.key(), x5.key(), x2.key(), x4.key(), x3.key(), x1.key(), x0.key()};
                        std::unique_ptr<double[]> data = rdm4->get_block(i0hash);
                        for (int i0 = 0; i0 != x0.size(); ++i0) {
                          for (int i1 = 0; i1 != x1.size(); ++i1) {
                            for (int i2 = 0; i2 != x2.size(); ++i2) {
                              for (int i5 = 0; i5 != x5.size(); ++i5) {
                                for (int i6 = 0; i6 != x6.size(); ++i6) {
                                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                                    for (int i3 = 0; i3 != x3.size(); ++i3) {
                                      for (int i4 = 0; i4 != x4.size(); ++i4) {
                                        odata[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i2+x2.size()*(i1+x1.size()*(i0)))))]
                                          += (1.0) * data[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))] * fdata[i4+x4.size()*(i3)];
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                  Gamma2->put_block(ohash, odata);
                }
              }
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
      Gamma2  = t[0];
      rdm2    = t[1];
      rdm3    = t[2];
      rdm4    = t[3];
      f1      = t[4];
    };
    ~Task4() {};
};

template <typename T>
class Task5 : public Task<T> {  // associated with gamma
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > Gamma3;
    std::shared_ptr<Tensor<T> > rdm1;
    std::shared_ptr<Tensor<T> > rdm2;
    std::shared_ptr<Tensor<T> > rdm3;

    void compute_() {
      for (auto& x0 : active_) {
        for (auto& x1 : active_) {
          for (auto& x2 : active_) {
            for (auto& x3 : active_) {
              for (auto& x4 : active_) {
                for (auto& x5 : active_) {
                  std::vector<size_t> ohash = {x5.key(), x4.key(), x3.key(), x2.key(), x1.key(), x0.key()};
                  std::unique_ptr<double[]> odata = Gamma3->move_block(ohash);
                  {
                    if (x1 == x4) {
                      std::vector<size_t> i0hash = {x5.key(), x0.key(), x3.key(), x2.key()};
                      std::unique_ptr<double[]> data = rdm2->get_block(i0hash);
                      for (int i0 = 0; i0 != x0.size(); ++i0) {
                        for (int i2 = 0; i2 != x2.size(); ++i2) {
                          for (int i3 = 0; i3 != x3.size(); ++i3) {
                            for (int i4 = 0; i4 != x4.size(); ++i4) {
                              for (int i5 = 0; i5 != x5.size(); ++i5) {
                                odata[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i4+x1.size()*(i0)))))]
                                  += (-1.0) * data[i5+x5.size()*(i0+x0.size()*(i3+x3.size()*(i2)))];
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                  {
                    if (x1 == x2 && x3 == x4) {
                      std::vector<size_t> i0hash = {x5.key(), x0.key()};
                      std::unique_ptr<double[]> data = rdm1->get_block(i0hash);
                      for (int i0 = 0; i0 != x0.size(); ++i0) {
                        for (int i2 = 0; i2 != x2.size(); ++i2) {
                          for (int i4 = 0; i4 != x4.size(); ++i4) {
                            for (int i5 = 0; i5 != x5.size(); ++i5) {
                              odata[i5+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2+x2.size()*(i2+x1.size()*(i0)))))]
                                += (1.0) * data[i5+x5.size()*(i0)];
                            }
                          }
                        }
                      }
                    }
                  }
                  {
                    if (x1 == x2) {
                      std::vector<size_t> i0hash = {x5.key(), x4.key(), x3.key(), x0.key()};
                      std::unique_ptr<double[]> data = rdm2->get_block(i0hash);
                      for (int i0 = 0; i0 != x0.size(); ++i0) {
                        for (int i2 = 0; i2 != x2.size(); ++i2) {
                          for (int i3 = 0; i3 != x3.size(); ++i3) {
                            for (int i4 = 0; i4 != x4.size(); ++i4) {
                              for (int i5 = 0; i5 != x5.size(); ++i5) {
                                odata[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i2+x1.size()*(i0)))))]
                                  += (1.0) * data[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i0)))];
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                  {
                    if (x3 == x4) {
                      std::vector<size_t> i0hash = {x5.key(), x2.key(), x1.key(), x0.key()};
                      std::unique_ptr<double[]> data = rdm2->get_block(i0hash);
                      for (int i0 = 0; i0 != x0.size(); ++i0) {
                        for (int i1 = 0; i1 != x1.size(); ++i1) {
                          for (int i2 = 0; i2 != x2.size(); ++i2) {
                            for (int i4 = 0; i4 != x4.size(); ++i4) {
                              for (int i5 = 0; i5 != x5.size(); ++i5) {
                                odata[i5+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2+x2.size()*(i1+x1.size()*(i0)))))]
                                  += (1.0) * data[i5+x5.size()*(i2+x2.size()*(i1+x1.size()*(i0)))];
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                  {
                    std::vector<size_t> i0hash = {x5.key(), x4.key(), x3.key(), x2.key(), x1.key(), x0.key()};
                    std::unique_ptr<double[]> data = rdm3->get_block(i0hash);
                    sort_indices<0,1,2,3,4,5,1,1,1,1>(data, odata, x5.size(), x4.size(), x3.size(), x2.size(), x1.size(), x0.size());
                  }
                  Gamma3->put_block(ohash, odata);
                }
              }
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
      Gamma3  = t[0];
      rdm1    = t[1];
      rdm2    = t[2];
      rdm3    = t[3];
    };
    ~Task5() {};
};

template <typename T>
class Task6 : public Task<T> {  // associated with gamma
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > Gamma4;
    std::shared_ptr<Tensor<T> > rdm2;
    std::shared_ptr<Tensor<T> > rdm3;

    void compute_() {
      for (auto& x0 : active_) {
        for (auto& x1 : active_) {
          for (auto& x3 : active_) {
            for (auto& x4 : active_) {
              for (auto& x2 : active_) {
                for (auto& x5 : active_) {
                  std::vector<size_t> ohash = {x5.key(), x2.key(), x4.key(), x3.key(), x1.key(), x0.key()};
                  std::unique_ptr<double[]> odata = Gamma4->move_block(ohash);
                  {
                    if (x1 == x2) {
                      std::vector<size_t> i0hash = {x5.key(), x0.key(), x4.key(), x3.key()};
                      std::unique_ptr<double[]> data = rdm2->get_block(i0hash);
                      for (int i0 = 0; i0 != x0.size(); ++i0) {
                        for (int i3 = 0; i3 != x3.size(); ++i3) {
                          for (int i4 = 0; i4 != x4.size(); ++i4) {
                            for (int i2 = 0; i2 != x2.size(); ++i2) {
                              for (int i5 = 0; i5 != x5.size(); ++i5) {
                                odata[i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i3+x3.size()*(i2+x1.size()*(i0)))))]
                                  += (-1.0) * data[i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i3)))];
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                  {
                    if (x1 == x3) {
                      std::vector<size_t> i0hash = {x5.key(), x2.key(), x4.key(), x0.key()};
                      std::unique_ptr<double[]> data = rdm2->get_block(i0hash);
                      for (int i0 = 0; i0 != x0.size(); ++i0) {
                        for (int i3 = 0; i3 != x3.size(); ++i3) {
                          for (int i4 = 0; i4 != x4.size(); ++i4) {
                            for (int i2 = 0; i2 != x2.size(); ++i2) {
                              for (int i5 = 0; i5 != x5.size(); ++i5) {
                                odata[i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))]
                                  += (1.0) * data[i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i0)))];
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                  {
                    std::vector<size_t> i0hash = {x5.key(), x2.key(), x4.key(), x3.key(), x1.key(), x0.key()};
                    std::unique_ptr<double[]> data = rdm3->get_block(i0hash);
                    sort_indices<0,1,2,3,4,5,1,1,1,1>(data, odata, x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
                  }
                  Gamma4->put_block(ohash, odata);
                }
              }
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
      Gamma4  = t[0];
      rdm2    = t[1];
      rdm3    = t[2];
    };
    ~Task6() {};
};

template <typename T>
class Task7 : public Task<T> {
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
    Task7(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      r = t[0];
      I0 = t[1];
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
    Task8(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I0 = t[0];
      t2 = t[1];
      I1 = t[2];
    };
    ~Task8() {};
};

template <typename T>
class Task9 : public Task<T> {
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
    Task9(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I1 = t[0];
      Gamma0 = t[1];
    };
    ~Task9() {};
};

template <typename T>
class Task10 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I0;
    std::shared_ptr<Tensor<T> > t2;
    std::shared_ptr<Tensor<T> > I17;

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
                  std::unique_ptr<double[]> i1data = I17->get_block(i1hash);
                  std::unique_ptr<double[]> i1data_sorted(new double[I17->get_size(i1hash)]);
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
    Task10(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I0 = t[0];
      t2 = t[1];
      I17 = t[2];
    };
    ~Task10() {};
};

template <typename T>
class Task11 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I17;
    std::shared_ptr<Tensor<T> > Gamma6;
    double e0_;

    void compute_() {
      for (auto& x1 : active_) {
        for (auto& x2 : active_) {
          for (auto& x0 : active_) {
            for (auto& x3 : active_) {
              std::vector<size_t> ohash = {x3.key(), x0.key(), x2.key(), x1.key()};
              std::unique_ptr<double[]> odata = I17->move_block(ohash);
              {
                std::vector<size_t> i0hash = {x3.key(), x0.key(), x2.key(), x1.key()};
                std::unique_ptr<double[]> i0data = Gamma6->get_block(i0hash);
                dscal_(x3.size()*x0.size()*x2.size()*x1.size(), -e0_, i0data.get(), 1);
                sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
              }
              I17->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task11(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i, const double e) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I17 = t[0];
      Gamma6 = t[1];
      e0_ = e; 
    };
    ~Task11() {};
};

template <typename T>
class Task12 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I0;
    std::shared_ptr<Tensor<T> > v2;
    std::shared_ptr<Tensor<T> > I19;

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
                  std::unique_ptr<double[]> i1data = I19->get_block(i1hash);
                  std::unique_ptr<double[]> i1data_sorted(new double[I19->get_size(i1hash)]);
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
    Task12(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I0 = t[0];
      v2 = t[1];
      I19 = t[2];
    };
    ~Task12() {};
};

template <typename T>
class Task13 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I19;
    std::shared_ptr<Tensor<T> > Gamma6;

    void compute_() {
      for (auto& x1 : active_) {
        for (auto& x2 : active_) {
          for (auto& x0 : active_) {
            for (auto& x3 : active_) {
              std::vector<size_t> ohash = {x3.key(), x0.key(), x2.key(), x1.key()};
              std::unique_ptr<double[]> odata = I19->move_block(ohash);
              {
                std::vector<size_t> i0hash = {x3.key(), x0.key(), x2.key(), x1.key()};
                std::unique_ptr<double[]> i0data = Gamma6->get_block(i0hash);
                sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
              }
              I19->put_block(ohash, odata);
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
      I19 = t[0];
      Gamma6 = t[1];
    };
    ~Task13() {};
};

template <typename T>
class Task14 : public Task<T> {
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
    Task14(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      r = t[0];
      I2 = t[1];
    };
    ~Task14() {};
};

template <typename T>
class Task15 : public Task<T> {
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
    Task15(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I2 = t[0];
      f1 = t[1];
      I3 = t[2];
    };
    ~Task15() {};
};

template <typename T>
class Task16 : public Task<T> {
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
    Task16(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I3 = t[0];
      t2 = t[1];
      I4 = t[2];
    };
    ~Task16() {};
};

template <typename T>
class Task17 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I4;
    std::shared_ptr<Tensor<T> > Gamma6;

    void compute_() {
      for (auto& x1 : active_) {
        for (auto& x2 : active_) {
          for (auto& x0 : active_) {
            for (auto& x3 : active_) {
              std::vector<size_t> ohash = {x3.key(), x0.key(), x2.key(), x1.key()};
              std::unique_ptr<double[]> odata = I4->move_block(ohash);
              {
                std::vector<size_t> i0hash = {x3.key(), x0.key(), x2.key(), x1.key()};
                std::unique_ptr<double[]> i0data = Gamma6->get_block(i0hash);
                sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
              }
              I4->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task17(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I4 = t[0];
      Gamma6 = t[1];
    };
    ~Task17() {};
};

template <typename T>
class Task18 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I2;
    std::shared_ptr<Tensor<T> > f1;
    std::shared_ptr<Tensor<T> > I14;

    void compute_() {
      for (auto& a2 : virt_) {
        for (auto& a1 : virt_) {
          for (auto& x1 : active_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), x1.key(), a1.key(), a2.key()};
              std::unique_ptr<double[]> odata = I2->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I2->get_size(ohash)]);
              std::fill(odata_sorted.get(), odata_sorted.get()+I2->get_size(ohash), 0.0);

              for (auto& x2 : active_) {
                std::vector<size_t> i0hash = {x2.key(), a2.key()};
                std::unique_ptr<double[]> i0data = f1->get_block(i0hash);
                std::unique_ptr<double[]> i0data_sorted(new double[f1->get_size(i0hash)]);
                sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), a2.size());

                std::vector<size_t> i1hash = {x0.key(), x2.key(), x1.key(), a1.key()};
                std::unique_ptr<double[]> i1data = I14->get_block(i1hash);
                std::unique_ptr<double[]> i1data_sorted(new double[I14->get_size(i1hash)]);
                sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), x2.size(), x1.size(), a1.size());

                dgemm_("T", "N", a2.size(), x0.size()*x1.size()*a1.size(), x2.size(),
                       1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
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
    Task18(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I2 = t[0];
      f1 = t[1];
      I14 = t[2];
    };
    ~Task18() {};
};

template <typename T>
class Task19 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I14;
    std::shared_ptr<Tensor<T> > t2;
    std::shared_ptr<Tensor<T> > I15;

    void compute_() {
      for (auto& a1 : virt_) {
        for (auto& x1 : active_) {
          for (auto& x2 : active_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), x2.key(), x1.key(), a1.key()};
              std::unique_ptr<double[]> odata = I14->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I14->get_size(ohash)]);
              std::fill(odata_sorted.get(), odata_sorted.get()+I14->get_size(ohash), 0.0);

              for (auto& x5 : active_) {
                for (auto& x4 : active_) {
                  for (auto& x3 : active_) {
                    std::vector<size_t> i0hash = {x5.key(), x4.key(), x3.key(), a1.key()};
                    std::unique_ptr<double[]> i0data = t2->get_block(i0hash);
                    std::unique_ptr<double[]> i0data_sorted(new double[t2->get_size(i0hash)]);
                    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), a1.size());

                    std::vector<size_t> i1hash = {x5.key(), x4.key(), x3.key(), x0.key(), x2.key(), x1.key()};
                    std::unique_ptr<double[]> i1data = I15->get_block(i1hash);
                    std::unique_ptr<double[]> i1data_sorted(new double[I15->get_size(i1hash)]);
                    sort_indices<0,1,2,3,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());

                    dgemm_("T", "N", a1.size(), x0.size()*x2.size()*x1.size(), x5.size()*x4.size()*x3.size(),
                           1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
                           1.0, odata_sorted, a1.size());
                  }
                }
              }

              sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x2.size(), x1.size());
              I14->put_block(ohash, odata);
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
      I14 = t[0];
      t2 = t[1];
      I15 = t[2];
    };
    ~Task19() {};
};

template <typename T>
class Task20 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I15;
    std::shared_ptr<Tensor<T> > Gamma5;

    void compute_() {
      for (auto& x1 : active_) {
        for (auto& x2 : active_) {
          for (auto& x0 : active_) {
            for (auto& x3 : active_) {
              for (auto& x4 : active_) {
                for (auto& x5 : active_) {
                  std::vector<size_t> ohash = {x5.key(), x4.key(), x3.key(), x0.key(), x2.key(), x1.key()};
                  std::unique_ptr<double[]> odata = I15->move_block(ohash);
                  {
                    std::vector<size_t> i0hash = {x5.key(), x4.key(), x3.key(), x0.key(), x2.key(), x1.key()};
                    std::unique_ptr<double[]> i0data = Gamma5->get_block(i0hash);
                    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
                  }
                  I15->put_block(ohash, odata);
                }
              }
            }
          }
        }
      }
    };

  public:
    Task20(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I15 = t[0];
      Gamma5 = t[1];
    };
    ~Task20() {};
};

template <typename T>
class Task21 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > r;
    std::shared_ptr<Tensor<T> > I5;

    void compute_() {
      for (auto& x1 : active_) {
        for (auto& x0 : active_) {
          for (auto& a1 : virt_) {
            for (auto& x2 : active_) {
              std::vector<size_t> ohash = {x2.key(), a1.key(), x0.key(), x1.key()};
              std::unique_ptr<double[]> odata = r->move_block(ohash);
              {
                std::vector<size_t> i0hash = {x2.key(), x1.key(), x0.key(), a1.key()};
                std::unique_ptr<double[]> i0data = I5->get_block(i0hash);
                sort_indices<0,3,2,1,1,1,1,1>(i0data, odata, x2.size(), x1.size(), x0.size(), a1.size());
              }
              r->put_block(ohash, odata);
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
      r = t[0];
      I5 = t[1];
    };
    ~Task21() {};
};

template <typename T>
class Task22 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I5;
    std::shared_ptr<Tensor<T> > t2;
    std::shared_ptr<Tensor<T> > I6;

    void compute_() {
      for (auto& a1 : virt_) {
        for (auto& x0 : active_) {
          for (auto& x1 : active_) {
            for (auto& x2 : active_) {
              std::vector<size_t> ohash = {x2.key(), x1.key(), x0.key(), a1.key()};
              std::unique_ptr<double[]> odata = I5->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I5->get_size(ohash)]);
              std::fill(odata_sorted.get(), odata_sorted.get()+I5->get_size(ohash), 0.0);

              for (auto& x7 : active_) {
                for (auto& x6 : active_) {
                  for (auto& x5 : active_) {
                    std::vector<size_t> i0hash = {x7.key(), x6.key(), x5.key(), a1.key()};
                    std::unique_ptr<double[]> i0data = t2->get_block(i0hash);
                    std::unique_ptr<double[]> i0data_sorted(new double[t2->get_size(i0hash)]);
                    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x5.size(), a1.size());

                    std::vector<size_t> i1hash = {x7.key(), x6.key(), x5.key(), x2.key(), x1.key(), x0.key()};
                    std::unique_ptr<double[]> i1data = I6->get_block(i1hash);
                    std::unique_ptr<double[]> i1data_sorted(new double[I6->get_size(i1hash)]);
                    sort_indices<0,1,2,3,4,5,0,1,1,1>(i1data, i1data_sorted, x7.size(), x6.size(), x5.size(), x2.size(), x1.size(), x0.size());

                    dgemm_("T", "N", a1.size(), x2.size()*x1.size()*x0.size(), x7.size()*x6.size()*x5.size(),
                           1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
                           1.0, odata_sorted, a1.size());
                  }
                }
              }

              sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x2.size(), x1.size(), x0.size());
              I5->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task22(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I5 = t[0];
      t2 = t[1];
      I6 = t[2];
    };
    ~Task22() {};
};

template <typename T>
class Task23 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I6;
    std::shared_ptr<Tensor<T> > Gamma2;

    void compute_() {
      for (auto& x0 : active_) {
        for (auto& x1 : active_) {
          for (auto& x2 : active_) {
            for (auto& x5 : active_) {
              for (auto& x6 : active_) {
                for (auto& x7 : active_) {
                  std::vector<size_t> ohash = {x7.key(), x6.key(), x5.key(), x2.key(), x1.key(), x0.key()};
                  std::unique_ptr<double[]> odata = I6->move_block(ohash);
                  {
                    std::vector<size_t> i0hash = {x7.key(), x6.key(), x5.key(), x2.key(), x1.key(), x0.key()};
                    std::unique_ptr<double[]> i0data = Gamma2->get_block(i0hash);
                    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x7.size(), x6.size(), x5.size(), x2.size(), x1.size(), x0.size());
                  }
                  I6->put_block(ohash, odata);
                }
              }
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
      I6 = t[0];
      Gamma2 = t[1];
    };
    ~Task23() {};
};

template <typename T>
class Task24 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I5;
    std::shared_ptr<Tensor<T> > f1;
    std::shared_ptr<Tensor<T> > I8;

    void compute_() {
      for (auto& a1 : virt_) {
        for (auto& x0 : active_) {
          for (auto& x1 : active_) {
            for (auto& x2 : active_) {
              std::vector<size_t> ohash = {x2.key(), x1.key(), x0.key(), a1.key()};
              std::unique_ptr<double[]> odata = I5->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I5->get_size(ohash)]);
              std::fill(odata_sorted.get(), odata_sorted.get()+I5->get_size(ohash), 0.0);

              for (auto& a2 : virt_) {
                std::vector<size_t> i0hash = {a2.key(), a1.key()};
                std::unique_ptr<double[]> i0data = f1->get_block(i0hash);
                std::unique_ptr<double[]> i0data_sorted(new double[f1->get_size(i0hash)]);
                sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a2.size(), a1.size());

                std::vector<size_t> i1hash = {x2.key(), x1.key(), x0.key(), a2.key()};
                std::unique_ptr<double[]> i1data = I8->get_block(i1hash);
                std::unique_ptr<double[]> i1data_sorted(new double[I8->get_size(i1hash)]);
                sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), x0.size(), a2.size());

                dgemm_("T", "N", a1.size(), x2.size()*x1.size()*x0.size(), a2.size(),
                       1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
                       1.0, odata_sorted, a1.size());
              }

              sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x2.size(), x1.size(), x0.size());
              I5->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task24(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I5 = t[0];
      f1 = t[1];
      I8 = t[2];
    };
    ~Task24() {};
};

template <typename T>
class Task25 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I8;
    std::shared_ptr<Tensor<T> > t2;
    std::shared_ptr<Tensor<T> > I9;

    void compute_() {
      for (auto& a2 : virt_) {
        for (auto& x0 : active_) {
          for (auto& x1 : active_) {
            for (auto& x2 : active_) {
              std::vector<size_t> ohash = {x2.key(), x1.key(), x0.key(), a2.key()};
              std::unique_ptr<double[]> odata = I8->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I8->get_size(ohash)]);
              std::fill(odata_sorted.get(), odata_sorted.get()+I8->get_size(ohash), 0.0);

              for (auto& x5 : active_) {
                for (auto& x4 : active_) {
                  for (auto& x3 : active_) {
                    std::vector<size_t> i0hash = {x5.key(), x4.key(), x3.key(), a2.key()};
                    std::unique_ptr<double[]> i0data = t2->get_block(i0hash);
                    std::unique_ptr<double[]> i0data_sorted(new double[t2->get_size(i0hash)]);
                    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), a2.size());

                    std::vector<size_t> i1hash = {x5.key(), x4.key(), x3.key(), x2.key(), x1.key(), x0.key()};
                    std::unique_ptr<double[]> i1data = I9->get_block(i1hash);
                    std::unique_ptr<double[]> i1data_sorted(new double[I9->get_size(i1hash)]);
                    sort_indices<0,1,2,3,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x3.size(), x2.size(), x1.size(), x0.size());

                    dgemm_("T", "N", a2.size(), x2.size()*x1.size()*x0.size(), x5.size()*x4.size()*x3.size(),
                           1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
                           1.0, odata_sorted, a2.size());
                  }
                }
              }

              sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), x2.size(), x1.size(), x0.size());
              I8->put_block(ohash, odata);
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
      I8 = t[0];
      t2 = t[1];
      I9 = t[2];
    };
    ~Task25() {};
};

template <typename T>
class Task26 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I9;
    std::shared_ptr<Tensor<T> > Gamma3;

    void compute_() {
      for (auto& x0 : active_) {
        for (auto& x1 : active_) {
          for (auto& x2 : active_) {
            for (auto& x3 : active_) {
              for (auto& x4 : active_) {
                for (auto& x5 : active_) {
                  std::vector<size_t> ohash = {x5.key(), x4.key(), x3.key(), x2.key(), x1.key(), x0.key()};
                  std::unique_ptr<double[]> odata = I9->move_block(ohash);
                  {
                    std::vector<size_t> i0hash = {x5.key(), x4.key(), x3.key(), x2.key(), x1.key(), x0.key()};
                    std::unique_ptr<double[]> i0data = Gamma3->get_block(i0hash);
                    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x4.size(), x3.size(), x2.size(), x1.size(), x0.size());
                  }
                  I9->put_block(ohash, odata);
                }
              }
            }
          }
        }
      }
    };

  public:
    Task26(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I9 = t[0];
      Gamma3 = t[1];
    };
    ~Task26() {};
};

template <typename T>
class Task27 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I5;
    std::shared_ptr<Tensor<T> > f1;
    std::shared_ptr<Tensor<T> > I11;

    void compute_() {
      for (auto& a1 : virt_) {
        for (auto& x0 : active_) {
          for (auto& x1 : active_) {
            for (auto& x2 : active_) {
              std::vector<size_t> ohash = {x2.key(), x1.key(), x0.key(), a1.key()};
              std::unique_ptr<double[]> odata = I5->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I5->get_size(ohash)]);
              std::fill(odata_sorted.get(), odata_sorted.get()+I5->get_size(ohash), 0.0);

              for (auto& a2 : virt_) {
                for (auto& x3 : active_) {
                  std::vector<size_t> i0hash = {a2.key(), x3.key()};
                  std::unique_ptr<double[]> i0data = f1->get_block(i0hash);
                  std::unique_ptr<double[]> i0data_sorted(new double[f1->get_size(i0hash)]);
                  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a2.size(), x3.size());

                  std::vector<size_t> i1hash = {x2.key(), x3.key(), x1.key(), x0.key(), a1.key(), a2.key()};
                  std::unique_ptr<double[]> i1data = I11->get_block(i1hash);
                  std::unique_ptr<double[]> i1data_sorted(new double[I11->get_size(i1hash)]);
                  sort_indices<5,1,0,2,3,4,0,1,1,1>(i1data, i1data_sorted, x2.size(), x3.size(), x1.size(), x0.size(), a1.size(), a2.size());

                  dgemm_("T", "N", 1, x2.size()*x1.size()*x0.size()*a1.size(), x3.size()*a2.size(),
                         1.0, i0data_sorted, x3.size()*a2.size(), i1data_sorted, x3.size()*a2.size(),
                         1.0, odata_sorted, 1);
                }
              }

              sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), a1.size());
              I5->put_block(ohash, odata);
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
      I5 = t[0];
      f1 = t[1];
      I11 = t[2];
    };
    ~Task27() {};
};

template <typename T>
class Task28 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I11;
    std::shared_ptr<Tensor<T> > t2;
    std::shared_ptr<Tensor<T> > I12;

    void compute_() {
      for (auto& a2 : virt_) {
        for (auto& a1 : virt_) {
          for (auto& x0 : active_) {
            for (auto& x1 : active_) {
              for (auto& x3 : active_) {
                for (auto& x2 : active_) {
                  std::vector<size_t> ohash = {x2.key(), x3.key(), x1.key(), x0.key(), a1.key(), a2.key()};
                  std::unique_ptr<double[]> odata = I11->move_block(ohash);
                  std::unique_ptr<double[]> odata_sorted(new double[I11->get_size(ohash)]);
                  std::fill(odata_sorted.get(), odata_sorted.get()+I11->get_size(ohash), 0.0);

                  for (auto& x5 : active_) {
                    for (auto& x4 : active_) {
                      std::vector<size_t> i0hash = {x5.key(), a1.key(), x4.key(), a2.key()};
                      std::unique_ptr<double[]> i0data = t2->get_block(i0hash);
                      std::unique_ptr<double[]> i0data_sorted(new double[t2->get_size(i0hash)]);
                      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a2.size());

                      std::vector<size_t> i1hash = {x5.key(), x2.key(), x4.key(), x3.key(), x1.key(), x0.key()};
                      std::unique_ptr<double[]> i1data = I12->get_block(i1hash);
                      std::unique_ptr<double[]> i1data_sorted(new double[I12->get_size(i1hash)]);
                      sort_indices<0,2,1,3,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());

                      dgemm_("T", "N", a1.size()*a2.size(), x2.size()*x3.size()*x1.size()*x0.size(), x5.size()*x4.size(),
                             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
                             1.0, odata_sorted, a1.size()*a2.size());
                    }
                  }

                  sort_indices<2,3,4,5,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x2.size(), x3.size(), x1.size(), x0.size());
                  I11->put_block(ohash, odata);
                }
              }
            }
          }
        }
      }
    };

  public:
    Task28(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I11 = t[0];
      t2 = t[1];
      I12 = t[2];
    };
    ~Task28() {};
};

template <typename T>
class Task29 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I12;
    std::shared_ptr<Tensor<T> > Gamma4;

    void compute_() {
      for (auto& x0 : active_) {
        for (auto& x1 : active_) {
          for (auto& x3 : active_) {
            for (auto& x4 : active_) {
              for (auto& x2 : active_) {
                for (auto& x5 : active_) {
                  std::vector<size_t> ohash = {x5.key(), x2.key(), x4.key(), x3.key(), x1.key(), x0.key()};
                  std::unique_ptr<double[]> odata = I12->move_block(ohash);
                  {
                    std::vector<size_t> i0hash = {x5.key(), x2.key(), x4.key(), x3.key(), x1.key(), x0.key()};
                    std::unique_ptr<double[]> i0data = Gamma4->get_block(i0hash);
                    sort_indices<0,1,2,3,4,5,1,1,2,1>(i0data, odata, x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
                  }
                  I12->put_block(ohash, odata);
                }
              }
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
      I12 = t[0];
      Gamma4 = t[1];
    };
    ~Task29() {};
};

template <typename T>
class Task30 : public EnergyTask<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I20;
    std::shared_ptr<Tensor<T> > t2;
    std::shared_ptr<Tensor<T> > I21;

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
              std::unique_ptr<double[]> i1data = I21->get_block(i1hash);
              std::unique_ptr<double[]> i1data_sorted(new double[I21->get_size(i1hash)]);
              sort_indices<3,1,2,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), a1.size(), a2.size());

              this->energy_ += ddot_(x0.size()*x1.size()*a1.size()*a2.size(), i0data_sorted, 1, i1data_sorted, 1);
            }
          }
        }
      }

    };

  public:
    Task30(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : EnergyTask<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I20 = t[0];
      t2 = t[1];
      I21 = t[2];
    };
    ~Task30() {};
};

template <typename T>
class Task31 : public EnergyTask<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I21;
    std::shared_ptr<Tensor<T> > v2;
    std::shared_ptr<Tensor<T> > I22;

    void compute_() {
      this->energy_ = 0.0;
      for (auto& a2 : virt_) {
        for (auto& a1 : virt_) {
          for (auto& x1 : active_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), x1.key(), a1.key(), a2.key()};
              std::unique_ptr<double[]> odata = I21->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I21->get_size(ohash)]);
              std::fill(odata_sorted.get(), odata_sorted.get()+I21->get_size(ohash), 0.0);

              for (auto& x3 : active_) {
                for (auto& x2 : active_) {
                  std::vector<size_t> i0hash = {x3.key(), a1.key(), x2.key(), a2.key()};
                  std::unique_ptr<double[]> i0data = v2->get_block(i0hash);
                  std::unique_ptr<double[]> i0data_sorted(new double[v2->get_size(i0hash)]);
                  sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

                  std::vector<size_t> i1hash = {x3.key(), x0.key(), x2.key(), x1.key()};
                  std::unique_ptr<double[]> i1data = I22->get_block(i1hash);
                  std::unique_ptr<double[]> i1data_sorted(new double[I22->get_size(i1hash)]);
                  sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

                  dgemm_("T", "N", a1.size()*a2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
                         1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
                         1.0, odata_sorted, a1.size()*a2.size());
                }
              }

              sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
              I21->put_block(ohash, odata);
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
      I21 = t[0];
      v2 = t[1];
      I22 = t[2];
    };
    ~Task31() {};
};

template <typename T>
class Task32 : public EnergyTask<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I22;
    std::shared_ptr<Tensor<T> > Gamma6;

    void compute_() {
      this->energy_ = 0.0;
      for (auto& x1 : active_) {
        for (auto& x2 : active_) {
          for (auto& x0 : active_) {
            for (auto& x3 : active_) {
              std::vector<size_t> ohash = {x3.key(), x0.key(), x2.key(), x1.key()};
              std::unique_ptr<double[]> odata = I22->move_block(ohash);
              {
                std::vector<size_t> i0hash = {x3.key(), x0.key(), x2.key(), x1.key()};
                std::unique_ptr<double[]> i0data = Gamma6->get_block(i0hash);
                sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
              }
              I22->put_block(ohash, odata);
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
      I22 = t[0];
      Gamma6 = t[1];
    };
    ~Task32() {};
};

template <typename T>
class Task33 : public EnergyTask<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I21;
    std::shared_ptr<Tensor<T> > r;
    std::shared_ptr<Tensor<T> > I25;

    void compute_() {
      this->energy_ = 0.0;
      for (auto& a2 : virt_) {
        for (auto& a1 : virt_) {
          for (auto& x1 : active_) {
            for (auto& x0 : active_) {
              std::vector<size_t> ohash = {x0.key(), x1.key(), a1.key(), a2.key()};
              std::unique_ptr<double[]> odata = I21->move_block(ohash);
              std::unique_ptr<double[]> odata_sorted(new double[I21->get_size(ohash)]);
              std::fill(odata_sorted.get(), odata_sorted.get()+I21->get_size(ohash), 0.0);

              for (auto& x3 : active_) {
                for (auto& x2 : active_) {
                  std::vector<size_t> i0hash = {x3.key(), a1.key(), x2.key(), a2.key()};
                  std::unique_ptr<double[]> i0data = r->get_block(i0hash);
                  std::unique_ptr<double[]> i0data_sorted(new double[r->get_size(i0hash)]);
                  sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

                  std::vector<size_t> i1hash = {x3.key(), x0.key(), x2.key(), x1.key()};
                  std::unique_ptr<double[]> i1data = I25->get_block(i1hash);
                  std::unique_ptr<double[]> i1data_sorted(new double[I25->get_size(i1hash)]);
                  sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

                  dgemm_("T", "N", a1.size()*a2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
                         1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
                         1.0, odata_sorted, a1.size()*a2.size());
                }
              }

              sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
              I21->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task33(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : EnergyTask<T>() {
      closed_ = i[0];
      active_ = i[1];
      virt_   = i[2];
      I21 = t[0];
      r = t[1];
      I25 = t[2];
    };
    ~Task33() {};
};

template <typename T>
class Task34 : public EnergyTask<T> {
  protected:
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I25;
    std::shared_ptr<Tensor<T> > Gamma6;

    void compute_() {
      this->energy_ = 0.0;
      for (auto& x1 : active_) {
        for (auto& x2 : active_) {
          for (auto& x0 : active_) {
            for (auto& x3 : active_) {
              std::vector<size_t> ohash = {x3.key(), x0.key(), x2.key(), x1.key()};
              std::unique_ptr<double[]> odata = I25->move_block(ohash);
              {
                std::vector<size_t> i0hash = {x3.key(), x0.key(), x2.key(), x1.key()};
                std::unique_ptr<double[]> i0data = Gamma6->get_block(i0hash);
                sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
              }
              I25->put_block(ohash, odata);
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
      I25 = t[0];
      Gamma6 = t[1];
    };
    ~Task34() {};
};


}
}
}
#endif

