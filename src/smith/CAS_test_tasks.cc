//
// BAGEL - Parallel electron correlation program.
// Filename: CAS_test_tasks.cc
// Copyright (C) 2014 Shiozaki group
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


#include <src/smith/CAS_test_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CAS_test;

void Task1::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x5 = b(2);
  const Index x6 = b(3);
  const Index x0 = b(4);
  const Index x7 = b(5);
  const Index x3 = b(6);
  const Index x4 = b(7);
  // tensor label: Gamma0
  std::unique_ptr<double[]> odata = out()->move_block(x7, x0, x6, x5, x2, x1);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(3)->get_block(x4, x3);
  if (x2 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x6, x1, x4, x3);
    for (int i3 = 0; i3 != x3.size(); ++i3) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i5+x2.size()*(i1)))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i1+x1.size()*(i4+x4.size()*(i3)))))] * fdata[i4+x4.size()*(i3)];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x5 && x2 == x3) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i6 = 0; i6 != x6.size(); ++i6) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i3+x2.size()*(i1)))))]
                  += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i1)))] * fdata[i5+x4.size()*(i3)];
              }
            }
          }
        }
      }
    }
  }
  if (x2 == x3) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x6, x5, x4, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i3+x2.size()*(i1)))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i1)))))] * fdata[i4+x4.size()*(i3)];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x6, x3, x2, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i2+x2.size()*(i1)))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i3+x3.size()*(i2+x2.size()*(i1)))))] * fdata[i5+x4.size()*(i3)];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x0, x6, x5, x4, x3, x2, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i2+x2.size()*(i1)))))]
                      += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i1)))))))] * fdata[i4+x4.size()*(i3)];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  out()->put_block(odata, x7, x0, x6, x5, x2, x1);
}


void Task2::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index x0 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma1
  std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x3, x2, x1);
  {
    if (x2 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i3+x3.size()*(i3+x2.size()*(i1)))))]
                  += (1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i1)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x0, x4, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x5, x0, x4, x3, x2, x1);
}


void Task3::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  const Index x4 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma4
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, x3, x0, x2, x1);
  {
    if (x2 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x1, x3, x0);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i0+x0.size()*(i4+x2.size()*(i1)))))]
                  += (1.0) * i0data[i5+x5.size()*(i1+x1.size()*(i3+x3.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i4+x3.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))]
                  += (1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i2+x2.size()*(i1)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x5, x4, x3, x0, x2, x1);
}


void Task4::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  // tensor label: Gamma5
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task5::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index x5 = b(4);
  const Index x6 = b(5);
  const Index x0 = b(6);
  const Index x7 = b(7);
  // tensor label: Gamma13
  std::unique_ptr<double[]> odata = out()->move_block(x7, x0, x6, x5, x4, x3, x2, x1);
  {
    if (x2 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x6, x1, x4, x3);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i5+x2.size()*(i1)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i1+x1.size()*(i4+x4.size()*(i3)))))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x4 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i3+x2.size()*(i1)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i1)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x6, x5, x4, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i3+x2.size()*(i1)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i1)))))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x6, x3, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i2+x2.size()*(i1)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i3+x3.size()*(i2+x2.size()*(i1)))))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x0, x6, x5, x4, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,6,7,1,1,1,1>(i0data, odata, x7.size(), x0.size(), x6.size(), x5.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x7, x0, x6, x5, x4, x3, x2, x1);
}


void Task6::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x7 = b(1);
  const Index x0 = b(2);
  const Index x6 = b(3);
  const Index x5 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);
  const Index x4 = b(7);
  const Index x3 = b(8);
  // tensor label: Gamma17
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x7, x0, x6, x5, x2, x1);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(3)->get_block(x4, x3);
  if (x2 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x7, x0, x6, x1, x4, x3);
    for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
      for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
        for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
          for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
            for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
              for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
                for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                  for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                    odata[ici0+ci0.size()*(ix7+x7.size()*(ix0+x0.size()*(ix6+x6.size()*(ix5+x5.size()*(ix5+x2.size()*(ix1))))))]
                      += (1.0) * i0data[ici0+ci0.size()*(ix7+x7.size()*(ix0+x0.size()*(ix6+x6.size()*(ix1+x1.size()*(ix4+x4.size()*(ix3))))))] * fdata[ix4+x4.size()*(ix3)];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x2 == x3 && x4 == x5) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x7, x0, x6, x1);
    for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
      for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
        for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
          for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
            for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix7+x7.size()*(ix0+x0.size()*(ix6+x6.size()*(ix5+x5.size()*(ix3+x2.size()*(ix1))))))]
                    += (1.0) * i0data[ici0+ci0.size()*(ix7+x7.size()*(ix0+x0.size()*(ix6+x6.size()*(ix1))))] * fdata[ix5+x4.size()*(ix3)];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x2 == x3) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x7, x0, x6, x5, x4, x1);
    for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
      for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
        for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
          for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
            for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
              for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
                for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
                  for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                    odata[ici0+ci0.size()*(ix7+x7.size()*(ix0+x0.size()*(ix6+x6.size()*(ix5+x5.size()*(ix3+x2.size()*(ix1))))))]
                      += (1.0) * i0data[ici0+ci0.size()*(ix7+x7.size()*(ix0+x0.size()*(ix6+x6.size()*(ix5+x5.size()*(ix4+x4.size()*(ix1))))))] * fdata[ix4+x4.size()*(ix3)];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x7, x0, x6, x3, x2, x1);
    for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
            for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
              for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
                for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                  for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                    odata[ici0+ci0.size()*(ix7+x7.size()*(ix0+x0.size()*(ix6+x6.size()*(ix5+x5.size()*(ix2+x2.size()*(ix1))))))]
                      += (1.0) * i0data[ici0+ci0.size()*(ix7+x7.size()*(ix0+x0.size()*(ix6+x6.size()*(ix3+x3.size()*(ix2+x2.size()*(ix1))))))] * fdata[ix5+x4.size()*(ix3)];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(ci0, x7, x0, x6, x5, x2, x1);
    for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
          for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
            for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
              for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix7+x7.size()*(ix0+x0.size()*(ix6+x6.size()*(ix5+x5.size()*(ix2+x2.size()*(ix1))))))]
                    += (1.0) * i0data[ici0+ci0.size()*(ix7+x7.size()*(ix0+x0.size()*(ix6+x6.size()*(ix5+x5.size()*(ix2+x2.size()*(ix1))))))];
                }
              }
            }
          }
        }
      }
    }
  }
  out()->put_block(odata, ci0, x7, x0, x6, x5, x2, x1);
}


void Task7::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);
  // tensor label: Gamma18
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x4, x3, x2, x1);
  {
    if (x2 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x1);
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix4+x4.size()*(ix3+x3.size()*(ix3+x2.size()*(ix1))))))]
                    += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix4+x4.size()*(ix1))))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x0, x4, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,6,1,1,1,1>(i0data, odata, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x0, x4, x3, x2, x1);
}


void Task8::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x0 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);
  // tensor label: Gamma24
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x3, x0, x2, x1);
  {
    if (x2 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x1, x3, x0);
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix3+x3.size()*(ix0+x0.size()*(ix4+x2.size()*(ix1))))))]
                    += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix1+x1.size()*(ix3+x3.size()*(ix0))))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x2, x1);
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
          for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix4+x3.size()*(ix0+x0.size()*(ix2+x2.size()*(ix1))))))]
                    += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix2+x2.size()*(ix1))))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x4, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,5,6,1,1,1,1>(i0data, odata, ci0.size(), x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x4, x3, x0, x2, x1);
}


void Task9::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  // tensor label: Gamma27
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}


void Task10::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x2 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x1 = b(5);
  const Index x0 = b(6);
  // tensor label: Gamma26
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x2, x4, x3, x1, x0);
  {
    if (x1 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x3);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix2+x2.size()*(ix4+x4.size()*(ix3+x3.size()*(ix2+x1.size()*(ix0))))))]
                    += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix4+x4.size()*(ix3))))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x2, x4, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix2+x2.size()*(ix4+x4.size()*(ix3+x3.size()*(ix3+x1.size()*(ix0))))))]
                    += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix2+x2.size()*(ix4+x4.size()*(ix0))))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x2, x4, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,6,1,1,1,1>(i0data, odata, ci0.size(), x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x5, x2, x4, x3, x1, x0);
}


void Task11::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);

  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(x1, x2, x0, a1);
  {
    // tensor label: I0
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x2, x1, a1);
    sort_indices<2,1,0,3,1,1,1,1>(i0data, odata, x0.size(), x2.size(), x1.size(), a1.size());
  }
  out()->put_block(odata, x1, x2, x0, a1);
}


void Task12::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);

  // tensor label: I0
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a1), 0.0);

  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a1, x6, x5);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a1, x6, x5)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x7.size(), a1.size(), x6.size(), x5.size());

        // tensor label: I1
        std::unique_ptr<double[]> i1data = in(1)->get_block(x7, x0, x6, x5, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x7, x0, x6, x5, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i1data, i1data_sorted, x7.size(), x0.size(), x6.size(), x5.size(), x2.size(), x1.size());

        dgemm_("T", "N", a1.size(), x0.size()*x2.size()*x1.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }

  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a1);
}


void Task13::Task_local::compute() {
  const Index x7 = b(0);
  const Index x0 = b(1);
  const Index x6 = b(2);
  const Index x5 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);

  // tensor label: I1
  std::unique_ptr<double[]> odata = out()->move_block(x7, x0, x6, x5, x2, x1);
  {
    // tensor label: Gamma0
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x5, x2, x1);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x7.size(), x0.size(), x6.size(), x5.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x7, x0, x6, x5, x2, x1);
}


void Task14::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);

  // tensor label: I0
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a1), 0.0);

  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a2, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a2.size(), a1.size());

    // tensor label: I3
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x2, x1, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x2, x1, a2)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x0.size(), x2.size(), x1.size(), a2.size());

    dgemm_("T", "N", a1.size(), x0.size()*x2.size()*x1.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, a1.size());
  }

  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a1);
}


void Task15::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);

  // tensor label: I3
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a2), 0.0);

  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, x4, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), x4.size(), x3.size());

        // tensor label: I4
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());

        dgemm_("T", "N", a2.size(), x0.size()*x2.size()*x1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }

  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a2);
}


void Task16::Task_local::compute() {
  const Index x5 = b(0);
  const Index x0 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);

  // tensor label: I4
  std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x3, x2, x1);
  {
    // tensor label: Gamma1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x5, x0, x4, x3, x2, x1);
}


void Task17::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);

  // tensor label: I0
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a1), 0.0);

  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), x3.size());

        // tensor label: I6
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());

        dgemm_("T", "N", a1.size(), x0.size()*x2.size()*x1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }

  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a1);
}


void Task18::Task_local::compute() {
  const Index x5 = b(0);
  const Index x0 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);

  // tensor label: I6
  std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x3, x2, x1);
  {
    // tensor label: Gamma1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
    dscal_(x5.size()*x0.size()*x4.size()*x3.size()*x2.size()*x1.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x5, x0, x4, x3, x2, x1);
}


void Task19::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);

  // tensor label: I0
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a1), 0.0);

  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), x3.size());

        // tensor label: I8
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());

        dgemm_("T", "N", a1.size(), x0.size()*x2.size()*x1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }

  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a1);
}


void Task20::Task_local::compute() {
  const Index x5 = b(0);
  const Index x0 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);

  // tensor label: I8
  std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x3, x2, x1);
  {
    // tensor label: Gamma1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x5, x0, x4, x3, x2, x1);
}


void Task21::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);

  // tensor label: I0
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a1), 0.0);

  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, a1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), a1.size());

        // tensor label: I10
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x3, x0, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x3, x0, x2, x1)]);
        sort_indices<0,1,2,3,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());

        dgemm_("T", "N", a1.size(), x0.size()*x2.size()*x1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }

  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a1);
}


void Task22::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);

  // tensor label: I10
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, x3, x0, x2, x1);
  {
    // tensor label: Gamma4
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x5, x4, x3, x0, x2, x1);
}


void Task23::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);

  // tensor label: I0
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a1), 0.0);

  for (auto& x3 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size());

    // tensor label: I12
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

    dgemm_("T", "N", a1.size(), x0.size()*x2.size()*x1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a1.size());
  }

  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a1);
}


void Task24::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);

  // tensor label: I12
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma5
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task25::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a1 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);

  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
  sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());

  // tensor label: I14
  std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x2, x1, a1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x2, x1, a1)]);
  sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x2.size(), x1.size(), a1.size());

  energy_ += ddot_(x0.size()*x2.size()*x1.size()*a1.size(), i0data_sorted, 1, i1data_sorted, 1);
}


void Task26::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);

  // tensor label: I14
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a1), 0.0);

  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a1, x6, x5);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a1, x6, x5)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x7.size(), a1.size(), x6.size(), x5.size());

        // tensor label: I15
        std::unique_ptr<double[]> i1data = in(1)->get_block(x7, x0, x6, x5, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x7, x0, x6, x5, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i1data, i1data_sorted, x7.size(), x0.size(), x6.size(), x5.size(), x2.size(), x1.size());

        dgemm_("T", "N", a1.size(), x0.size()*x2.size()*x1.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }

  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a1);
}


void Task27::Task_local::compute() {
  energy_ = 0.0;
  const Index x7 = b(0);
  const Index x0 = b(1);
  const Index x6 = b(2);
  const Index x5 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);

  // tensor label: I15
  std::unique_ptr<double[]> odata = out()->move_block(x7, x0, x6, x5, x2, x1);
  {
    // tensor label: Gamma0
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x5, x2, x1);
    sort_indices<0,1,2,3,4,5,1,1,1,4>(i0data, odata, x7.size(), x0.size(), x6.size(), x5.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x7, x0, x6, x5, x2, x1);
}


void Task28::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);

  // tensor label: I14
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a1), 0.0);

  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a2, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a2.size(), a1.size());

    // tensor label: I18
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x2, x1, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x2, x1, a2)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x0.size(), x2.size(), x1.size(), a2.size());

    dgemm_("T", "N", a1.size(), x0.size()*x2.size()*x1.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, a1.size());
  }

  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a1);
}


void Task29::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);

  // tensor label: I18
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a2), 0.0);

  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, x4, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), x4.size(), x3.size());

        // tensor label: I19
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());

        dgemm_("T", "N", a2.size(), x0.size()*x2.size()*x1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }

  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a2);
}


void Task30::Task_local::compute() {
  energy_ = 0.0;
  const Index x5 = b(0);
  const Index x0 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);

  // tensor label: I19
  std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x3, x2, x1);
  {
    // tensor label: Gamma1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,1,1,1,4>(i0data, odata, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x5, x0, x4, x3, x2, x1);
}


void Task31::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);

  // tensor label: I14
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a1), 0.0);

  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), x3.size());

        // tensor label: I22
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());

        dgemm_("T", "N", a1.size(), x0.size()*x2.size()*x1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }

  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a1);
}


void Task32::Task_local::compute() {
  energy_ = 0.0;
  const Index x5 = b(0);
  const Index x0 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);

  // tensor label: I22
  std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x3, x2, x1);
  {
    // tensor label: Gamma1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
    dscal_(x5.size()*x0.size()*x4.size()*x3.size()*x2.size()*x1.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,4,5,1,1,-1,4>(i0data, odata, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x5, x0, x4, x3, x2, x1);
}


void Task33::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);

  // tensor label: I14
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a1), 0.0);

  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), x3.size());

        // tensor label: I25
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());

        dgemm_("T", "N", a1.size(), x0.size()*x2.size()*x1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }

  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a1);
}


void Task34::Task_local::compute() {
  energy_ = 0.0;
  const Index x5 = b(0);
  const Index x0 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);

  // tensor label: I25
  std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x3, x2, x1);
  {
    // tensor label: Gamma1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,1,1,1,2>(i0data, odata, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x5, x0, x4, x3, x2, x1);
}


void Task35::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);

  // tensor label: I14
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a1), 0.0);

  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, a1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), a1.size());

        // tensor label: I28
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x3, x0, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x3, x0, x2, x1)]);
        sort_indices<0,1,2,3,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());

        dgemm_("T", "N", a1.size(), x0.size()*x2.size()*x1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }

  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a1);
}


void Task36::Task_local::compute() {
  energy_ = 0.0;
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);

  // tensor label: I28
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, x3, x0, x2, x1);
  {
    // tensor label: Gamma4
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,5,1,1,1,2>(i0data, odata, x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x5, x4, x3, x0, x2, x1);
}


void Task37::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);

  // tensor label: I14
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a1), 0.0);

  for (auto& x3 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size());

    // tensor label: I31
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

    dgemm_("T", "N", a1.size(), x0.size()*x2.size()*x1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a1.size());
  }

  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a1);
}


void Task38::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);

  // tensor label: I31
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma5
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task39::Task_local::compute(){
  correction_ = 0.0;
  const Index x0 = b(0);
  const Index a1 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);

  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
  sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());

  // tensor label: I33
  std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x2, x1, a1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x2, x1, a1)]);
  sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x2.size(), x1.size(), a1.size());

  correction_ += ddot_(x0.size()*x2.size()*x1.size()*a1.size(), i0data_sorted, 1, i1data_sorted, 1);
}


void Task40::Task_local::compute(){
  correction_ = 0.0;
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);

  // tensor label: I33
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a1), 0.0);

  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), x3.size());

        // tensor label: I34
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());

        dgemm_("T", "N", a1.size(), x0.size()*x2.size()*x1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }

  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a1);
}


void Task41::Task_local::compute(){
  correction_ = 0.0;
  const Index x5 = b(0);
  const Index x0 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);

  // tensor label: I34
  std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x3, x2, x1);
  {
    // tensor label: Gamma1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,1,1,1,4>(i0data, odata, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x5, x0, x4, x3, x2, x1);
}


void Task43::Task_local::compute() {
  const Index x3 = b(0);
  const Index x4 = b(1);

  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(x3, x4);
  {
    // tensor label: I35
    std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x3);
    sort_indices<1,0,1,1,1,1>(i0data, odata, x4.size(), x3.size());
  }
  out()->put_block(odata, x3, x4);
}


void Task44::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);

  // tensor label: I35
  std::unique_ptr<double[]> odata = out()->move_block(x4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x3), 0.0);

  for (auto& x0 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& x1 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());

          // tensor label: I36
          std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x4, x3, x2, x1, a1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x4, x3, x2, x1, a1)]);
          sort_indices<0,5,4,3,1,2,0,1,1,1>(i1data, i1data_sorted, x0.size(), x4.size(), x3.size(), x2.size(), x1.size(), a1.size());

          dgemm_("T", "N", 1, x4.size()*x3.size(), x0.size()*x2.size()*x1.size()*a1.size(),
                 1.0, i0data_sorted, x0.size()*x2.size()*x1.size()*a1.size(), i1data_sorted, x0.size()*x2.size()*x1.size()*a1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }

  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x4.size(), x3.size());
  out()->put_block(odata, x4, x3);
}


void Task45::Task_local::compute() {
  const Index x0 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  const Index a1 = b(5);

  // tensor label: I36
  std::unique_ptr<double[]> odata = out()->move_block(x0, x4, x3, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x4, x3, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x4, x3, x2, x1, a1), 0.0);

  for (auto& x5 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x7 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a1, x6, x5);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a1, x6, x5)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x7.size(), a1.size(), x6.size(), x5.size());

        // tensor label: I37
        std::unique_ptr<double[]> i1data = in(1)->get_block(x7, x0, x6, x5, x4, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x7, x0, x6, x5, x4, x3, x2, x1)]);
        sort_indices<3,2,0,1,4,5,6,7,0,1,1,1>(i1data, i1data_sorted, x7.size(), x0.size(), x6.size(), x5.size(), x4.size(), x3.size(), x2.size(), x1.size());

        dgemm_("T", "N", a1.size(), x0.size()*x4.size()*x3.size()*x2.size()*x1.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }

  sort_indices<1,2,3,4,5,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x4, x3, x2, x1, a1);
}


void Task46::Task_local::compute() {
  const Index x7 = b(0);
  const Index x0 = b(1);
  const Index x6 = b(2);
  const Index x5 = b(3);
  const Index x4 = b(4);
  const Index x3 = b(5);
  const Index x2 = b(6);
  const Index x1 = b(7);

  // tensor label: I37
  std::unique_ptr<double[]> odata = out()->move_block(x7, x0, x6, x5, x4, x3, x2, x1);
  {
    // tensor label: Gamma13
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x5, x4, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,6,7,1,1,1,4>(i0data, odata, x7.size(), x0.size(), x6.size(), x5.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x7, x0, x6, x5, x4, x3, x2, x1);
}


void Task47::Task_local::compute() {
  const Index a1 = b(0);
  const Index a2 = b(1);

  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(a1, a2);
  {
    // tensor label: I38
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, a1);
    sort_indices<1,0,1,1,1,1>(i0data, odata, a2.size(), a1.size());
  }
  out()->put_block(odata, a1, a2);
}


void Task48::Task_local::compute() {
  const Index a2 = b(0);
  const Index a1 = b(1);

  // tensor label: I38
  std::unique_ptr<double[]> odata = out()->move_block(a2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, a1), 0.0);

  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());

        // tensor label: I39
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x2, x1, a2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x2, x1, a2)]);
        sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), x2.size(), x1.size(), a2.size());

        dgemm_("T", "N", a1.size(), a2.size(), x0.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x0.size()*x2.size()*x1.size(), i1data_sorted, x0.size()*x2.size()*x1.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }

  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size());
  out()->put_block(odata, a2, a1);
}


void Task49::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);

  // tensor label: I39
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a2), 0.0);

  for (auto& x3 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, x4, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, x4, x3)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), x4.size(), x3.size());

        // tensor label: I40
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<3,2,0,1,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());

        dgemm_("T", "N", a2.size(), x0.size()*x2.size()*x1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }

  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a2);
}


void Task50::Task_local::compute() {
  const Index x5 = b(0);
  const Index x0 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);

  // tensor label: I40
  std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x3, x2, x1);
  {
    // tensor label: Gamma1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,1,1,1,4>(i0data, odata, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x5, x0, x4, x3, x2, x1);
}


void Task52::Task_local::compute() {
  const Index x0 = b(0);
  const Index a1 = b(1);

  // tensor label: den1
  std::unique_ptr<double[]> odata = out()->move_block(x0, a1);
  {
    // tensor label: I41
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x0.size(), a1.size());
  }
  out()->put_block(odata, x0, a1);
}


void Task53::Task_local::compute() {
  const Index x0 = b(0);
  const Index a1 = b(1);

  // tensor label: I41
  std::unique_ptr<double[]> odata = out()->move_block(x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a1), 0.0);

  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, x1)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), x1.size());

        // tensor label: I42
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
        sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

        dgemm_("T", "N", a1.size(), x0.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }

  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size());
  out()->put_block(odata, x0, a1);
}


void Task54::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);

  // tensor label: I42
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma5
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}


void Task56::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);

  // tensor label: Den1
  std::unique_ptr<double[]> odata = out()->move_block(x1, x2, x0, a1);
  {
    // tensor label: I43
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x2, x1, a1);
    sort_indices<2,1,0,3,1,1,1,1>(i0data, odata, x0.size(), x2.size(), x1.size(), a1.size());
  }
  out()->put_block(odata, x1, x2, x0, a1);
}


void Task57::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);

  // tensor label: I43
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a1), 0.0);

  for (auto& x3 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, x3)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), x3.size());

        // tensor label: I44
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<3,2,0,1,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());

        dgemm_("T", "N", a1.size(), x0.size()*x2.size()*x1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }

  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a1);
}


void Task58::Task_local::compute() {
  const Index x5 = b(0);
  const Index x0 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);

  // tensor label: I44
  std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x3, x2, x1);
  {
    // tensor label: Gamma1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,1,1,1,2>(i0data, odata, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x5, x0, x4, x3, x2, x1);
}


void Task60::Task_local::compute() {
  const Index ci0 = b(0);

  // tensor label: deci
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  {
    // tensor label: I45
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0);
    sort_indices<0,1,1,1,1>(i0data, odata, ci0.size());
  }
  out()->put_block(odata, ci0);
}


void Task61::Task_local::compute() {
  const Index ci0 = b(0);

  // tensor label: I45
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);

  for (auto& x0 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& x1 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());

          // tensor label: I46
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x2, x1, a1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x2, x1, a1)]);
          sort_indices<1,4,3,2,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x2.size(), x1.size(), a1.size());

          dgemm_("T", "N", 1, ci0.size(), x0.size()*x2.size()*x1.size()*a1.size(),
                 1.0, i0data_sorted, x0.size()*x2.size()*x1.size()*a1.size(), i1data_sorted, x0.size()*x2.size()*x1.size()*a1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }

  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}


void Task62::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a1 = b(4);

  // tensor label: I46
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x2, x1, a1), 0.0);

  for (auto& x5 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x7 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a1, x6, x5);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a1, x6, x5)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x7.size(), a1.size(), x6.size(), x5.size());

        // tensor label: I47
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x7, x0, x6, x5, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x7, x0, x6, x5, x2, x1)]);
        sort_indices<4,3,1,0,2,5,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x7.size(), x0.size(), x6.size(), x5.size(), x2.size(), x1.size());

        dgemm_("T", "N", a1.size(), ci0.size()*x0.size()*x2.size()*x1.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }

  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x0, x2, x1, a1);
}


void Task63::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x7 = b(1);
  const Index x0 = b(2);
  const Index x6 = b(3);
  const Index x5 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);

  // tensor label: I47
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x7, x0, x6, x5, x2, x1);
  {
    // tensor label: Gamma17
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x7, x0, x6, x5, x2, x1);
    sort_indices<0,1,2,3,4,5,6,1,1,1,4>(i0data, odata, ci0.size(), x7.size(), x0.size(), x6.size(), x5.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x7, x0, x6, x5, x2, x1);
}


void Task64::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a1 = b(4);

  // tensor label: I46
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x2, x1, a1), 0.0);

  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a2, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a2.size(), a1.size());

    // tensor label: I50
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x2, x1, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x2, x1, a2)]);
    sort_indices<4,0,1,2,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x2.size(), x1.size(), a2.size());

    dgemm_("T", "N", a1.size(), ci0.size()*x0.size()*x2.size()*x1.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, a1.size());
  }

  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x0, x2, x1, a1);
}


void Task65::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a2 = b(4);

  // tensor label: I50
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x2, x1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x2, x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x2, x1, a2), 0.0);

  for (auto& x3 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, x4, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, x4, x3)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), x4.size(), x3.size());

        // tensor label: I51
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x0, x4, x3, x2, x1)]);
        sort_indices<4,3,1,0,2,5,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());

        dgemm_("T", "N", a2.size(), ci0.size()*x0.size()*x2.size()*x1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }

  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a2.size(), ci0.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x0, x2, x1, a2);
}


void Task66::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);

  // tensor label: I51
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x4, x3, x2, x1);
  {
    // tensor label: Gamma18
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,6,1,1,1,4>(i0data, odata, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x0, x4, x3, x2, x1);
}


void Task67::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a1 = b(4);

  // tensor label: I46
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x2, x1, a1), 0.0);

  for (auto& x3 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, x3)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), x3.size());

        // tensor label: I61
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x0, x4, x3, x2, x1)]);
        sort_indices<4,3,1,0,2,5,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());

        dgemm_("T", "N", a1.size(), ci0.size()*x0.size()*x2.size()*x1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }

  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x0, x2, x1, a1);
}


void Task68::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);

  // tensor label: I61
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x4, x3, x2, x1);
  {
    // tensor label: Gamma18
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x3, x2, x1);
    dscal_(ci0.size()*x5.size()*x0.size()*x4.size()*x3.size()*x2.size()*x1.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,4,5,6,1,1,-1,4>(i0data, odata, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x0, x4, x3, x2, x1);
}


void Task69::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a1 = b(4);

  // tensor label: I46
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x2, x1, a1), 0.0);

  for (auto& x3 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, x3)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), x3.size());

        // tensor label: I67
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x0, x4, x3, x2, x1)]);
        sort_indices<4,3,1,0,2,5,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());

        dgemm_("T", "N", a1.size(), ci0.size()*x0.size()*x2.size()*x1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }

  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x0, x2, x1, a1);
}


void Task70::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);

  // tensor label: I67
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x4, x3, x2, x1);
  {
    // tensor label: Gamma18
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,6,1,1,1,2>(i0data, odata, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x0, x4, x3, x2, x1);
}


void Task71::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a1 = b(4);

  // tensor label: I46
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x2, x1, a1), 0.0);

  for (auto& x3 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, a1)]);
        sort_indices<2,1,0,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), a1.size());

        // tensor label: I70
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x4, x3, x0, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x4, x3, x0, x2, x1)]);
        sort_indices<3,2,1,0,4,5,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());

        dgemm_("T", "N", a1.size(), ci0.size()*x0.size()*x2.size()*x1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }

  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x0, x2, x1, a1);
}


void Task72::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x0 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);

  // tensor label: I70
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x3, x0, x2, x1);
  {
    // tensor label: Gamma24
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x4, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,5,6,1,1,1,2>(i0data, odata, ci0.size(), x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x4, x3, x0, x2, x1);
}


void Task73::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a1 = b(4);

  // tensor label: I46
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x2, x1, a1), 0.0);

  for (auto& x3 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size());

    // tensor label: I79
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
    sort_indices<1,0,2,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());

    dgemm_("T", "N", a1.size(), ci0.size()*x0.size()*x2.size()*x1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a1.size());
  }

  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x0, x2, x1, a1);
}


void Task74::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);

  // tensor label: I79
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma27
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}


void Task75::Task_local::compute() {
  const Index ci0 = b(0);

  // tensor label: I45
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);

  for (auto& x7 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& x6 : *range_[1]) {
        for (auto& x5 : *range_[1]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a1, x6, x5);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a1, x6, x5)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), a1.size(), x6.size(), x5.size());

          // tensor label: I53
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x7, x6, x5, a1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x7, x6, x5, a1)]);
          sort_indices<1,4,2,3,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x7.size(), x6.size(), x5.size(), a1.size());

          dgemm_("T", "N", 1, ci0.size(), x7.size()*x6.size()*x5.size()*a1.size(),
                 1.0, i0data_sorted, x7.size()*x6.size()*x5.size()*a1.size(), i1data_sorted, x7.size()*x6.size()*x5.size()*a1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }

  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}


void Task76::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x7 = b(1);
  const Index x6 = b(2);
  const Index x5 = b(3);
  const Index a1 = b(4);

  // tensor label: I53
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x7, x6, x5, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x7, x6, x5, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x7, x6, x5, a1), 0.0);

  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());

        // tensor label: I54
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x7, x0, x6, x5, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x7, x0, x6, x5, x2, x1)]);
        sort_indices<5,6,2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x7.size(), x0.size(), x6.size(), x5.size(), x2.size(), x1.size());

        dgemm_("T", "N", a1.size(), ci0.size()*x7.size()*x6.size()*x5.size(), x0.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x0.size()*x2.size()*x1.size(), i1data_sorted, x0.size()*x2.size()*x1.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }

  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x7.size(), x6.size(), x5.size());
  out()->put_block(odata, ci0, x7, x6, x5, a1);
}


void Task77::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x7 = b(1);
  const Index x0 = b(2);
  const Index x6 = b(3);
  const Index x5 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);

  // tensor label: I54
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x7, x0, x6, x5, x2, x1);
  {
    // tensor label: Gamma17
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x7, x0, x6, x5, x2, x1);
    sort_indices<0,1,2,3,4,5,6,1,1,1,4>(i0data, odata, ci0.size(), x7.size(), x0.size(), x6.size(), x5.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x7, x0, x6, x5, x2, x1);
}


void Task78::Task_local::compute() {
  const Index ci0 = b(0);

  // tensor label: I45
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);

  for (auto& x5 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& x4 : *range_[1]) {
        for (auto& x3 : *range_[1]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, x3)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), x3.size());

          // tensor label: I56
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x4, x3, a1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x4, x3, a1)]);
          sort_indices<1,4,2,3,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x4.size(), x3.size(), a1.size());

          dgemm_("T", "N", 1, ci0.size(), x5.size()*x4.size()*x3.size()*a1.size(),
                 1.0, i0data_sorted, x5.size()*x4.size()*x3.size()*a1.size(), i1data_sorted, x5.size()*x4.size()*x3.size()*a1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }

  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}


void Task79::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index a1 = b(4);

  // tensor label: I56
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x5, x4, x3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x5, x4, x3, a1), 0.0);

  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a2, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a2.size(), a1.size());

    // tensor label: I57
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x4, x3, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x4, x3, a2)]);
    sort_indices<4,0,1,2,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x4.size(), x3.size(), a2.size());

    dgemm_("T", "N", a1.size(), ci0.size()*x5.size()*x4.size()*x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, a1.size());
  }

  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x5.size(), x4.size(), x3.size());
  out()->put_block(odata, ci0, x5, x4, x3, a1);
}


void Task80::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index a2 = b(4);

  // tensor label: I57
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x5, x4, x3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x5, x4, x3, a2), 0.0);

  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a2, x1, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a2, x1, x2)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), a2.size(), x1.size(), x2.size());

        // tensor label: I58
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x0, x4, x3, x2, x1)]);
        sort_indices<5,6,2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());

        dgemm_("T", "N", a2.size(), ci0.size()*x5.size()*x4.size()*x3.size(), x0.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x0.size()*x2.size()*x1.size(), i1data_sorted, x0.size()*x2.size()*x1.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }

  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a2.size(), ci0.size(), x5.size(), x4.size(), x3.size());
  out()->put_block(odata, ci0, x5, x4, x3, a2);
}


void Task81::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);

  // tensor label: I58
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x4, x3, x2, x1);
  {
    // tensor label: Gamma18
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,6,1,1,1,4>(i0data, odata, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x0, x4, x3, x2, x1);
}


void Task82::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index a1 = b(4);

  // tensor label: I56
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x5, x4, x3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x5, x4, x3, a1), 0.0);

  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());

        // tensor label: I64
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x0, x4, x3, x2, x1)]);
        sort_indices<5,6,2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());

        dgemm_("T", "N", a1.size(), ci0.size()*x5.size()*x4.size()*x3.size(), x0.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x0.size()*x2.size()*x1.size(), i1data_sorted, x0.size()*x2.size()*x1.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }

  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x5.size(), x4.size(), x3.size());
  out()->put_block(odata, ci0, x5, x4, x3, a1);
}


void Task83::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);

  // tensor label: I64
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x4, x3, x2, x1);
  {
    // tensor label: Gamma18
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x3, x2, x1);
    dscal_(ci0.size()*x5.size()*x0.size()*x4.size()*x3.size()*x2.size()*x1.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,4,5,6,1,1,-1,4>(i0data, odata, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x0, x4, x3, x2, x1);
}


void Task84::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index a1 = b(4);

  // tensor label: I56
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x5, x4, x3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x5, x4, x3, a1), 0.0);

  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());

        // tensor label: I73
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x0, x4, x3, x2, x1)]);
        sort_indices<5,6,2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());

        dgemm_("T", "N", a1.size(), ci0.size()*x5.size()*x4.size()*x3.size(), x0.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x0.size()*x2.size()*x1.size(), i1data_sorted, x0.size()*x2.size()*x1.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }

  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x5.size(), x4.size(), x3.size());
  out()->put_block(odata, ci0, x5, x4, x3, a1);
}


void Task85::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);

  // tensor label: I73
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x4, x3, x2, x1);
  {
    // tensor label: Gamma18
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,6,1,1,1,2>(i0data, odata, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x0, x4, x3, x2, x1);
}


void Task86::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index a1 = b(4);

  // tensor label: I56
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x5, x4, x3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x5, x4, x3, a1), 0.0);

  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, x2, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1, x2, a1)]);
        sort_indices<2,1,0,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size(), x2.size(), a1.size());

        // tensor label: I76
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x2, x4, x3, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x2, x4, x3, x1, x0)]);
        sort_indices<2,5,6,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());

        dgemm_("T", "N", a1.size(), ci0.size()*x5.size()*x4.size()*x3.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }

  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x5.size(), x4.size(), x3.size());
  out()->put_block(odata, ci0, x5, x4, x3, a1);
}


void Task87::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x2 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x1 = b(5);
  const Index x0 = b(6);

  // tensor label: I76
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x2, x4, x3, x1, x0);
  {
    // tensor label: Gamma26
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x2, x4, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,6,1,1,1,2>(i0data, odata, ci0.size(), x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x5, x2, x4, x3, x1, x0);
}


void Task88::Task_local::compute() {
  const Index ci0 = b(0);

  // tensor label: I45
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);

  for (auto& x3 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& x2 : *range_[1]) {
        for (auto& x1 : *range_[1]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, x1);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, x1)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), x1.size());

          // tensor label: I81
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, a1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, a1)]);
          sort_indices<1,4,2,3,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), a1.size());

          dgemm_("T", "N", 1, ci0.size(), x3.size()*x2.size()*x1.size()*a1.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x1.size()*a1.size(), i1data_sorted, x3.size()*x2.size()*x1.size()*a1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }

  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}


void Task89::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a1 = b(4);

  // tensor label: I81
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, x1, a1), 0.0);

  for (auto& x0 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size());

    // tensor label: I82
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());

    dgemm_("T", "N", a1.size(), ci0.size()*x3.size()*x2.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a1.size());
  }

  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x3.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x3, x2, x1, a1);
}


void Task90::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);

  // tensor label: I82
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma27
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}


