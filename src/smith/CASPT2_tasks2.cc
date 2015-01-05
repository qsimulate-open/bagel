//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks2.cc
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


#include <src/smith/CASPT2_tasks2.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task50::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  const Index x2 = b(4);
  // tensor label: Gamma294
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x1, x0, x2);
  {
    if (x0 == x1) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
              odata[ici0+ci0.size()*(ix3+x3.size()*(ix1+x1.size()*(ix1+x0.size()*(ix2))))]
                += (1.0) * i0data[ici0+ci0.size()*(ix3+x3.size()*(ix2))];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x1);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
              odata[ici0+ci0.size()*(ix3+x3.size()*(ix1+x1.size()*(ix2+x0.size()*(ix2))))]
                += (-2.0) * i0data[ici0+ci0.size()*(ix3+x3.size()*(ix1))];
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x3, x1, x0, x2);
    sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x1.size(), x0.size(), x2.size());
  }
  out()->put_block(odata, ci0, x3, x1, x0, x2);
}

void Task51::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x1 = b(3);
  const Index x3 = b(4);
  const Index x2 = b(5);
  const Index x0 = b(6);
  // tensor label: Gamma300
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x1, x3, x2, x0);
  {
    if (x2 == x4 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix3+x1.size()*(ix3+x3.size()*(ix4+x2.size()*(ix0))))))]
                  += (-2.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x4, x2, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix3+x1.size()*(ix3+x3.size()*(ix2+x2.size()*(ix0))))))]
                    += (-2.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix2+x2.size()*(ix0))))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix4+x1.size()*(ix3+x3.size()*(ix3+x2.size()*(ix0))))))]
                  += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x3, x2, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix4+x1.size()*(ix3+x3.size()*(ix2+x2.size()*(ix0))))))]
                    += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix3+x3.size()*(ix2+x2.size()*(ix0))))];
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
      std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x4, x1, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix1+x1.size()*(ix3+x3.size()*(ix3+x2.size()*(ix0))))))]
                    += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix1+x1.size()*(ix0))))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x0, x1, x3);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix1+x1.size()*(ix3+x3.size()*(ix4+x2.size()*(ix0))))))]
                    += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix1+x1.size()*(ix3))))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(ci0, x5, x4, x1, x3, x2, x0);
    sort_indices<0,1,2,3,4,5,6,1,1,1,1>(i0data, odata, ci0.size(), x5.size(), x4.size(), x1.size(), x3.size(), x2.size(), x0.size());
  }
  out()->put_block(odata, ci0, x5, x4, x1, x3, x2, x0);
}

void Task52::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index x0 = b(4);
  // tensor label: Gamma301
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x3, x2, x0);
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x2, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
              odata[ici0+ci0.size()*(ix3+x1.size()*(ix3+x3.size()*(ix2+x2.size()*(ix0))))]
                += (-2.0) * i0data[ici0+ci0.size()*(ix2+x2.size()*(ix0))];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
            for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
              odata[ici0+ci0.size()*(ix1+x1.size()*(ix3+x3.size()*(ix3+x2.size()*(ix0))))]
                += (1.0) * i0data[ici0+ci0.size()*(ix1+x1.size()*(ix0))];
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x1, x3, x2, x0);
    sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x1.size(), x3.size(), x2.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x3, x2, x0);
}

void Task53::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  const Index x4 = b(4);
  const Index x3 = b(5);
  const Index x2 = b(6);
  // tensor label: Gamma303
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x1, x4);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(3)->get_block(x3, x2);
  if (x1 == x4) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x0, x3, x2);
    for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
      for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
        for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
          for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix4+x1.size()*(ix4))))]
                  += (2.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix3+x3.size()*(ix2))))] * fdata[ix3+x3.size()*(ix2)];
              }
            }
          }
        }
      }
    }
  }
  if (x3 == x4 && x1 == x2) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
        for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
          for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
            for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
              odata[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix2+x1.size()*(ix4))))]
                += (2.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0))] * fdata[ix4+x3.size()*(ix2)];
            }
          }
        }
      }
    }
  }
  if (x1 == x2) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x0, x3, x4);
    for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
      for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
        for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
          for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
            for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix2+x1.size()*(ix4))))]
                  += (-1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix3+x3.size()*(ix4))))] * fdata[ix3+x3.size()*(ix2)];
              }
            }
          }
        }
      }
    }
  }
  if (x3 == x4) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x0, x1, x2);
    for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
          for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix1+x1.size()*(ix4))))]
                  += (-1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix1+x1.size()*(ix2))))] * fdata[ix4+x3.size()*(ix2)];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(ci0, x5, x0, x1, x4, x3, x2);
    for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
      for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
        for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
          for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
            for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix1+x1.size()*(ix4))))]
                    += (-1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix1+x1.size()*(ix4+x4.size()*(ix3+x3.size()*(ix2))))))] * fdata[ix3+x3.size()*(ix2)];
                }
              }
            }
          }
        }
      }
    }
  }
  out()->put_block(odata, ci0, x5, x0, x1, x4);
}

void Task54::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  const Index x2 = b(4);
  // tensor label: Gamma304
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x1, x2);
  {
    if (x1 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
              odata[ici0+ci0.size()*(ix3+x3.size()*(ix0+x0.size()*(ix2+x1.size()*(ix2))))]
                += (2.0) * i0data[ici0+ci0.size()*(ix3+x3.size()*(ix0))];
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x3, x0, x1, x2);
    sort_indices<0,1,2,3,4,1,1,-1,1>(i0data, odata, ci0.size(), x3.size(), x0.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, ci0, x3, x0, x1, x2);
}

void Task55::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  const Index x3 = b(5);
  const Index x2 = b(6);
  // tensor label: Gamma306
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x1, x0);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(3)->get_block(x3, x2);
  if (x1 == x4) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x0, x3, x2);
    for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
      for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
        for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
          for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix4+x1.size()*(ix0))))]
                  += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix3+x3.size()*(ix2))))] * fdata[ix3+x3.size()*(ix2)];
              }
            }
          }
        }
      }
    }
  }
  if (x3 == x4 && x1 == x2) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
        for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
          for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
            for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
              odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix2+x1.size()*(ix0))))]
                += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0))] * fdata[ix4+x3.size()*(ix2)];
            }
          }
        }
      }
    }
  }
  if (x1 == x2) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x4, x3, x0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
        for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
          for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
            for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix2+x1.size()*(ix0))))]
                  += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix3+x3.size()*(ix0))))] * fdata[ix3+x3.size()*(ix2)];
              }
            }
          }
        }
      }
    }
  }
  if (x3 == x4) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x2, x1, x0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
          for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix1+x1.size()*(ix0))))]
                  += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix2+x2.size()*(ix1+x1.size()*(ix0))))] * fdata[ix4+x3.size()*(ix2)];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(ci0, x5, x4, x3, x2, x1, x0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix1+x1.size()*(ix0))))]
                    += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix3+x3.size()*(ix2+x2.size()*(ix1+x1.size()*(ix0))))))] * fdata[ix3+x3.size()*(ix2)];
                }
              }
            }
          }
        }
      }
    }
  }
  out()->put_block(odata, ci0, x5, x4, x1, x0);
}

void Task56::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: Gamma307
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    if (x1 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
              odata[ici0+ci0.size()*(ix3+x3.size()*(ix2+x2.size()*(ix2+x1.size()*(ix0))))]
                += (1.0) * i0data[ici0+ci0.size()*(ix3+x3.size()*(ix0))];
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task57::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x1 = b(5);
  const Index x2 = b(6);
  // tensor label: Gamma309
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x4, x3, x1, x2);
  {
    if (x1 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x3);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix4+x4.size()*(ix3+x3.size()*(ix2+x1.size()*(ix2))))))]
                    += (2.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix4+x4.size()*(ix3))))];
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
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x2);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix4+x4.size()*(ix3+x3.size()*(ix3+x1.size()*(ix2))))))]
                    += (-1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix4+x4.size()*(ix2))))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x0, x4, x3, x1, x2);
    sort_indices<0,1,2,3,4,5,6,1,1,-1,1>(i0data, odata, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, ci0, x5, x0, x4, x3, x1, x2);
}

void Task58::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: Gamma310
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task59::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x2 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x1 = b(5);
  const Index x0 = b(6);
  // tensor label: Gamma323
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

void Task60::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  const Index x4 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);
  // tensor label: Gamma328
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x3, x4, x2, x1);
  {
    if (x2 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x3, x1);
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix3+x3.size()*(ix4+x4.size()*(ix4+x2.size()*(ix1))))))]
                    += (-1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix3+x3.size()*(ix1))))];
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
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix4+x3.size()*(ix4+x4.size()*(ix2+x2.size()*(ix1))))))]
                    += (2.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix2+x2.size()*(ix1))))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x0, x3, x4, x2, x1);
    sort_indices<0,1,2,3,4,5,6,1,1,-1,1>(i0data, odata, ci0.size(), x5.size(), x0.size(), x3.size(), x4.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x0, x3, x4, x2, x1);
}

void Task61::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x0 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);
  // tensor label: Gamma329
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

void Task62::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x7 = b(1);
  const Index x0 = b(2);
  const Index x6 = b(3);
  const Index x5 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);
  const Index x4 = b(7);
  const Index x3 = b(8);
  // tensor label: Gamma330
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
  if (x4 == x5 && x2 == x3) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x7, x0, x6, x1);
    for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
      for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
        for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
          for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
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

void Task63::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);
  // tensor label: Gamma331
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

void Task64::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  // tensor label: Gamma332
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}

void Task65::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: Gamma341
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(1)->get_block(x1, x0);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
          odata[ici0]
            += (1.0) * i0data[ici0+ci0.size()*(ix1+x1.size()*(ix0))] * fdata[ix1+x1.size()*(ix0)];
        }
      }
    }
  }
  out()->put_block(odata, ci0);
}

void Task66::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  // tensor label: Gamma353
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(1)->get_block(x2, x1);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
              odata[ici0+ci0.size()*(ix3+x3.size()*(ix0))]
                += (1.0) * i0data[ici0+ci0.size()*(ix3+x3.size()*(ix0+x0.size()*(ix2+x2.size()*(ix1))))] * fdata[ix2+x2.size()*(ix1)];
            }
          }
        }
      }
    }
  }
  out()->put_block(odata, ci0, x3, x0);
}

void Task67::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x1 = b(4);
  const Index x3 = b(5);
  const Index x2 = b(6);
  // tensor label: Gamma364
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x4, x1);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(1)->get_block(x3, x2);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x1, x3, x2);
    for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
      for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
        for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix4+x4.size()*(ix1))))]
                    += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix4+x4.size()*(ix1+x1.size()*(ix3+x3.size()*(ix2))))))] * fdata[ix3+x3.size()*(ix2)];
                }
              }
            }
          }
        }
      }
    }
  }
  out()->put_block(odata, ci0, x5, x0, x4, x1);
}

void Task68::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x2 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x1 = b(5);
  const Index x0 = b(6);
  // tensor label: Gamma379
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x2, x5, x4, x3, x1, x0);
  {
    if (x2 == x5 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x4, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix5+x2.size()*(ix5+x5.size()*(ix4+x4.size()*(ix3+x3.size()*(ix3+x1.size()*(ix0))))))]
                  += (2.0) * i0data[ici0+ci0.size()*(ix4+x4.size()*(ix0))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x4, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix3+x2.size()*(ix5+x5.size()*(ix4+x4.size()*(ix3+x3.size()*(ix5+x1.size()*(ix0))))))]
                  += (-1.0) * i0data[ici0+ci0.size()*(ix4+x4.size()*(ix0))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x4, x3, x2, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix2+x2.size()*(ix5+x5.size()*(ix4+x4.size()*(ix3+x3.size()*(ix5+x1.size()*(ix0))))))]
                    += (-1.0) * i0data[ici0+ci0.size()*(ix4+x4.size()*(ix3+x3.size()*(ix2+x2.size()*(ix0))))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x4, x3, x1, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x2.size()*(ix5+x5.size()*(ix4+x4.size()*(ix3+x3.size()*(ix1+x1.size()*(ix0))))))]
                    += (2.0) * i0data[ici0+ci0.size()*(ix4+x4.size()*(ix3+x3.size()*(ix1+x1.size()*(ix0))))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x2, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
            for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix2+x2.size()*(ix5+x5.size()*(ix5+x4.size()*(ix3+x3.size()*(ix3+x1.size()*(ix0))))))]
                  += (-1.0) * i0data[ici0+ci0.size()*(ix2+x2.size()*(ix0))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x2, x5, x4, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix2+x2.size()*(ix5+x5.size()*(ix4+x4.size()*(ix3+x3.size()*(ix3+x1.size()*(ix0))))))]
                    += (-1.0) * i0data[ici0+ci0.size()*(ix2+x2.size()*(ix5+x5.size()*(ix4+x4.size()*(ix0))))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5 && x2 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix3+x2.size()*(ix5+x5.size()*(ix5+x4.size()*(ix3+x3.size()*(ix1+x1.size()*(ix0))))))]
                  += (2.0) * i0data[ici0+ci0.size()*(ix1+x1.size()*(ix0))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x4, x5, x1, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix3+x2.size()*(ix5+x5.size()*(ix4+x4.size()*(ix3+x3.size()*(ix1+x1.size()*(ix0))))))]
                    += (-1.0) * i0data[ici0+ci0.size()*(ix4+x4.size()*(ix5+x5.size()*(ix1+x1.size()*(ix0))))];
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
      std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x2, x3, x1, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix2+x2.size()*(ix5+x5.size()*(ix5+x4.size()*(ix3+x3.size()*(ix1+x1.size()*(ix0))))))]
                    += (-1.0) * i0data[ici0+ci0.size()*(ix2+x2.size()*(ix3+x3.size()*(ix1+x1.size()*(ix0))))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(ci0, x2, x5, x4, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,6,1,1,-1,1>(i0data, odata, ci0.size(), x2.size(), x5.size(), x4.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x2, x5, x4, x3, x1, x0);
}

void Task69::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  const Index x3 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);
  // tensor label: Gamma397
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x0, x3, x2, x1);
  {
    if (x2 == x3 && x0 == x1) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x4);
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix1+x0.size()*(ix3+x3.size()*(ix3+x2.size()*(ix1))))))]
                  += (2.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix4))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4 && x0 == x1) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x3);
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix1+x0.size()*(ix3+x3.size()*(ix4+x2.size()*(ix1))))))]
                  += (-1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix3))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x1) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x4, x2, x3);
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix1+x0.size()*(ix3+x3.size()*(ix2+x2.size()*(ix1))))))]
                    += (-1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix2+x2.size()*(ix3))))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x1);
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix3+x0.size()*(ix3+x3.size()*(ix4+x2.size()*(ix1))))))]
                  += (2.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix1))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x4, x2, x1);
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix3+x0.size()*(ix3+x3.size()*(ix2+x2.size()*(ix1))))))]
                    += (2.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix2+x2.size()*(ix1))))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x1);
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix4+x0.size()*(ix3+x3.size()*(ix3+x2.size()*(ix1))))))]
                  += (-1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix1))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x3, x2, x1);
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix4+x0.size()*(ix3+x3.size()*(ix2+x2.size()*(ix1))))))]
                    += (-1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix3+x3.size()*(ix2+x2.size()*(ix1))))];
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
      std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x4, x0, x1);
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix0+x0.size()*(ix3+x3.size()*(ix3+x2.size()*(ix1))))))]
                    += (-1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix0+x0.size()*(ix1))))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x1, x0, x3);
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix0+x0.size()*(ix3+x3.size()*(ix4+x2.size()*(ix1))))))]
                    += (-1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix1+x1.size()*(ix0+x0.size()*(ix3))))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(ci0, x5, x4, x0, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,6,1,1,-1,1>(i0data, odata, ci0.size(), x5.size(), x4.size(), x0.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x4, x0, x3, x2, x1);
}

void Task71::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index x0 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, c1, x0);
  {
    // tensor label: I0
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c2, x0, x1);
    sort_indices<1,3,0,2,1,1,1,1>(i0data, odata, c1.size(), c2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, c2, x1, c1, x0);
}

void Task72::Task_local::compute() {
  const Index c1 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I0
  std::unique_ptr<double[]> odata = out()->move_block(c1, c2, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c2, x0, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: Gamma0
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x5, x1, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x5, x1, x4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), x5.size(), x1.size(), x4.size());
      // tensor label: I1
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x5, c2, x4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x5, c2, x4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), x5.size(), c2.size(), x4.size());
      dgemm_("T", "N", x0.size()*x1.size(), c1.size()*c2.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c1.size(), c2.size());
  out()->put_block(odata, c1, c2, x0, x1);
}

void Task73::Task_local::compute() {
  const Index c1 = b(0);
  const Index x5 = b(1);
  const Index c2 = b(2);
  const Index x4 = b(3);
  // tensor label: I1
  std::unique_ptr<double[]> odata = out()->move_block(c1, x5, c2, x4);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x5, c2, x4);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, c1.size(), x5.size(), c2.size(), x4.size());
  }
  out()->put_block(odata, c1, x5, c2, x4);
}

void Task74::Task_local::compute() {
  const Index c1 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I0
  std::unique_ptr<double[]> odata = out()->move_block(c1, c2, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c2, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma94
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3, x1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x3, x1, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size(), x1.size(), x2.size());
      // tensor label: I277
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x3, c2, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x3, c2, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), x3.size(), c2.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), c1.size()*c2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c1.size(), c2.size());
  out()->put_block(odata, c1, c2, x0, x1);
}

void Task75::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index c2 = b(2);
  const Index x2 = b(3);
  // tensor label: I277
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, c2, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x3, c2, x2);
    dscal_(c1.size()*x3.size()*c2.size()*x2.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c1.size(), x3.size(), c2.size(), x2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x3, c2, x2);
    sort_indices<0,1,2,3,1,1,2,1>(i1data, odata, c1.size(), x3.size(), c2.size(), x2.size());
  }
  out()->put_block(odata, c1, x3, c2, x2);
}

void Task76::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index x0 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, c1, x0);
  {
    // tensor label: I2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c1, x0, x1);
    sort_indices<0,3,1,2,1,1,1,1>(i0data, odata, c2.size(), c1.size(), x0.size(), x1.size());
  }
  {
    // tensor label: I2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c2, x1, x0);
    sort_indices<1,2,0,3,1,1,1,1>(i0data, odata, c1.size(), c2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, c2, x1, c1, x0);
}

void Task77::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I2
  std::unique_ptr<double[]> odata = out()->move_block(c2, c1, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma94
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3, x1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x3, x1, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size(), x1.size(), x2.size());
      // tensor label: I3
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c1, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c1, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), c1.size(), x3.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), c2.size()*c1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c2.size(), c1.size());
  out()->put_block(odata, c2, c1, x0, x1);
}

void Task78::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I3
  std::unique_ptr<double[]> odata = out()->move_block(c2, c1, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c1, x3, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x3, c3, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x3, c3, x2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x3.size(), c3.size(), x2.size());
    // tensor label: I4
    std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c2.size(), c3.size());
    dgemm_("T", "N", c1.size()*x3.size()*x2.size(), c2.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, c1.size()*x3.size()*x2.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), x3.size(), x2.size(), c2.size());
  out()->put_block(odata, c2, c1, x3, x2);
}

void Task79::Task_local::compute() {
  const Index c2 = b(0);
  const Index c3 = b(1);
  // tensor label: I4
  std::unique_ptr<double[]> odata = out()->move_block(c2, c3);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c3);
    sort_indices<0,1,1,1,-2,1>(i0data, odata, c2.size(), c3.size());
  }
  out()->put_block(odata, c2, c3);
}

void Task80::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I2
  std::unique_ptr<double[]> odata = out()->move_block(c2, c1, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), x2.size());
    // tensor label: I6
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x0, x1, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x0, x1, x2)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), x0.size(), x1.size(), x2.size());
    dgemm_("T", "N", c2.size(), c1.size()*x0.size()*x1.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, c2.size());
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, c2.size(), c1.size(), x0.size(), x1.size());
  out()->put_block(odata, c2, c1, x0, x1);
}

void Task81::Task_local::compute() {
  const Index c1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: I6
  std::unique_ptr<double[]> odata = out()->move_block(c1, x0, x1, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x0, x1, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x0, x1, x2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x0, x3, x1, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x0, x3, x1, x2)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x0.size(), x3.size(), x1.size(), x2.size());
        // tensor label: I7
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, c1, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, c1, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), c1.size(), x3.size());
        dgemm_("T", "N", x0.size()*x1.size()*x2.size(), c1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x0.size()*x1.size()*x2.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), c1.size());
  out()->put_block(odata, c1, x0, x1, x2);
}

void Task82::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I7
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, c1, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c1, x3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x5.size(), x4.size(), c1.size(), x3.size());
  }
  out()->put_block(odata, x5, x4, c1, x3);
}

void Task83::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I2
  std::unique_ptr<double[]> odata = out()->move_block(c2, c1, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma3
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x3, x0, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x3, x0, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), x3.size(), x0.size(), x2.size());
      // tensor label: I9
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c1, c2, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c1, c2, x3)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), c1.size(), c2.size(), x3.size());
      dgemm_("T", "N", x1.size()*x0.size(), c1.size()*c2.size(), x2.size()*x3.size(),
             1.0, i0data_sorted, x2.size()*x3.size(), i1data_sorted, x2.size()*x3.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), c2.size());
  out()->put_block(odata, c2, c1, x0, x1);
}

void Task84::Task_local::compute() {
  const Index x2 = b(0);
  const Index c1 = b(1);
  const Index c2 = b(2);
  const Index x3 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata = out()->move_block(x2, c1, c2, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, c1, c2, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, c1, c2, x3), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, c2, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, c2, x3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), c2.size(), x3.size());
    // tensor label: I10
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x2)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), x2.size());
    dgemm_("T", "N", c1.size()*c2.size()*x3.size(), x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, c1.size()*c2.size()*x3.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), c2.size(), x3.size(), x2.size());
  out()->put_block(odata, x2, c1, c2, x3);
}

void Task85::Task_local::compute() {
  const Index a3 = b(0);
  const Index x2 = b(1);
  // tensor label: I10
  std::unique_ptr<double[]> odata = out()->move_block(a3, x2);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x2);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, a3.size(), x2.size());
  }
  out()->put_block(odata, a3, x2);
}

void Task86::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c1, x2, x0, x1);
  {
    // tensor label: I11
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x2, x1, x0);
    sort_indices<0,1,3,2,1,1,1,1>(i0data, odata, c1.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, c1, x2, x0, x1);
}

void Task87::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I11
  std::unique_ptr<double[]> odata = out()->move_block(c1, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: Gamma4
        std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x5, x3, x4, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x5, x3, x4, x1, x0)]);
        sort_indices<1,2,3,0,4,5,0,1,1,1>(i0data, i0data_sorted, x2.size(), x5.size(), x3.size(), x4.size(), x1.size(), x0.size());
        // tensor label: I12
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c1, x5, x4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c1, x5, x4)]);
        sort_indices<2,0,3,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), c1.size(), x5.size(), x4.size());
        dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x3.size()*x5.size()*x4.size(),
               1.0, i0data_sorted, x3.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x5.size()*x4.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->put_block(odata, c1, x2, x1, x0);
}

void Task88::Task_local::compute() {
  const Index x3 = b(0);
  const Index c1 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  // tensor label: I12
  std::unique_ptr<double[]> odata = out()->move_block(x3, c1, x5, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, c1, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, c1, x5, x4), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x5, c2, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x5, c2, x4)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x5.size(), c2.size(), x4.size());
    // tensor label: I13
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c2)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x3.size(), c2.size());
    dgemm_("T", "N", c1.size()*x5.size()*x4.size(), x3.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, c1.size()*x5.size()*x4.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), x5.size(), x4.size(), x3.size());
  out()->put_block(odata, x3, c1, x5, x4);
}

void Task89::Task_local::compute() {
  const Index x3 = b(0);
  const Index c2 = b(1);
  // tensor label: I13
  std::unique_ptr<double[]> odata = out()->move_block(x3, c2);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, c2);
    sort_indices<0,1,1,1,2,1>(i0data, odata, x3.size(), c2.size());
  }
  out()->put_block(odata, x3, c2);
}

void Task90::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I11
  std::unique_ptr<double[]> odata = out()->move_block(c1, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: Gamma5
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, x2, x5, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, x2, x5, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x2.size(), x5.size(), x1.size(), x0.size());
        // tensor label: I15
        std::unique_ptr<double[]> i1data = in(1)->get_block(x7, x6, c1, x5);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x7, x6, c1, x5)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x7.size(), x6.size(), c1.size(), x5.size());
        dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->put_block(odata, c1, x2, x1, x0);
}

void Task91::Task_local::compute() {
  const Index x7 = b(0);
  const Index x6 = b(1);
  const Index c1 = b(2);
  const Index x5 = b(3);
  // tensor label: I15
  std::unique_ptr<double[]> odata = out()->move_block(x7, x6, c1, x5);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, c1, x5);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x7.size(), x6.size(), c1.size(), x5.size());
  }
  out()->put_block(odata, x7, x6, c1, x5);
}

void Task92::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I11
  std::unique_ptr<double[]> odata = out()->move_block(c1, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma6
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x2, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x2, x3, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
        // tensor label: I17
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x5, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x5, x4, x3)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), x5.size(), x4.size(), x3.size());
        dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->put_block(odata, c1, x2, x1, x0);
}

void Task93::Task_local::compute() {
  const Index c1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I17
  std::unique_ptr<double[]> odata = out()->move_block(c1, x5, x4, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c1, x3);
    dscal_(x5.size()*x4.size()*c1.size()*x3.size(), e0_, i0data.get(), 1);
    sort_indices<2,0,1,3,1,1,-1,1>(i0data, odata, x5.size(), x4.size(), c1.size(), x3.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, c1, x3);
    sort_indices<2,0,1,3,1,1,1,1>(i1data, odata, x5.size(), x4.size(), c1.size(), x3.size());
  }
  out()->put_block(odata, c1, x5, x4, x3);
}

void Task94::Task_local::compute() {
  const Index c1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I17
  std::unique_ptr<double[]> odata = out()->move_block(c1, x5, x4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x5, x4, x3), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c2, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, c2, x3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), c2.size(), x3.size());
    // tensor label: I18
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c2)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), c2.size());
    dgemm_("T", "N", x5.size()*x4.size()*x3.size(), c1.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, x5.size()*x4.size()*x3.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), c1.size());
  out()->put_block(odata, c1, x5, x4, x3);
}

void Task95::Task_local::compute() {
  const Index c1 = b(0);
  const Index c2 = b(1);
  // tensor label: I18
  std::unique_ptr<double[]> odata = out()->move_block(c1, c2);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c2);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, c1.size(), c2.size());
  }
  out()->put_block(odata, c1, c2);
}

void Task96::Task_local::compute() {
  const Index c1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I17
  std::unique_ptr<double[]> odata = out()->move_block(c1, x5, x4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x5, x4, x3), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x5, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x5, x4)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x5.size(), x4.size());
    // tensor label: I30
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x3)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a2.size(), x3.size());
    dgemm_("T", "N", c1.size()*x5.size()*x4.size(), x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, c1.size()*x5.size()*x4.size());
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), x5.size(), x4.size(), x3.size());
  out()->put_block(odata, c1, x5, x4, x3);
}

void Task97::Task_local::compute() {
  const Index a2 = b(0);
  const Index x3 = b(1);
  // tensor label: I30
  std::unique_ptr<double[]> odata = out()->move_block(a2, x3);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, x3);
    sort_indices<0,1,1,1,1,1>(i0data, odata, a2.size(), x3.size());
  }
  out()->put_block(odata, a2, x3);
}

void Task98::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I11
  std::unique_ptr<double[]> odata = out()->move_block(c1, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: Gamma7
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x3, x1, x0)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
    // tensor label: I20
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), x3.size());
    dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, x2.size()*x1.size()*x0.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->put_block(odata, c1, x2, x1, x0);
}

void Task99::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  // tensor label: I20
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x3);
    sort_indices<0,1,1,1,2,1>(i0data, odata, c1.size(), x3.size());
  }
  out()->put_block(odata, c1, x3);
}

