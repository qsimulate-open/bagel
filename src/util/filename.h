//
// Author : Toru Shiozaki
// Date   : August 2009
//

#ifndef __src_util_filename_h
#define __src_util_filename_h

#include <string>

class Filename {
  private:

  public:
    Filename();
    ~Filename();

    const std::string filename_next() const;
};

#endif
