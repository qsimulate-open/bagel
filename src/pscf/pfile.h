//
// BAGEL - Parallel electron correlation program.
// Filename: pfile.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#ifndef __src_util_pfile_h
#define __src_util_pfile_h

#include <stddef.h>
#include <fstream>
#include <string>
#include <complex>
#include <sstream>
#include <algorithm>
#include <src/pscf/filename.h>
#include <src/pscf/cache.h>

namespace bagel {

template<class T>
class PFile {
  protected:
    std::shared_ptr<std::fstream> file_;
    long filesize_;
    std::string filename_;

    const int K_;

  public:
    PFile(const long, const int, const bool late_init = false);
    ~PFile();

    void append(const long, const T*);
    void add_block(const long, const long, const T*);
    void put_block(const long, const long, const T*);
    void get_block(const long, const long, T*) const;
    void clear();
    void reopen_with_inout();

    const std::string filename() const { return filename_; };
    int K() const { return K_; };
};


template<class T>
PFile<T>::PFile(const long fsize, const int k, const bool late_init) : filesize_(fsize), K_(k) {

  {
    Filename tmpf;
    filename_ = tmpf.filename_next();
  }

  std::shared_ptr<std::fstream> tmp(new std::fstream(filename_.c_str(), std::ios::out | std::ios::trunc | std::ios::binary));
  file_ = tmp;

  if (!late_init) {
    const T czero = static_cast<T>(0.0);
    T* work = (T*) work_char;
    std::fill(work, work + cachesize, czero);

    long remaining = filesize_;
    while (remaining > 0L) {
      const size_t writesize = std::min(cachesize, remaining) * sizeof(T);
      file_->write((const char*)work, writesize);
      remaining -= cachesize;
    }

    // reopen with in/out/binary
    reopen_with_inout();
  }

};


#include <unistd.h>
template<class T>
PFile<T>::~PFile() {
  unlink(filename_.c_str());
};


template<class T>
void PFile<T>::reopen_with_inout() {
  file_->close();
  file_->open(filename_.c_str(), std::ios::in | std::ios::out | std::ios::binary);
};


template<class T>
void PFile<T>::append(const long length, const T* data) {
  file_->clear();
  file_->seekp(0, std::ios::end);
  long remaining = length;
  long current = 0L;

  while (remaining > 0L) {
    const long writesize = std::min(remaining, cachesize) * sizeof(T);
    file_->write((const char*)(data + current), writesize);

    remaining -= cachesize;
    current += cachesize;
  }
};


template<class T>
void PFile<T>::add_block(const long position, const long length, const T* data) {
  long remaining = length;
  long current = 0L;
  T* work = (T*) work_char;

  while (remaining > 0L) {
    const size_t rsize = std::min(remaining, cachesize);
    const size_t readsize = rsize * sizeof(T);

    file_->clear();
    file_->seekg((position + current) * sizeof(T));
    file_->read((char*)work, readsize);
    for (size_t i = 0; i != rsize; ++i) work[i] += data[current + i];

    file_->clear();
    file_->seekp((position + current) * sizeof(T));
    file_->write((const char*)work, readsize);

    remaining -= cachesize;
    current += cachesize;
  }
};


template<class T>
void PFile<T>::get_block(const long position, const long length, T* data) const {
  long remaining = length;
  long current = 0L;

  while (remaining > 0L) {
    const long readsize = std::min(remaining, cachesize) * sizeof(T);

    file_->clear();
    file_->seekg((position + current) * sizeof(T));
    file_->read((char*)(data + current), readsize);

    remaining -= cachesize;
    current += cachesize;
  }
};


template<class T>
void PFile<T>::put_block(const long position, const long length, const T* data) {
  long remaining = length;
  long current = 0L;

  while (remaining > 0L) {
    const long writesize = std::min(remaining, cachesize) * sizeof(T);
    file_->clear();
    file_->seekp((position + current) * sizeof(T));
    file_->write((const char*)(data + current), writesize);

    remaining -= cachesize;
    current += cachesize;
  }
};


template<class T>
void PFile<T>::clear() {
  long remaining = filesize_;
  long current = 0L;
  T zero = static_cast<T>(0.0);
  T* work = (T*) work_char;
  std::fill(work, work + cachesize, zero);

  while (remaining > 0L) {
    const long writesize = std::min(remaining, cachesize) * sizeof(T);
    file_->clear();
    file_->seekp(current * sizeof(T));
    file_->write((const char*)work, writesize);

    remaining -= cachesize;
    current += cachesize;
  }
};

}

#endif

