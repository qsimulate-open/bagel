//
// Author : Toru Shiozaki
// Date   : August 2009
//

#ifndef __src_util_pfile_h
#define __src_util_pfile_h

#include <fstream>
#include <string>
#include <complex>
#include <sstream>
#include <algorithm>
#include <src/util/filename.h>
#include <src/util/cache.h>

template<class T>
class PFile {
  protected:
    boost::shared_ptr<std::fstream> file_;
    long filesize_;
    std::string filename_;

    const int K_;

  public:
    PFile(const long, const int, const bool late_init = false); 
    ~PFile(); 

    void append(const long, const T*);
    void add_block(const long, const long, const T*);
    void put_block(const long, const long, const T*);
    void get_block(const long, const long, T*);
    void clear();
    void reopen_with_inout();

    const std::string filename() const { return filename_; };
    const int K() const { return K_; };
};


template<class T>
PFile<T>::PFile(const long fsize, const int k, const bool late_init) : filesize_(fsize), K_(k) {

  {
    Filename tmpf;
    filename_ = tmpf.filename_next();
  }

  boost::shared_ptr<std::fstream> tmp(new std::fstream(filename_.c_str(), std::ios::binary | std::ios::out));
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


template<class T>
PFile<T>::~PFile() {
  unlink(filename_.c_str());
};


template<class T>
void PFile<T>::reopen_with_inout() {
  file_->close();
  file_->open(filename_.c_str(), std::ios::binary | std::ios::in | std::ios::out);
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
    const int rsize = std::min(remaining, cachesize);
    const size_t readsize = rsize * sizeof(T);

    file_->clear();
    file_->seekg((position + current) * sizeof(T));
    file_->read((char*)work, readsize); 
    for (int i = 0; i != rsize; ++i) work[i] += data[current + i];

    file_->clear();
    file_->seekp((position + current) * sizeof(T));
    file_->write((const char*)work, readsize);

    remaining -= cachesize;
    current += cachesize;
  } 
};


template<class T>
void PFile<T>::get_block(const long position, const long length, T* data) {
  long remaining = length;
  long current = 0L;

  while (remaining > 0L) {
    const long readsize = std::min(remaining, cachesize) * sizeof(T);

    file_->clear();
    file_->seekg((position + current) * sizeof(T));
    file_->read((char*)(&data[current]), readsize); 

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

#endif

