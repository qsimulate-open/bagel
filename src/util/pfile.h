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
    std::fstream* file_;
    long filesize_;
    std::string filename_;

    const int K_;

  public:
    PFile(const long, const int); 
    ~PFile(); 

    void add_block(const long, const long, const T*);
    void put_block(const long, const long, const T*);
    void get_block(const long, const long, T*);
    void clear();

    const std::string filename() const { return filename_; };
    const int K() const { return K_; };
};


template<class T>
PFile<T>::PFile(const long fsize, const int k) : filesize_(fsize), K_(k) {

  {
    Filename tmpf;
    filename_ = tmpf.filename_next();
  }

  const T czero = static_cast<T>(0.0);
  T* work = (T*) work_char;
  std::fill(work, work + cachesize, czero);

  file_ = new std::fstream(filename_.c_str(), std::ios::binary | std::ios::out);
  long remaining = filesize_;
  while (remaining > 0L) {
    const size_t writesize = std::min(cachesize, remaining) * sizeof(T);
    file_->write((const char*)work, writesize);
    remaining -= cachesize;
  }
  file_->close();

  // reopen with in/out/binary
  file_->open(filename_.c_str(), std::ios::binary | std::ios::in | std::ios::out);

};


template<class T>
PFile<T>::~PFile() {
  delete file_;
  unlink(filename_.c_str());
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

