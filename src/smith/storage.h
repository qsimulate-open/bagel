//
// Author : Toru Shiozaki
// Date   : Feb 2012
//
//
// this is a base class of storage object
// ... It can be distributed memory, distributed disk, etc.
// ... All of them should be derived from here with the same interface.

#ifndef __SRC_SMITH_STORAGE_H
#define __SRC_SMITH_STORAGE_H

#include <map>
#include <memory>
#include <tuple>
#include <stdexcept>

namespace SMITH {

class Storage_base {
  protected:
    // length of this storage
    size_t length_;
    // this relates hash keys, offsets, and block lengths (in this order).
    std::map<size_t, std::pair<size_t, size_t> > hashtable_;
    

  public:
    Storage_base(const std::map<size_t, size_t>& size) {
      length_ = 0lu;
      for (auto i = size.begin(); i != size.end(); ++i) {
        auto j = hashtable_.insert(std::make_pair(i->first, std::make_pair(length_, i->second)));
        if (!j.second) throw std::logic_error("duplicated hash keys in Storage::Storage"); 
        length_ += i->second; 
      }
    };
    ~Storage_base() {};

    // functions that return protected members
    size_t length() const { return length_; };

    // get, put, and add a block from the storage and returns unique_ptr<double[]>, which is local
    virtual std::unique_ptr<double[]> get_block(const size_t& key) const = 0;
    virtual void put_block(const size_t& key, const std::unique_ptr<double[]>& dat) = 0;
    virtual void add_block(const size_t& key, const std::unique_ptr<double[]>& dat) = 0;
};

class Storage_Incore : public Storage_base {
  protected:
    std::unique_ptr<double[]> data_;

  public:
    Storage_Incore(const std::map<size_t, size_t>& size) : Storage_base(size) {
      std::unique_ptr<double[]> tmp(new double[length()]);
      data_ = std::move(tmp);
    };
    ~Storage_Incore() {};

    std::unique_ptr<double[]> get_block(const size_t& key) const;
    void put_block(const size_t& key, const std::unique_ptr<double[]>& dat);
    void add_block(const size_t& key, const std::unique_ptr<double[]>& dat);

};

}

#endif
