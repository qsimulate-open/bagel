//
// Author : Toru Shiozaki
// Date   : Feb 2012
//

#include <src/util/f77.h>
#include <src/smith/storage.h> 
#include <stdexcept>
#include <iostream>
#include <algorithm>

using namespace SMITH;
using namespace std;


Storage_Incore::Storage_Incore(const map<size_t, size_t>& size) : Storage_base(size) {
  cout << "creating a field of " << length() << endl;
  for (auto i = size.begin(); i != size.end(); ++i) {
    unique_ptr<double[]> tmp(new double[i->second]); 
    data_.push_back(move(tmp));
  }
}


unique_ptr<double[]> Storage_Incore::get_block(const size_t& key) const {
  // first find a key
  auto hash = hashtable_.find(key);
  if (hash == hashtable_.end())
    throw logic_error("a key was not found in Storage::get_block(const size_t&)");

  // then create a memory 
  const size_t blocksize = hash->second.second;
  const size_t blocknum  = hash->second.first;
  unique_ptr<double[]> buf(new double[blocksize]);

  // then copy...
  dcopy_(blocksize, data_.at(blocknum), 1, buf, 1); 

  return move(buf);
}


unique_ptr<double[]> Storage_Incore::move_block(const size_t& key) {
  // first find a key
  auto hash = hashtable_.find(key);
  if (hash == hashtable_.end())
    throw logic_error("a key was not found in Storage::get_block(const size_t&)");
  const size_t blocknum  = hash->second.first;
  return move(data_.at(blocknum));
}


void Storage_Incore::put_block(const size_t& key, unique_ptr<double[]>& dat) {
  auto hash = hashtable_.find(key);
  if (hash == hashtable_.end())
    throw logic_error("a key was not found in Storage::put_block(const size_t&)");
  const size_t blocknum  = hash->second.first;
  data_[blocknum] = move(dat);
}


void Storage_Incore::add_block(const size_t& key, const unique_ptr<double[]>& dat) {
  auto hash = hashtable_.find(key);
  if (hash == hashtable_.end())
    throw logic_error("a key was not found in Storage::put_block(const size_t&)");

  const size_t blocksize = hash->second.second;
  const size_t blocknum  = hash->second.first;

  daxpy_(blocksize, 1.0, dat, 1, data_.at(blocknum), 1); 
}


void Storage_Incore::zero() {
  for (auto i = hashtable_.begin(); i != hashtable_.end(); ++i) {
    const size_t bn = i->second.first;
    const size_t ln = i->second.second;
    fill(data_[bn].get(), data_[bn].get()+ln, 0.0);
  }
}
