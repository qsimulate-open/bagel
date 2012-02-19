//
// Author : Toru Shiozaki
// Date   : Feb 2012
//

#include <src/util/f77.h>
#include <src/smith/storage.h> 
#include <stdexcept>
#include <iostream>

using namespace SMITH;
using namespace std;


Storage_Incore::Storage_Incore(const map<size_t, size_t>& size) : Storage_base(size) {
  cout << "creating a field of " << length() << endl;
  unique_ptr<double[]> tmp(new double[length()]);
#if 0
  // for debug - initializing tensors with some values
  fill(&tmp[0], &tmp[length()], 1.0);
#endif
  data_ = move(tmp);
}


unique_ptr<double[]> Storage_Incore::get_block(const size_t& key) const {
  // first find a key
  auto hash = hashtable_.find(key);
  if (hash == hashtable_.end())
    throw logic_error("a key was not found in Storage::get_block(const size_t&)");

  // then create a memory 
  const size_t blocksize = hash->second.second;
  const size_t blockoff  = hash->second.first;
  unique_ptr<double[]> buf(new double[blocksize]);

  if (blockoff + blocksize > length_)
    throw logic_error("going beyond boundary in Storage::get_block(const size_t&)");

  // then copy...
  dcopy_(blocksize, &data_[blockoff], 1, &buf[0], 1); 

  return move(buf);
}


void Storage_Incore::put_block(const size_t& key, const unique_ptr<double[]>& dat) {
  auto hash = hashtable_.find(key);
  if (hash == hashtable_.end())
    throw logic_error("a key was not found in Storage::put_block(const size_t&)");

  const size_t blocksize = hash->second.second;
  const size_t blockoff  = hash->second.first;

  if (blockoff + blocksize > length_)
    throw logic_error("going beyond boundary in Storage::put_block(const size_t&)");

  dcopy_(blocksize, &dat[0], 1, &data_[blockoff], 1); 
}


void Storage_Incore::add_block(const size_t& key, const unique_ptr<double[]>& dat) {
  auto hash = hashtable_.find(key);
  if (hash == hashtable_.end())
    throw logic_error("a key was not found in Storage::put_block(const size_t&)");

  const size_t blocksize = hash->second.second;
  const size_t blockoff  = hash->second.first;

  if (blockoff + blocksize > length_)
    throw logic_error("going beyond boundary in Storage::put_block(const size_t&)");

  daxpy_(blocksize, 1.0, &dat[0], 1, &data_[blockoff], 1); 
}

