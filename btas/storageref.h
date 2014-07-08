/*
 * storageref.h
 *
 *  Created on: Dec 27, 2013
 *      Author: evaleev
 */

#ifndef BTAS_STORAGEREF_H_
#define BTAS_STORAGEREF_H_

#include <btas/storage_traits.h>

namespace btas {

  /// Wraps pointer to \c _Storage in _Storage-like interface
  /// If \c _Storage is const, will provide read-access only
  template <typename _Storage>
  class StorageRef {
      struct Enabler {};
    public:
      using storage_type = _Storage;
      using nonconst_storage_type = typename std::remove_const<storage_type>::type;
      using value_type = typename storage_type::value_type;
      using size_type = typename storage_type::size_type;
      using iterator = typename std::conditional<std::is_const<storage_type>::value,
                                        typename storage_type::const_iterator,
                                        typename storage_type::iterator>::type;
      using const_iterator = typename storage_type::const_iterator;

      StorageRef() : begin_(), end_() { end_ = begin_; }

      ~StorageRef() {}

      template <typename S = _Storage>
      StorageRef(S& stor)
                 : begin_(stor.begin()), end_(stor.end())
                 { }

      // TODO dangerous hack
      template <typename S = _Storage>
      StorageRef(const S& stor)
                 : begin_(const_cast<S*>(&stor)->begin()), end_(const_cast<S*>(&stor)->end())
                 { }

      template <typename Iter1, typename Iter2>
      StorageRef(Iter1 b, Iter2 e) : begin_(b), end_(e) {}

      StorageRef(const StorageRef& other) : begin_(other.begin_), end_(other.end_) {}

      StorageRef& operator=(const StorageRef& other) {
        begin_ = other.begin_;
        end_ = other.end_;
        return *this;
      }

      template <typename S = _Storage>
      typename std::enable_if<not std::is_const<S>::value,StorageRef&>::type
      operator=(nonconst_storage_type& stor) {
        begin_ = std::begin(stor);
        end_ = std::end(stor);
        return *this;
      }

      template <typename S = _Storage>
      typename std::enable_if<std::is_const<S>::value,StorageRef&>::type
      operator=(const nonconst_storage_type& stor) {
        begin_ = stor.cbegin();
        end_ = stor.cend();
        return *this;
      }

      // no point in move functionality

      template <typename S = _Storage>
      typename std::enable_if<not std::is_const<S>::value,value_type&>::type
      operator[](size_t i) {
        return *(begin_ + i);
      }
      const value_type& operator[](size_t i) const {
        return *(begin_ + i);
      }

      template <typename S = _Storage>
      typename std::enable_if<not std::is_const<S>::value,iterator>::type
      begin() {
        return begin_;
      }
      const_iterator begin() const {
        return cbegin();
      }
      const_iterator cbegin() const {
        return begin_;
      }

      template <typename S = _Storage>
      typename std::enable_if<not std::is_const<S>::value,iterator>::type
      end() {
        return end_;
      }
      const_iterator end() const {
        return cend();
      }
      const_iterator cend() const {
        return end_;
      }

      size_type size() const {
        return cend() - cbegin();
      }

    private:
      iterator begin_; // begin
      iterator end_;   // end
  }; // StorageRef

}


#endif /* BTAS_STORAGEREF_H_ */
