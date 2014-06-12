/*
 * tensorview_iterator.h
 *
 *  Created on: Dec 28, 2013
 *      Author: evaleev
 */

#ifndef BTAS_TENSORVIEW_ITERATOR_H_
#define BTAS_TENSORVIEW_ITERATOR_H_

#include <cstddef>

#include <btas/storageref.h>

namespace btas {

  /// Iterates over elements of \c Storage using ordinal values of indices in \c Range

  template <typename Range, typename Storage>
  class TensorViewIterator : public std::iterator<typename std::conditional<std::is_const<Storage>::value,
                                                                            std::forward_iterator_tag,
                                                                            std::output_iterator_tag>::type,
                                                  typename std::conditional<std::is_const<Storage>::value,
                                                  const typename StorageRef<Storage>::value_type,
                                                        typename StorageRef<Storage>::value_type>::type>
  {
      struct Enabler {};

    public:
      typedef Storage storage_type;
      typedef StorageRef<storage_type> storageref_type;
      typedef StorageRef<typename std::remove_const<storage_type>::type> ncstorageref_type;
      typedef std::iterator<typename std::conditional<std::is_const<Storage>::value,
          std::forward_iterator_tag,
          std::output_iterator_tag>::type,
          typename std::conditional<std::is_const<Storage>::value,
          const typename StorageRef<Storage>::value_type,
          typename StorageRef<Storage>::value_type>::type> base_type;
      using typename base_type::value_type;
      using typename base_type::pointer;
      using typename base_type::reference;
      using typename base_type::difference_type;
      using typename base_type::iterator_category;

    private:
      typedef typename Range::ordinal_subiterator subiterator;
      typedef typename Range::ordinal_iterator iterator;
      typedef typename iterator::value_type ordinal_type;
      typedef typename Range::index_type index_type;

    public:
      /// Default constructor
      TensorViewIterator() {}
      /// Destructor
      ~TensorViewIterator() {}

      TensorViewIterator(const typename Range::iterator& index_iter,
                         Storage& storage) :
        iter_(subiterator(std::make_pair(*index_iter,index_iter.range()->ordinal(*index_iter)),index_iter.range())),
        storageref_(storage) {}

      TensorViewIterator(const typename Range::iterator& index_iter,
                         const storageref_type& storage) :
        iter_(subiterator(std::make_pair(*index_iter,index_iter.range()->ordinal(*index_iter)),index_iter.range())),
        storageref_(storage) {}

      template <typename S = Storage>
      TensorViewIterator(const typename Range::iterator& index_iter,
                         const ncstorageref_type& storage,
                         typename std::enable_if<std::is_const<S>::value>::type* = 0) :
        iter_(subiterator(std::make_pair(*index_iter,index_iter.range()->ordinal(*index_iter)),index_iter.range())),
        // standard const_cast cannot "map" const into nontrivial structures, have to reinterpret here
        storageref_(reinterpret_cast<const storageref_type&>(storage)) {}

      TensorViewIterator(const typename Range::iterator& index_iter,
                         const ordinal_type& ord,
                         Storage& storage) :
        iter_(subiterator(std::make_pair(*index_iter,ord),index_iter.range())),
        storageref_(storage) {}


      TensorViewIterator(const iterator& iter,
                         Storage& storage) :
        iter_(iter), storageref_(storage) {}

      TensorViewIterator(iterator&& iter,
                         Storage& storage) :
        iter_(iter), storageref_(storage) {}

      TensorViewIterator& operator++() {
        ++iter_;
        return *this;
      }

      const reference operator*() const {
        return *(storageref_.cbegin() + *iter_);
      }

      //template <class = typename std::enable_if<not std::is_const<storage_type>::value,Enabler>::type>
      template <typename S = Storage>
      typename std::enable_if<not std::is_const<S>::value,reference>::type
      operator*() {
        return *(storageref_.begin() + *iter_);
      }

      const index_type& index() const {
        return first(*iter_.base());
      }

      template <typename R, typename S>
      friend bool operator==(const TensorViewIterator<R,S>&, const TensorViewIterator<R,S>&);

    private:
      iterator iter_;
      storageref_type storageref_;
  };

  template <typename Range, typename Storage>
  inline bool operator==(const TensorViewIterator<Range,Storage>& i1,
                         const TensorViewIterator<Range,Storage>& i2) {
    return i1.iter_ == i2.iter_;
  }

  template <typename Range, typename Storage>
  inline bool operator!=(const TensorViewIterator<Range,Storage>& i1,
                         const TensorViewIterator<Range,Storage>& i2) {
    return not (i1 == i2);
  }

}


#endif /* BTAS_TENSORVIEW_ITERATOR_H_ */
