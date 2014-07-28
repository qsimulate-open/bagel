/*
 * storage_traits.h
 *
 *  Created on: Dec 27, 2013
 *      Author: evaleev
 */

#ifndef BTAS_STORAGE_TRAITS_H_
#define BTAS_STORAGE_TRAITS_H_

#include <btas/index_traits.h>

namespace btas {

  template <typename _Storage>
  struct storage_traits;

  template <typename _T>
  struct storage_traits<_T*> {
      typedef typename std::remove_const<_T>::type value_type;
      typedef _T* pointer;
      typedef typename std::add_const<pointer>::type const_pointer;
      typedef          value_type& reference;
      typedef    const value_type& const_reference;
      typedef size_t size_type;
      typedef ptrdiff_t difference_type;

      typedef pointer iterator;
      typedef const_pointer const_iterator;
  };

  template <typename _T>
  struct storage_traits<_T* const> {
      typedef typename std::remove_const<_T>::type value_type;
      typedef _T* pointer;
      typedef typename std::add_const<pointer>::type const_pointer;
      typedef          value_type& reference;
      typedef    const value_type& const_reference;
      typedef size_t size_type;
      typedef ptrdiff_t difference_type;

      typedef pointer iterator;
      typedef const_pointer const_iterator;
  };

  template <typename _Storage>
  struct storage_traits {
      typedef typename _Storage::value_type value_type;
      typedef typename _Storage::reference reference;
      typedef typename _Storage::const_reference const_reference;
      typedef typename _Storage::size_type size_type;
      typedef typename _Storage::difference_type difference_type;

      typedef typename _Storage::iterator iterator;
      typedef typename _Storage::const_iterator const_iterator;
  };

  /// test _Storage conforms the TWG.Storage concept
  /// in addition to Storage, check extent() member and extent_type
  template<class _Storage>
  class is_storage {
  public:
     static constexpr const bool
     value = has_begin<_Storage>::value &
             has_end<_Storage>::value;
  };


}


#endif /* BTAS_STORAGE_TRAITS_H_ */
