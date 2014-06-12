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
  class storage_traits;

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
