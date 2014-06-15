/*
 * defaults.h
 *
 *  Created on: Dec 19, 2013
 *      Author: evaleev
 */

#ifndef BTAS_DEFAULTS_H_
#define BTAS_DEFAULTS_H_

//
//  Default types
//

#include <btas/range.h>
#include <btas/varray/varray.h>

namespace btas {

namespace DEFAULT {

/// default storage class
template<typename _T>
using storage = btas::varray<_T>;

/// default range type
using range = btas::Range;

} // namespace DEFAULT

}


#endif /* BTAS_DEFAULTS_H_ */
