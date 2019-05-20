/*
 * ordinal.h
 *
 *  Created on: Dec 26, 2013
 *      Author: evaleev
 */

#ifndef BTAS_ORDINAL_H_
#define BTAS_ORDINAL_H_

#include <btas/types.h>
#include <btas/array_adaptor.h>
#include <btas/index_traits.h>
#include <btas/varray/varray.h>
#include <boost/assert.hpp>

namespace btas {

  /// BoxOrdinal is an implementation detail of BoxRange.
  /// It maps the index to its ordinal value. It also knows whether
  /// the map is contiguous (i.e. whether adjacent indices have adjacent ordinal
  /// values).
  template <CBLAS_ORDER _Order = CblasRowMajor,
            typename _Index = btas::varray<long>,
            class = typename std::enable_if<btas::is_index<_Index>::value>
           >
  class BoxOrdinal {
    public:
      typedef _Index index_type;
      const static CBLAS_ORDER order = _Order;
      typedef int64_t value_type;
      typedef typename btas::replace_value_type<_Index,value_type>::type stride_type;    ///< stride type

      template <CBLAS_ORDER _O,
                typename _I,
                class _X
               >
      friend class BoxOrdinal;

      BoxOrdinal() {
        BOOST_ASSERT((contiguous_ = false) || true); // workaround for Boost serialization
                                                     // it breaks Debug builds when reading uninitialized bools
      }

      template <typename Index1,
                typename Index2,
                class = typename std::enable_if<btas::is_index<Index1>::value && btas::is_index<Index2>::value>::type
               >
      BoxOrdinal(const Index1& lobound,
                 const Index2& upbound) {
          init(lobound, upbound);
      }

      template <typename Index1,
                typename Index2,
                typename Weight,
                class = typename std::enable_if<btas::is_index<Index1>::value &&
                                                btas::is_index<Index2>::value &&
                                                btas::is_index<Weight>::value>::type
               >
      BoxOrdinal(const Index1& lobound,
                 const Index2& upbound,
                 const Weight& stride) {
          init(lobound, upbound, stride);
      }

      BoxOrdinal(stride_type&& stride,
                 value_type&& offset,
                 bool cont) : stride_(stride), offset_(offset), contiguous_(cont) {
      }


      BoxOrdinal(const BoxOrdinal& other) :
        stride_(other.stride_),
        offset_(other.offset_),
        contiguous_ (other.contiguous_) {
      }

      template <CBLAS_ORDER _O,
                typename _I,
                class = typename std::enable_if<btas::is_index<_I>::value>
               >
      BoxOrdinal(const BoxOrdinal<_O,_I>& other) {
          auto n = other.rank();
          stride_ = array_adaptor<stride_type>::construct(n);
          std::copy(std::begin(other.stride_), std::end(other.stride_),
                    std::begin(stride_));
          offset_ = other.offset_;
          contiguous_ = other.contiguous_;
      }

      ~BoxOrdinal() {}

      std::size_t rank() const {
        using btas::rank;
        return rank(stride_);
      }

      const stride_type& stride() const {
        return stride_;
      }

      value_type offset() const {
        return offset_;
      }

      bool contiguous() const {
        return contiguous_;
      }

      BoxOrdinal& operator=(const BoxOrdinal&) = default;
      BoxOrdinal& operator=(BoxOrdinal&&) = default;

      template <typename Index>
      typename std::enable_if<btas::is_index<Index>::value, value_type>::type
      operator()(const Index& index) const {
        value_type o = 0;
        const auto end = this->rank();
        for(auto i = 0ul; i != end; ++i)
          o += *(std::begin(index) + i) * *(std::begin(this->stride_) + i);

        return o - offset_;
      }

      /// Does ordinal value belong to this ordinal range?
      template <typename I>
      typename std::enable_if<std::is_integral<I>::value, bool>::type
      includes(const I& ord) const {
        assert(false); // "BoxOrdinal::includes() is not not yet implemented"
      }

    private:

      template <typename Index1,
                typename Index2,
                class = typename std::enable_if<btas::is_index<Index1>::value && btas::is_index<Index2>::value>::type
               >
      void init(const Index1& lobound,
                const Index2& upbound) {
        using btas::rank;
        auto n = rank(lobound);
        if (n == 0) return;

        value_type volume = 1;
        offset_ = 0;
        stride_ = array_adaptor<stride_type>::construct(n);

        // Compute range data
        if (order == CblasRowMajor) {
          for(int i = n - 1; i >= 0; --i) {
            stride_[i] = volume;
            auto li = *(std::begin(lobound) + i);
            auto ui = *(std::begin(upbound) + i);
            offset_ += li * volume;
            volume *= (ui - li);
          }
        }
        else {
          for(auto i = 0; i != n; ++i) {
            stride_[i] = volume;
            auto li = *(std::begin(lobound) + i);
            auto ui = *(std::begin(upbound) + i);
            offset_ += li * volume;
            volume *= (ui - li);
          }
        }
        contiguous_ = true;
      }

      /// upbound only needed to check contiguousness
      template <typename Index1,
                typename Index2,
                typename Weight,
                class = typename std::enable_if<btas::is_index<Index1>::value &&
                                                btas::is_index<Index2>::value &&
                                                btas::is_index<Weight>::value>::type
               >
      void init(const Index1& lobound,
                const Index2& upbound,
                const Weight& stride) {
        using btas::rank;
        auto n = rank(lobound);
        if (n == 0) return;

        value_type volume = 1;
        offset_ = 0;
        stride_ = array_adaptor<stride_type>::construct(n);
        std::copy(std::begin(stride), std::end(stride), std::begin(stride_));

        // Compute offset and check whether contiguous
        contiguous_ = true;
        stride_type tmpstride = array_adaptor<stride_type>::construct(n);
        if (order == CblasRowMajor) {
          for(int i = n - 1; i >= 0; --i) {
            tmpstride[i] = volume;
            contiguous_ &= (tmpstride[i] == stride_[i]);
            auto li = *(std::begin(lobound) + i);
            auto ui = *(std::begin(upbound) + i);
            offset_ += li * stride_[i];
            volume *= (ui - li);
          }
        }
        else {
          for(auto i = 0; i != n; ++i) {
            tmpstride[i] = volume;
            contiguous_ &= (tmpstride[i] == stride_[i]);
            auto li = *(std::begin(lobound) + i);
            auto ui = *(std::begin(upbound) + i);
            offset_ += li * stride_[i];
            volume *= (ui - li);
          }
        }
      }

      friend class boost::serialization::access;
      template<class Archive>
      void serialize(Archive& ar, const unsigned int version) {
        ar & BOOST_SERIALIZATION_NVP(stride_) & BOOST_SERIALIZATION_NVP(offset_) & BOOST_SERIALIZATION_NVP(contiguous_);
      }

      stride_type stride_; //!< stride of each dimension (stride in the language of NumPy)
      value_type offset_; //!< lobound . stride so that easy to compute ordinal: ordinal(index) = index . stride - offset
      bool contiguous_; //!< whether index iterator traverses a contiguous sequence of ordinals
  };

  /// Permutes BoxOrdinal

  /// permutes the dimensions using permutation \c p = {p[0], p[1], ... }; for example, if \c stride() initially returned
  /// {s[0], s[1], ... }, after this call \c stride() will return {s[p[0]], s[p[1]], ...}.
  /// \param perm an array specifying permutation of the dimensions
  template <CBLAS_ORDER _Order,
            typename _Index,
            typename AxisPermutation,
            class = typename std::enable_if<btas::is_index<AxisPermutation>::value>::type>
  BoxOrdinal<_Order, _Index> permute(const BoxOrdinal<_Order, _Index>& ord,
                                     const AxisPermutation& perm)
  {
    const auto rank = ord.rank();
    auto st = ord.stride();

    typedef typename BoxOrdinal<_Order, _Index>::stride_type stride_type;
    stride_type stride;
    stride = array_adaptor<stride_type>::construct(rank);

    std::for_each(std::begin(perm), std::end(perm), [&](const typename AxisPermutation::value_type& i){
      const auto pi = *(std::begin(perm) + i);
      *(std::begin(stride)+i) = *(std::begin(st) + pi);
    });

    return BoxOrdinal<_Order, _Index>(std::move(stride), ord.offset(), ord.contiguous());
  }

  /// Range output operator

  /// \param os The output stream that will be used to print \c r
  /// \param r The range to be printed
  /// \return A reference to the output stream
  template <CBLAS_ORDER _Order,
            typename _Index>
  inline std::ostream& operator<<(std::ostream& os, const BoxOrdinal<_Order,_Index>& ord) {
    array_adaptor<typename BoxOrdinal<_Order,_Index>::stride_type>::print(ord.stride(), os);
    return os;
  }

} // namespace btas

#endif /* BTAS_ORDINAL_H_ */
