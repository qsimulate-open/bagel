/*
 * tensorview.h
 *
 *  Created on: Dec 28, 2013
 *      Author: evaleev
 */

#ifndef BTAS_TENSORVIEW_H_
#define BTAS_TENSORVIEW_H_

#include <btas/tensorview_iterator.h>

namespace btas {


  /// View (aka generalized slice) of a tensor

  /**
      @tparam _T apparent element type, TensorView will present tensor elements as values of this type
      @tparam _Range Range type
      @tparam _Storage Storage type
  */
  template<typename _T,
           class _Range = btas::DEFAULT::range,
           class _Storage = btas::DEFAULT::storage<_T>
           >
  class TensorView {

    public:

      /// value type
      typedef _T value_type;

      /// type of Range
      typedef _Range range_type;

      /// type of index
      typedef typename _Range::index_type index_type;

      /// type of underlying data storage
      typedef _Storage storage_type;

      // for convenience
      typedef typename std::remove_const<storage_type>::type nonconst_storage_type;

      /// type of data storage reference
      typedef StorageRef<storage_type> storageref_type;

      /// size type
      typedef typename storageref_type::size_type size_type;

      /// element iterator
      typedef TensorViewIterator<range_type, storage_type> iterator;

      /// element iterator
      typedef TensorViewIterator<range_type, const storage_type> const_iterator;

    private:
      struct Enabler {};

    public:

      /// default constructor
      TensorView () { }

      /// destructor
      ~TensorView () { }

      /// construct from \c range and \c storage
      template <typename S = _Storage>
      explicit
      TensorView (const range_type& range, S& storage)
      : range_(range), storageref_(storage)
      {
      }

      /// construct from \c range and \c storage
      template <typename S = _Storage>
      explicit
      TensorView (const range_type& range, const S& storage)
      : range_(range), storageref_(*const_cast<S*>(&storage))
      {
      }

      /// move-construct from \c range and \c storage
      explicit
      TensorView (range_type&& range, storage_type&& storage) :
      range_(range), storageref_(storage)
      {
      }

      /// conversion from Tensor
      template<class _Tensor, class = typename std::enable_if<is_boxtensor<_Tensor>::value>::type>
      explicit
      TensorView (const _Tensor& x)
      : range_ (x.range()),
      // TODO this can be optimized to bitewise copy if x::value_type and my value_type are equal, and storage is linear
        storageref_(x.storage())
      {
      }

      /// conversion from Tensor
      template<class _Tensor, class = typename std::enable_if<is_boxtensor<_Tensor>::value>::type>
      explicit
      TensorView (_Tensor& x)
      : range_ (x.range()),
      // TODO this can be optimized to bitewise copy if x::value_type and my value_type are equal, and storage is linear
        storageref_(x.storage())
      {
      }


      /// copy constructor
      TensorView (const TensorView& x)
      : range_ (x.range()), storageref_(x.storageref_)
      {
      }

      /// copy assignment
      TensorView&
      operator= (const TensorView& x)
      {
        range_ = x.range_;
        storageref_ = x.storageref_;
        return *this;
      }

      /// move constructor
      explicit
      TensorView (TensorView&& x)
      {
        std::swap(range_, x.range_);
        std::swap(storageref_, x.storageref_);
      }

      /// move assignment operator
      TensorView&
      operator= (TensorView&& x)
      {
        std::swap(range_, x.range_);
        std::swap(storageref_, x.storageref_);
        return *this;
      }

      /// number of indices (tensor rank)
      size_type
      rank () const
      {
        return range_.rank();
      }

      /// \return number of elements
      size_type
      size () const
      {
        return range_.area();
      }

      /// \return range object
      const range_type&
      range() const
      {
        return range_;
      }

      /// \param d dimension
      /// \return subrange for dimension \d
      const Range1d<typename index_type::value_type>
      range(size_t d) const
      {
        return range_.range(d);
      }

      /// \return range's extent object
      typename range_type::extent_type
      extent() const
      {
        return range_.extent();
      }

      /// \return extent of range along dimension \c d
      typename range_type::extent_type::value_type
      extent(size_t d) const
      {
        return range_.extent(d);
      }

      /// \return storage object
      const storageref_type&
      storage() const
      {
        return storageref_;
      }

      /// \return storage object
      storageref_type&
      storage()
      {
        return storageref_;
      }

      /// test whether TensorView is empty
      bool
      empty() const
      {
        return range_.area() == 0;
      }

      /// \return const iterator begin
      const_iterator
      begin() const
      {
        return cbegin();
      }

      /// \return const iterator end
      const_iterator
      end() const
      {
        return cend();
      }

      /// \return const iterator begin, even if this is not itself const
      const_iterator
      cbegin() const
      {
        return const_iterator(range().begin(), storage());
      }

      /// \return const iterator end, even if this is not itself const
      const_iterator
      cend() const
      {
        return const_iterator(range().end(), storage());
      }

      /// \return iterator begin
      template <typename S = _Storage>
      typename std::enable_if<not std::is_const<S>::value,iterator>::type
      begin()
      {
        return iterator(range().begin(), storage());
      }

      /// \return iterator end
      template <typename S = _Storage>
      typename std::enable_if<not std::is_const<S>::value,iterator>::type
      end()
      {
        return iterator(range().end(), storage());
      }

      /// \return element without range check
      template<typename index0, typename... _args>
      typename std::enable_if<std::is_integral<index0>::value, const value_type&>::type
      operator() (const index0& first, const _args&... rest) const
      {
        typedef typename common_signed_type<index0, typename index_type::value_type>::type ctype;
        auto indexv = {static_cast<ctype>(first), static_cast<ctype>(rest)...};
        index_type index = array_adaptor<index_type>::construct(indexv.size());
        std::copy(std::begin(indexv), std::end(indexv), std::begin(index));
        return storageref_[ range_.ordinal(index) ];
      }

      /// \return element without range check (rank() == general)
      template <typename Index>
      typename std::enable_if<is_index<Index>::value, const value_type&>::type
      operator() (const Index& index) const
      {
        return storageref_[range_.ordinal(index)];
      }

      /// access element without range check
      template<typename index0, typename... _args>
      typename std::enable_if<std::is_integral<index0>::value, value_type&>::type
      operator() (const index0& first, const _args&... rest)
      {
        typedef typename common_signed_type<index0, typename index_type::value_type>::type ctype;
        auto indexv = {static_cast<ctype>(first), static_cast<ctype>(rest)...};
        index_type index = array_adaptor<index_type>::construct(indexv.size());
        std::copy(std::begin(indexv), std::end(indexv), std::begin(index));
        return storageref_[ range_.ordinal(index) ];
      }

      /// access element without range check (rank() == general)
      template <typename Index>
      typename std::enable_if<is_index<Index>::value, value_type&>::type
      operator() (const Index& index)
      {
        return storageref_[range_.ordinal(index)];
      }

      /// \return element without range check
      template<typename index0, typename... _args>
      typename std::enable_if<std::is_integral<index0>::value, const value_type&>::type
      at (const index0& first, const _args&... rest) const
      {
        typedef typename common_signed_type<index0, typename index_type::value_type>::type ctype;
        auto indexv = {static_cast<ctype>(first), static_cast<ctype>(rest)...};
        index_type index = array_adaptor<index_type>::construct(indexv.size());
        std::copy(std::begin(indexv), std::end(indexv), std::begin(index));
        assert( range_.includes(index) );
        return storageref_[ range_.ordinal(index) ];
      }

      /// \return element without range check (rank() == general)
      template <typename Index>
      typename std::enable_if<is_index<Index>::value, const value_type&>::type
      at (const Index& index) const
      {
        assert( range_.includes(index) );
        return storageref_[ range_.ordinal(index) ];
      }

      /// access element without range check
      template<typename index0, typename... _args>
      typename std::enable_if<std::is_integral<index0>::value, value_type&>::type
      at (const index0& first, const _args&... rest)
      {
        typedef typename common_signed_type<index0, typename index_type::value_type>::type ctype;
        auto indexv = {static_cast<ctype>(first), static_cast<ctype>(rest)...};
        index_type index = array_adaptor<index_type>::construct(indexv.size());
        std::copy(std::begin(indexv), std::end(indexv), std::begin(index));
        assert( range_.includes(index) );
        return storageref_[ range_.ordinal(index) ];
      }

      /// access element without range check (rank() == general)
      template <typename Index>
      typename std::enable_if<is_index<Index>::value, value_type&>::type
      at (const Index& index)
      {
        assert( range_.includes(index) );
        return storageref_[ range_.ordinal(index) ];
      }

      /// swap this and x
      void
      swap (TensorView& x)
      {
        std::swap(range_, x.range_);
        std::swap(storageref_, x.storageref_);
      }

      //  ========== Finished Public Interface and Its Reference Implementations ==========

      //
      //  Here come Non-Standard members (to be discussed)
      //
#if 0
      /// addition assignment
      TensorView&
      operator+= (const TensorView& x)
      {
        assert( std::equal(range_.begin(), range_.end(), x.range_.begin()) );
        std::transform(storageref_.begin(), storageref_.end(), x.storageref_.begin(), storageref_.begin(), std::plus<value_type>());
        return *this;
      }

      /// addition of tensors
      TensorView
      operator+ (const TensorView& x) const
      {
        TensorView y(*this); y += x;
        return y; /* automatically called move semantics */
      }

      /// subtraction assignment
      TensorView&
      operator-= (const TensorView& x)
      {
        assert(
            std::equal(range_.begin(), range_.end(), x.range_.begin()));
        std::transform(storageref_.begin(), storageref_.end(), x.storageref_.begin(), storageref_.begin(), std::minus<value_type>());
        return *this;
      }

      /// subtraction of tensors
      TensorView
      operator- (const TensorView& x) const
      {
        TensorView y(*this); y -= x;
        return y; /* automatically called move semantics */
      }

      /// fill all elements by val
      void
      fill (const value_type& val)
      {
        std::fill(storageref_.begin(), storageref_.end(), val);
      }

      /// generate all elements by gen()
      template<class Generator>
      void
      generate (Generator gen)
      {
          std::generate(storageref_.begin(), storageref_.end(), gen);
      }
#endif

    private:

      range_type range_;///< range object
      storageref_type storageref_;///< dataref

  }; // end of TensorView

  /// TensorConstView is a read-only variant of TensorView
  template <typename _T,
            class _Range   = btas::DEFAULT::range,
            class _Storage = btas::DEFAULT::storage<_T>
           >
  using TensorConstView = TensorView<_T, _Range, const _Storage>;

  template <typename _T, typename _Range, typename _Storage>
  auto cbegin(const btas::TensorView<_T, _Range, _Storage>& x) -> decltype(x.cbegin()) {
    return x.cbegin();
  }
  template <typename _T, typename _Range, typename _Storage>
  auto cend(const btas::TensorView<_T, _Range, _Storage>& x) -> decltype(x.cbegin()) {
    return x.cend();
  }

  /// maps TensorView -> Range
  template <typename _T, typename _Range, typename _Storage>
  auto
  range (const btas::TensorView<_T, _Range, _Storage>& t) -> decltype(t.range()) {
    return t.range();
  }

  /// maps TensorView -> Range extent
  template <typename _T, typename _Range, typename _Storage>
  auto
  extent (const btas::TensorView<_T, _Range, _Storage>& t) -> decltype(t.range().extent()) {
    return t.range().extent();
  }

  /// TensorView stream output operator

  /// prints TensorView in row-major form. To be implemented elsewhere using slices.
  /// \param os The output stream that will be used to print \c t
  /// \param t The TensorView to be printed
  /// \return A reference to the output stream
  template <typename _T, typename _Range, typename _Storage>
  std::ostream& operator<<(std::ostream& os, const btas::TensorView<_T, _Range, _Storage>& t) {
    os << "TensorView:\n  Range: " << t.range() << std::endl;
    return os;
  }

} // namespace btas


#endif /* TENSORVIEW_H_ */
