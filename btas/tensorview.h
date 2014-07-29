/*
 * tensorview.h
 *
 *  Created on: Dec 28, 2013
 *      Author: evaleev
 */

#ifndef BTAS_TENSORVIEW_H_
#define BTAS_TENSORVIEW_H_

#include <functional>

#include <btas/storage_traits.h>
#include <btas/util/sequence_adaptor.h>
#include <btas/tensorview_iterator.h>
#include <btas/defaults.h>
#include <btas/util/functional.h>

namespace btas {

  enum TensorViewPolicy_ConstnessPolicy {
    TensorViewPolicy_RuntimeConst = 1,
    TensorViewPolicy_CompiletimeConst = 0
  };

  /// TensorViewPolicy configures behavior of certain features of TensorView
  /// \tparam RuntimeConst: if true, constness of data access is checked at runtime. This involves
  ///      extra space overhead (enough to store a boolean readwrite flag). Non-const data access members
  ///      will also check whether readwrite is set using assert (hence runtime overhead can be eliminated after
  ///      testing. This feature is needed if you want to use a single TensorView<T,Range,Storage> type
  ///      for mutable (non-const) and immutable (const) views. The default value is false, which requires use
  ///      of TensorView<T,Range,const Storage>, aka TensorConstView<T.Range,Storage>, for immutable views.
  template <TensorViewPolicy_ConstnessPolicy ConstnessPolicy = TensorViewPolicy_CompiletimeConst>
  struct TensorViewPolicy {
      /// true if constness tracked at runtime
      static constexpr bool runtimeconst = (ConstnessPolicy == TensorViewPolicy_RuntimeConst);
  };

  /// View (aka generalized slice) of a tensor

  /**
      @tparam _T apparent element type, TensorView will present tensor elements as values of this type
      @tparam _Range Range type
      @tparam _Storage Storage type
  */
  template<typename _T,
           class _Range = btas::DEFAULT::range,
           class _Storage = btas::DEFAULT::storage<_T>,
           class _Policy = btas::TensorViewPolicy<>
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

      /// type of data storage reference
      typedef std::reference_wrapper<storage_type> storageref_type;

      /// size type
      typedef typename storage_traits<storage_type>::size_type size_type;

      /// element iterator
      typedef TensorViewIterator<range_type, storage_type> iterator;

      /// element iterator
      typedef TensorViewIterator<range_type, const storage_type> const_iterator;

    private:
      struct Enabler {};

    public:

      /// destructor
      ~TensorView () { }

      /// move-construct from \c range and \c storageref ; write access must be passed explicitly if \c _Policy requires
      template<class Policy = _Policy, class = typename std::enable_if<not Policy::runtimeconst>::type>
      explicit
      TensorView (range_type&& range,
                  storageref_type&& storageref,
                  bool can_write = not _Policy::runtimeconst ? not std::is_const<storage_type>::value : false) :
      range_(range), storageref_(storageref), can_write_(can_write)
      {
      }

      /// conversion from const Tensor into TensorConstView
      template<class _Tensor,
               class Storage = _Storage,
               class = typename std::enable_if<is_boxtensor<_Tensor>::value &&
                                               std::is_const<Storage>::value>::type
              >
      TensorView (const _Tensor& x)
      : range_ (x.range()),
        storageref_(std::cref(x.storage())),
        can_write_(false)
      {
      }

      /// conversion from const Tensor to non-const View only possible if \c Policy::runtimeconst is \c true
      template<class _Tensor,
               class Storage = _Storage,
               class Policy = _Policy,
               class = typename std::enable_if<is_boxtensor<_Tensor>::value &&
                                               not std::is_const<Storage>::value &&
                                               Policy::runtimeconst>::type
              >
      TensorView (const _Tensor& x)
      : range_ (x.range()),
        storageref_(std::ref(const_cast<storage_type&>(x.storage()))),
        can_write_(false)
      {
      }

      /// conversion from non-const Tensor
      template<class _Tensor,
               class Storage = _Storage,
               class = typename std::enable_if<is_boxtensor<_Tensor>::value &&
                                               std::is_same<typename _Tensor::storage_type,Storage>::value>::type>
      TensorView (_Tensor& x)
      : range_ (x.range()),
        storageref_(std::ref(x.storage())),
        can_write_(true)
      {
      }

      /// conversion from non-const TensorView
      template<class __T,
               class __Range,
               class __Storage,
               class __Policy,
               class = typename std::enable_if<not std::is_const<__Storage>::value>::type>
      TensorView (TensorView<__T,__Range,__Storage,__Policy>& x)
      : range_ (x.range()),
        storageref_(std::ref(x.storage())),
        can_write_(_Policy::runtimeconst ? bool(x.can_write_) : not std::is_const<storage_type>::value)
      {
      }

      /// standard copy constructor
      TensorView (const TensorView& x) :
        range_ (x.range_),
        storageref_(x.storageref_),
        can_write_(false)
      {
      }

      /// copy assignment
      TensorView&
      operator= (const TensorView& x)
      {
        range_ = x.range_;
        storageref_ = x.storageref_;
        can_write_ = x.can_write_;
        return *this;
      }

      /// move constructor
      TensorView (TensorView&& x) : range_(), storageref_(x.storageref_), can_write_(x.can_write_)
      {
        std::swap(range_, x.range_);
      }

      /// move assignment operator
      TensorView&
      operator= (TensorView&& x)
      {
        std::swap(range_, x.range_);
        std::swap(storageref_, x.storageref_);
        std::swap(can_write_, x.can_write_);
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
      const storage_type&
      storage() const
      {
        return storageref_.get();
      }

      /// \return storage object
      storage_type&
      storage()
      {
        assert_writable();
        return storageref_.get();
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

      /// \return begin iterator
      iterator
      begin()
      {
        assert_writable();
        return iterator(range().begin(), storage());
      }

      /// \return const end iterator
      const_iterator
      end() const
      {
        return cend();
      }

      /// \return const end iterator
      iterator
      end()
      {
        assert_writable();
        return iterator(range().end(), storageref_);
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

      /// Immutable access to an element without range check.

      /// Available when \c value_type == \c storage_type::value_type
      /// \return const reference to the element indexed by {\c first, \c rest}
      template<typename index0, typename... _args>
      typename std::enable_if<std::is_integral<index0>::value &&
                               std::is_same<value_type,typename storage_type::value_type>::value,
                              const value_type&
                             >::type
      operator() (const index0& first, const _args&... rest) const
      {
        typedef typename common_signed_type<index0, typename index_type::value_type>::type ctype;
        auto indexv = {static_cast<ctype>(first), static_cast<ctype>(rest)...};
        index_type index = array_adaptor<index_type>::construct(indexv.size());
        std::copy(std::begin(indexv), std::end(indexv), std::begin(index));
        return storageref_.get()[ range_.ordinal(index) ];
      }

      /// Immutable access to an element without range check.

      /// Available when \c value_type == \c storage_type::value_type
      /// \return const reference to the element indexed by \c index
      template <typename Index>
      typename std::enable_if<is_index<Index>::value &&
                               std::is_same<value_type,typename storage_type::value_type>::value,
                              const value_type&
                             >::type
      operator() (const Index& index) const
      {
        return storageref_.get()[range_.ordinal(index)];
      }

      /// Mutable access to an element without range check.

      /// Available when \c value_type == \c storage_type::value_type
      /// \return reference to the element indexed by {\c first, \c rest}
      template<typename index0, typename... _args>
      typename std::enable_if<std::is_integral<index0>::value &&
                               std::is_same<value_type,typename storage_type::value_type>::value,
                              value_type&
                             >::type
      operator() (const index0& first, const _args&... rest)
      {
        assert_writable();
        typedef typename common_signed_type<index0, typename index_type::value_type>::type ctype;
        auto indexv = {static_cast<ctype>(first), static_cast<ctype>(rest)...};
        index_type index = array_adaptor<index_type>::construct(indexv.size());
        std::copy(std::begin(indexv), std::end(indexv), std::begin(index));
        return storageref_.get()[ range_.ordinal(index) ];
      }

      /// Mutable access to an element without range check (rank() == general)

      /// Available when \c value_type == \c storag_type::value_type
      /// \return reference to the element indexed by \c index
      template <typename Index>
      typename std::enable_if<is_index<Index>::value &&
                               std::is_same<value_type,typename storage_type::value_type>::value,
                              value_type&
                             >::type
      operator() (const Index& index)
      {
        assert_writable();
        return storageref_.get()[range_.ordinal(index)];
      }

      /// Immutable access to an element without range check.

      /// Available when \c value_type != \c storage_type::value_type
      /// \return value of the element indexed by {\c first, \c rest}
      template<typename index0, typename... _args>
      typename std::enable_if<std::is_integral<index0>::value &&
                              not std::is_same<value_type,typename storage_type::value_type>::value,
                              value_type>::type
      operator() (const index0& first, const _args&... rest) const
      {
        typedef typename common_signed_type<index0, typename index_type::value_type>::type ctype;
        auto indexv = {static_cast<ctype>(first), static_cast<ctype>(rest)...};
        index_type index = array_adaptor<index_type>::construct(indexv.size());
        std::copy(std::begin(indexv), std::end(indexv), std::begin(index));
        return storageref_.get()[ range_.ordinal(index) ];
      }

      /// Immutable access to an element without range check (rank() == general)

      /// Available when \c value_type != \c storage_type::value_type
      /// \return value of the element indexed by \c index
      template <typename Index>
      typename std::enable_if<is_index<Index>::value &&
                               not std::is_same<value_type,typename storage_type::value_type>::value,
                              value_type
                             >::type
      operator() (const Index& index) const
      {
        return storageref_.get()[range_.ordinal(index)];
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
        return storageref_.get()[ range_.ordinal(index) ];
      }

      /// \return element without range check (rank() == general)
      template <typename Index>
      typename std::enable_if<is_index<Index>::value, const value_type&>::type
      at (const Index& index) const
      {
        assert( range_.includes(index) );
        return storageref_.get()[ range_.ordinal(index) ];
      }

      /// access element without range check
      template<typename index0, typename... _args>
      typename std::enable_if<std::is_integral<index0>::value, value_type&>::type
      at (const index0& first, const _args&... rest)
      {
        assert_writable();
        typedef typename common_signed_type<index0, typename index_type::value_type>::type ctype;
        auto indexv = {static_cast<ctype>(first), static_cast<ctype>(rest)...};
        index_type index = array_adaptor<index_type>::construct(indexv.size());
        std::copy(std::begin(indexv), std::end(indexv), std::begin(index));
        assert( range_.includes(index) );
        return storageref_.get()[ range_.ordinal(index) ];
      }

      /// access element without range check (rank() == general)
      template <typename Index>
      typename std::enable_if<is_index<Index>::value, value_type&>::type
      at (const Index& index)
      {
        assert_writable();
        assert( range_.includes(index) );
        return storageref_.get()[ range_.ordinal(index) ];
      }

      /// swap this and x
      void
      swap (TensorView& x)
      {
        std::swap(range_, x.range_);
        std::swap(storageref_, x.storageref_);
        std::swap(can_write_, x.can_write_);
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
      typedef typename std::conditional<_Policy::runtimeconst,
                                        bool,
                                        btas::detail::bool_type<not std::is_const<storage_type>::value>
                                       >::type writable_type;
      writable_type can_write_;

      /// use this in non-const members to assert writability if Policy calls for runtime const check
      void assert_writable() const {
        if (_Policy::runtimeconst)
          assert(can_write_ == true);
      }

      /// construct from \c range and \c storage; pass \c can_write explicitly if needed
      explicit TensorView (range_type&& range, storage_type& storage,
                           bool can_write = not _Policy::runtimeconst ? not std::is_const<storage_type>::value : false) :
          range_(std::move(range)), storageref_(std::ref(storage)), can_write_(can_write) {
      }

      template <typename T,
                typename Range,
                typename Storage,
                typename Policy>
      friend TensorView<T,
                 Range,
                 typename std::conditional<std::is_same<T,typename Storage::value_type>::value,
                                           Storage,
                                           typename std::add_const<Storage>::type
                                          >::type,
                 Policy>
      __make_view(Range&& range, Storage& storage,
                  Policy,
                  bool can_write);
      template <typename T,
                typename Range,
                typename Storage,
                typename Policy>
      friend TensorView<T, Range, const Storage, Policy> __make_cview(Range&& range, const Storage& storage, Policy);

      template <class __T,
                class __Range,
                class __Storage,
                class __Policy>
      friend class TensorView;
  }; // end of TensorView

  /// TensorConstView is a read-only variant of TensorView
  template <typename _T,
            class _Range   = btas::DEFAULT::range,
            class _Storage = btas::DEFAULT::storage<_T>,
            class _Policy  = btas::TensorViewPolicy<>
           >
  using TensorConstView = TensorView<_T, _Range, const _Storage, _Policy>;

  /// TensorRWView is a variant of TensorView with runtime write access check
  template <typename _T,
            class _Range   = btas::DEFAULT::range,
            class _Storage = btas::DEFAULT::storage<_T>,
            class _Policy  = btas::TensorViewPolicy<TensorViewPolicy_RuntimeConst>
           >
  using TensorRWView = TensorView<_T, _Range, typename std::remove_const<_Storage>::type, _Policy>;


  /// Helper function (friendly to TensorView) that constructs a view with an explicitly-specified element type of the view. Useful if need to
  /// view a tensor of floats as a tensor of complex floats.
  /// \tparam T the element type of the resulting view
  /// \tparam Range the range type
  /// \tparam Storage the storage type
  /// \tparam Policy the TensorViewPolicy type
  /// \param range the range object defining the view
  /// \param storage the storage object that will be viewed into
  /// \return TensorView into \c storage using \c range and policy \c Policy
  /// \attention use __make_cview if you must force a const view; this will provide const view, however, if \c storage is a const reference.
  template <typename T,
            typename Range,
            typename Storage,
            typename Policy>
  TensorView<T,
             Range,
             typename std::conditional<std::is_same<T,typename Storage::value_type>::value,
                                       Storage,
                                       typename std::add_const<Storage>::type
                                      >::type,
             Policy>
  __make_view(Range&& range, Storage& storage,
              Policy = Policy(),
              bool can_write = not Policy::runtimeconst
                               ? (not std::is_const<Storage>::value && std::is_same<T,typename Storage::value_type>::value)
                               : false)
  {
    typedef  TensorView<T,
        Range,
        typename std::conditional<std::is_same<T,typename Storage::value_type>::value,
                                  Storage,
                                  typename std::add_const<Storage>::type
                                 >::type,
        Policy> result_type;
    return result_type(std::move(range), storage, can_write);
  }

  /// Helper function (friendly to TensorView) that constructs a view, with an explicitly-specified element type of the view. Useful if need to
  /// view a tensor of floats as a tensor of complex floats. \sa TensorConstView
  /// \tparam T the element type of the resulting view
  /// \tparam Range the range type
  /// \tparam Storage the storage type
  /// \param range the range object defining the view
  /// \param storage the storage object that will be viewed into
  /// \return TensorView into \c storage using \c range and policy \c Policy
  template <typename T,
            typename Range,
            typename Storage,
            typename Policy>
  TensorView<T, Range, const Storage, Policy>
  __make_cview(Range&& range, const Storage& storage, Policy = Policy())
  {
    return TensorView<T, Range, const Storage, Policy>(std::move(range), storage, false);
  }

  /// Helper function that constructs TensorView.
  /// \tparam Range the range type
  /// \tparam Storage the storage type
  /// \tparam Policy the TensorViewPolicy type; if the Policy requires additional runtime parameters use __make_view instead
  /// \param range the range object defining the view
  /// \param storage the storage object that will be viewed into
  /// \return TensorView into \c storage using \c range, with policy \c Policy
  /// \attention use make_cview if you must force a const view; this will provide const view, however, if \c storage is a const reference.
  template <typename Range,
            typename Storage,
            typename Policy = TensorViewPolicy<TensorViewPolicy_CompiletimeConst>,
            class = typename std::enable_if<not std::is_reference<Range>::value>::type>
  TensorView<typename Storage::value_type, Range, Storage, Policy>
  make_view(const Range& range, Storage& storage, Policy = Policy())
  {
    return make_view<typename Storage::value_type, Range, Storage, Policy>(range, storage);
  }

  /// Helper function that constructs TensorView.
  /// \tparam Range the range type
  /// \tparam Storage the storage type
  /// \tparam Policy the TensorViewPolicy type
  /// \param range the range object defining the view
  /// \param storage the storage object that will be viewed into
  /// \return TensorView into \c storage using \c range, with policy \c Policy
  /// \attention use make_cview if you must force a const view; this will provide const view, however, if \c storage is a const reference.
  template <typename Range,
            typename Storage,
            typename Policy = TensorViewPolicy<TensorViewPolicy_CompiletimeConst>,
            class = typename std::enable_if<not std::is_reference<Range>::value>::type>
  TensorView<typename Storage::value_type, Range, Storage, Policy>
  make_view(Range&& range, Storage& storage, Policy = Policy())
  {
    return make_view<typename Storage::value_type, Range, Storage, Policy>(range, storage);
  }


  /// Helper function that constructs TensorView, with an explicitly-specified element type of the view. Useful if need to
  /// view a tensor of floats as a tensor of complex floats.
  /// \tparam T the element type of the resulting view
  /// \tparam Range the range type
  /// \tparam Storage the storage type
  /// \tparam Policy the TensorViewPolicy type
  /// \param range the range object defining the view
  /// \param storage the storage object that will be viewed into
  /// \return TensorView into \c storage using \c range, with policy \c Policy
  /// \attention use make_cview if you must force a const view; this will provide const view, however, if \c storage is a const reference.
  template <typename T,
            typename Range,
            typename Storage,
            typename Policy = TensorViewPolicy<TensorViewPolicy_CompiletimeConst>,
            class = typename std::enable_if<not std::is_reference<Range>::value>::type>
  auto
  make_view(const Range& range, Storage& storage, Policy = Policy()) -> decltype(__make_view<T, Range, Storage, Policy>(Range(range), storage))
  {
    return __make_view<T, Range, Storage, Policy>(Range(range), storage);
  }

  /// Helper function that constructs TensorView, with an explicitly-specified element type of the view. Useful if need to
  /// view a tensor of floats as a tensor of complex floats.
  /// \tparam T the element type of the resulting view
  /// \tparam Range the range type
  /// \tparam Storage the storage type
  /// \tparam Policy the TensorViewPolicy type
  /// \param range the range object defining the view
  /// \param storage the storage object that will be viewed into
  /// \return TensorView into \c storage using \c range, with policy \c Policy
  /// \attention use make_cview if you must force a const view; this will provide const view, however, if \c storage is a const reference.
  template <typename T,
            typename Range,
            typename Storage,
            typename Policy = TensorViewPolicy<TensorViewPolicy_CompiletimeConst>,
            class = typename std::enable_if<not std::is_reference<Range>::value>::type>
  auto
  make_view(Range&& range, Storage& storage, Policy = Policy()) -> decltype(__make_view<T, Range, Storage, Policy>(range, storage))
  {
    return __make_view<T, Range, Storage, Policy>(range, storage);
  }

  /// Helper function that constructs a full TensorView of a Tensor.
  /// \tparam Tensor the tensor type
  /// \param tensor the Tensor object
  /// \return TensorView, a full view of the \c tensor
  /// \attention use make_cview if you must force a const view; this will provide const view, however, if \c tensor is a const reference.
  /// \note Provided for completeness.
  template <typename Tensor,
            typename Policy = TensorViewPolicy<TensorViewPolicy_CompiletimeConst>,
            class = typename std::enable_if<is_boxtensor<Tensor>::value>::type>
  TensorView<typename Tensor::value_type,
             typename Tensor::range_type,
             typename Tensor::storage_type,
             Policy>
  make_view(Tensor& tensor, Policy = Policy())
  {
    return TensorView<typename Tensor::value_type,
                      typename Tensor::range_type,
                      typename Tensor::storage_type,
                      Policy>(tensor);
  }

  /// Helper function that constructs a full TensorView of a Tensor,
  /// with an explicitly-specified element type of the view. Useful if need to
  /// view a tensor of floats as a tensor of complex floats.
  /// \tparam T the element type of the resulting view
  /// \tparam Tensor the tensor type
  /// \param tensor the Tensor object
  /// \return TensorView, a full view of the \c tensor
  /// \attention use make_cview if you must force a const view; this will provide const view, however, if \c tensor is a const reference.
  /// \note Provided for completeness.
  template <typename T, typename Tensor,
            typename Policy = TensorViewPolicy<TensorViewPolicy_CompiletimeConst>,
            class = typename std::enable_if<is_boxtensor<Tensor>::value>::type>
  TensorView<T,
             typename Tensor::range_type,
             typename std::conditional<std::is_same<T,typename Tensor::storage_type::value_type>::value,
                                       typename Tensor::storage_type,
                                       typename std::add_const<typename Tensor::storage_type>::type
                                      >::type,
             Policy>
  make_view(Tensor& tensor, Policy = Policy())
  {
      typedef   TensorView<T,
          typename Tensor::range_type,
          typename std::conditional<std::is_same<T,typename Tensor::storage_type::value_type>::value,
                                    typename Tensor::storage_type,
                                    typename std::add_const<typename Tensor::storage_type>::type
                                   >::type,
          Policy> result_type;
    return result_type(tensor);
  }

  /// Helper function that constructs a constant TensorView. \sa TensorConstView
  /// \tparam Range the range type
  /// \tparam Storage the storage type
  /// \param range the range object defining the view
  /// \param storage the storage object that will be viewed into
  /// \return TensorView into \c storage using \c range
  template <typename Range,
            typename Storage,
            typename Policy = btas::TensorViewPolicy<TensorViewPolicy_CompiletimeConst>,
            class = typename std::enable_if<not std::is_reference<Range>::value>::type>
  TensorView<typename std::add_const<typename Storage::value_type>::type, Range, const Storage, Policy>
  make_cview(const Range& range, const Storage& storage, Policy = Policy())
  {
    return make_cview<typename std::add_const<typename Storage::value_type>::type, Range, Storage, Policy>(range, storage);
  }

  /// Helper function that constructs a constant TensorView, with an explicitly-specified element type of the view. Useful if need to
  /// view a tensor of floats as a tensor of complex floats. \sa TensorConstView
  /// \tparam T the element type of the resulting view
  /// \tparam Range the range type
  /// \tparam Storage the storage type
  /// \param range the range object defining the view
  /// \param storage the storage object that will be viewed into
  /// \return TensorView into \c storage using \c range
  template <typename T,
            typename Range,
            typename Storage,
            typename Policy = TensorViewPolicy<TensorViewPolicy_CompiletimeConst>,
            class = typename std::enable_if<not std::is_reference<Range>::value>::type>
  TensorView<typename std::add_const<T>::type, Range, const Storage, Policy>
  make_cview(const Range& range, const Storage& storage, Policy = Policy())
  {
    return __make_cview<typename std::add_const<T>::type, Range, const Storage, Policy>(Range(range), storage);
  }

  /// Helper function that constructs a full constant TensorView of a Tensor.
  /// \tparam Tensor the tensor type
  /// \param tensor the Tensor object
  /// \return TensorView, a full view of the \c tensor
  /// \note Provided for completeness.
  template <typename Tensor,
            typename Policy = TensorViewPolicy<TensorViewPolicy_CompiletimeConst>,
            class = typename std::enable_if<is_boxtensor<Tensor>::value>::type>
  TensorView<typename std::add_const<typename Tensor::value_type>::type,
             typename Tensor::range_type,
             const typename Tensor::storage_type,
             Policy>
  make_cview(const Tensor& tensor)
  {
    return TensorView<typename std::add_const<typename Tensor::value_type>::type,
                      typename Tensor::range_type,
                      const typename Tensor::storage_type,
                      Policy>(tensor);
  }

  /// Helper function that constructs a full constant TensorView of a Tensor,
  /// with an explicitly-specified element type of the view. Useful if need to
  /// view a tensor of floats as a tensor of complex floats.
  /// \tparam T the element type of the resulting view
  /// \tparam Tensor the tensor type
  /// \param tensor the Tensor object
  /// \return TensorView, a full view of the \c tensor
  /// \note Provided for completeness.
  template <typename T, typename Tensor,
            typename Policy = TensorViewPolicy<TensorViewPolicy_CompiletimeConst>,
            class = typename std::enable_if<is_boxtensor<Tensor>::value>::type>
  TensorView<typename std::add_const<T>::type,
             typename Tensor::range_type,
             const typename Tensor::storage_type,
             Policy>
  make_cview(const Tensor& tensor)
  {
    return TensorView<typename std::add_const<T>::type,
                      typename Tensor::range_type,
                      const typename Tensor::storage_type,
                      Policy>(tensor);
  }

  /// Helper function that constructs writable TensorView.
  /// \tparam Range the range type
  /// \tparam Storage the storage type
  /// \param range the range object defining the view
  /// \param storage the storage object that will be viewed into
  /// \return TensorView into \c storage using \c range
  /// \attention use make_cview if you must force a const view; this will provide const view, however, if \c storage is a const reference.
  template <typename Range,
            typename Storage,
            class = typename std::enable_if<not std::is_reference<Range>::value>::type>
  TensorRWView<typename Storage::value_type, Range, Storage>
  make_rwview(const Range& range,
              Storage& storage,
              bool can_write = not std::is_const<Storage>::value)
  {
    // enforce mutability
    can_write = can_write && (not std::is_const<Storage>::value);
    return make_rwview(Range(range), storage, can_write);
  }

  /// Helper function that constructs writable TensorView.
  /// \tparam Range the range type
  /// \tparam Storage the storage type
  /// \param range the range object defining the view
  /// \param storage the storage object that will be viewed into
  /// \return TensorView into \c storage using \c range
  /// \attention use make_cview if you must force a const view; this will provide const view, however, if \c storage is a const reference.
  template <typename Range,
            typename Storage,
            class = typename std::enable_if<not std::is_reference<Range>::value>::type>
  TensorRWView<typename Storage::value_type, Range, Storage>
  make_rwview(Range&& range,
              Storage& storage,
              bool can_write = not std::is_const<Storage>::value)
  {
    // enforce mutability
    can_write = can_write && (not std::is_const<Storage>::value);
    return make_rwview<typename Storage::value_type, Range, Storage>(std::move(range), storage, can_write);
  }

  /// Helper function that constructs writable TensorView, with an explicitly-specified element type of the view. Useful if need to
  /// view a tensor of floats as a tensor of complex floats.
  /// \tparam T the element type of the resulting view
  /// \tparam Range the range type
  /// \tparam Storage the storage type
  /// \param range the range object defining the view
  /// \param storage the storage object that will be viewed into
  /// \return TensorView into \c storage using \c range
  /// \attention use make_cview if you must force a const view; this will provide const view, however, if \c storage is a const reference.
  template <typename T,
            typename Range,
            typename Storage,
            class = typename std::enable_if<not std::is_reference<Range>::value>::type>
  TensorRWView<T, Range, Storage>
  make_rwview(const Range& range, Storage& storage,
              bool can_write = not std::is_const<Storage>::value &&
                               std::is_same<T,typename Storage::value_type>::value)
  {
    // enforce mutability
    can_write = can_write && (not std::is_const<Storage>::value &&
                              std::is_same<T,typename Storage::value_type>::value);
    return make_rwview(Range(range),
                       storage,
                       can_write);
  }

  /// Helper function that constructs writable TensorView, with an explicitly-specified element type of the view. Useful if need to
  /// view a tensor of floats as a tensor of complex floats.
  /// \tparam T the element type of the resulting view
  /// \tparam Range the range type
  /// \tparam Storage the storage type
  /// \param range the range object defining the view
  /// \param storage the storage object that will be viewed into
  /// \return TensorView into \c storage using \c range
  /// \attention use make_cview if you must force a const view; this will provide const view, however, if \c storage is a const reference.
  template <typename T,
            typename Range,
            typename Storage,
            class = typename std::enable_if<not std::is_reference<Range>::value>::type>
  TensorRWView<T, Range, Storage>
  make_rwview(Range&& range, Storage& storage,
              bool can_write = not std::is_const<Storage>::value &&
                               std::is_same<T,typename Storage::value_type>::value)
  {
    // enforce mutability
    can_write = can_write && (not std::is_const<Storage>::value &&
                              std::is_same<T,typename Storage::value_type>::value);
    return __make_view<T,
                       Range,
                       typename std::remove_const<Storage>::type,
                       TensorViewPolicy<TensorViewPolicy_RuntimeConst> >(std::move(range),
                          const_cast<typename std::remove_const<Storage>::type&>(storage),
                          TensorViewPolicy<TensorViewPolicy_RuntimeConst>(),
                          can_write);
  }

  /// Helper function that constructs a full writable TensorView of a Tensor.
  /// \tparam Tensor the tensor type
  /// \param tensor the Tensor object
  /// \return TensorView, a full view of the \c tensor
  /// \note Provided for completeness.
  template <typename Tensor, class = typename std::enable_if<is_boxtensor<Tensor>::value>::type>
  TensorRWView<typename Tensor::value_type,
               typename Tensor::range_type,
               typename Tensor::storage_type>
  make_rwview(Tensor& tensor,
              bool can_write = not std::is_const<Tensor>::value &&
                               not std::is_const<typename Tensor::storage_type>::value)
  {
      // enforce mutability
      can_write = can_write && (not std::is_const<Tensor>::value && not std::is_const<typename Tensor::storage_type>::value);
      return make_rwview(tensor.range(), tensor.storage(), can_write);
  }

  /// Helper function that constructs a full writable TensorView of a Tensor,
  /// with an explicitly-specified element type of the view. Useful if need to
  /// view a tensor of floats as a tensor of complex floats.
  /// \tparam T the element type of the resulting view
  /// \tparam Tensor the tensor type
  /// \param tensor the Tensor object
  /// \return TensorView, a full view of the \c tensor
  /// \note Provided for completeness.
  template <typename T, typename Tensor, class = typename std::enable_if<is_boxtensor<Tensor>::value>::type>
  TensorRWView<T,
               typename Tensor::range_type,
               typename Tensor::storage_type>
  make_rwview(Tensor& tensor,
              bool can_write = not std::is_const<Tensor>::value &&
                               not std::is_const<typename Tensor::storage_type>::value &&
                               std::is_same<T,typename Tensor::storage_type::value_type>::value)
  {
      // enforce mutability
      can_write = can_write &&
                  (not std::is_const<Tensor>::value &&
                   not std::is_const<typename Tensor::storage_type>::value &&
                   std::is_same<T,typename Tensor::storage_type::value_type>::value);
      return make_rwview(tensor.range(), tensor.storage(), can_write);
  }

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

  /// TensorMap views a sequence of values as a Tensor
  template <typename _T,
            class _Range = btas::DEFAULT::range>
  using TensorMap = TensorView<_T, _Range, btas::infinite_sequence_adaptor<_T*>>;
  /// TensorConstMap const-views a sequence of values as a Tensor
  template <typename _T,
            class _Range = btas::DEFAULT::range>
  using TensorConstMap = TensorView<_T, _Range, const btas::infinite_sequence_adaptor<const _T*>>;

  /// Helper function that constructs TensorMap.
  /// \tparam T the element type returned by the view
  /// \tparam Range the range type
  /// \param range the range object defining the view
  /// \return TensorView into \c storage using \c range
  /// \attention use make_cmap if you must force a const view; this will provide const view, however, if \c storage is a const reference.
  template <typename T,
            typename Range>
  TensorMap<T, Range>
  make_map(T* data, Range&& range)
  {
    return TensorMap<T, Range>(std::move(range),
                               std::ref(btas::infinite_sequence_adaptor<T*>(data)));
  }

  /// Helper function that constructs TensorConstMap.
  /// \tparam T the element type returned by the view
  /// \tparam Range the range type
  /// \param range the range object defining the view
  /// \return TensorView into \c storage using \c range
  /// \attention use make_cmap if you must force a const view; this will provide const view, however, if \c storage is a const reference.
  template <typename T,
            typename Range>
  TensorConstMap<T, Range>
  make_map(const T* data, Range&& range)
  {
    return TensorConstMap<T, Range>(std::move(range),
                                    std::cref(btas::infinite_sequence_adaptor<const T*>(data)));
  }

  /// Helper function that constructs TensorConstMap.
  /// \tparam Range the range type
  /// \param range the range object defining the view
  /// \return TensorView into \c storage using \c range
  /// \attention use make_cmap if you must force a const view; this will provide const view, however, if \c storage is a const reference.
  template <typename T,
            typename Range>
  TensorConstMap<typename std::remove_const<T>::type, Range>
  make_cmap(T* data, Range&& range)
  {
    typedef typename std::remove_const<T>::type value_type;
    typedef TensorConstMap<value_type, Range> result_type;
    return result_type(std::move(range),
                       std::cref(btas::infinite_sequence_adaptor<const T*>(const_cast<const T*>(data))));
  }

} // namespace btas


#endif /* TENSORVIEW_H_ */
