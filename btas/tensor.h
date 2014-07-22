#ifndef __BTAS_TENSOR_H
#define __BTAS_TENSOR_H 1

#include <cassert>
#include <algorithm>
#include <functional>
#include <type_traits>
#include <vector>

#include <btas/types.h>
#include <btas/defaults.h>
#include <btas/tensor_traits.h>
#include <btas/tensorview.h>
#include <btas/array_adaptor.h>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>

namespace btas {

  /** BTAS implementation of "dense" tensor class that models \ref labelTWGTensor "TWG.BoxTensor" concept
      @tparam _T element type, Tensor contains values of this type
      @tparam _Range Range type, models \ref labelTWGRange "TWG.Range" concept
      @tparam _Storage Storage type, models \ref labelTWGStorage "TWG.Storage" concept
  */
  template<typename _T,
           class _Range = btas::DEFAULT::range,
           class _Storage = btas::DEFAULT::storage<_T>,
           class = typename std::enable_if<std::is_same<_T, typename _Storage::value_type>::value>::type
          >
  class Tensor {

    public:

      /// type of underlying data storage
      typedef _Storage storage_type;

      /// type of Range
      typedef _Range range_type;

      /// type of index
      typedef typename _Range::index_type index_type;

      ///\name Container requirements (c++std:[container.requirements.general]).
      ///@{

      /// value type
      typedef _T value_type;

      /// lvalue type of _T
      typedef value_type& reference;

      /// const lvalue type of _T
      typedef const value_type& const_reference;

      /// element iterator
      typedef typename storage_type::iterator iterator;

      /// constant element iterator
      typedef typename storage_type::const_iterator const_iterator;

      /// size type
      typedef typename storage_type::size_type size_type;

      ///@}

    private:
      struct Enabler {};

    public:

      /// default constructor
      Tensor () { }

      /// destructor
      ~Tensor () { }

      /// constructor with index extent
      template<typename... _args>
      explicit
      Tensor (const size_type& first, const _args&... rest) :
      range_(range_type(first, rest...))
      {
        // TODO make this disablable in all constructors
        //assert(range_.ordinal(range_.lobound()) == 0);
        array_adaptor<storage_type>::resize(storage_, range_.area());
      }

      /// construct from \c range, allocate data, but not initialized
      template <typename Range>
      explicit
      Tensor (const Range& range, typename std::enable_if<btas::is_boxrange<Range>::value>::type* = 0) :
      range_(range.lobound(), range.upbound())
      {
        array_adaptor<storage_type>::resize(storage_, range_.area());
      }

      /// construct from \c range object, set all elements to \c v
      template <typename Range>
      explicit
      Tensor (const Range& range,
              value_type v,
              typename std::enable_if<btas::is_boxrange<Range>::value>::type* = 0) :
              range_(range.lobound(), range.upbound())
      {
        array_adaptor<storage_type>::resize(storage_, range_.area());
        std::fill(begin(), end(), v);
      }

      /// construct from \c range and \c storage
      template <typename Range, typename Storage>
      explicit
      Tensor (const Range& range,
              const Storage& storage,
              typename std::enable_if<btas::is_boxrange<Range>::value &
                                      not std::is_same<Range,range_type>::value &
                                      not std::is_same<Storage,storage_type>::value
                                     >::type* = 0) :
      range_(range.lobound(), range.upbound()), storage_(storage)
      {
        if (storage_.size() != range_.area())
          array_adaptor<storage_type>::resize(storage_, range_.area());
      }

      /// copy-copy-construct from \c range and \c storage
      explicit
      Tensor (const range_type& range, const storage_type& storage) :
      range_(range.ordinal(*range.begin()) == 0 ? range : range_type(range.lobound(), range.upbound())),
      storage_(storage)
      {
        if (storage_.size() != range_.area())
          array_adaptor<storage_type>::resize(storage_, range_.area());
      }

      /// copy-move-construct from \c range and \c storage
      explicit
      Tensor (const range_type& range, storage_type&& storage) :
      range_(range.ordinal(*range.begin()) == 0 ? range : range_type(range.lobound(), range.upbound())),
      storage_(storage)
      {
        if (storage_.size() != range_.area())
          array_adaptor<storage_type>::resize(storage_, range_.area());
      }

      /// move-construct from \c range and \c storage
      explicit
      Tensor (range_type&& range, storage_type&& storage) :
      range_(range.ordinal(*range.begin()) == 0 ? range : range_type(range.lobound(), range.upbound())),
      storage_(storage)
      {
        if (storage_.size() != range_.area())
          array_adaptor<storage_type>::resize(storage_, range_.area());
      }

      /// Construct an evaluated tensor

      /// This constructor will allocate memory for \c range.area() elements. Each element
      /// will be initialized as:
      /// \code
      /// for(int i =Range An input Range type.
      /// \tparam InIter An input iterator type.
      /// \tparam Op A unary operation type
      /// \param range the input range type
      /// \param first An input iterator for the argument
      /// \param op The unary operation to be applied to the argument data
      template <typename Range, typename InIter, typename Op>
      explicit
      Tensor (const Range& range, InIter it, const Op& op,
              typename std::enable_if<btas::is_boxrange<Range>::value>::type* = 0) :
              range_(range.lobound(), range.upbound())
      {
        auto size = range_.area();
        array_adaptor<storage_type>::resize(storage_, size);
        std::transform(it, it+size, begin(), op);
      }

      /// copy constructor
      /// It will accept Tensors and TensorViews
      template<class _Tensor, class = typename std::enable_if<is_boxtensor<_Tensor>::value>::type>
      Tensor (const _Tensor& x)
        :
        range_ (x.range().lobound(), x.range().upbound()),
        storage_(x.cbegin(),x.cend())
      {
      }

      /// copy constructor
      explicit
      Tensor (const Tensor& x)
      : range_ (x.range()), storage_(x.storage_)
      {
      }

      /// move constructor
      Tensor (Tensor&& x)
      : range_ (x.range()), storage_(x.storage_)
      {
      }

      /// copy assignment operator
      template<class _Tensor, class = typename std::enable_if<is_boxtensor<_Tensor>::value>::type>
      Tensor&
      operator= (const _Tensor& x)
      {
          range_ = range_type(x.range().lobound(), x.range().upbound());
          array_adaptor<storage_type>::resize(storage_, range_.area());
          std::copy(std::begin(x), std::end(x), std::begin(storage_));
          return *this;
      }

      /// copy assignment operator
      template<class _Tensor, class = typename std::enable_if<is_boxtensor<_Tensor>::value>::type>
      Tensor&
      operator= (_Tensor&& x)
      {
          range_ = range_type(x.range().lobound(), x.range().upbound());
          storage_ = x.storage();
          return *this;
      }

      /// copy assignment
      Tensor&
      operator= (const Tensor& x)
      {
        range_ = x.range_;
        storage_ = x.storage_;
        return *this;
      }

      /// move assignment operator
      Tensor&
      operator= (Tensor&& x)
      {
        std::swap(range_, x.range_);
        std::swap(storage_, x.storage_);
        return *this;
      }

      /// number of indices (tensor rank)
      size_type
      rank () const
      {
        return range_.rank();
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
        return storage_;
      }

      /// \return storage object
      storage_type&
      storage()
      {
        return storage_;
      }


      ///\name Container requirements (c++std:[container.requirements.general]).
      ///@{

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
        return std::begin(const_cast<const storage_type&>(storage_));
      }

      /// \return const iterator end, even if this is not itself const
      const_iterator
      cend() const
      {
        return std::end(const_cast<const storage_type&>(storage_));
      }

      /// \return iterator begin
      iterator
      begin()
      {
        return std::begin(storage_);
      }

      /// \return iterator end
      iterator
      end()
      {
        return std::end(storage_);
      }

      /// \return number of elements
      size_type
      size () const
      {
        return range_.area();
      }

      /// \return maximum number of elements that can be be contained Tensor
      size_type
      max_size () const
      {
        return std::numeric_limits<size_type>::max();
      }

      /// test whether Tensor is empty
      bool
      empty() const
      {
        return range_.area() == 0;
      }

      /// swap this and x
      void
      swap (Tensor& x)
      {
        std::swap(range_, x.range_);
        std::swap(storage_, x.storage_);
      }

      ///@} // container requirements

      /// @name Element accessors without range check
      /// @{

      /// accesses element using its index, given as a pack of integers
      template<typename index0, typename... _args>
      typename std::enable_if<std::is_integral<index0>::value, const_reference>::type
      operator() (const index0& first, const _args&... rest) const
      {
        typedef typename common_signed_type<index0, typename index_type::value_type>::type ctype;
        auto indexv = {static_cast<ctype>(first), static_cast<ctype>(rest)...};
        index_type index = array_adaptor<index_type>::construct(indexv.size());
        std::copy(std::begin(indexv), std::end(indexv), std::begin(index));
        return storage_[ range_.ordinal(index) ];
      }

      template <typename Index>
      typename std::enable_if<is_index<Index>::value, const_reference>::type
      operator() (const Index& index) const
      {
        return storage_[range_.ordinal(index)];
      }

      template <typename Index>
      typename std::enable_if<is_index<Index>::value, const_reference>::type
      operator[] (const Index& index) const
      {
        return storage_[range_.ordinal(index)];
      }

      /// accesses element using its ordinal value
      /// \param indexord ordinal value of the index
      template <typename IndexOrdinal>
      typename std::enable_if<std::is_integral<IndexOrdinal>::value, const_reference>::type
      operator[] (const IndexOrdinal& indexord) const
      {
        return storage_[indexord];
      }

      template<typename index0, typename... _args>
      typename std::enable_if<std::is_integral<index0>::value, reference>::type
      operator() (const index0& first, const _args&... rest)
      {
        typedef typename common_signed_type<index0, typename index_type::value_type>::type ctype;
        auto indexv = {static_cast<ctype>(first), static_cast<ctype>(rest)...};
        index_type index = array_adaptor<index_type>::construct(indexv.size());
        std::copy(std::begin(indexv), std::end(indexv), std::begin(index));
        return storage_[ range_.ordinal(index) ];
      }

      template <typename Index>
      typename std::enable_if<is_index<Index>::value, reference>::type
      operator() (const Index& index)
      {
        return storage_[range_.ordinal(index)];
      }

      template <typename Index>
      typename std::enable_if<is_index<Index>::value, reference>::type
      operator[] (const Index& index)
      {
        return storage_[range_.ordinal(index)];
      }

      /// accesses element using its ordinal value
      /// \param indexord ordinal value of the index
      template <typename IndexOrdinal>
      typename std::enable_if<std::is_integral<IndexOrdinal>::value, reference>::type
      operator[] (const IndexOrdinal& indexord)
      {
        return storage_[indexord];
      }

      ///@} // element accessors with range check

      /// @name Element accessors with range check
      /// @{

      /// accesses element using its index, given as a pack of integers
      template<typename index0, typename... _args>
      typename std::enable_if<std::is_integral<index0>::value, const_reference>::type
      at (const index0& first, const _args&... rest) const
      {
        typedef typename common_signed_type<index0, typename index_type::value_type>::type ctype;
        auto indexv = {static_cast<ctype>(first), static_cast<ctype>(rest)...};
        index_type index = array_adaptor<index_type>::construct(indexv.size());
        std::copy(std::begin(indexv), std::end(indexv), std::begin(index));
        assert( range_.includes(index) );
        return storage_[ range_.ordinal(index) ];
      }

      template <typename Index>
      typename std::enable_if<is_index<Index>::value, const_reference>::type
      at (const Index& index) const
      {
        assert( range_.includes(index) );
        return storage_[ range_.ordinal(index) ];
      }

      /// accesses element using its ordinal value
      /// \param indexord ordinal value of the index
      template <typename IndexOrdinal>
      typename std::enable_if<std::is_integral<IndexOrdinal>::value, const_reference>::type
      at (const IndexOrdinal& indexord) const
      {
        assert( range_.includes(indexord) );
        return storage_[ indexord ];
      }

      /// accesses element using its index, given as a pack of integers
      template<typename index0, typename... _args>
      typename std::enable_if<std::is_integral<index0>::value, reference>::type
      at (const index0& first, const _args&... rest)
      {
        typedef typename common_signed_type<index0, typename index_type::value_type>::type ctype;
        auto indexv = {static_cast<ctype>(first), static_cast<ctype>(rest)...};
        index_type index = array_adaptor<index_type>::construct(indexv.size());
        std::copy(std::begin(indexv), std::end(indexv), std::begin(index));
        assert( range_.includes(index) );
        return storage_[ range_.ordinal(index) ];
      }

      template <typename Index>
      typename std::enable_if<is_index<Index>::value, reference>::type
      at (const Index& index)
      {
        assert( range_.includes(index) );
        return storage_[ range_.ordinal(index) ];
      }

      /// accesses element using its ordinal value
      /// \param indexord ordinal value of the index
      template <typename IndexOrdinal>
      typename std::enable_if<std::is_integral<IndexOrdinal>::value, reference>::type
      at (const IndexOrdinal& indexord)
      {
        assert( range_.includes(indexord) );
        return storage_[ indexord ];
      }

      ///@} // element accessors with range check

      /// resize array with range object
      template <typename Range>
      void
      resize (const Range& range, typename std::enable_if<is_boxrange<Range>::value,Enabler>::type = Enabler())
      {
        range_ = range_type(range.lobound(),range.upbound());
        array_adaptor<storage_type>::resize(storage_, range_.area());
      }

      /// resize array with extent object
      template <typename Extent>
      void
      resize (const Extent& extent, typename std::enable_if<is_index<Extent>::value &&
                                                            not is_boxrange<Extent>::value,
                                                            Enabler>::type = Enabler())
      {
        range_ = range_type(extent);
        array_adaptor<storage_type>::resize(storage_, range_.area());
      }

      /// clear all members
      void
      clear()
      {
        range_ = range_type();
        storage_ = storage_type();
      }

      //  ========== Finished Public Interface and Its Reference Implementations ==========

      //
      //  Here come Non-Standard members (to be discussed)
      //

      /// Constructs a Tensor slice defined by a subrange for each dimension
      template <typename U>
      TensorView<value_type, range_type, const storage_type>
      slice(std::initializer_list<Range1d<U>> range1s) const
      {
        return TensorView<value_type, range_type, const storage_type>{this->range().slice(range1s), this->storage()};
      }

      /// addition assignment
      Tensor&
      operator+= (const Tensor& x)
      {
        assert( std::equal(range_.begin(), range_.end(), x.range_.begin()) );
        std::transform(std::begin(storage_), std::end(storage_), std::begin(x.storage_), std::begin(storage_), std::plus<value_type>());
        return *this;
      }

      /// addition of tensors
      Tensor
      operator+ (const Tensor& x) const
      {
        Tensor y(*this); y += x;
        return y; /* automatically called move semantics */
      }

      /// subtraction assignment
      Tensor&
      operator-= (const Tensor& x)
      {
        assert(
            std::equal(range_.begin(), range_.end(), x.range_.begin()));
        std::transform(std::begin(storage_), std::end(storage_), std::begin(x.storage_), std::begin(storage_), std::minus<value_type>());
        return *this;
      }

      /// subtraction of tensors
      Tensor
      operator- (const Tensor& x) const
      {
        Tensor y(*this); y -= x;
        return y; /* automatically called move semantics */
      }

      /// \return bare const pointer to the first element of data_
      /// this enables to call BLAS functions
      const value_type*
      data () const
      {
        return storage_.data();
      }

      /// \return bare pointer to the first element of data_
      /// this enables to call BLAS functions
      value_type*
      data()
      {
        return storage_.data();
      }

      /// fill all elements by val
      void
      fill (const value_type& val)
      {
        std::fill(std::begin(storage_), std::end(storage_), val);
      }

      /// generate all elements by gen()
      template<class Generator>
      void
      generate (Generator gen)
      {
          std::generate(std::begin(storage_), std::end(storage_), gen);
      }

    private:

      range_type range_;///< range object
      storage_type storage_;///< data

  }; // end of Tensor

  template <class _Tensor, class = typename std::enable_if<btas::is_boxtensor<_Tensor>::value>::type>
  auto cbegin(const _Tensor& x) -> decltype(x.cbegin()) {
    return x.cbegin();
  }
  template <class _Tensor, class = typename std::enable_if<btas::is_boxtensor<_Tensor>::value>::type>
  auto cend(const _Tensor& x) -> decltype(x.cbegin()) {
    return x.cend();
  }

  /// maps Tensor -> Range
  template <class _Tensor, class = typename std::enable_if<btas::is_boxtensor<_Tensor>::value>::type>
  auto
  range (const _Tensor& t) -> decltype(t.range()) {
    return t.range();
  }

  /// maps Tensor -> Range extent
  template <class _Tensor, class = typename std::enable_if<btas::is_boxtensor<_Tensor>::value>::type>
  auto
  extent (const _Tensor& t) -> decltype(t.range().extent()) {
    return t.range().extent();
  }

  /// maps Tensor -> Range rank
  template <class _Tensor, class = typename std::enable_if<btas::is_boxtensor<_Tensor>::value>::type>
  auto
  rank (const _Tensor& t) -> decltype(t.rank()) {
    return t.rank();
  }

  /// Tensor stream output operator

  /// prints Tensor in row-major form. To be implemented elsewhere using slices.
  /// \param os The output stream that will be used to print \c t
  /// \param t The Tensor to be printed
  /// \return A reference to the output stream
  template <class _Tensor, class = typename std::enable_if<btas::is_boxtensor<_Tensor>::value>::type>
  std::ostream& operator<<(std::ostream& os, const _Tensor& t) {
    os << "Tensor:\n  Range: " << t.range() << std::endl;
    return os;
  }

  /// Tensor comparison operator

  template <class _Tensor1, class _Tensor2,
            class = typename std::enable_if<btas::is_boxtensor<_Tensor1>::value>::type,
            class = typename std::enable_if<btas::is_boxtensor<_Tensor2>::value>::type >
  bool operator==(const _Tensor1& t1, const _Tensor2& t2) {
      if (t1.range().order == t2.range().order &&
          t1.range().ordinal().contiguous() &&
          t2.range().ordinal().contiguous()) // plain Tensor
        return congruent(t1.range(), t2.range()) && std::equal(std::cbegin(t1.storage()),
                                                               std::cend(t1.storage()),
                                                               std::cbegin(t2.storage()));
      else { // not plain, or different orders
        auto cong = congruent(t1.range(), t2.range());
        if (not cong)
          return false;
        typedef TensorView<typename _Tensor1::value_type, typename _Tensor1::range_type, const typename _Tensor1::storage_type>  cview1;
        typedef TensorView<typename _Tensor2::value_type, typename _Tensor2::range_type, const typename _Tensor2::storage_type>  cview2;
        cview1 vt1(t1);
        cview2 vt2(t2);
        return std::equal(vt1.cbegin(), vt1.cend(), vt2.cbegin());
      }
  }

} // namespace btas

namespace boost {
  namespace serialization {

    /// boost serialization
    template<class Archive,
             typename _T,
             class _Range,
             class _Storage>
    void serialize(Archive& ar, btas::Tensor<_T,_Range,_Storage>& t, const unsigned int version) {
      boost::serialization::split_free(ar, t, version);
    }
    template<class Archive,
             typename _T,
             class _Range,
             class _Storage>
    void save(Archive& ar, const btas::Tensor<_T,_Range,_Storage>& t, const unsigned int version) {
      ar << t.range() << t.storage();
    }
    template<class Archive,
             typename _T,
             class _Range,
             class _Storage>
    void load(Archive& ar, btas::Tensor<_T,_Range,_Storage>& t, const unsigned int version) {
      _Range range;
      _Storage storage;
      ar >> range >> storage;
      t = btas::Tensor<_T,_Range,_Storage>(range, storage);
    }

  } // namespace serialization
} // namespace boost

#endif // __BTAS_TENSOR_H
