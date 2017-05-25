#ifndef __BTAS_VARRAY_H
#define __BTAS_VARRAY_H 1

#include <algorithm>
#include <cassert>
#include <btas/serialization.h>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/collection_size_type.hpp>
#if BOOST_VERSION >= 106400
#include <boost/serialization/array_wrapper.hpp>
#endif

namespace btas {

/// variable size array class without capacity info
template <typename _T,
          typename _Allocator = std::allocator<_T> >
class varray : private _Allocator {
public:

   typedef std::allocator_traits<_Allocator> allocator_traits;            ///< Allocator traits
   typedef typename allocator_traits::allocator_type allocator_type;      ///< Allocator type

   typedef typename allocator_traits::value_type value_type;
   typedef          value_type& reference;
   typedef    const value_type& const_reference;
   typedef typename allocator_traits::pointer pointer;
   typedef typename allocator_traits::const_pointer const_pointer;
   typedef typename allocator_traits::difference_type difference_type;
   typedef typename allocator_traits::size_type size_type;

   typedef pointer iterator;
   typedef const_pointer const_iterator;
   typedef iterator reverse_iterator;
   typedef const_iterator const_reverse_iterator;

private:

   struct _M_impl {
     pointer _M_start;
     pointer _M_finish;

     _M_impl() : _M_start(nullptr), _M_finish(nullptr) { }
     _M_impl(pointer s, pointer f) : _M_start(s), _M_finish(f) { }
     _M_impl(_M_impl&& other) {
       _M_start = other._M_start;
       _M_finish = other._M_finish;
       other._M_start = nullptr;
       other._M_finish = nullptr;
     }

     size_type size() const { return _M_finish - _M_start; }
     bool empty() const { return _M_start == _M_finish; }

     pointer begin() { return _M_start; }
     pointer end() { return _M_finish; }
     const_pointer cbegin() const { return const_cast<const_pointer>(_M_start); }
     const_pointer cend() const { return const_cast<const_pointer>(_M_finish); }

     pointer rbegin() { return _M_finish-1; }
     pointer rend() { return _M_start-1; }
     const_pointer crbegin() const { return const_cast<const_pointer>(_M_finish-1); }
     const_pointer crend() const { return const_cast<const_pointer>(_M_start-1); }

     reference front() { return *begin(); }
     reference back() { return *rbegin(); }
     const_reference front() const { return *cbegin(); }
     const_reference back() const { return *crbegin(); }

     reference operator[](size_type i) {
       return _M_start[i];
     }
     const_reference operator[](size_type i) const {
       return const_cast<const_reference>(_M_start[i]);
     }
     reference at(size_type i) {
       assert(i < size());
       return _M_start[i];
     }
     const_reference at(size_type i) const {
       assert(i < size());
       return const_cast<const_reference>(_M_start[i]);
     }

     pointer data() {
       return _M_start;
     }
     const_pointer data() const {
       return const_cast<const_pointer>(_M_start);
     }

     void swap (_M_impl& other)
     {
       _M_start = other._M_start;
       _M_finish = other._M_finish;
       other._M_start = nullptr;
       other._M_finish = nullptr;
     }
   };
   _M_impl data_;

   allocator_type& alloc() { return static_cast<allocator_type&>(*this); }
   const allocator_type& alloc() const { return static_cast<const allocator_type&>(*this); }

   friend class boost::serialization::access;

public:

   varray () : allocator_type()
   { }

   explicit
   varray (const allocator_type& a) : allocator_type(a)
   { }

   ~varray ()
   {
     deallocate();
   }

   explicit
   varray (size_type n, const allocator_type& a = allocator_type()) : allocator_type(a)
   {
     if (n > 0) { // this ensures that if n == 0, pointers are null
       allocate(n);
       construct(n);
     }
   }

   varray (size_type n, const_reference val,
           const allocator_type& a = allocator_type()) : allocator_type(a)
   {
     if (n > 0) {
       allocate(n);
       construct(n, val);
     }
   }

   template <class InputIterator>
   varray (InputIterator first, InputIterator last)
   {
     const auto n = std::distance(first, last);
     if (n > 0) {
       allocate(n);
       construct(first, last);
     }
   }

   varray (const varray& x) : allocator_type(x)
   {
     const auto n = x.size();
     if (n > 0) {
       allocate(n);
       construct(x.cbegin(), x.cend());
     }
   }

   varray (const varray& x, const allocator_type& a) : allocator_type(a)
   {
     const auto n = x.size();
     if (n > 0) {
       allocate(n);
       construct(x.cbegin(), x.cend());
     }
   }

   varray (varray&& x)
   : allocator_type(std::move(static_cast<allocator_type&&>(x))), data_(std::move(x.data_))
   {
   }

   template <typename U, class = typename std::enable_if< std::is_convertible<U, value_type>::value >::type >
   varray (std::initializer_list<U> il)
   {
      size_type n = il.size();
      if (n > 0) {
        allocate(n);
        construct(il.begin(), il.end());
      }
   }

   varray& operator= (const varray& x) {
     const auto n = x.size();
     if (n != data_.size()) {
       deallocate();
       if (n > 0)
         allocate(n);
     }

     if (n > 0) {
       construct(x.cbegin(), x.cend());
     }
     return *this;
   }

   varray& operator= (varray&& x)
   {
      swap (x); // if something in this object, it will be destructed by x
      return *this;
   }

   template <typename U, class = typename std::enable_if< std::is_convertible<U, value_type>::value >::type >
   varray& operator= (std::initializer_list<U> il)
   {
       const auto n = il.size();
       if (n != data_.size()) {
         deallocate();
         if (n > 0)
           allocate(n);
       }

       if (n > 0) {
         construct(il.cbegin(), il.cend());
       }
       return *this;
   }

   iterator begin () noexcept
   {
      return data_.begin();
   }

   const_iterator begin () const noexcept
   {
      return cbegin();
   }

   const_iterator cbegin () const noexcept
   {
      return data_.cbegin();
   }

   iterator end () noexcept
   {
      return data_.end();
   }

   const_iterator end () const noexcept
   {
      return cend();
   }

   const_iterator cend () const noexcept
   {
      return data_.cend();
   }

   reverse_iterator rbegin () noexcept
   {
      return data_.rbegin();
   }

   const_reverse_iterator rbegin () const noexcept
   {
      return data_.rbegin();
   }

   reverse_iterator rend () noexcept
   {
      return data_.rend();
   }

   const_reverse_iterator rend () const noexcept
   {
      return data_.rend();
   }

   size_type size () const noexcept
   { return data_.size(); }

   void resize (size_type n)
   {
     if (size() != n) {
       if (!empty()) {
         deallocate();
       }
       if (n > 0) {
         allocate(n);
       }
     }
   }

   void resize (size_type n, const value_type& val)
   {
     resize(n);
     construct(n, val);
   }

   bool empty () const noexcept
   { return data_.empty(); }

   reference operator [] (size_type n)
   { return data_[n]; }

   const_reference operator [] (size_type n) const
   { return data_[n]; }

   reference at (size_type n)
   { return data_.at(n); }

   const_reference at (size_type n) const
   { return data_.at(n); }

   reference front ()
   { return data_.front(); }

   const_reference front () const
   { return data_.front(); }

   reference back ()
   { return data_.back(); }

   const_reference back () const
   { return data_.back(); }

   value_type* data () noexcept
   { return data_.data(); }

   const value_type* data () const noexcept
   { return data_.data(); }

   void swap (varray& x)
   { data_.swap(x.data_); }

   void clear ()
   {
      if (!empty()) {
        deallocate();
      }
   }

  private:

   void allocate(size_type n) {
     assert(n <= allocator_traits::max_size(alloc()));
     data_._M_start = allocator_traits::allocate(alloc(), n);
     data_._M_finish = data_._M_start + n;
   }

   void deallocate() {
     if (!data_.empty())
       allocator_traits::deallocate(alloc(), data_._M_start, data_.size());
     data_._M_start = data_._M_finish = nullptr;
   }

   void
   construct(size_type n)
   {
     auto ptr = data_._M_start;
       do
       {
           allocator_traits::construct(alloc(), ptr);
           ++ptr;
           --n;
       } while (n > 0);
   }

   void
   construct(size_type n, const_reference x)
   {
     auto ptr = data_._M_start;
     do
     {
       allocator_traits::construct(alloc(), ptr, x);
       ++ptr;
       --n;
     } while (n > 0);
   }

   template <typename InputIterator>
   void
   construct(const InputIterator& begin, const InputIterator& end) {
     auto ptr = data_._M_start;
     for(auto i = begin; i != end; ++i) {
       allocator_traits::construct(alloc(), ptr, *i);
       ++ptr;
     }
   }
};

template <typename T>
inline bool operator== (const btas::varray<T>& a,
                        const btas::varray<T>& b) {
  return std::equal(a.begin(), a.end(), b.begin());
}

template <typename T>
inline bool operator!= (const btas::varray<T>& a,
                        const btas::varray<T>& b) {
  return not (a == b);
}

} // namespace btas

namespace boost {
  namespace serialization {

  /// boost serialization for varray
  template<class Archive, typename T>
  void serialize (Archive& ar, btas::varray<T>& x, const unsigned int version)
  {
      boost::serialization::split_free(ar, x, version);
  }
  template<class Archive, typename T>
  void save (Archive& ar, const btas::varray<T>& x, const unsigned int version)
  {
      const boost::serialization::collection_size_type count(x.size());
      ar << BOOST_SERIALIZATION_NVP(count);
      if (count != 0)
        ar << boost::serialization::make_array(x.data(), count);
  }
  template<class Archive, typename T>
  void load (Archive& ar, btas::varray<T>& x, const unsigned int version)
  {
      boost::serialization::collection_size_type count;
      ar >> BOOST_SERIALIZATION_NVP(count);
      x.resize(count);
      if (count != 0)
        ar >> boost::serialization::make_array(x.data(), count);
  }

  } // namespace serialization
} // namespace boost

template <typename T>
inline bool operator== (const btas::varray<T>& a,
                        const btas::varray<T>& b) {
  return std::equal(a.begin(), a.end(), b.begin());
}

template <typename T>
inline bool operator!= (const btas::varray<T>& a,
                        const btas::varray<T>& b) {
  return not (a == b);
}

#endif // __BTAS_VARRAY_H
