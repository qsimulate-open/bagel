#ifndef __BTAS_UTIL_SEQUENCEADAPTOR_H_
#define __BTAS_UTIL_SEQUENCEADAPTOR_H_

namespace btas {

  /// infinite_sequence_adaptor represents pointer \c ptr as a \c ptr[0] , \c ptr[1] , .. sequence.
  /// Because the sequence is infinite, several attributes of sequence are not supported (end, cend, size, resize)
  /// \tparam _Ptr pointer type
  template <typename _Ptr,
            class = typename std::enable_if<std::is_pointer<_Ptr>::value>::type >
  class infinite_sequence_adaptor {
    public:
      typedef typename std::remove_pointer<_Ptr>::type value_type;
      typedef _Ptr pointer;
      typedef typename std::add_const<pointer>::type const_pointer;
      typedef          value_type& reference;
      typedef    const value_type& const_reference;
      typedef size_t size_type;
      typedef ptrdiff_t difference_type;

      typedef pointer iterator;
      typedef const_pointer const_iterator;

      infinite_sequence_adaptor() : ptr_(nullptr) {}
      infinite_sequence_adaptor(pointer ptr) : ptr_(ptr) {}
      infinite_sequence_adaptor(const infinite_sequence_adaptor&) = default;
      infinite_sequence_adaptor(infinite_sequence_adaptor&&) = default;
      ~infinite_sequence_adaptor() {}

      const_iterator cbegin() const {
        return const_cast<const_iterator>(ptr_);
      }
      const_iterator begin() const {
        return cbegin();
      }
      iterator begin() {
        return ptr_;
      }

      const_reference operator[](size_type i) const {
        return ptr_[i];
      }
      reference operator[](size_type i) {
        return ptr_[i];
      }


    private:
      pointer ptr_;
  }; // class infinite_sequence_adaptor

  template <typename _T>
  infinite_sequence_adaptor<_T*>
  make_infinite_sequence_adaptor(_T* data) {
    return infinite_sequence_adaptor<_T*>(data);
  }

} // namespace btas

#include <functional>

namespace std {

  // re-implement reference_wrapper<infinite_sequence_adaptor<Ptr>> to act like Ptr
  template <typename _T>
  class reference_wrapper<btas::infinite_sequence_adaptor<_T*>> : private btas::infinite_sequence_adaptor<_T*> {
    public:

      // types
      typedef btas::infinite_sequence_adaptor<_T*> type;

      // construct/copy/destroy
      reference_wrapper(type& x) noexcept : type(x) {}
      // DO bind to temps
      reference_wrapper(type&& x) noexcept : type(x) {}
      reference_wrapper(const reference_wrapper<type>& x) noexcept : type(x) {}

      // assignment
      reference_wrapper& operator=(const reference_wrapper<type>& x) noexcept {
        static_cast<type&>(*this) = static_cast<type&>(x);
      }

      // access
      operator type& () noexcept { return *this; }
      type& get() noexcept { return *this; }
      const type& get() const noexcept { return *this; }

  };

  // re-implement reference_wrapper<const infinite_sequence_adaptor<Ptr>> to act like const Ptr
  template <typename _T>
  class reference_wrapper<const btas::infinite_sequence_adaptor<_T*>> : private btas::infinite_sequence_adaptor<_T*> {
    public:

      // types
      typedef const btas::infinite_sequence_adaptor<_T*> ctype;
      typedef btas::infinite_sequence_adaptor<_T*> nctype;

      // construct/copy/destroy
      reference_wrapper(ctype& x) noexcept : ctype(x) {}
      // DO bind to temps
      reference_wrapper(nctype&& x) noexcept : ctype(x) {}
      reference_wrapper(const reference_wrapper<ctype>& x) noexcept : ctype(x) {}

      // assignment
      reference_wrapper& operator=(const reference_wrapper<ctype>& x) noexcept {
        static_cast<ctype&>(*this) = static_cast<ctype&>(x);
      }

      // access
      operator ctype& () const noexcept { return *this; }
      ctype& get() const noexcept { return *this; }

  };

  template <class _T>
  reference_wrapper<btas::infinite_sequence_adaptor<_T*>> ref(btas::infinite_sequence_adaptor<_T*>&& t) {
      return reference_wrapper<btas::infinite_sequence_adaptor<_T*>>(t);
  }
  template <class _T>
  reference_wrapper<const btas::infinite_sequence_adaptor<_T*>> cref(const btas::infinite_sequence_adaptor<_T*>&& t) {
      return reference_wrapper<const btas::infinite_sequence_adaptor<_T*>>(t);
  }

}

#endif /* __BTAS_UTIL_SEQUENCEADAPTOR_H_ */
