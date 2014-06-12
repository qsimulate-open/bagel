/*
 * allocators.h
 *
 *  Created on: Jan 4, 2014
 *      Author: evaleev
 */

#ifndef BTAS_VARRAY_ALLOCATORS_H_
#define BTAS_VARRAY_ALLOCATORS_H_

#include <cstddef>
#include <algorithm>
#include <memory>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <iomanip>

namespace btas {

  struct stack_arena {
      size_t size;            //!< in bytes
      char* const buffer;     //!< buffer begin
      char* current;          //!< ptr to the free space in buffer

      template<typename T>
      stack_arena(T* b, size_t s) :
        size(s),
        buffer(reinterpret_cast<char*>(b)),
        current(reinterpret_cast<char*>(b))
        {
      }

      void increment(const size_t n) { current += n; }
      void decrement(const size_t n) { current -= n; }
  };


  /// This is a very simple allocator implementation that uses an externally managed memory stack.
  /// It's mostly for demonstration purposes.

  /// stack_allocator is a first-in-last-out allocator,
  /// i.e. deallocation of the memory must happen in the opposite order of allocation.
  template<typename T>
  class stack_allocator {
    private:
      std::shared_ptr<stack_arena> arena_;

      size_t size() const {
        return arena_->size / sizeof(T);
      }
      T* buffer() const {
        return reinterpret_cast<T*>(arena_->buffer);
      }
      T* current() const {
        return reinterpret_cast<T*>(arena_->current);
      }

      void increment(const size_t n) {
        arena_->increment(n * sizeof(T));
      }
      void decrement(const size_t n) {
        arena_->decrement(n * sizeof(T));
      }

    public:
      // do not allocate here
      stack_allocator(std::shared_ptr<stack_arena> a) :
          arena_(a) {
      }

      stack_allocator(const stack_allocator<T>& o) :
          arena_(o.arena_) {
      }
      stack_allocator& operator=(const stack_allocator<T>& o) {
        arena_ = o.arena_;
        return *this;
      }
      struct rebind {
          typedef stack_allocator<T> other;
      };

      typedef T* pointer;
      typedef T& reference;
      typedef const T* const_pointer;
      typedef const T& const_reference;
      typedef T value_type;
      typedef size_t size_type;
      typedef ptrdiff_t difference_type;

      pointer allocate(size_type n, const void* = 0) {
        //std::cout << "Allocating " << std::setw(6) << sizeof(T)*n << " bytes. "
        //          << "Available memory: "
        //          << std::setw(6) << sizeof(T) * (size() - std::distance(buffer(), current())) << " bytes."<< std::endl;
        pointer out = current();
        increment(n);
        if (std::distance(buffer(), current()) > size())
          throw std::runtime_error("preallocated memory exhausted");
        return out;
      }

      pointer address(reference x) const {
        return &x;
      }
      const_pointer address(const_reference x) const {
        return &x;
      }

      void deallocate(pointer p, size_type n) {
        assert(p == current() - n);
        decrement(n);
      }

      void construct(pointer p, const T& val) {
        *p = val;
      }
      void destroy(pointer p) {
      }

      size_type max_size() const {
        return size();
      }
  };

}


#endif /* BTAS_VARRAY_ALLOCATORS_H_ */
