#ifndef __BTAS_UTIL_OPTIONAL_PTR_H
#define __BTAS_UTIL_OPTIONAL_PTR_H 1

#include <memory>

namespace btas {

/**

  optional_ptr<T> functions either as a raw unmanaged pointer
  or as a smart pointer of type managed_ptr depending on whether
  it is initialized through the set_external method or
  the set_managed method


*/
template <typename T, typename managed_ptr = std::unique_ptr<T>>
class optional_ptr
    {
    public:

    using ptr = T*;

    optional_ptr() : p_(nullptr) { }

    optional_ptr(optional_ptr&& other)
        :
        p_(other.p_),
        up_(std::move(other.up_))
        { }

    T&
    operator*() const { return *p_; }

    ptr
    operator->() const { return p_; }

    void
    set_managed(ptr new_p)
        {
        up_ = std::move(managed_ptr(new_p));
        p_ = up_.get();
        }

    void
    set_external(ptr ext_p)
        {
        p_ = ext_p;
        up_.reset();
        }

    private:
    ptr p_;
    managed_ptr up_;
    };

} // namespace btas

#endif // __BTAS_UTIL_OPTIONAL_PTR_H
