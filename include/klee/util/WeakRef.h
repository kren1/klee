//===-- Ref.h ---------------------------------------------------*- C++ -*-===//
//
//                     The KLEE Symbolic Virtual Machine
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//

/**
 * @file Ref.h
 * @brief Implements smart-pointer weak_ref<> used by KLEE.
 *
 * ## Basic usage:
 *
 * Add the following to your struct/class to enable weak_ref<> pointer usage
 * @code{.cpp}
 *
 * struct MyStruct{
 *   ...
 *   /// @brief Required by klee::ref-managed objects
 *   class ReferenceCounter _refCount;
 *   ...
 * }
 * @endcode
 *
 */

#ifndef KLEE_WEAKREF_H
#define KLEE_WEAKREF_H

#include "Ref.h"
namespace klee {
template<class T>
class weak_ref {
  T *ptr;

public:
  // default constructor: create a NULL weak_reference
  weak_ref() : ptr(nullptr) {}
  ~weak_ref () { dec (); }

private:
  void inc() const {
    if (ptr)
      ++ptr->_refCount.weakrefCount;
  }

  void dec() const {
    if (ptr && --ptr->_refCount.weakrefCount == 0 && ptr->_refCount.refCount == 0)
      delete ptr;
  }

public:
  template<class U> friend class weak_ref;

  // constructor from pointer
  weak_ref(T *p) : ptr(p) {
    inc();
  }

  // normal copy constructor
  weak_ref(const weak_ref<T> &r) : ptr(r.ptr) {
    inc();
  }

  // conversion constructor
  template<class U>
  weak_ref (const weak_ref<U> &r) : ptr(r.ptr) {
    inc();
  }

  // normal move constructor
  weak_ref(weak_ref<T> &&r) noexcept : ptr(std::move(r.ptr)) { r.ptr = nullptr; }

  // conversion move constructor
  template <class U> weak_ref(weak_ref<U> &&r) noexcept : ptr(std::move(r.ptr)) {
    r.ptr = nullptr;
  }

  // pointer operations
  T *get () const {
    return ptr;
  }

  /* The copy assignment operator must also explicitly be defined,
   * despite a redundant template. */
  weak_ref<T> &operator= (const weak_ref<T> &r) {
    r.inc();
    // Create a copy of the pointer as the
    // referenced object might get destroyed by the following dec(),
    // like in the following example:
    // ````````````````````````
    //    Expr {
    //        ref<Expr> next;
    //    }
    //
    //    ref<Expr> root;
    //    root = root->next;
    // ````````````````````````
    T *saved_ptr = r.ptr;
    dec();
    ptr = saved_ptr;

    return *this;
  }

  template<class U> weak_ref<T> &operator= (const weak_ref<U> &r) {
    r.inc();
    // Create a copy of the pointer as the currently
    // referenced object might get destroyed by the following dec(),
    // like in the following example:
    // ````````````````````````
    //    Expr {
    //        weak_ref<Expr> next;
    //    }
    //
    //    weak_ref<Expr> root;
    //    root = root->next;
    // ````````````````````````

    U *saved_ptr = r.ptr;
    dec();
    ptr = saved_ptr;

    return *this;
  }

  // Move assignment operator
  weak_ref<T> &operator=(weak_ref<T> &&r) noexcept {
    using std::swap;
    swap(ptr, r.ptr);
    return *this;
  }

  // Move assignment operator
  template <class U> weak_ref<T> &operator=(weak_ref<U> &&r) {
    // swap
    auto *tmp = cast_or_null<T>(r.ptr);
    r.ptr = cast_or_null<U>(ptr);
    ptr = tmp;
    return *this;
  }

  T& operator*() const {
    return *ptr;
  }

  T* operator->() const {
    return ptr;
  }

  bool isNull() const { return ptr == nullptr; }

  // assumes non-null arguments
  int compare(const weak_ref &rhs) const {
    assert(!isNull() && !rhs.isNull() && "Invalid call to compare()");
    return get()->compare(*rhs.get());
  }

  // assumes non-null arguments
  bool operator<(const weak_ref &rhs) const { return compare(rhs)<0; }
  bool operator==(const weak_ref &rhs) const { return compare(rhs)==0; }
  bool operator!=(const weak_ref &rhs) const { return compare(rhs)!=0; }
};

template<class T>
inline llvm::raw_ostream &operator<<(llvm::raw_ostream &os, const weak_ref<T> &e) {
  os << *e;
  return os;
}

template<class T>
inline std::stringstream &operator<<(std::stringstream &os, const weak_ref<T> &e) {
  os << *e;
  return os;
}

} // end namespace klee

namespace llvm {
  // simplify_type implementation for weak_ref<>, which allows dyn_cast from on a
  // weak_ref<> to apply to the wrapper type. Conceptually the result of such a
  // dyn_cast should probably be a weak_ref of the casted type, but that breaks the
  // idiom of initializing a variable to the result of a dyn_cast inside an if
  // condition, or we would have to implement operator(bool) for weak_ref<> with
  // isNull semantics, which doesn't seem like a good idea.
template<typename T>
struct simplify_type<const ::klee::weak_ref<T> > {
  using SimpleType = T *;
  static SimpleType getSimplifiedValue(const ::klee::weak_ref<T> &Ref) {
    return Ref.get();
  }
};

template<typename T>
struct simplify_type< ::klee::weak_ref<T> >
  : public simplify_type<const ::klee::weak_ref<T> > {};
}  // namespace llvm

#endif /* KLEE_REF_H */
