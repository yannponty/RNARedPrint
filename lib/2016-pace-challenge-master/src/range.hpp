#pragma once
#include <cstddef>  // for size_t
#include <utility>

template <class Iter>
class IterWrap {
 private:
  Iter s_, e_;

 public:
  IterWrap(Iter s, Iter e) : s_(s), e_(e) {}
  Iter begin() const { return s_; }
  Iter end() const { return e_; }
};
template <class Iter>
inline IterWrap<Iter> as_range(Iter s, Iter e) {
  return {s, e};
}

template <class Iter>
inline IterWrap<Iter> as_range(std::pair<Iter, Iter> pair) {
  return {pair.first, pair.second};
}

template <class T, size_t incr>
class RangeIterator {
 private:
  T t_;

 public:
  RangeIterator(T t) : t_(t) {}
  T operator*() const { return t_; }
  RangeIterator<T, incr> &operator++() {
    t_ += incr;
    return *this;
  }
  bool operator==(const RangeIterator<T, incr> &other) const {
    return t_ == other.t_;
  }
  bool operator!=(const RangeIterator<T, incr> &other) const {
    return t_ != other.t_;
  }
};

template <class T, size_t incr = 1>
using Range = IterWrap<RangeIterator<T, incr>>;

inline Range<size_t, 1> range(size_t i) { return Range<size_t, 1>(0, i); }
