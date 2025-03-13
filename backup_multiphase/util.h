#ifndef OPT_UTIL_H
#define OPT_UTIL_H

#include <array>
#include <memory>
#include <vector>
#include <limits>

// TODO: remove me? -> only use case for long double precision 
// may be numerically ill conditioned operations with high order Radau schemes 
typedef double gNumber;
typedef std::vector<gNumber> gVector;

template <std::size_t SIZE>
using gArray = std::array<gNumber, SIZE>;

// update me?
const gNumber PLUS_INFINITY = std::numeric_limits<gNumber>::infinity();
const gNumber MINUS_INFINITY = -std::numeric_limits<gNumber>::infinity();

#endif  // OPT_UTIL_H
