#pragma once

/** Defines always_assert macro that runs
 *  regardless of DEBUG or NDEBUG setting.
 */

#include <cstdio>   // for fprintf, stderr
#include <cstdlib>  // for std::abort

#define always_assert(X)                                               \
  do {                                                                 \
    if (!(((X)))) {                                                    \
      fprintf(stderr, "%s:%d in %s: assertion failed: %s\n", __FILE__, \
              __LINE__, __func__, #X);                                 \
      std::abort();                                                    \
    }                                                                  \
  } while (0)
