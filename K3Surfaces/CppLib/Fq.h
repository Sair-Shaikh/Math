#include <ctype.h>
#include <assert.h>
#include <iostream>
#include <tuple>
#include "constants.h"


int initialized_ff_bitsize() {
  return FINITEFIELDBITSIZE;
}

ff2k_t ff2k_mult(ff2k_t x, ff2k_t y) {
  // Stock multiplication algorithm.
  unsigned a = x, b = y, ab = 0;

  while (b != 0) {
    // b = const + higher-deg part
    // if const = 1 then ret += a
    if (b & 1) ab ^= a;
    // kill const, divide higher-deg part by t
    b >>= 1;
    // multiply a by t
    a <<= 1;
    if (a & q) a ^= p;
  }

  return ab;
}

// For evaluating complicated polynomials it can help to delay reduction.
static inline uint64_t _ff2k_mult_delay_reduction(unsigned a, unsigned b) {
  return ff2k_mult(a, b); 
}

static inline ff2k_t _ff2k_reduce(uint64_t a) {
  return a;
}

inline uint64_t _ff2k_square_delay_reduction(ff2k_t a) {
  return ff2k_mult(a, a);
}

inline ff2k_t ff2k_square(ff2k_t a) {
  return ff2k_mult(a, a);
}

// NOTE: No division by zero check.
inline ff2k_t ff2k_inv(ff2k_t a) {
  const int n = FINITEFIELDBITSIZE;
  assert(a!=0);
  // Use Lagrange to invert. that is, compute a^(2^n-2).
  // Note the binary expansion of 2^n-2 = 111...110.
  ff2k_t inva = 1;
  ff2k_t sqac = ff2k_square(a);
  for (int i=1; i<n; i++) {
    inva = ff2k_mult(inva, sqac);
    sqac = ff2k_square(sqac);
  }
  return inva;
}

ff2k_t ff2k_divi(ff2k_t a, ff2k_t b) {
  return ff2k_mult(a, ff2k_inv(b));
}
ff2k_t ff2k_sqrt(ff2k_t a) {
  const int n = FINITEFIELDBITSIZE;
  ff2k_t s = a;
  for (int i=1; i<=n-1; i++) {
    s = ff2k_square(s); // Compute a^(2^(n-1)).
  }
  return s;
}



