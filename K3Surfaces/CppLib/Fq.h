#include <ctype.h>
#include <assert.h>
#include <iostream>
#include "constants.h"
#include <cstdint>

#ifndef N
#define N 15
#endif

const unsigned q = 1 << N;
const unsigned FINITEFIELDBITSIZE = N;
const unsigned p = polynomials[N];

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

// #endif /////////// End ARM chip check

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

////////////////////////////////////////////////////////////////////////////////
//
// Arf invariants
//
////////////////////////////////////////////////////////////////////////////////

// In the special case b=1, we save a few operations.
int Arf_invariant_b_equals_1(unsigned a, unsigned c) {
  const unsigned  n = FINITEFIELDBITSIZE;
  const unsigned* trace_basis = trace_bases[n];
  const unsigned* pretrace_basis = pretrace_bases[n];

  unsigned acbb  = ff2k_mult(a, c);
  unsigned pivot = 1 << (n-1);
  
  for (int i = 0; i < n; i++) {
    if (acbb & pivot)
      acbb ^= trace_basis[i];
    pivot >>= 1;
  }

  // The trace is zero if and only if what remains is zero.
  return (int)(acbb != 0);
}

// Arf_invariant(unsigned a, unsigned b, unsigned c)
// 
// Given a quadratic polynomial ax^2 + bx + c with b != 0, determine the
// Arf invariant (Equal to 0 or 1).
//
int Arf_invariant(unsigned a, unsigned b, unsigned c) {
  
  unsigned binv2 = ff2k_square(ff2k_inv(b));
  unsigned acbb  = ff2k_mult(ff2k_mult(a, c), binv2);

  return Arf_invariant_b_equals_1(1, acbb);
}

int Arf_invariant_mu2_b_equals_1(unsigned a, unsigned c) {
  return Arf_invariant_b_equals_1(a, c) == 1 ? -1 : 1;
}

int Arf_invariant_mu2(unsigned X, unsigned Y, unsigned Z) {
  return Arf_invariant(X, Y, Z) == 1 ? -1 : 1;
}


////////////////////////////////////////////////////////////////////////////////
//
// Polynomial utilities.
//
////////////////////////////////////////////////////////////////////////////////

void make_monic_nodiv(unsigned* f, int d) {
  // Use a division-free algorithm to make f monic via the transform
  // f(x) -> a^(d-1) f(x/a). It is assumed the leading coefficient of f is non-zero and
  // that d = degree(f).

  unsigned powa = 1;
  for (int i=1; i<=d; i++) {
    f[d-i] = ff2k_mult(f[d-i], powa);
    powa = ff2k_mult(powa, f[d]);
  }
  f[d] = 1;
  return;
}

int gcd_degree(unsigned* A, unsigned* B, int dA, int dB) {

  int degA = dA; // Note: degA must be accurate. degB is just a bound.
  int degB = dB;

  // Possibly correct the degree of B.
  while (degB >= 0) {
    if (B[degB] != 0) break;
    degB--;
  }
  
  /////////////////
  // Do Euclid's algorithm.
  
  while (degB > 0) {

    // Main Euclid loop.
    for (int i=0; i <=degA - degB; i++) {
      int lt = A[degA-i];

      // Scale up A. (Instead of making B monic.)
      // Note we deliberately avoid touching the leading coefficient, since
      // we still need this value for the reductions.
      for (int j=0; j<degA-i; j++) {
        A[j] = ff2k_mult(A[j], B[degB]);
      }

      // Perform the reduction.
      for (int j=1; j <= degB; j++) {
        A[degA-i-j] ^= ff2k_mult(B[degB-j], lt);
      }
      A[degA-i] = 0; // Discard the leading coefficient, which is theoretically killed.
    }
    
    // Swap
    unsigned* C = A;
    A = B;
    B = C;
    degA = degB;
    degB = degB - 1;

    while (degB >= 0) {
      if ( B[degB] != 0) break;
      degB--;
    }
  }

  // Final count.
  return (degB < 0) ? degA : 0;
}


////////////////////////////////////////////////////////////////////////////////
//
// Polynomial root counting.
//
////////////////////////////////////////////////////////////////////////////////

int check_gcd(unsigned* poly) {


  // The key trick is to quickly compute the GCD with x^q-x, then check the degree.
  int deg = 3;

  // Declare f. Note this mutates poly.
  const unsigned n = FINITEFIELDBITSIZE;
  unsigned *f = poly;
  
  // q = 2^n.
  // Compute x^q mod f(x). 
  unsigned powx[deg];
  unsigned sqpowx[2*deg-1];
 
  // Initialize
  for (int ll=0; ll <= deg-1; ll++) {
    powx[ll] = 0;
  }
  for (int ll=0; ll <= 2*deg-2; ll++) {
    sqpowx[ll] = 0;
  }          
  powx[1] = 1;

  // Square-Reduce loop
  for (int i=0; i < n; i++) {
    // Square
    for (int j=0; j<=deg-1; j++) {
      sqpowx[2*j] = ff2k_square(powx[j]);
    }

    // Reduce
    for (int j = deg-2; j>=0; j--) {
      if (sqpowx[j+deg] != 0) {
        for (int ll=0; ll < deg; ll++) {
          sqpowx[j+ll] ^= ff2k_mult(sqpowx[j+deg], f[ll]);
        }
      }
    }

    // Repeat
    for (int ll=0; ll <= deg-1; ll++) {
      powx[ll] = sqpowx[ll];
    }
    for (int ll=0; ll <= 2*deg-2; ll++) {
      sqpowx[ll] = 0;
    }          
  }

  // Subtract x
  powx[1] ^= 1;

  // Compute the GCD degree normally from this point.
  return gcd_degree(f, powx, deg, deg-1);
}

int count_quadratic_roots(unsigned A, unsigned B, unsigned C) {
  // Return an array {n, r1, r2}. The first number indicates the number
  // of distinct roots. The remaining values are the actual roots.
  //
  // It is assumed that f is genuinely a quadratic.
  const int n = FINITEFIELDBITSIZE;
  
  if (A == 0) {
    if (B == 0) {
      if (C == 0) {
        return q;
      } else {
        return 0;
      }
    } else {
      return 1;
    }  
  } else if (B == 0) { // The quadratic has a double root.
      return 1;
  } else {
      int arf = Arf_invariant(A, B, C);
      return (arf ^ 1) << 1;
  }
  return 0;
}

int count_cubic_roots(unsigned A, unsigned B, unsigned C, unsigned D) {


    // unsigned count = 0;
    // for (int i = 0; i < q; i++) {
    //     if ((ff2k_mult(A, ff2k_mult(i, ff2k_square(i))) ^ ff2k_mult(B, ff2k_square(i)) ^ ff2k_mult(C, i) ^ D) == 0) {
    //         count++;
    //     }
    // }
    // return count;
    
    if (A == 0) { // In this case, we have an extra solution at infinity
        // Equation is a quadratic
        return 1+count_quadratic_roots(B, C, D);
    } else {

        ff2k_t P, Q, R;
        // Check if already depressed:
        if (A == 1 & B == 0) {
            P = C;
            Q = D;
        } else {
            // Make depressed
            unsigned B2 = ff2k_divi(B,A);
            P = ff2k_mult(B2, B2) ^ ff2k_divi(C, A);
            Q = ff2k_mult(B2, ff2k_divi(C, A)) ^ ff2k_divi(D,A);            
        }

        // unsigned count = 0;
        // for (int i = 0; i < q; i++) {
        //     if ((ff2k_mult(i, ff2k_square(i)) ^ ff2k_mult(P, i) ^ Q) == 0) {
        //         count++;
        //     }
        // }
        // return count;
        

        if (P == 0 & Q == 0) {
            return 1;
        } else if (Q == 0) {
            return 2;
        } else {

            // Check deg of gcd with x^q - x
            unsigned coeffs[4] = {Q, P, 0, 1};
            unsigned *f = coeffs;
            
            return check_gcd(f);
        }
    } 

    return 0;
}




////////////////////////////////////////////////////////////////////////////////
//
// Printing.
//
////////////////////////////////////////////////////////////////////////////////

// REMARK: Only works for nonnegative numbers.
std::string int128_to_string(__int128 num) {
    std::string str;
    do {
        int digit = num % 10;
        str = std::to_string(digit) + str;
        num = (num - digit) / 10;
    } while (num != 0);
    return str;
}