#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <tuple>
#include "tableio.h"

std::tuple<unsigned**, unsigned**, unsigned**, unsigned**, unsigned*, unsigned*> generate_all_tables(unsigned, unsigned);

int main() {

  unsigned** quadratic_roots;
  unsigned** depressed_cubic_roots;
  unsigned** mult;
  unsigned** divi;
  unsigned* orbit_rep;
  unsigned* orbit_size;

  for (int N = 1; N < 12; N++) {

    const unsigned q = 1 << N;
    const unsigned p = polynomials[N];

    std::tie(quadratic_roots, depressed_cubic_roots, mult, divi, orbit_rep, orbit_size) = generate_all_tables(q, p);

    // Write data to file.
    std::string qq = std::to_string(q);

    printf("deb2: %d \n", depressed_cubic_roots[0][1]);


    write_table(quadratic_roots, q, q, "quadratic_roots_" + qq);
    write_table(depressed_cubic_roots, q, q, "depressed_cubic_roots_" + qq);
    write_table(mult, q, q, "mult_" + qq);
    write_table(divi, q, q, "divi_" + qq);
    write_table(orbit_rep, q, "orbit_rep_" + qq);
    write_table(orbit_size, q, "orbit_size_" + qq);

  }

  return 0;
}

std::tuple<unsigned**, unsigned**, unsigned**, unsigned**, unsigned*, unsigned*> generate_all_tables(unsigned q, unsigned p) {

  unsigned **mult, **divi;
  unsigned **quadratic_roots;
  unsigned **depressed_cubic_roots;


  // Generate multiplication and division tables
  mult = new unsigned*[q];
  divi = new unsigned*[q];
  for (unsigned i = 0; i < q; i++) {
    mult[i] = new unsigned[q];
    divi[i] = new unsigned[q];
  }

  // Fill multiplication and division tables
  for (unsigned i = 0; i < q; i++) {
    for (unsigned j = 0; j <= i; j++) {
      
      // main multiplication algorithm
      unsigned a = i, b = j, ij = 0;
      while (b != 0) {
        // b = const + higher-deg part
        if (b & 1) ij ^= a;
        // kill const, divide higher-deg part by x
        b >>= 1;
        // multiply a by x
        a <<= 1;
        if (a & q) a ^= p;
            }
        mult[i][j] = ij;
        mult[j][i] = ij;
        divi[ij][i] = j;
        divi[ij][j] = i;
    }
  }

  // Initialize table for number of roots of monic quadratic 
  quadratic_roots = new unsigned*[q];
  for (int i = 0; i < q; i++) {
    quadratic_roots[i] = new unsigned[q];
    for (int j = 0; j < q; j++) {
      quadratic_roots[i][j] = 0;
    }
  }
  // Fill the quadratics table
  for (int i = 0; i < q; i++) {
    for (int j = i; j < q; j++) {
      // x^2 + ax + b = (x+i)(x+j)
      unsigned a = i ^ j, b = mult[i][j];
      quadratic_roots[a][b] = 1;
      if (j > i)
        quadratic_roots[a][b] = 2;
    }
  }


  // Initialize table for number of roots of depressed cubic 
  depressed_cubic_roots = new unsigned*[q];
  for (int i = 0; i < q; i++) {
    depressed_cubic_roots[i] = new unsigned[q];
    for (int j = 0; j < q; j++) {
      depressed_cubic_roots[i][j] = 0;
    }
  }
  
  // Fill the cubics table
  for (int a = 0; a < q; a++) {
    for (int b = 0; b < q; b++) {
      // (x^2 + ax + b)(x+a) = x^3 + sx + t
      unsigned s = mult[a][a] ^ b, t = mult[a][b];
      if (b == 0) { // if a is already a root of x^2 + ax + b
        depressed_cubic_roots[s][t] = quadratic_roots[a][b];
      } else {

        depressed_cubic_roots[s][t] = 1+quadratic_roots[a][b];
      }
    }
  }

  // Frobenius orbits
  unsigned* orbit_size = new unsigned[q]; // only valid if i == orbit_rep[i]
  unsigned* orbit_rep = new unsigned[q];
  for (unsigned i = 0; i < q; i++) {
    orbit_rep[i] = q; // means null
  }

  // Fill the Frobenius Orbits
  for (unsigned i = 0; i < q; i++) {
    if (orbit_rep[i] != q) continue; // already filled
    unsigned j = i;
    int size = 1;
    while (1) {
      orbit_rep[j] = i;
      j = mult[j][j];
      if (j == i) {
        orbit_size[i] = size;
        break;
      }
      size++;
    } 
  }

  return std::make_tuple(quadratic_roots, depressed_cubic_roots,
                         mult, divi, orbit_rep, orbit_size);

}
