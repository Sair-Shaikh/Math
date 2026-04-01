#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <tuple>
#include "tableio.h"
#include "Fq.h"


std::tuple<unsigned*, unsigned*, unsigned*> generate_bigFq_tables(unsigned, unsigned);

int main() {
    unsigned *art_sch_sols;
    unsigned *orbit_rep;
    unsigned *orbit_size;

    std::tie(orbit_rep, orbit_size, art_sch_sols) = generate_bigFq_tables(q, p);

    // Write data to file.
    std::string qq = std::to_string(q);
    write_table(art_sch_sols, q, "art_sch_sols_" + qq);
    write_table(orbit_rep, q, "orbit_rep_" + qq);
    write_table(orbit_size, q, "orbit_size_" + qq);

    return 0;

  }


std::tuple<unsigned*, unsigned*, unsigned*> generate_bigFq_tables(unsigned q, unsigned p) {

  // Solutions to Artin-Schreier Quadratics
  unsigned *art_sch_sols = new unsigned[q];
  for (unsigned i = 0; i < q; i++) {
      art_sch_sols[i] = ff2k_square(i) ^ i;  
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
      j = ff2k_mult(j, j);
      if (j == i) {
        orbit_size[i] = size;
        break;
      }
      size++;
    } 
  }

  return std::make_tuple(orbit_rep, orbit_size, art_sch_sols);

}