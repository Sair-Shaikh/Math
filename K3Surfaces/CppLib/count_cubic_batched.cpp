#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <assert.h>
#include <set>
#include "tableio.h"

#ifndef EXT_COEFFS
#include "batched_coeffs_1.h"
#endif

#ifndef BATCH_SIZE
#define BATCH_SIZE 100
#endif

unsigned** mult;
unsigned** divi;
unsigned** quadratic_roots;
unsigned** depressed_cubic_roots;

// function prototypes
int* contribution_of_fibre_over_P2_point(unsigned, unsigned, unsigned, unsigned, int);

int main(int argc, char **argv) {

    for (int N = 1; N < 12; N ++ ) {
        
        const unsigned q = 1 << N;

        std::string qq = std::to_string(q); 

        mult = read_table(q, q, "mult_" + qq);
        divi = read_table(q, q, "divi_" + qq);
        quadratic_roots = read_table(q, q, "quadratic_roots_" + qq);
        depressed_cubic_roots = read_table(q, q, "depressed_cubic_roots_" + qq);
        unsigned* orbit_size = read_table(q, "orbit_size_" + qq);
        unsigned* orbit_rep = read_table(q, "orbit_rep_" + qq);
        
        int* count = new int[BATCH_SIZE];
        int* contributions = new int[BATCH_SIZE];

        contributions = contribution_of_fibre_over_P2_point(q, 0, 0, 1, 2);
        for (int i = 0; i < BATCH_SIZE; i++) {
            count[i] = contributions[i];
        }

        // The contribution over the hyperplane at infinity.
        for (unsigned y_2 = 0; y_2 < q; y_2++) {
            if (orbit_rep[y_2] == y_2) {
                contributions = contribution_of_fibre_over_P2_point(q, 0, 1, y_2, 2);
                for (int i = 0; i < BATCH_SIZE; i++) {
                    count[i] += contributions[i]*orbit_size[y_2];
                }
                
            }
        }

        // The contribution from the A2 part. 
        for (unsigned y_1 = 0; y_1 < q; y_1++) {
            if (orbit_rep[y_1] == y_1) {
                for (unsigned y_2 = 0; y_2 < q; y_2++) {
                    contributions = contribution_of_fibre_over_P2_point(q, 1, y_1, y_2, 1);
                    for (int i = 0; i < BATCH_SIZE; i++) {
                        count[i] += contributions[i]*orbit_size[y_1];
                    }
                }
            }
        }

        for (int i = 0; i < BATCH_SIZE; i++) {
            printf("%d\n", count[i]);
        }                    
    }

        return 0;

}

// Count points on a quadratic
int* contribution_of_fibre_over_P2_point(unsigned q, unsigned y_0, unsigned y_1, unsigned y_2, int type) {
    unsigned* As = new unsigned[BATCH_SIZE], *Bs = new unsigned[BATCH_SIZE], *Cs = new unsigned[BATCH_SIZE], *Ds = new unsigned[BATCH_SIZE];
    int* results = new int[BATCH_SIZE];
    if (type == 1) {
        BATCHED_ABCD;
    } else if (type == 2) {
        BATCHED_ABCD2;
    } else {
        const char type_error[] = "Fiber Type needs to be 1 or 2.\n";
        std::cerr << type_error << std::endl;
        return new int[BATCH_SIZE];
    }

    for (int i = 0; i < BATCH_SIZE; i++) {
        unsigned A = As[i], B = Bs[i], C = Cs[i], D = Ds[i];
        if (A == 0) { // In this case, we have an extra solution at infinity
            // Equation is a quadratic
            if (B == 0) {
                if (C == 0) {
                    if (D == 0) {
                        results[i] = 1+q;
                    } else {
                        results[i] = 1;
                    }
                } else {
                    results[i] = 1+1;
                }
            } else {
                ff2k_t L = divi[C][B];
                ff2k_t M = divi[D][B];
                
                results[i] = 1+quadratic_roots[L][M];
            }
        } else {
            // Check if already depressed:
            if (A == 1 & B == 0) {
                results[i] = depressed_cubic_roots[C][D];
            } else {
                // Make depressed
                unsigned B2 = divi[B][A];
                ff2k_t P = mult[B2][B2]^divi[C][A];
                ff2k_t Q = mult[B2][divi[C][A]]^divi[D][A];

                results[i] = depressed_cubic_roots[P][Q];
            }
        }
    }
    return results;

}
