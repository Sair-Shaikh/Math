#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <assert.h>
#include <set>
#include "tableio.h"
#include "monomials_bigq.h"
#include "Fq.h"

#ifndef EXT_COEFFS
#include "batched_coeffs_quad_fast_test.h"
#endif

#ifndef BATCH_SIZE
#define BATCH_SIZE 10
#endif

// function prototypes
int* contribution_of_fibre_over_P2_point(unsigned, unsigned, unsigned, unsigned, std::array<std::array<unsigned,30>,BATCH_SIZE>& As, std::array<std::array<unsigned,30>,BATCH_SIZE>& Bs, std::array<std::array<unsigned,30>,BATCH_SIZE>& Cs);

int main(int argc, char **argv) {

    std::string qq = std::to_string(q); 
    
    uint64_t* count = new uint64_t[BATCH_SIZE];
    int* contributions = new int[BATCH_SIZE];

    unsigned* orbit_size = read_table(q, "orbit_size_" + qq);
    unsigned* orbit_rep = read_table(q, "orbit_rep_" + qq);

    std::array<std::array<unsigned,30>, BATCH_SIZE> As, Bs, Cs;
    BATCHED_ABC2;
    delete[] contributions;
    contributions = contribution_of_fibre_over_P2_point(q, 0, 0, 1, As, Bs, Cs);
    for (int i = 0; i < BATCH_SIZE; i++) {
        count[i] = contributions[i];
    }


    // The contribution over the hyperplane at infinity.
    for (unsigned y_2 = 0; y_2 < q; y_2++) {
        if (orbit_rep[y_2] == y_2) {
            delete[] contributions;
            contributions = contribution_of_fibre_over_P2_point(q, 0, 1, y_2, As, Bs, Cs);
            for (int i = 0; i < BATCH_SIZE; i++) {
                count[i] += contributions[i]*orbit_size[y_2];
            }
            
        }
    }


    BATCHED_ABC;
    // The contribution from the A2 part. 
    for (unsigned y_1 = 0; y_1 < q; y_1++) {
        if (orbit_rep[y_1] == y_1) {
            for (unsigned y_2 = 0; y_2 < q; y_2++) {
                delete[] contributions;
                contributions = contribution_of_fibre_over_P2_point(q, 1, y_1, y_2, As, Bs, Cs);
                for (int i = 0; i < BATCH_SIZE; i++) {
                    count[i] += contributions[i]*orbit_size[y_1];
                }
            }
        }
    }

    for (int i = 0; i < BATCH_SIZE; i++) {
        printf("%llu\n", count[i]);
    }       
    
    delete[] contributions;
    delete[] count;
    delete[] orbit_size;
    delete[] orbit_rep;

    return 0;
}

unsigned to_coeff(std::array<unsigned, 30> Amons, unsigned *mons) {

    unsigned A = 0;

    for (int idx : Amons) {
        if (idx == 0) {
            return A;          // stop immediately
        }

        A ^= mons[idx - 1];     // mons index is 1-based
    }

    return A;
}


int* contribution_of_fibre_over_P2_point(unsigned q, unsigned y_0, unsigned y_1, unsigned y_2, std::array<std::array<unsigned,30>,BATCH_SIZE>& As, std::array<std::array<unsigned,30>,BATCH_SIZE>& Bs, std::array<std::array<unsigned,30>,BATCH_SIZE>& Cs) {
    unsigned *mons = new unsigned[84];
    MONS;
    int *results = new int[BATCH_SIZE];

    for (int i = 0; i < BATCH_SIZE; i++) {
        std::array<unsigned, 30>  Amons = As[i], Bmons = Bs[i], Cmons = Cs[i];
        unsigned A = to_coeff(Amons, mons), B = to_coeff(Bmons, mons), C = to_coeff(Cmons, mons);
        results[i] = count_quadratic_roots(A, B, C);
    }

    delete[] mons;

    return results;
}