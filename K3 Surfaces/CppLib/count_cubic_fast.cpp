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
#include "batched_coeffs_cubic_fast.h"
#endif

#ifndef BATCH_SIZE
#define BATCH_SIZE 10
#endif

// function prototypes
int* contribution_of_fibre_over_P2_point(unsigned, unsigned, unsigned, unsigned, unsigned __int128*, unsigned __int128*, unsigned __int128*, unsigned __int128*);

int main(int argc, char **argv) {

    std::string qq = std::to_string(q); 
    
    uint64_t* count = new uint64_t[BATCH_SIZE];
    int* contributions = new int[BATCH_SIZE];

    unsigned* orbit_size = read_table(q, "orbit_size_" + qq);
    unsigned* orbit_rep = read_table(q, "orbit_rep_" + qq);

    unsigned __int128 *As = new unsigned __int128[BATCH_SIZE], *Bs = new unsigned __int128[BATCH_SIZE], *Cs = new unsigned __int128[BATCH_SIZE], *Ds = new unsigned __int128[BATCH_SIZE];
    BATCHED_ABC2;
    delete[] contributions;
    contributions = contribution_of_fibre_over_P2_point(q, 0, 0, 1, As, Bs, Cs, Ds);
    for (int i = 0; i < BATCH_SIZE; i++) {
        count[i] = contributions[i];
        // printf("ct[%o]=%llu\n", i, count[i]);
    }


    // The contribution over the hyperplane at infinity.
    for (unsigned y_2 = 0; y_2 < q; y_2++) {
        if (orbit_rep[y_2] == y_2) {
            delete[] contributions;
            contributions = contribution_of_fibre_over_P2_point(q, 0, 1, y_2, As, Bs, Cs, Ds);
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
                contributions = contribution_of_fibre_over_P2_point(q, 1, y_1, y_2, As, Bs, Cs, Ds);
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

unsigned to_coeff(unsigned __int128 num, unsigned *mons) {

    unsigned A = 0;                   
    int i = 0;
    do {
        num = num >> 1;
        if (num & 1) {
            A ^= mons[i];
            // printf("i=%d, A=%d\n", i, A);

        }
        i++;
    } while (num != 0);

    return A;
}


int* contribution_of_fibre_over_P2_point(unsigned q, unsigned y_0, unsigned y_1, unsigned y_2, unsigned __int128 *As, unsigned __int128 *Bs, unsigned __int128 *Cs, unsigned __int128 *Ds) {
    unsigned *mons = new unsigned[84];
    MONS;
    int *results = new int[BATCH_SIZE];

    // for (int i = 0; i < 84; i ++) {
    //     printf("i=%d, mon=%d\n", i, mons[i]);
    // }


    for (int i = 0; i < BATCH_SIZE; i++) {
        unsigned __int128 Aint = As[i], Bint = Bs[i], Cint = Cs[i], Dint = Ds[i];
        unsigned A = to_coeff(Aint, mons), B = to_coeff(Bint, mons), C = to_coeff(Cint, mons), D = to_coeff(Dint, mons);
        // printf("A=%d\n", A);
        // printf("B=%d\n", B);
        // printf("C=%d\n", C);
        // printf("D=%d\n", D);


        results[i] = count_cubic_roots(A, B, C, D);
    }

    delete[] mons;

    return results;
}