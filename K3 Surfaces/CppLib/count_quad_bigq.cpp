#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <assert.h>
#include <set>
#include "tableio.h"
#include "Fq.h"

#ifndef EXT_COEFFS
#include "batched_coeffs_quad_bigq.h"
#endif

#ifndef BATCH_SIZE
#define BATCH_SIZE 10
#endif

// function prototypes
int* contribution_of_fibre_over_P2_point(unsigned, unsigned, unsigned, unsigned, int);

int main(int argc, char **argv) {

    std::string qq = std::to_string(q); 
    
    uint64_t* count = new uint64_t[BATCH_SIZE];
    int* contributions = new int[BATCH_SIZE];

    unsigned* orbit_size = read_table(q, "orbit_size_" + qq);
    unsigned* orbit_rep = read_table(q, "orbit_rep_" + qq);

    delete[] contributions;
    contributions = contribution_of_fibre_over_P2_point(q, 0, 0, 1, 2);
    for (int i = 0; i < BATCH_SIZE; i++) {
        count[i] = contributions[i];
    }

    // The contribution over the hyperplane at infinity.
    for (unsigned y_2 = 0; y_2 < q; y_2++) {
        if (orbit_rep[y_2] == y_2) {
            delete[] contributions;
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
                delete[] contributions;
                contributions = contribution_of_fibre_over_P2_point(q, 1, y_1, y_2, 1);
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

// Count points on a quadratic
int* contribution_of_fibre_over_P2_point(unsigned q, unsigned y_0, unsigned y_1, unsigned y_2, int type) {
    unsigned* As = new unsigned[BATCH_SIZE], *Bs = new unsigned[BATCH_SIZE], *Cs = new unsigned[BATCH_SIZE];
    int* results = new int[BATCH_SIZE];
    if (type == 1) {
        BATCHED_ABC;
    } else if (type == 2) {
        BATCHED_ABC2;
    } else {
        const char type_error[] = "Fiber Type needs to be 1 or 2.\n";
        std::cerr << type_error << std::endl;
        return new int[BATCH_SIZE];
    }

    for (int i = 0; i < BATCH_SIZE; i++) {
        unsigned A = As[i], B = Bs[i], C = Cs[i];
        results[i] = count_quadratic_roots(A, B, C);
    }

    delete[] As;
    delete[] Bs;
    delete[] Cs;


    return results;
}