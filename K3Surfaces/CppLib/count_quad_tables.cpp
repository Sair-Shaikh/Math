#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <assert.h>
#include <set>
#include "tables.h"

#ifndef EXT_COEFFS
#include "coeffs_quad.h"
#endif

unsigned NULL_Fq_elt = q;
unsigned** mult;
unsigned** divi;
unsigned** quadratic_roots;


// function prototypes
// int contribution_of_fibre_over_P2_point(unsigned, unsigned, unsigned, const std::set<ff2k_t>&, int);
int contribution_of_fibre_over_P2_point(unsigned, unsigned, unsigned, int);

int main(int argc, char **argv) {

    std::string qq = std::to_string(q); 

    mult = read_table(q, q, "mult_" + qq);
    divi = read_table(q, q, "divi_" + qq);
    quadratic_roots = read_table(q, q, "quadratic_roots_" + qq);
    unsigned* orbit_size = read_table(q, "orbit_size_" + qq);
    unsigned* orbit_rep = read_table(q, "orbit_rep_" + qq);

    // unsigned* ASSols = read_table(q, "ASSols_" + qq);
    // std::set<ff2k_t> S(ASSols, ASSols + q);
    
    int count = 0;
    // count += contribution_of_fibre_over_P2_point(0, 0, 1, S, 2);
    count += contribution_of_fibre_over_P2_point(0, 0, 1, 2);


    // The contribution over the hyperplane at infinity.
    for (unsigned y_2 = 0; y_2 < q; y_2++)
        if (orbit_rep[y_2] == y_2) {
            // count += contribution_of_fibre_over_P2_point(0, 1, y_2, S, 2) * orbit_size[y_2];
            count += contribution_of_fibre_over_P2_point(0, 1, y_2, 2) * orbit_size[y_2];
        }


    // The contribution from the A2 part. 
    for (unsigned y_1 = 0; y_1 < q; y_1++)
        if (orbit_rep[y_1] == y_1)
            for (unsigned y_2 = 0; y_2 < q; y_2++) {
                // count += contribution_of_fibre_over_P2_point(1, y_1, y_2, S, 1)*orbit_size[y_1];
                count += contribution_of_fibre_over_P2_point(1, y_1, y_2, 1)*orbit_size[y_1];
            }

                
    printf("%d\n", count);

    return 0;
}

// Count points on a quadratic
int contribution_of_fibre_over_P2_point(unsigned y_0, unsigned y_1, unsigned y_2, int type) {
    unsigned A, B, C;
    if (type == 1) {
        ABC;
    } else if (type == 2) {
        ABC2;
    } else {
        const char type_error[] = "Fiber Type needs to be 1 or 2.\n";
        std::cerr << type_error << std::endl;
        return 1;
    }

    if (A == 0) {
        if (B == 0) {
            if (C == 0) {
                return (q + 1);
            } else {
                return 1;
            }
        } else {
            return 2;
        }
    } else {
        ff2k_t L = divi[B][A];
        ff2k_t M = divi[C][A];
        return quadratic_roots[L][M];

        // ff2k_t L = ff2k_divi(B, A);
        // ff2k_t M = ff2k_divi(C, A);
        // if (L == 0) {
        //     return 1;
        // } else if (S.count(ff2k_divi(M, ff2k_square(L))) > 0) {
        //     return 2;
        // } else {
        //     return 0;
        // }
    }

    // default return (shouldn't reach here)
    return 0;

}
