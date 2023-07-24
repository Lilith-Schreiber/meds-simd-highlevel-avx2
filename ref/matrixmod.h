#ifndef MATRIXMOD_H
#define MATRIXMOD_H

#include <stdint.h>
#include <stdio.h>

#include "params.h"
#include "vec.h"

#define pmod_mat_entry(M, M_r, M_c, r, c) M[(M_c) * (r) + (c)]

#define pmod_mat_set_entry(M, M_r, M_c, r, c, v) (M[(M_c) * (r) + (c)] = v)

#define pmod_mat_t GFq_t

#define align16 __attribute__((aligned(16)))
#define align32 __attribute__((aligned(32)))
#define align64 __attribute__((aligned(64)))

void pmod_mat_print(pmod_mat_t *M, int M_r, int M_c);
void pmod_mat_fprint(FILE *stream, pmod_mat_t *M, int M_r, int M_c);


void pmod_mat_mul(pmod_mat_t *C, int C_r, int C_c, pmod_mat_t *A, int A_r,
                  int A_c, pmod_mat_t *B, int B_r, int B_c);
void pmod_mat_mul_vec(pmod_vec_t *C, int C_r, int C_c, pmod_vec_t *A, int A_r,
                  int A_c, pmod_vec_t *B, int B_r, int B_c);

pmod_vec_mask_t pmod_mat_syst_ct_vec(pmod_vec_t *M, int M_r, int M_c);
void pmod_mat_back_substitution_ct_vec(pmod_vec_t *M, int M_r, int M_c);

int pmod_mat_syst_ct(pmod_mat_t *M, int M_r, int M_c);
int pmod_mat_row_echelon_ct(pmod_mat_t *M, int M_r, int M_c);
int pmod_mat_back_substitution_ct(pmod_mat_t *M, int M_r, int M_c);

int pmod_mat_inv(pmod_mat_t *B, pmod_mat_t *A, int A_r, int A_c);

GFq_t GF_inv(GFq_t val);

#endif
