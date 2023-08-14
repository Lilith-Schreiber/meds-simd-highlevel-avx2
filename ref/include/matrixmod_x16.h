#ifndef MATRIXMOD_X16_H
#define MATRIXMOD_X16_H

#include <immintrin.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "params.h"
#include "matrixmod.h"

#include "vec_x16.h"

pmod_mat_x16_t pmod_mat_entry_x16(uint16_t *M[], int M_r, int M_c, int r, int c);
void pmod_mat_set_entry_x16(uint16_t *M[], int M_r, int M_c, int r, int c,
                             pmod_mat_x16_t val, int num);

void pmod_mat_mul_x16(pmod_mat_x16_t *C, int C_r, int C_c, pmod_mat_x16_t *A, int A_r,
                      int A_c, pmod_mat_x16_t *B, int B_r, int B_c);
void pmod_mat_mask_mul_x16(pmod_mat_x16_t *C, int C_r, int C_c, pmod_mat_x16_t *A,
                           int A_r, int A_c, pmod_mat_x16_t *B, int B_r, int B_c,
                           int B_cm);

pmod_mat_mask_x16_t pmod_mat_syst_ct_x16(pmod_mat_x16_t *M, int M_r, int M_c);
pmod_mat_mask_x16_t pmod_mat_inv_x16(pmod_mat_x16_t *M_inv, pmod_mat_x16_t *M, int M_r,
                                 int M_c);

#endif

