#ifndef MATRIXMOD_W32_H
#define MATRIXMOD_W32_H

#include <immintrin.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "params.h"
#include "matrixmod.h"

#include "vec_w32.h"

pmod_mat_w32_t pmod_mat_entry_w32(uint16_t *M[], int M_r, int M_c, int r, int c);
void pmod_mat_set_entry_w32(uint16_t *M[], int M_r, int M_c, int r, int c,
                             pmod_mat_w32_t val, int num);

void pmod_mat_mul_w32(pmod_mat_w32_t *C, int C_r, int C_c, pmod_mat_w32_t *A, int A_r,
                      int A_c, pmod_mat_w32_t *B, int B_r, int B_c);
void pmod_mat_mask_mul_w32(pmod_mat_w32_t *C, int C_r, int C_c, pmod_mat_w32_t *A,
                           int A_r, int A_c, pmod_mat_w32_t *B, int B_r, int B_c,
                           int B_cm);

pmod_mat_mask_w32_t pmod_mat_syst_ct_w32(pmod_mat_w32_t *M, int M_r, int M_c);
pmod_mat_mask_w32_t pmod_mat_inv_w32(pmod_mat_w32_t *M_inv, pmod_mat_w32_t *M, int M_r,
                                 int M_c);

#endif

