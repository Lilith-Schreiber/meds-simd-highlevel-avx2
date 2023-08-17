#ifndef UTILS_W32_H
#define UTILS_W32_H

#include <stdint.h>
#include <string.h>

#include "params.h"
#include "vec_w32.h"
#include "matrixmod_vec.h"

pmod_mat_mask_w32_t solve_w32(pmod_mat_w32_t *A, pmod_mat_w32_t *B_inv, pmod_mat_w32_t *G0prime, pmod_mat_w32_t Amm);

void pi_w32(pmod_mat_w32_t *Go, pmod_mat_w32_t *A, pmod_mat_w32_t *B, pmod_mat_w32_t *G0);

#endif

