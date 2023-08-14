#ifndef UTILS_X16_H
#define UTILS_X16_H

#include <stdint.h>
#include <string.h>

#include "params.h"
#include "vec_x16.h"
#include "matrixmod_x16.h"

pmod_mat_mask_x16_t solve_x16(pmod_mat_x16_t *A, pmod_mat_x16_t *B_inv, pmod_mat_x16_t *G0prime, pmod_mat_x16_t Amm);

void pi_x16(pmod_mat_x16_t *Go, pmod_mat_x16_t *A, pmod_mat_x16_t *B, pmod_mat_x16_t *G0);

#endif

