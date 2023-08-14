#ifndef VEC_X16_H
#define VEC_X16_H

#include <immintrin.h>
#include <stdint.h>

#include "params.h"

#define align16 __attribute__((aligned(16)))
#define align32 __attribute__((aligned(32)))
#define align64 __attribute__((aligned(64)))

#define ADD_x16(a, b) _mm512_add_epi32(a, b)
#define ADD_M_x16(a, b, m) _mm512_mask_add_epi32(a, m, a, b)
#define SUB_x16(a, b) _mm512_sub_epi32(a, b)
#define SUB_M_x16(a, b, m) _mm512_mask_sub_epi32(a, m, a, b)

#define MULLO_x16(a, b) _mm512_mullo_epi32(a, b)

#define SRLI_x16(a, b) _mm512_srli_epi32(a, b)
#define SLLI_x16(a, b) _mm512_slli_epi32(a, b)

#define AND_x16(a, b) _mm512_and_epi32(a, b)
#define OR_x16(a, b) _mm512_or_epi32(a, b)
#define XOR_x16(a, b) _mm512_xor_epi32(a, b)

#define GT_x16(a, b) _mm512_cmpgt_epi32_mask(a, b)
#define GE_x16(a, b) _mm512_cmpge_epi32_mask(a, b)
#define LT_x16(a, b) _mm512_cmplt_epi32_mask(a, b)
#define LE_x16(a, b) _mm512_cmple_epi32_mask(a, b)
#define EQ_x16(a, b) _mm512_cmpeq_epi32_mask(a, b)
#define NEQ_x16(a, b) _mm512_cmpneq_epi32_mask(a, b)

#define SET1_x16(a) _mm512_set1_epi32(a)
#define LOAD_x16(a) _mm512_load_epi32(a)
#define STORE_x16(a, b) _mm512_store_epi32(a, b)

#define pmod_mat_x16_t __m512i
#define pmod_mat_mask_x16_t __mmask16

int pmod_mask_count_x16(pmod_mat_mask_x16_t mask);

uint32_t extract_vec_x16(__m512i x, int pos);
uint32_t extract_mask_x16(__mmask16 x, int pos);

__m512i GF_reduc_x16(const __m512i u);
__m512i GF_mod_x16(const __m512i u);
pmod_mat_x16_t GF_inv_x16(pmod_mat_x16_t x);

#endif
