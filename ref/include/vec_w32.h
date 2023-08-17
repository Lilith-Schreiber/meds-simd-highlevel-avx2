#ifndef VEC_W32_H
#define VEC_W32_H

#include <immintrin.h>
#include <stdint.h>

#include "params.h"

#define aligned __attribute__((aligned(64)))

#define ADD_w32(a, b) _mm512_add_epi32(a, b)
#define ADD_M_w32(a, b, m) _mm512_mask_add_epi32(a, m, a, b)
#define SUB_w32(a, b) _mm512_sub_epi32(a, b)
#define SUB_M_w32(a, b, m) _mm512_mask_sub_epi32(a, m, a, b)

#define MULLO_w32(a, b) _mm512_mullo_epi32(a, b)

#define SRLI_w32(a, b) _mm512_srli_epi32(a, b)
#define SLLI_w32(a, b) _mm512_slli_epi32(a, b)

#define AND_w32(a, b) _mm512_and_epi32(a, b)
#define OR_w32(a, b) _mm512_or_epi32(a, b)
#define XOR_w32(a, b) _mm512_xor_epi32(a, b)

#define GT_w32(a, b) _mm512_cmpgt_epi32_mask(a, b)
#define GE_w32(a, b) _mm512_cmpge_epi32_mask(a, b)
#define LT_w32(a, b) _mm512_cmplt_epi32_mask(a, b)
#define LE_w32(a, b) _mm512_cmple_epi32_mask(a, b)
#define EQ_w32(a, b) _mm512_cmpeq_epi32_mask(a, b)
#define NEQ_w32(a, b) _mm512_cmpneq_epi32_mask(a, b)

#define SET1_w32(a) _mm512_set1_epi32(a)
#define LOAD_w32(a) _mm512_load_epi32(a)
#define STORE_w32(a, b) _mm512_store_epi32(a, b)

#define pmod_mat_w32_t __m512i
#define pmod_mat_mask_w32_t __mmask16

int pmod_mask_count_w32(pmod_mat_mask_w32_t mask);

uint32_t extract_vec_w32(__m512i x, int pos);
uint32_t extract_mask_w32(__mmask16 x, int pos);

__m512i GF_reduc_w32(const __m512i u);
__m512i GF_mod_w32(const __m512i u);
pmod_mat_w32_t GF_inv_w32(pmod_mat_w32_t x);

#endif
