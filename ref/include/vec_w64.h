#ifndef VEC_W64_H
#define VEC_W64_H

#include <immintrin.h>
#include <stdint.h>

#include "params.h"

#define aligned __attribute__((aligned(64)))

#define ADD_w64(a, b) _mm512_add_epi64(a, b)
#define ADD_M_w64(a, b, m) _mm512_mask_add_epi64(a, m, a, b)
#define SUB_w64(a, b) _mm512_sub_epi64(a, b)
#define SUB_M_w64(a, b, m) _mm512_mask_sub_epi64(a, m, a, b)

#define MULLO_w64(a, b) _mm512_mullo_epi64(a, b)

#define SRLI_w64(a, b) _mm512_srli_epi64(a, b)
#define SLLI_w64(a, b) _mm512_slli_epi64(a, b)

#define AND_w64(a, b) _mm512_and_epi64(a, b)
#define OR_w64(a, b) _mm512_or_epi64(a, b)
#define XOR_w64(a, b) _mm512_xor_epi64(a, b)

#define GT_w64(a, b) _mm512_cmpgt_epi64_mask(a, b)
#define GE_w64(a, b) _mm512_cmpge_epi64_mask(a, b)
#define LT_w64(a, b) _mm512_cmplt_epi64_mask(a, b)
#define LE_w64(a, b) _mm512_cmple_epi64_mask(a, b)
#define EQ_w64(a, b) _mm512_cmpeq_epi64_mask(a, b)
#define NEQ_w64(a, b) _mm512_cmpneq_epi64_mask(a, b)

#define SET1_CT_w64(a) {a, a, a, a, a, a, a, a}
#define SET1_w64(a) _mm512_set1_epi64(a)

#define LOAD_w64(a) _mm512_load_epi64(a)
#define STORE_w64(a, b) _mm512_store_epi64(a, b)

#define pmod_mat_w64_t __m512i
#define pmod_mat_mask_w64_t __mmask8

uint64_t extract_vec_w64(pmod_mat_w64_t x, int pos);
uint32_t extract_mask_w64(pmod_mat_mask_w64_t x, int pos);

#endif
