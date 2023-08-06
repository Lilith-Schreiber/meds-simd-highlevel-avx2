#ifndef VEC_X8_H
#define VEC_X8_H

#include <immintrin.h>
#include <stdint.h>

#define align64 __attribute__((aligned(64)))

#define ADD_x8(a, b) _mm512_add_epi64(a, b)
#define ADD_M_x8(a, b, m) _mm512_mask_add_epi64(a, m, a, b)
#define SUB_x8(a, b) _mm512_sub_epi64(a, b)
#define SUB_M_x8(a, b, m) _mm512_mask_sub_epi64(a, m, a, b)

#define MULLO_x8(a, b) _mm512_mullo_epi64(a, b)

#define SRLI_x8(a, b) _mm512_srli_epi64(a, b)
#define SLLI_x8(a, b) _mm512_slli_epi64(a, b)

#define AND_x8(a, b) _mm512_and_epi64(a, b)
#define OR_x8(a, b) _mm512_or_epi64(a, b)
#define XOR_x8(a, b) _mm512_xor_epi64(a, b)

#define GT_x8(a, b) _mm512_cmpgt_epi64_mask(a, b)
#define GE_x8(a, b) _mm512_cmpge_epi64_mask(a, b)
#define LT_x8(a, b) _mm512_cmplt_epi64_mask(a, b)
#define LE_x8(a, b) _mm512_cmple_epi64_mask(a, b)
#define EQ_x8(a, b) _mm512_cmpeq_epi64_mask(a, b)
#define NEQ_x8(a, b) _mm512_cmpneq_epi64_mask(a, b)

#define SET1_CT_x8(a) {a, a, a, a, a, a, a, a}
#define SET1_x8(a) _mm512_set1_epi64(a)

#define LOAD_x8(a) _mm512_load_epi64(a)
#define STORE_x8(a, b) _mm512_store_epi64(a, b)

#define pmod_vec_x8_t __m512i
#define pmod_vec_mask_x8_t __mmask8

uint64_t extract_vec_x8(pmod_vec_x8_t x, int pos);
uint32_t extract_mask_x8(pmod_vec_mask_x8_t x, int pos);

#endif
