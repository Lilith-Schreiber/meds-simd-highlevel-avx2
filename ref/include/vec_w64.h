#ifndef VEC_W64_H
#define VEC_W64_H

#include <immintrin.h>
#include <stdint.h>

#include "params.h"

#define aligned __attribute__((aligned(64)))

#define pmod_mat_w64_t __m256i
#define pmod_mat_mask_w64_t __mmask8

#define ADD_w64(a, b) _mm256_add_epi64(a, b)
// #define ADD_M_w64(a, b, m) _mm_mask_add_epi64(a, m, a, b) // need to modify
#define SUB_w64(a, b) _mm256_sub_epi64(a, b)
// #define SUB_M_w64(a, b, m) _mm_mask_sub_epi64(a, m, a, b) // need to modify

pmod_mat_w64_t ADD_M_w64(pmod_mat_w64_t a, pmod_mat_w64_t b, pmod_mat_mask_w64_t m);
pmod_mat_w64_t SUB_M_w64(pmod_mat_w64_t a, pmod_mat_w64_t b, pmod_mat_mask_w64_t m);

#define MULLO_w64(a, b) _mm256_mullo_epi32(a, b)

#define SRLI_w64(a, b) _mm256_srli_epi64(a, b)
#define SLLI_w64(a, b) _mm256_slli_epi64(a, b)

#define AND_w64(a, b) _mm256_and_si256(a, b)
#define OR_w64(a, b) _mm256_or_si256(a, b)
#define XOR_w64(a, b) _mm256_xor_si256(a, b)

// #define GT_w64(a, b) _mm256_cmpgt_epi64_mask(a, b)
// #define GE_w64(a, b) _mm256_cmpge_epi64_mask(a, b)
// #define LT_w64(a, b) _mm256_cmplt_epi64_mask(a, b)
// #define LE_w64(a, b) _mm256_cmple_epi64_mask(a, b)
// #define EQ_w64(a, b) _mm256_cmpeq_epi64_mask(a, b)
// #define NEQ_w64(a, b) _mm256_cmpneq_epi64_mask(a, b)

__mmask8 GT_w64(__m256i a, __m256i b);
__mmask8 GE_w64(__m256i a, __m256i b);
__mmask8 LT_w64(__m256i a, __m256i b);
__mmask8 LE_w64(__m256i a, __m256i b);
__mmask8 EQ_w64(__m256i a, __m256i b);
__mmask8 NEQ_w64(__m256i a, __m256i b);

#define SET1_CT_w64(a) {a, a, a, a, a, a, a, a}
#define SET1_w64(a) _mm256_set1_epi64x(a)

#define LOAD_w64(a) _mm256_loadu_si256(a)
#define STORE_w64(a, b) _mm256_store_si256(a, b)

uint64_t extract_vec_w64(pmod_mat_w64_t x, int pos);
uint32_t extract_mask_w64(pmod_mat_mask_w64_t x, int pos);

#endif
