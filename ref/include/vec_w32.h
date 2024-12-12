#ifndef VEC_W32_H
#define VEC_W32_H

#include <immintrin.h>
#include <stdint.h>

#include "params.h"

#define aligned __attribute__((aligned(32)))

#define pmod_mat_w32_t __m256i
#define pmod_mat_mask_w32_t __mmask8

#define ADD_w32(a, b) _mm256_add_epi32(a, b)
// #define ADD_M_w32(a, b, m) _mm_mask_add_epi32(a, m, a, b)
#define SUB_w32(a, b) _mm256_sub_epi32(a, b)
// #define SUB_M_w32(a, b, m) _mm_mask_sub_epi32(a, m, a, b)

pmod_mat_w32_t ADD_M_w32(pmod_mat_w32_t a, pmod_mat_w32_t b, __mmask8 m);
pmod_mat_w32_t SUB_M_w32(pmod_mat_w32_t a, pmod_mat_w32_t b, __mmask8 m);

#define MULLO_w32(a, b) _mm256_mullo_epi32(a, b)

#define SRLI_w32(a, b) _mm256_srli_epi32(a, b)
#define SLLI_w32(a, b) _mm256_slli_epi32(a, b)

// #define AND_w32(a, b) _mm256_and_si256(a, b)
#define OR_w32(a, b) _mm256_or_si256(a, b) 
#define XOR_w32(a, b) _mm256_xor_si256(a, b)

pmod_mat_w32_t AND_w32(pmod_mat_w32_t a, pmod_mat_w32_t b);

// #define GT_w32(a, b) _mm256_cmp_ps(a, b, 30)
// #define GE_w32(a, b) _mm256_cmp_ps(a, b, 29)
// #define LT_w32(a, b) _mm256_cmp_ps(a, b, 1)
// #define LE_w32(a, b) _mm256_cmp_ps_(a, b, 2)
// #define EQ_w32(a, b) _mm256_cm_ps(a, b, 0)
// #define NEQ_w32(a, b) _mm256_cmp_ps(a, b, 28)

__mmask8 GT_w32(__m256i a, __m256i b);
__mmask8 GE_w32(__m256i a, __m256i b);
__mmask8 LT_w32(__m256i a, __m256i b);
__mmask8 LE_w32(__m256i a, __m256i b);
__mmask8 EQ_w32(__m256i a, __m256i b);
__mmask8 NEQ_w32(__m256i a, __m256i b);

/*
__m256d _mm256_cmp_ps (__m256d a, __m256d b, const int imm8);
GT --> 30
GE --> 29
LT --> 1
LE --> 2
EQ --> 0
NEQ --> 28
*/

#define SET1_w32(a) _mm256_set1_epi32(a)
#define LOAD_w32(a) _mm256_loadu_si256(a) // if there is something wrong with this function, add `&` for the parameter input
#define STORE_w32(a, b) _mm256_storeu_si256(a, b)

int pmod_mask_count_w32(pmod_mat_mask_w32_t mask);

uint32_t extract_vec_w32(__m256i x, int pos);
uint32_t extract_mask_w32(__mmask16 x, int pos);

pmod_mat_w32_t GF_reduc_w32(const __m256i u);
pmod_mat_w32_t GF_mod_w32(const __m256i u);
pmod_mat_w32_t GF_inv_w32(pmod_mat_w32_t x);

#endif
