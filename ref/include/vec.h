#ifndef VEC_H
#define VEC_H

#include <immintrin.h>
#include <stdint.h>

#include "params.h"

#define ADD(a, b) _mm512_add_epi32(a, b)
#define ADD_M(a, b, m) _mm512_mask_add_epi32(a, m, a, b)
#define SUB(a, b) _mm512_sub_epi32(a, b)
#define SUB_M(a, b, m) _mm512_mask_sub_epi32(a, m, a, b)

#define MULLO(a, b) _mm512_mullo_epi32(a, b)
#define SRLI(a, b) _mm512_srli_epi32(a, b)
#define AND(a, b) _mm512_and_epi32(a, b)

#define GT(a, b) _mm512_cmpgt_epi32_mask(a, b)
#define GE(a, b) _mm512_cmpge_epi32_mask(a, b)
#define EQ(a, b) _mm512_cmpeq_epi32_mask(a, b)
#define NEQ(a, b) _mm512_cmpneq_epi32_mask(a, b)
#define LT(a, b) _mm512_cmplt_epi32_mask(a, b)

#define SET1(a) _mm512_set1_epi32(a)
#define LOAD(a) _mm512_load_epi32(a)
#define STORE(a, b) _mm512_store_epi32(a, b)

#define pmod_vec_t __m512i
#define pmod_vec_mask_t __mmask16

pmod_vec_t pmod_mat_entry_vec(uint16_t *M[], int M_r, int M_c, int r, int c);
void pmod_mat_set_entry_vec(uint16_t *M[], int M_r, int M_c, int r, int c,
                             pmod_vec_t val);

int pmod_mask_count(pmod_vec_mask_t mask);

// uint32_t extract_vec(__m512i x);
uint32_t extract_vec(__m512i x, int pos);
uint32_t extract_mask(__mmask16 x, int pos);

void print_256_vec(__m256i a, __m256i mask);
void print_512_vec(__m512i a);

__m512i GF_reduc_vec(const __m512i u);
__m512i GF_mod_vec(const __m512i u);
pmod_vec_t GF_inv_vec(pmod_vec_t x);

#endif
