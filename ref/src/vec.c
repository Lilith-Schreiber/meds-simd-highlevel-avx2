#include "vec.h"

#include <immintrin.h>

#include "matrixmod.h"
#include "params.h"

pmod_vec_t pmod_mat_entry_vec(uint16_t *M[], int M_r, int M_c, int r, int c) {
  uint32_t M_buf[16] align64 = {0};
  for (int i = 0; i < 16; i++) {
    M_buf[i] = pmod_mat_entry(M[i], M_r, M_c, r, c);
  }

  return LOAD(M_buf);
}

void pmod_mat_set_entry_vec(uint16_t *M[], int M_r, int M_c, int r, int c,
                            pmod_vec_t val, int num) {
  uint32_t buf[16] align64 = {0};
  STORE((int *)buf, val);

  for (int i = 0; i < num; i++) {
    const pmod_mat_t entry = buf[i];
    pmod_mat_set_entry(M[i], M_r, M_c, r, c, entry);
  }
}

void print_256_vec(__m256i a, __m256i mask) {
  uint32_t a_buf[16 >> 1];
  _mm256_maskstore_epi32((int *)a_buf, mask, a);
  for (int i = 0; i < 8; i++) {
    printf("%d ", a_buf[i] >> 16);
    printf("%d ", a_buf[i] & 0xFFFF);
  }
  printf("\n");
}

void print_512_vec(__m512i a) {
  uint32_t a_buf[16];
  STORE((int *)a_buf, a);
  for (int i = 0; i < 16; i++) {
    printf("%6d", a_buf[i]);
  }
  printf("\n");
}

int pmod_mask_count(pmod_vec_mask_t mask) {
  return __builtin_popcount(((int)mask) & 0xffff);
}

uint32_t extract_vec(__m512i x, int pos) {  
  uint32_t buf[16] align64;
  STORE(buf, x);
  return buf[pos];
}

uint32_t extract_mask(__mmask16 x, int pos) { return ((uint32_t)x) & (1 << pos); }

__m512i GF_reduc_vec(const __m512i u) {
  const __m512i beta_m = SET1((1 << GFq_bits) - 1);
  const __m512i one = SET1(1);

  const __m512i d = SET1(MEDS_p);
  const __m512i v = SET1(MEDS_rep);

  const __m512i u1 = SRLI(u, GFq_bits);
  const __m512i u0 = AND(u, beta_m);

  const __m512i q = ADD(MULLO(v, u1), u);
  const __m512i q1 = AND(ADD(SRLI(q, GFq_bits), one), beta_m);
  const __m512i q0 = AND(q, beta_m);

  const __m512i q1d = MULLO(q1, d);
  __m512i r = AND(SUB(u0, q1d), beta_m);

  const __mmask16 rgtq0 = GT(r, q0);
  r = AND(ADD_M(r, d, rgtq0), beta_m);

  const __mmask16 rged = GE(r, d);
  r = AND(SUB_M(r, d, rged), beta_m);

  return r;
}

// __m512i GF_mod_vec(const __m512i u) {
//   // return u;
//   const __m512i beta2_m = SET1((1 << (GFq_bits << 1)) - 1);
//   const __m512i beta2_r = SET1(GFq_b2r);
//   const __m512i P = SET1(MEDS_p);
//
//   const __m512i u1 = SRLI(u, (GFq_bits << 1));
//   const __m512i u0 = AND(u, beta2_m);
//
//   const __m512i r1 = GF_reduc_vec(u1);
//   const __m512i r0 = GF_reduc_vec(u0);
//
//   const __m512i u2 = MULLO(r1, beta2_r);
//   const __m512i r2 = GF_reduc_vec(u2);
//
//   __m512i r = ADD(r0, r2);
//
//   const __mmask16 rgep = GE(r, P);
//   r = SUB_M(r, P, rgep);
//
//   return r;
// }

__m512i GF_mod_vec(const __m512i u) {
  // return u;

  const __m512i beta2_m = SET1((1 << (GFq_bits << 1)) - 1);
  const __m512i beta2_r = SET1(GFq_b2r);
  const __m512i P = SET1(MEDS_p);

  const __m512i u1 = SRLI(u, (GFq_bits << 1));
  const __m512i u0 = AND(u, beta2_m);

  const __m512i r1 = GF_reduc_vec(u1);
  const __m512i r0 = GF_reduc_vec(u0);

  const __m512i u2 = MULLO(r1, beta2_r);
  const __m512i r2 = GF_reduc_vec(u2);

  __m512i r = ADD(r0, r2);

  const __mmask16 rgep = GE(r, P);
  r = SUB_M(r, P, rgep);

  return r;
}

pmod_vec_t GF_inv_vec(pmod_vec_t x) {
  uint16_t exp = MEDS_p - 2;
  pmod_vec_t t = SET1(1);

  while (exp > 0) {
    if (exp & 1) t = GF_mod_vec(MULLO(t, x));
    x = GF_mod_vec(MULLO(x, x));
    exp >>= 1;
  }

  return t;
}
