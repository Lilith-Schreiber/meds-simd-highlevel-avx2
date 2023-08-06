#include "vec_x16.h"

int pmod_mask_count(pmod_vec_mask_t mask) {
  return __builtin_popcount(((int)mask) & 0xffff);
}

uint32_t extract_vec_x16(__m512i x, int pos) {  
  uint32_t buf[16] align64;
  STORE(buf, x);
  return buf[pos];
}

uint32_t extract_mask_x16(__mmask16 x, int pos) { return ((uint32_t)x) & (1 << pos); }

__m512i GF_reduc_x16(const __m512i u) {
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

__m512i GF_mod_x16(const __m512i u) {
  const __m512i beta2_m = SET1((1 << (GFq_bits << 1)) - 1);
  const __m512i beta2_r = SET1(GFq_b2r);
  const __m512i P = SET1(MEDS_p);

  const __m512i u1 = SRLI(u, (GFq_bits << 1));
  const __m512i u0 = AND(u, beta2_m);

  const __m512i r1 = GF_reduc_x16(u1);
  const __m512i r0 = GF_reduc_x16(u0);

  const __m512i u2 = MULLO(r1, beta2_r);
  const __m512i r2 = GF_reduc_x16(u2);

  __m512i r = ADD(r0, r2);

  const __mmask16 rgep = GE(r, P);
  r = SUB_M(r, P, rgep);

  return r;
}

pmod_mat_x16_t GF_inv_x16(pmod_mat_x16_t x) {
  uint16_t exp = MEDS_p - 2;
  pmod_mat_x16_t t = SET1(1);

  while (exp > 0) {
    if (exp & 1) t = GF_mod_x16(MULLO(t, x));
    x = GF_mod_x16(MULLO(x, x));
    exp >>= 1;
  }

  return t;
}
