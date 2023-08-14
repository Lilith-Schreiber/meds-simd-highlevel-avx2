#include "vec_x16.h"

int pmod_mask_count_x16(pmod_mat_mask_x16_t mask) {
  return __builtin_popcount(((int)mask) & 0xffff);
}

uint32_t extract_vec_x16(__m512i x, int pos) {  
  uint32_t buf[16] align64;
  STORE_x16(buf, x);
  return buf[pos];
}

uint32_t extract_mask_x16(__mmask16 x, int pos) { return ((uint32_t)x) & (1 << pos); }

__m512i GF_reduc_x16(const __m512i u) {
  const __m512i beta_m = SET1_x16((1 << GFq_bits) - 1);
  const __m512i one = SET1_x16(1);

  const __m512i d = SET1_x16(MEDS_p);
  const __m512i v = SET1_x16(MEDS_rep);

  const __m512i u1 = SRLI_x16(u, GFq_bits);
  const __m512i u0 = AND_x16(u, beta_m);

  const __m512i q = ADD_x16(MULLO_x16(v, u1), u);
  const __m512i q1 = AND_x16(ADD_x16(SRLI_x16(q, GFq_bits), one), beta_m);
  const __m512i q0 = AND_x16(q, beta_m);

  const __m512i q1d = MULLO_x16(q1, d);
  __m512i r = AND_x16(SUB_x16(u0, q1d), beta_m);

  const __mmask16 rgtq0 = GT_x16(r, q0);
  r = AND_x16(ADD_M_x16(r, d, rgtq0), beta_m);

  const __mmask16 rged = GE_x16(r, d);
  r = AND_x16(SUB_M_x16(r, d, rged), beta_m);

  return r;
}

__m512i GF_mod_x16(const __m512i u) {
  const __m512i beta2_m = SET1_x16((1 << (GFq_bits << 1)) - 1);
  const __m512i beta2_r = SET1_x16(GFq_b2r);
  const __m512i P = SET1_x16(MEDS_p);

  const __m512i u1 = SRLI_x16(u, (GFq_bits << 1));
  const __m512i u0 = AND_x16(u, beta2_m);

  const __m512i r1 = GF_reduc_x16(u1);
  const __m512i r0 = GF_reduc_x16(u0);

  const __m512i u2 = MULLO_x16(r1, beta2_r);
  const __m512i r2 = GF_reduc_x16(u2);

  __m512i r = ADD_x16(r0, r2);

  const __mmask16 rgep = GE_x16(r, P);
  r = SUB_M_x16(r, P, rgep);

  return r;
}

pmod_mat_x16_t GF_inv_x16(pmod_mat_x16_t x) {
  uint16_t exp = MEDS_p - 2;
  pmod_mat_x16_t t = SET1_x16(1);

  while (exp > 0) {
    if (exp & 1) t = GF_mod_x16(MULLO_x16(t, x));
    x = GF_mod_x16(MULLO_x16(x, x));
    exp >>= 1;
  }

  return t;
}
