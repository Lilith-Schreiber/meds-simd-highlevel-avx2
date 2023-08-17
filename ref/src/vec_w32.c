#include "vec_w32.h"

int pmod_mask_count_w32(pmod_mat_mask_w32_t mask) {
  return __builtin_popcount(((int)mask) & 0xffff);
}

uint32_t extract_vec_w32(__m512i x, int pos) {  
  uint32_t buf[16] aligned;
  STORE_w32(buf, x);
  return buf[pos];
}

uint32_t extract_mask_w32(__mmask16 x, int pos) { return ((uint32_t)x) & (1 << pos); }

__m512i GF_reduc_w32(const __m512i u) {
  const __m512i beta_m = SET1_w32((1 << GFq_bits) - 1);
  const __m512i one = SET1_w32(1);

  const __m512i d = SET1_w32(MEDS_p);
  const __m512i v = SET1_w32(MEDS_rep);

  const __m512i u1 = SRLI_w32(u, GFq_bits);
  const __m512i u0 = AND_w32(u, beta_m);

  const __m512i q = ADD_w32(MULLO_w32(v, u1), u);
  const __m512i q1 = AND_w32(ADD_w32(SRLI_w32(q, GFq_bits), one), beta_m);
  const __m512i q0 = AND_w32(q, beta_m);

  const __m512i q1d = MULLO_w32(q1, d);
  __m512i r = AND_w32(SUB_w32(u0, q1d), beta_m);

  const __mmask16 rgtq0 = GT_w32(r, q0);
  r = AND_w32(ADD_M_w32(r, d, rgtq0), beta_m);

  const __mmask16 rged = GE_w32(r, d);
  r = AND_w32(SUB_M_w32(r, d, rged), beta_m);

  return r;
}

__m512i GF_mod_w32(const __m512i u) {
  const __m512i beta2_m = SET1_w32((1 << (GFq_bits << 1)) - 1);
  const __m512i beta2_r = SET1_w32(GFq_b2r);
  const __m512i P = SET1_w32(MEDS_p);

  const __m512i u1 = SRLI_w32(u, (GFq_bits << 1));
  const __m512i u0 = AND_w32(u, beta2_m);

  const __m512i r1 = GF_reduc_w32(u1);
  const __m512i r0 = GF_reduc_w32(u0);

  const __m512i u2 = MULLO_w32(r1, beta2_r);
  const __m512i r2 = GF_reduc_w32(u2);

  __m512i r = ADD_w32(r0, r2);

  const __mmask16 rgep = GE_w32(r, P);
  r = SUB_M_w32(r, P, rgep);

  return r;
}

pmod_mat_w32_t GF_inv_w32(pmod_mat_w32_t x) {
  uint16_t exp = MEDS_p - 2;
  pmod_mat_w32_t t = SET1_w32(1);

  while (exp > 0) {
    if (exp & 1) t = GF_mod_w32(MULLO_w32(t, x));
    x = GF_mod_w32(MULLO_w32(x, x));
    exp >>= 1;
  }

  return t;
}
