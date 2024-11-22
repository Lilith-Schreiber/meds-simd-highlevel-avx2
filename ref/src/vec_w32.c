#include "vec_w32.h"

pmod_mat_w32_t ADD_M_w32(pmod_mat_w32_t a, pmod_mat_w32_t b, pmod_mat_mask_w32_t m) {
        pmod_mat_w32_t r;

        uint32_t* r2 = (uint32_t*)&r;
        uint32_t* a2 = (uint32_t*)&a;
        uint32_t* b2 = (uint32_t*)&b;

        for (int i = 0; i < 8; i++) {
                ((m >> i) & 1)? r2[i] = a2[i] + b2[i] : a2[i];
        }

        return r;
}

pmod_mat_w32_t SUB_M_w32(pmod_mat_w32_t a, pmod_mat_w32_t b, pmod_mat_mask_w32_t m) {
        pmod_mat_w32_t r;

        uint32_t* r2 = (uint32_t*)&r;
        uint32_t* a2 = (uint32_t*)&a;
        uint32_t* b2 = (uint32_t*)&b;

        for (int i = 0; i < 8; i++) {
                ((m >> i) & 1)? r2[i] = a2[i] - b2[i] : a2[i];
        }

        return r;
}

__mmask8 GT_w32(__m256i a, __m256i b) {
        __mmask8 r = 0x0;

        __m256 tmpA = _mm256_cvtepi32_ps(a);
        __m256 tmpB = _mm256_cvtepi32_ps(b);

        __m256 ret = _mm256_cmp_ps(tmpA, tmpB, 30);
        float* retf = (float*)&ret;

        for (int i = 0; i < 8; i++) {
                r = r << 1;

                if (ret[i] == 0)
                        r = r | 0;
                else
                        r = r | 1;
        }

        return r;
}

__mmask8 GE_w32(__m256i a, __m256i b) {
        __mmask8 r = 0x0;

        __m256 tmpA = _mm256_cvtepi32_ps(a);
        __m256 tmpB = _mm256_cvtepi32_ps(b);

        __m256 ret = _mm256_cmp_ps(tmpA, tmpB, 29);
        float* retf = (float*)&ret;

        for (int i = 0; i < 8; i++) {
                r = r << 1;

                if (ret[i] == 0)
                        r = r | 0;
                else
                        r = r | 1;
        }

        return r;
}

__mmask8 LT_w32(__m256i a, __m256i b) {
        __mmask8 r = 0x0;

        __m256 tmpA = _mm256_cvtepi32_ps(a);
        __m256 tmpB = _mm256_cvtepi32_ps(b);

        __m256 ret = _mm256_cmp_ps(tmpA, tmpB, 1);
        float* retf = (float*)&ret;

        for (int i = 0; i < 8; i++) {
                r = r << 1;

                if (ret[i] == 0)
                        r = r | 0;
                else
                        r = r | 1;
        }

        return r;
}

__mmask8 LE_w32(__m256i a, __m256i b) {
        __mmask8 r = 0x0;

        __m256 tmpA = _mm256_cvtepi32_ps(a);
        __m256 tmpB = _mm256_cvtepi32_ps(b);

        __m256 ret = _mm256_cmp_ps(tmpA, tmpB, 2);
        float* retf = (float*)&ret;

        for (int i = 0; i < 8; i++) {
                r = r << 1;

                if (ret[i] == 0)
                        r = r | 0;
                else
                        r = r | 1;
        }

        return r;
}

__mmask8 EQ_w32(__m256i a, __m256i b) {
        __mmask8 r = 0x0;

        __m256 tmpA = _mm256_cvtepi32_ps(a);
        __m256 tmpB = _mm256_cvtepi32_ps(b);

        __m256 ret = _mm256_cmp_ps(tmpA, tmpB, 0);
        float* retf = (float*)&ret;

        for (int i = 0; i < 8; i++) {
                r = r << 1;

                if (ret[i] == 0)
                        r = r | 0;
                else
                        r = r | 1;
        }

        return r;
}

__mmask8 NEQ_w32(__m256i a, __m256i b) {
        __mmask8 r = 0x0;

        __m256 tmpA = _mm256_cvtepi32_ps(a);
        __m256 tmpB = _mm256_cvtepi32_ps(b);

        __m256 ret = _mm256_cmp_ps(tmpA, tmpB, 28);
        float* retf = (float*)&ret;

        for (int i = 0; i < 8; i++) {
                r = r << 1;

                if (ret[i] == 0)
                        r = r | 0;
                else
                        r = r | 1;
        }

        return r;
}

int pmod_mask_count_w32(pmod_mat_mask_w32_t mask) {
  return __builtin_popcount(((int)mask) & 0xffff);
}

uint32_t extract_vec_w32(__m256i x, int pos) {
  uint32_t buf[16] aligned;
  STORE_w32(buf, x);
  return buf[pos];
}

uint32_t extract_mask_w32(__mmask8 x, int pos) { return ((uint32_t)x) & (1 << pos); }

pmod_mat_w32_t GF_reduc_w32(const __m256i u) {
  const __m256i beta_m = SET1_w32((1 << GFq_bits) - 1);
  const __m256i one = SET1_w32(1);

  const __m256i d = SET1_w32(MEDS_p);
  const __m256i v = SET1_w32(MEDS_rep);

  const __m256i u1 = SRLI_w32(u, GFq_bits);
  const __m256i u0 = AND_w32(u, beta_m);

  const __m256i q = ADD_w32(MULLO_w32(v, u1), u);
  const __m256i q1 = AND_w32(ADD_w32(SRLI_w32(q, GFq_bits), one), beta_m);
  const __m256i q0 = AND_w32(q, beta_m);

  const __m256i q1d = MULLO_w32(q1, d);
  __m256i r = AND_w32(SUB_w32(u0, q1d), beta_m);

  const __mmask8 rgtq0 = GT_w32(r, q0);
  r = AND_w32(ADD_M_w32(r, d, rgtq0), beta_m);

  const __mmask8 rged = GE_w32(r, d);
  r = AND_w32(SUB_M_w32(r, d, rged), beta_m);

  return r;
}

pmod_mat_w32_t GF_mod_w32(const __m256i u) {
  const __m256i beta2_m = SET1_w32((1 << (GFq_bits << 1)) - 1);
  const __m256i beta2_r = SET1_w32(GFq_b2r);
  const __m256i P = SET1_w32(MEDS_p);

  const __m256i u1 = SRLI_w32(u, (GFq_bits << 1));
  const __m256i u0 = AND_w32(u, beta2_m);

  const __m256i r1 = GF_reduc_w32(u1);
  const __m256i r0 = GF_reduc_w32(u0);

  const __m256i u2 = MULLO_w32(r1, beta2_r);
  const __m256i r2 = GF_reduc_w32(u2);

  pmod_mat_w32_t r = ADD_w32(r0, r2);

  const __mmask8 rgep = GE_w32(r, P);
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
