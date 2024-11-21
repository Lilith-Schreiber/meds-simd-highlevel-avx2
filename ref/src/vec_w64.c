#include "vec_w64.h"

pmod_mat_w64_t ADD_M_w64(pmod_mat_w64_t a, pmod_mat_w64_t b, pmod_mat_mask_w64_t m) {
        __m256i r;

        uint64_t* r2 = (uint64_t*)&r;
        uint64_t* a2 = (uint64_t*)&a;
        uint64_t* b2 = (uint64_t*)&b;

        for (int i = 0; i < 3; i++) {
                ((m >> i) & 1)? r2[i] = a2[i] + b2[i] : a2[i];
        }

        return r;
}

pmod_mat_w64_t SUB_M_w64(pmod_mat_w64_t a, pmod_mat_w64_t b, pmod_mat_mask_w64_t m) {
        __m256i r;

        uint64_t* r2 = (uint64_t*)&r;
        uint64_t* a2 = (uint64_t*)&a;
        uint64_t* b2 = (uint64_t*)&b;

        for (int i = 0; i < 3; i++) {
                ((m >> i) & 1)? r2[i] = a2[i] - b2[i] : a2[i];
        }

        return r;
}

__mmask8 GT_w64(__m256i a, __m256i b) {
        __mmask8 r = 0x0;

        __m256i ret = _mm256_cmpgt_epi64(a, b);
        int* retf = (int*)&ret;

        for (int i = 0; i < 3; i++) {
                if (retf[3-i] == 0) {
                        r = r << 1;
                        r = r | 0;
                        r = r << 1;
                        r = r | 0;
                }
                else {
                        r = r << 1;
                        r = r | 1;
                        r = r << 1;
                        r = r | 1;
                }
        }

        return r;
}

__mmask8 GE_w64(__m256i a, __m256i b) {
        __mmask8 r = 0x0;

        __m256i ret = _mm256_cmpgt_epi64(a, b);
        __m256i ret2 = _mm256_cmpeq_epi64(a, b);
        int* retf = (int*)&ret;
        int* retf2 = (int*)&ret2;

        for (int i = 0; i < 3; i++) {
                if (retf[3-i] == 0 && retf[3-i] == 0) {
                        r = r << 1;
                        r = r | 0;
                        r = r << 1;
                        r = r | 0;
                }
                else {
                        r = r << 1;
                        r = r | 1;
                        r = r << 1;
                        r = r | 1;
                }
        }

        return r;
}

__mmask8 LT_w64(__m256i a, __m256i b) {
        __mmask8 r = 0x0;

        __m256i ret = _mm256_cmpgt_epi64(a, b);
        __m256i ret2 = _mm256_cmpeq_epi64(a, b);
        int* retf = (int*)&ret;
        int* retf2 = (int*)&ret2;

        for (int i = 0; i < 3; i++) {
                if (retf[3-i] == 0 && retf[3-i] == 0) {
                        r = r << 1;
                        r = r | 1;
                        r = r << 1;
                        r = r | 1;
                }
                else {
                        r = r << 1;
                        r = r | 0;
                        r = r << 1;
                        r = r | 0;
                }
        }

        return r;
}

__mmask8 LE_w64(__m256i a, __m256i b) {
        __mmask8 r = 0x0;

        __m256i ret = _mm256_cmpgt_epi64(a, b);
        int* retf = (int*)&ret;

        for (int i = 0; i < 3; i++) {
                if (retf[3-i] == 0) {
                        r = r << 1;
                        r = r | 1;
                        r = r << 1;
                        r = r | 1;
                }
                else {
                        r = r << 1;
                        r = r | 0;
                        r = r << 1;
                        r = r | 0;
                }
        }

        return r;
}

__mmask8 EQ_w64(__m256i a, __m256i b) {
        __mmask8 r = 0x0;

        __m256i ret = _mm256_cmpeq_epi64(a, b);
        int* retf = (int*)&ret;

        for (int i = 0; i < 3; i++) {
                if (retf[3-i] == 0) {
                        r = r << 1;
                        r = r | 0;
                        r = r << 1;
                        r = r | 0;
                }
                else {
                        r = r << 1;
                        r = r | 1;
                        r = r << 1;
                        r = r | 1;
                }
        }

        return r;
}

__mmask8 NEQ_w64(__m256i a, __m256i b) {
        __mmask8 r = 0x0;

        __m256i ret = _mm256_cmpeq_epi64(a, b);
        int* retf = (int*)&ret;

        for (int i = 0; i < 3; i++) {
                if (retf[3-i] == 0) {
                        r = r << 1;
                        r = r | 1;
                        r = r << 1;
                        r = r | 1;
                }
                else {
                        r = r << 1;
                        r = r | 0;
                        r = r << 1;
                        r = r | 0;
                }
        }

        return r;
}

uint64_t extract_vec_w64(__m512i x, int pos) {  
  uint64_t buf[8] aligned;
  STORE_w64(buf, x);
  return buf[pos];
}

uint32_t extract_mask_w64(__mmask8 x, int pos) { return ((uint32_t)x) & (1 << pos); }

