#include "matrixmod.h"

#include <immintrin.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "params.h"
#include "vec.h"

void pmod_mat_print(pmod_mat_t *M, int M_r, int M_c) {
  pmod_mat_fprint(stdout, M, M_r, M_c);
}

void pmod_mat_fprint(FILE *stream, pmod_mat_t *M, int M_r, int M_c) {
  for (int r = 0; r < M_r; r++) {
    fprintf(stream, "[");
    for (int c = 0; c < M_c - 1; c++)
      fprintf(stream, "%4i ", pmod_mat_entry(M, M_r, M_c, r, c));
    fprintf(stream, "%4i", pmod_mat_entry(M, M_r, M_c, r, M_c - 1));
    fprintf(stream, "]\n");
  }
}

void pmod_mat_mul(pmod_mat_t *C, int C_r, int C_c, pmod_mat_t *A, int A_r,
                  int A_c, pmod_mat_t *B, int B_r, int B_c) {
  GFq_t tmp[C_r * C_c];

  for (int c = 0; c < C_c; c++)
    for (int r = 0; r < C_r; r++) {
      uint64_t val = 0;

      for (int i = 0; i < A_r; i++) {
        val = (val + (uint64_t)pmod_mat_entry(A, A_r, A_c, r, i) *
                         (uint64_t)pmod_mat_entry(B, B_r, B_c, i, c));
      }

      tmp[r * C_c + c] = val % MEDS_p;
    }

  for (int c = 0; c < C_c; c++)
    for (int r = 0; r < C_r; r++)
      pmod_mat_set_entry(C, C_r, C_c, r, c, tmp[r * C_c + c]);
}

void pmod_mat_mask_mul_vec(pmod_vec_t *C, int C_r, int C_c, pmod_vec_t *A,
                           int A_r, int A_c, pmod_vec_t *B, int B_r, int B_c,
                           int B_cm) {
  pmod_vec_t tmp[C_r * C_c];
  const pmod_vec_t zero = SET1(0);

  for (int c = B_cm; c < C_c; c++)
    for (int r = 0; r < C_r; r++) {
      pmod_vec_t val = zero;

      for (int i = 0; i < A_r; i++) {
        val = ADD(val, MULLO(pmod_mat_entry(A, A_r, A_c, r, i),
                             pmod_mat_entry(B, B_r, B_c, i, c)));
      }
      tmp[r * C_c + c] = GF_mod_vec(val);
    }

  for (int c = B_cm; c < C_c; c++)
    for (int r = 0; r < C_r; r++)
      pmod_mat_set_entry(C, C_r, C_c, r, c, tmp[r * C_c + c]);
}

void pmod_mat_mul_vec(pmod_vec_t *C, int C_r, int C_c, pmod_vec_t *A, int A_r,
                      int A_c, pmod_vec_t *B, int B_r, int B_c) {
  pmod_vec_t tmp[C_r * C_c];
  const pmod_vec_t zero = SET1(0);

  for (int c = 0; c < C_c; c++)
    for (int r = 0; r < C_r; r++) {
      pmod_vec_t val = zero;

      for (int i = 0; i < A_r; i++) {
        val = ADD(val, MULLO(pmod_mat_entry(A, A_r, A_c, r, i),
                             pmod_mat_entry(B, B_r, B_c, i, c)));
      }
      tmp[r * C_c + c] = GF_mod_vec(val);
    }

  for (int c = 0; c < C_c; c++)
    for (int r = 0; r < C_r; r++)
      pmod_mat_set_entry(C, C_r, C_c, r, c, tmp[r * C_c + c]);
}

int pmod_mat_syst_ct(pmod_mat_t *M, int M_r, int M_c) {
  if (pmod_mat_row_echelon_ct(M, M_r, M_c) < 0) return -1;

  return pmod_mat_back_substitution_ct(M, M_r, M_c);
}

// void pmod_mat_back_substitution_ct_vec(pmod_vec_t *M, int M_r, int M_c) {
//   const pmod_vec_t zero = SET1(0);
//   const pmod_vec_t p = SET1(MEDS_p);
//
//   for (int r = M_r - 1; r >= 0; r--)
//     for (int r2 = 0; r2 < r; r2++) {
//       pmod_vec_t factor = pmod_mat_entry(M, M_r, M_c, r2, r);
//
//       pmod_vec_t tmp0 = pmod_mat_entry(M, M_r, M_c, r, r);
//       pmod_vec_t tmp1 = pmod_mat_entry(M, M_r, M_c, r2, r);
//
//       pmod_vec_t val = GF_reduc_vec(MULLO(tmp0, factor));
//
//       val = SUB(tmp1, val);
//       // val = SUB(tmp1, factor);
//       val = ADD_M(val, p, LT(val, zero));
//
//       pmod_mat_set_entry(M, M_r, M_c, r2, r, val);
//
//       for (int c = M_r; c < M_c; c++) {
//         pmod_vec_t tmp0 = pmod_mat_entry(M, M_r, M_c, r, c);
//         pmod_vec_t tmp1 = pmod_mat_entry(M, M_r, M_c, r2, c);
//
//         pmod_vec_t val = GF_reduc_vec(MULLO(tmp0, factor));
//
//         val = SUB(tmp1, val);
//         val = ADD_M(val, p, LT(val, zero));
//
//         pmod_mat_set_entry(M, M_r, M_c, r2, c, val);
//       }
//     }
// }

pmod_vec_mask_t pmod_mat_syst_ct_vec(pmod_vec_t *M, int M_r, int M_c) {
  // const pmod_vec_t zero = SET1(0);
  // const pmod_vec_t p = SET1(MEDS_p);
  //
  // __mmask16 valid = 0xFFFF;
  //
  // for (int r = 0; r < M_r; r++) {
  //   for (int r2 = r + 1; r2 < M_r; r2++) {
  //     pmod_vec_t Mrr = pmod_mat_entry(M, M_r, M_c, r, r);
  //     __mmask16 Mrr_eq_zero = EQ(Mrr, zero);
  //
  //     for (int c = r; c < M_c; c++) {
  //       pmod_vec_t val = pmod_mat_entry(M, M_r, M_c, r2, c);
  //
  //       pmod_vec_t Mrc = pmod_mat_entry(M, M_r, M_c, r, c);
  //       pmod_vec_t res = ADD_M(Mrc, val, Mrr_eq_zero);
  //       res = SUB_M(res, p, GE(res, p));
  //       pmod_mat_set_entry(M, M_r, M_c, r, c, res);
  //     }
  //   }
  //
  //   pmod_vec_t val = pmod_mat_entry(M, M_r, M_c, r, r);
  //   __mmask16 val_neq_zero = NEQ(val, zero);
  //   valid = valid & val_neq_zero;
  //
  //   val = GF_inv_vec(val);
  //   for (int c = r; c < M_c; c++) {
  //     pmod_vec_t Mrc = pmod_mat_entry(M, M_r, M_c, r, c);
  //     pmod_vec_t res = GF_reduc_vec(MULLO(Mrc, val));
  //     pmod_mat_set_entry(M, M_r, M_c, r, c, res);
  //   }
  //
  //   for (int r2 = r + 1; r2 < M_r; r2++) {
  //     pmod_vec_t factor = pmod_mat_entry(M, M_r, M_c, r2, r);
  //
  //     for (int c = r; c < M_c; c++) {
  //       pmod_vec_t tmp0 = pmod_mat_entry(M, M_r, M_c, r, c);
  //       pmod_vec_t tmp1 = pmod_mat_entry(M, M_r, M_c, r2, c);
  //       pmod_vec_t val = GF_reduc_vec(MULLO(tmp0, factor));
  //       val = SUB(tmp1, val);
  //       val = ADD_M(val, p, LT(val, zero));
  //       pmod_mat_set_entry(M, M_r, M_c, r2, c, val);
  //     }
  //   }
  // }
  //
  // pmod_mat_back_substitution_ct_vec(M, M_r, M_c);

  pmod_vec_t M_inv[M_r * M_r];
  pmod_vec_mask_t valid = pmod_mat_inv_vec(M_inv, M, M_r, M_c);
  pmod_mat_mask_mul_vec(M, M_r, M_c, M_inv, M_r, M_r, M, M_r, M_c, M_r);

  return valid;
}

int pmod_mat_row_echelon_ct(pmod_mat_t *M, int M_r, int M_c) {
  for (int r = 0; r < M_r; r++) {
    // swap
    for (int r2 = r + 1; r2 < M_r; r2++) {
      uint64_t Mrr = pmod_mat_entry(M, M_r, M_c, r, r);

      for (int c = r; c < M_c; c++) {
        uint64_t val = pmod_mat_entry(M, M_r, M_c, r2, c);

        uint64_t Mrc = pmod_mat_entry(M, M_r, M_c, r, c);

        pmod_mat_set_entry(M, M_r, M_c, r, c,
                           (Mrc + val * (Mrr == 0)) % MEDS_p);
      }
    }

    uint64_t val = pmod_mat_entry(M, M_r, M_c, r, r);

    if (val == 0) return -1;

    val = GF_inv(val);

    // normalize
    for (int c = r; c < M_c; c++) {
      uint64_t tmp =
          ((uint64_t)pmod_mat_entry(M, M_r, M_c, r, c) * val) % MEDS_p;
      pmod_mat_set_entry(M, M_r, M_c, r, c, tmp);
    }

    // eliminate
    for (int r2 = r + 1; r2 < M_r; r2++) {
      uint64_t factor = pmod_mat_entry(M, M_r, M_c, r2, r);

      for (int c = r; c < M_c; c++) {
        uint64_t tmp0 = pmod_mat_entry(M, M_r, M_c, r, c);
        uint64_t tmp1 = pmod_mat_entry(M, M_r, M_c, r2, c);

        int64_t val = (tmp0 * factor) % MEDS_p;

        val = tmp1 - val;

        val += MEDS_p * (val < 0);

        pmod_mat_set_entry(M, M_r, M_c, r2, c, val);
      }
    }
  }

  return 0;
}

int pmod_mat_back_substitution_ct(pmod_mat_t *M, int M_r, int M_c) {
  // back substitution
  for (int r = M_r - 1; r >= 0; r--)
    for (int r2 = 0; r2 < r; r2++) {
      uint64_t factor = pmod_mat_entry(M, M_r, M_c, r2, r);

      uint64_t tmp0 = pmod_mat_entry(M, M_r, M_c, r, r);
      uint64_t tmp1 = pmod_mat_entry(M, M_r, M_c, r2, r);

      int64_t val = (tmp0 * factor) % MEDS_p;

      val = tmp1 - val;

      val += MEDS_p * (val < 0);

      pmod_mat_set_entry(M, M_r, M_c, r2, r, val);

      for (int c = M_r; c < M_c; c++) {
        uint64_t tmp0 = pmod_mat_entry(M, M_r, M_c, r, c);
        uint64_t tmp1 = pmod_mat_entry(M, M_r, M_c, r2, c);

        int val = (tmp0 * factor) % MEDS_p;

        val = tmp1 - val;

        val += MEDS_p * (val < 0);

        pmod_mat_set_entry(M, M_r, M_c, r2, c, val);
      }
    }

  return 0;
}

GFq_t GF_inv(GFq_t val) {
  // if (MEDS_p == 8191) {
  if (0) {
    // Use optimal addition chain...
    uint64_t tmp_0 = val;
    uint64_t tmp_1 = (tmp_0 * tmp_0) % MEDS_p;
    uint64_t tmp_2 = (tmp_1 * tmp_0) % MEDS_p;
    uint64_t tmp_3 = (tmp_2 * tmp_1) % MEDS_p;
    uint64_t tmp_4 = (tmp_3 * tmp_3) % MEDS_p;
    uint64_t tmp_5 = (tmp_4 * tmp_3) % MEDS_p;
    uint64_t tmp_6 = (tmp_5 * tmp_5) % MEDS_p;
    uint64_t tmp_7 = (tmp_6 * tmp_6) % MEDS_p;
    uint64_t tmp_8 = (tmp_7 * tmp_7) % MEDS_p;
    uint64_t tmp_9 = (tmp_8 * tmp_8) % MEDS_p;
    uint64_t tmp_10 = (tmp_9 * tmp_5) % MEDS_p;
    uint64_t tmp_11 = (tmp_10 * tmp_10) % MEDS_p;
    uint64_t tmp_12 = (tmp_11 * tmp_11) % MEDS_p;
    uint64_t tmp_13 = (tmp_12 * tmp_2) % MEDS_p;
    uint64_t tmp_14 = (tmp_13 * tmp_13) % MEDS_p;
    uint64_t tmp_15 = (tmp_14 * tmp_14) % MEDS_p;
    uint64_t tmp_16 = (tmp_15 * tmp_15) % MEDS_p;
    uint64_t tmp_17 = (tmp_16 * tmp_3) % MEDS_p;

    return tmp_17;
  } else {
    uint64_t exponent = MEDS_p - 2;
    uint64_t t = 1;

    while (exponent > 0) {
      if ((exponent & 1) != 0) t = (t * (uint64_t)val) % MEDS_p;

      val = ((uint64_t)val * (uint64_t)val) % MEDS_p;

      exponent >>= 1;
    }

    return t;
  }
}

pmod_vec_mask_t pmod_mat_inv_vec(pmod_vec_t *M_inv, pmod_vec_t *M, int M_r, int M_c) {
  const pmod_vec_t zero = SET1(0);
  const pmod_vec_t one = SET1(1);
  const pmod_vec_t p = SET1(MEDS_p);

  pmod_vec_t M_prime[M_r * M_r];
  for (int r = 0; r < M_r; r++)
    memcpy(M_prime + r * M_r, M + r * M_c, sizeof(pmod_vec_t) * M_r);

  for (int i = 0; i < M_r; i++)
    for (int j = 0; j < M_r; j++)
      pmod_mat_set_entry(M_inv, M_r, M_r, i, j, zero);
  for (int i = 0; i < M_r; i++) pmod_mat_set_entry(M_inv, M_r, M_r, i, i, one);

  __mmask16 valid = 0xFFFF;

  for (int r = 0; r < M_r; r++) {
    for (int r2 = r + 1; r2 < M_r; r2++) {
      pmod_vec_t Mrr = pmod_mat_entry(M_prime, M_r, M_r, r, r);
      __mmask16 Mrr_eq_zero = EQ(Mrr, zero);

      for (int c = r; c < M_r; c++) {
        pmod_vec_t val = pmod_mat_entry(M_prime, M_r, M_r, r2, c);
        pmod_vec_t Mrc = pmod_mat_entry(M_prime, M_r, M_r, r, c);
        pmod_vec_t res = ADD_M(Mrc, val, Mrr_eq_zero);
        res = SUB_M(res, p, GE(res, p));
        pmod_mat_set_entry(M_prime, M_r, M_r, r, c, res);
      }
      for (int c = 0; c < M_r; c++) {
        pmod_vec_t val = pmod_mat_entry(M_inv, M_r, M_r, r2, c);
        pmod_vec_t Mrc = pmod_mat_entry(M_inv, M_r, M_r, r, c);
        pmod_vec_t res = ADD_M(Mrc, val, Mrr_eq_zero);
        res = SUB_M(res, p, GE(res, p));
        pmod_mat_set_entry(M_inv, M_r, M_r, r, c, res);
      }
    }

    pmod_vec_t val = pmod_mat_entry(M_prime, M_r, M_r, r, r);
    __mmask16 val_neq_zero = NEQ(val, zero);
    valid = valid & val_neq_zero;

    val = GF_inv_vec(val);
    for (int c = r; c < M_r; c++) {
      pmod_vec_t Mrc = pmod_mat_entry(M_prime, M_r, M_r, r, c);
      pmod_vec_t res = GF_reduc_vec(MULLO(Mrc, val));
      pmod_mat_set_entry(M_prime, M_r, M_r, r, c, res);
    }
    for (int c = 0; c < M_r; c++) {
      pmod_vec_t Mrc = pmod_mat_entry(M_inv, M_r, M_r, r, c);
      pmod_vec_t res = GF_reduc_vec(MULLO(Mrc, val));
      pmod_mat_set_entry(M_inv, M_r, M_r, r, c, res);
    }

    for (int r2 = r + 1; r2 < M_r; r2++) {
      pmod_vec_t factor = pmod_mat_entry(M_prime, M_r, M_r, r2, r);

      for (int c = r; c < M_r; c++) {
        pmod_vec_t tmp0 = pmod_mat_entry(M_prime, M_r, M_r, r, c);
        pmod_vec_t tmp1 = pmod_mat_entry(M_prime, M_r, M_r, r2, c);
        pmod_vec_t val = GF_reduc_vec(MULLO(tmp0, factor));
        val = SUB(tmp1, val);
        val = ADD_M(val, p, LT(val, zero));
        pmod_mat_set_entry(M_prime, M_r, M_r, r2, c, val);
      }
      for (int c = 0; c < M_r; c++) {
        pmod_vec_t tmp0 = pmod_mat_entry(M_inv, M_r, M_r, r, c);
        pmod_vec_t tmp1 = pmod_mat_entry(M_inv, M_r, M_r, r2, c);
        val = GF_reduc_vec(MULLO(tmp0, factor));
        val = SUB(tmp1, val);
        val = ADD_M(val, p, LT(val, zero));
        pmod_mat_set_entry(M_inv, M_r, M_r, r2, c, val);
      }
    }
  }

  // back substitution
  for (int r = M_r - 1; r >= 0; r--)
    for (int r2 = 0; r2 < r; r2++) {
      pmod_vec_t factor = pmod_mat_entry(M_prime, M_r, M_r, r2, r);

      pmod_vec_t tmp0 = pmod_mat_entry(M_prime, M_r, M_r, r, r);
      pmod_vec_t tmp1 = pmod_mat_entry(M_prime, M_r, M_r, r2, r);
      pmod_vec_t val = GF_reduc_vec(MULLO(tmp0, factor));
      val = SUB(tmp1, val);
      val = ADD_M(val, p, LT(val, zero));
      pmod_mat_set_entry(M_prime, M_r, M_r, r2, r, val);
      for (int c = 0; c < M_r; c++) {
        tmp0 = pmod_mat_entry(M_inv, M_r, M_r, r, c);
        tmp1 = pmod_mat_entry(M_inv, M_r, M_r, r2, c);
        val = GF_reduc_vec(MULLO(tmp0, factor));
        val = SUB(tmp1, val);
        val = ADD_M(val, p, LT(val, zero));
        pmod_mat_set_entry(M_inv, M_r, M_r, r2, c, val);
      }
    }

  return valid;
}

int pmod_mat_inv(pmod_mat_t *B, pmod_mat_t *A, int A_r, int A_c) {
  pmod_mat_t M[A_r * A_c * 2];

  for (int r = 0; r < A_r; r++) {
    memcpy(&M[r * A_c * 2], &A[r * A_c], A_c * sizeof(GFq_t));

    for (int c = 0; c < A_c; c++)
      pmod_mat_set_entry(M, A_r, A_c * 2, r, A_c + c, r == c ? 1 : 0);
  }

  int ret = pmod_mat_syst_ct(M, A_r, A_c * 2);

  if ((ret == 0) && B)
    for (int r = 0; r < A_r; r++)
      memcpy(&B[r * A_c], &M[r * A_c * 2 + A_c], A_c * sizeof(GFq_t));

  return ret;
}

// int pmod_mat_inv_vec(pmod_vec_t *B, pmod_vec_t *A, int A_r, int A_c) {
//   pmod_mat_t M[A_r * A_c * 2];
//
//   for (int r = 0; r < A_r; r++) {
//     memcpy(&M[r * A_c * 2], &A[r * A_c], A_c * sizeof(GFq_t));
//
//     for (int c = 0; c < A_c; c++)
//       pmod_mat_set_entry(M, A_r, A_c * 2, r, A_c + c, r == c ? 1 : 0);
//   }
//
//   int ret = pmod_mat_syst_ct(M, A_r, A_c * 2);
//
//   if ((ret == 0) && B)
//     for (int r = 0; r < A_r; r++)
//       memcpy(&B[r * A_c], &M[r * A_c * 2 + A_c], A_c * sizeof(GFq_t));
//
//   return ret;
// }
