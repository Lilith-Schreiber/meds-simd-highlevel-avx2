#include "matrixmod_x16.h"

pmod_mat_x16_t pmod_mat_entry_x16(uint16_t *M[], int M_r, int M_c, int r, int c) {
  uint32_t M_buf[16] align64 = {0};
  for (int i = 0; i < 16; i++) {
    M_buf[i] = pmod_mat_entry(M[i], M_r, M_c, r, c);
  }

  return LOAD(M_buf);
}

void pmod_mat_set_entry_x16(uint16_t *M[], int M_r, int M_c, int r, int c,
                            pmod_mat_x16_t val, int num) {
  uint32_t buf[16] align64 = {0};
  STORE((int *)buf, val);

  for (int i = 0; i < num; i++) {
    const pmod_mat_t entry = buf[i];
    pmod_mat_set_entry(M[i], M_r, M_c, r, c, entry);
  }
}

void pmod_mat_mul_x16(pmod_mat_x16_t *C, int C_r, int C_c, pmod_mat_x16_t *A, int A_r,
                      int A_c, pmod_mat_x16_t *B, int B_r, int B_c) {
  pmod_mat_x16_t tmp[C_r * C_c];
  const pmod_mat_x16_t zero = SET1(0);

  for (int c = 0; c < C_c; c++)
    for (int r = 0; r < C_r; r++) {
      pmod_mat_x16_t val = zero;

      for (int i = 0; i < A_r; i++) {
        val = ADD(val, MULLO(pmod_mat_entry(A, A_r, A_c, r, i),
                             pmod_mat_entry(B, B_r, B_c, i, c)));
      }
      tmp[r * C_c + c] = GF_mod_x16(val);
    }

  for (int c = 0; c < C_c; c++)
    for (int r = 0; r < C_r; r++)
      pmod_mat_set_entry(C, C_r, C_c, r, c, tmp[r * C_c + c]);
}

void pmod_mat_mask_mul_x16(pmod_mat_x16_t *C, int C_r, int C_c, pmod_mat_x16_t *A,
                           int A_r, int A_c, pmod_mat_x16_t *B, int B_r, int B_c,
                           int B_cm) {
  pmod_mat_x16_t tmp[C_r * C_c];
  const pmod_mat_x16_t zero = SET1(0);

  for (int c = B_cm; c < C_c; c++)
    for (int r = 0; r < C_r; r++) {
      pmod_mat_x16_t val = zero;

      for (int i = 0; i < A_r; i++) {
        val = ADD(val, MULLO(pmod_mat_entry(A, A_r, A_c, r, i),
                             pmod_mat_entry(B, B_r, B_c, i, c)));
      }
      tmp[r * C_c + c] = GF_mod_x16(val);
    }

  for (int c = B_cm; c < C_c; c++)
    for (int r = 0; r < C_r; r++)
      pmod_mat_set_entry(C, C_r, C_c, r, c, tmp[r * C_c + c]);
}

pmod_vec_mask_t pmod_mat_syst_ct_x16(pmod_mat_x16_t *M, int M_r, int M_c) {
  pmod_mat_x16_t M_inv[M_r * M_r];
  pmod_vec_mask_t valid = pmod_mat_inv_x16(M_inv, M, M_r, M_c);
  pmod_mat_mask_mul_x16(M, M_r, M_c, M_inv, M_r, M_r, M, M_r, M_c, M_r);

  return valid;
}

pmod_vec_mask_t pmod_mat_inv_x16(pmod_mat_x16_t *M_inv, pmod_mat_x16_t *M, int M_r, int M_c) {
  const pmod_mat_x16_t zero = SET1(0);
  const pmod_mat_x16_t one = SET1(1);
  const pmod_mat_x16_t p = SET1(MEDS_p);

  pmod_mat_x16_t M_prime[M_r * M_r];
  for (int r = 0; r < M_r; r++)
    memcpy(M_prime + r * M_r, M + r * M_c, sizeof(pmod_mat_x16_t) * M_r);

  for (int i = 0; i < M_r; i++)
    for (int j = 0; j < M_r; j++)
      pmod_mat_set_entry(M_inv, M_r, M_r, i, j, zero);
  for (int i = 0; i < M_r; i++) pmod_mat_set_entry(M_inv, M_r, M_r, i, i, one);

  __mmask16 valid = 0xFFFF;

  for (int r = 0; r < M_r; r++) {
    for (int r2 = r + 1; r2 < M_r; r2++) {
      pmod_mat_x16_t Mrr = pmod_mat_entry(M_prime, M_r, M_r, r, r);
      __mmask16 Mrr_eq_zero = EQ(Mrr, zero);

      for (int c = r; c < M_r; c++) {
        pmod_mat_x16_t val = pmod_mat_entry(M_prime, M_r, M_r, r2, c);
        pmod_mat_x16_t Mrc = pmod_mat_entry(M_prime, M_r, M_r, r, c);
        pmod_mat_x16_t res = ADD_M(Mrc, val, Mrr_eq_zero);
        res = SUB_M(res, p, GE(res, p));
        pmod_mat_set_entry(M_prime, M_r, M_r, r, c, res);
      }
      for (int c = 0; c < M_r; c++) {
        pmod_mat_x16_t val = pmod_mat_entry(M_inv, M_r, M_r, r2, c);
        pmod_mat_x16_t Mrc = pmod_mat_entry(M_inv, M_r, M_r, r, c);
        pmod_mat_x16_t res = ADD_M(Mrc, val, Mrr_eq_zero);
        res = SUB_M(res, p, GE(res, p));
        pmod_mat_set_entry(M_inv, M_r, M_r, r, c, res);
      }
    }

    pmod_mat_x16_t val = pmod_mat_entry(M_prime, M_r, M_r, r, r);
    __mmask16 val_neq_zero = NEQ(val, zero);
    valid = valid & val_neq_zero;

    val = GF_inv_x16(val);
    for (int c = r; c < M_r; c++) {
      pmod_mat_x16_t Mrc = pmod_mat_entry(M_prime, M_r, M_r, r, c);
      pmod_mat_x16_t res = GF_reduc_x16(MULLO(Mrc, val));
      pmod_mat_set_entry(M_prime, M_r, M_r, r, c, res);
    }
    for (int c = 0; c < M_r; c++) {
      pmod_mat_x16_t Mrc = pmod_mat_entry(M_inv, M_r, M_r, r, c);
      pmod_mat_x16_t res = GF_reduc_x16(MULLO(Mrc, val));
      pmod_mat_set_entry(M_inv, M_r, M_r, r, c, res);
    }

    for (int r2 = r + 1; r2 < M_r; r2++) {
      pmod_mat_x16_t factor = pmod_mat_entry(M_prime, M_r, M_r, r2, r);

      for (int c = r; c < M_r; c++) {
        pmod_mat_x16_t tmp0 = pmod_mat_entry(M_prime, M_r, M_r, r, c);
        pmod_mat_x16_t tmp1 = pmod_mat_entry(M_prime, M_r, M_r, r2, c);
        pmod_mat_x16_t val = GF_reduc_x16(MULLO(tmp0, factor));
        val = SUB(tmp1, val);
        val = ADD_M(val, p, LT(val, zero));
        pmod_mat_set_entry(M_prime, M_r, M_r, r2, c, val);
      }
      for (int c = 0; c < M_r; c++) {
        pmod_mat_x16_t tmp0 = pmod_mat_entry(M_inv, M_r, M_r, r, c);
        pmod_mat_x16_t tmp1 = pmod_mat_entry(M_inv, M_r, M_r, r2, c);
        val = GF_reduc_x16(MULLO(tmp0, factor));
        val = SUB(tmp1, val);
        val = ADD_M(val, p, LT(val, zero));
        pmod_mat_set_entry(M_inv, M_r, M_r, r2, c, val);
      }
    }
  }

  // back substitution
  for (int r = M_r - 1; r >= 0; r--)
    for (int r2 = 0; r2 < r; r2++) {
      pmod_mat_x16_t factor = pmod_mat_entry(M_prime, M_r, M_r, r2, r);

      pmod_mat_x16_t tmp0 = pmod_mat_entry(M_prime, M_r, M_r, r, r);
      pmod_mat_x16_t tmp1 = pmod_mat_entry(M_prime, M_r, M_r, r2, r);
      pmod_mat_x16_t val = GF_reduc_x16(MULLO(tmp0, factor));
      val = SUB(tmp1, val);
      val = ADD_M(val, p, LT(val, zero));
      pmod_mat_set_entry(M_prime, M_r, M_r, r2, r, val);
      for (int c = 0; c < M_r; c++) {
        tmp0 = pmod_mat_entry(M_inv, M_r, M_r, r, c);
        tmp1 = pmod_mat_entry(M_inv, M_r, M_r, r2, c);
        val = GF_reduc_x16(MULLO(tmp0, factor));
        val = SUB(tmp1, val);
        val = ADD_M(val, p, LT(val, zero));
        pmod_mat_set_entry(M_inv, M_r, M_r, r2, c, val);
      }
    }

  return valid;
}
