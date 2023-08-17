#include "util_vec.h"

pmod_mat_mask_w32_t solve_w32(pmod_mat_w32_t *A, pmod_mat_w32_t *B_inv, pmod_mat_w32_t *G0prime,
                          pmod_mat_w32_t Amm) {
  const pmod_mat_w32_t p = SET1_w32(MEDS_p);
  const pmod_mat_w32_t zero = SET1_w32(0);
  const pmod_mat_w32_t one = SET1_w32(1);

  pmod_mat_w32_t P0prime0[MEDS_m * MEDS_n];
  pmod_mat_w32_t P0prime1[MEDS_m * MEDS_n];

  for (int i = 0; i < MEDS_m * MEDS_n; i++) {
    P0prime0[i] = G0prime[i];
    P0prime1[i] = G0prime[i + MEDS_m * MEDS_n];
  }

  pmod_mat_w32_t N[MEDS_n * MEDS_m];

  for (int i = 0; i < MEDS_m; i++)
    for (int j = 0; j < MEDS_n; j++)
      N[j * MEDS_m + i] = GF_mod_w32(SUB_w32(p, P0prime0[i * MEDS_n + j]));

  pmod_mat_w32_t M[MEDS_n * (MEDS_m + MEDS_m + 2)] = {0};

  for (int i = 0; i < MEDS_m; i++)
    for (int j = 0; j < MEDS_n; j++)
      M[j * (MEDS_m + MEDS_m + 2) + i] =
          GF_mod_w32(SUB_w32(p, P0prime1[i * MEDS_n + j]));

  for (int i = 0; i < MEDS_m; i++)
    for (int j = 0; j < MEDS_n; j++)
      M[j * (MEDS_m + MEDS_m + 2) + i + MEDS_n] = P0prime0[i * MEDS_n + j];

  for (int j = 0; j < MEDS_n; j++)
    M[j * (MEDS_m + MEDS_m + 2) + MEDS_m + MEDS_n] =
        GF_mod_w32(MULLO_w32(P0prime0[(MEDS_m - 1) * MEDS_n + j], SUB_w32(p, Amm)));

  for (int j = 0; j < MEDS_n; j++)
    M[j * (MEDS_m + MEDS_m + 2) + MEDS_m + MEDS_n + 1] =
        GF_mod_w32(MULLO_w32(P0prime1[(MEDS_m - 1) * MEDS_n + j], Amm));

  pmod_mat_mask_w32_t valid =
      pmod_mat_syst_ct_w32(M, MEDS_n - 1, MEDS_m + MEDS_m + 2);

  // eliminate last row
  for (int r = 0; r < MEDS_n - 1; r++) {
    pmod_mat_w32_t factor =
        pmod_mat_entry(M, MEDS_n, MEDS_m + MEDS_m + 2, MEDS_n - 1, r);

    // ignore last column
    for (int c = MEDS_n - 1; c < MEDS_m + MEDS_m + 1; c++) {
      pmod_mat_w32_t tmp0 =
          pmod_mat_entry(M, MEDS_n, MEDS_m + MEDS_m + 2, MEDS_n - 1, c);
      pmod_mat_w32_t tmp1 = pmod_mat_entry(M, MEDS_n, MEDS_m + MEDS_m + 2, r, c);

      pmod_mat_w32_t val = GF_mod_w32(MULLO_w32(tmp1, factor));

      val = SUB_w32(tmp0, val);

      val = ADD_M_w32(val, p, LT_w32(val, zero));

      pmod_mat_set_entry(M, MEDS_n, MEDS_m + MEDS_m + 2, MEDS_n - 1, c, val);
    }

    pmod_mat_set_entry(M, MEDS_n, MEDS_m + MEDS_m + 2, MEDS_n - 1, r, zero);
  }

  // normalize last row
  {
    pmod_mat_w32_t val =
        pmod_mat_entry(M, MEDS_n, MEDS_m + MEDS_m + 2, MEDS_n - 1, MEDS_n - 1);

    valid = valid & NEQ_w32(val, zero);

    val = GF_inv_w32(val);

    // ignore last column
    for (int c = MEDS_n; c < MEDS_m + MEDS_m + 1; c++) {
      pmod_mat_w32_t tmp =
          pmod_mat_entry(M, MEDS_n, MEDS_m + MEDS_m + 2, MEDS_n - 1, c);

      tmp = GF_mod_w32(MULLO_w32(tmp, val));

      pmod_mat_set_entry(M, MEDS_n, MEDS_m + MEDS_m + 2, MEDS_n - 1, c, tmp);
    }
  }

  pmod_mat_set_entry(M, MEDS_n, MEDS_m + MEDS_m + 2, MEDS_n - 1, MEDS_n - 1,
                     one);

  M[MEDS_n * (MEDS_m + MEDS_m + 2) - 1] = zero;

  // back substitute
  for (int r = 0; r < MEDS_n - 1; r++) {
    pmod_mat_w32_t factor =
        pmod_mat_entry(M, MEDS_n, MEDS_m + MEDS_m + 2, r, MEDS_n - 1);

    // ignore last column
    for (int c = MEDS_n; c < MEDS_m + MEDS_m + 1; c++) {
      pmod_mat_w32_t tmp0 =
          pmod_mat_entry(M, MEDS_n, MEDS_m + MEDS_m + 2, MEDS_n - 1, c);
      pmod_mat_w32_t tmp1 = pmod_mat_entry(M, MEDS_n, MEDS_m + MEDS_m + 2, r, c);

      pmod_mat_w32_t val = GF_mod_w32(MULLO_w32(tmp0, factor));

      val = SUB_w32(tmp1, val);

      val = ADD_M_w32(val, p, LT_w32(val, zero));

      pmod_mat_set_entry(M, MEDS_n, MEDS_m + MEDS_m + 2, r, c, val);
    }

    pmod_mat_set_entry(M, M_r, MEDS_m + MEDS_m + 2, r, MEDS_n - 1, zero);
  }

  pmod_mat_w32_t sol[MEDS_n * MEDS_n + MEDS_m * MEDS_m] = {0};

  sol[MEDS_n * MEDS_n + MEDS_m * MEDS_m - 1] = Amm;

  for (int i = 0; i < MEDS_n - 1; i++)
    sol[MEDS_n * MEDS_n + MEDS_m * MEDS_m - MEDS_n + i] =
        M[(i + 1) * (MEDS_m + MEDS_m + 2) - 1];

  for (int i = 0; i < MEDS_n; i++)
    sol[MEDS_n * MEDS_n + MEDS_m * MEDS_m - 2 * MEDS_n + i] =
        M[(i + 1) * (MEDS_m + MEDS_m + 2) - 2];

  for (int i = 0; i < MEDS_n; i++)
    sol[MEDS_n * MEDS_n - MEDS_n + i] =
        GF_mod_w32(MULLO_w32(P0prime0[(MEDS_m - 1) * MEDS_n + i], Amm));

  // incomplete blocks:

  for (int i = 0; i < MEDS_n; i++)
    for (int j = 0; j < MEDS_n - 1; j++) {
      pmod_mat_w32_t tmp =
          SUB_w32(ADD_w32(sol[MEDS_n * MEDS_n + MEDS_m * MEDS_m - 2 * MEDS_n + i], p),
              GF_mod_w32(
                  MULLO_w32(M[i * (MEDS_m + MEDS_m + 2) + MEDS_n + MEDS_n - 2 - j],
                        sol[MEDS_n * MEDS_n + MEDS_m * MEDS_m - 2 - j])));
      sol[MEDS_n * MEDS_n + MEDS_m * MEDS_m - 2 * MEDS_n + i] = GF_mod_w32(tmp);
    }

  for (int i = 0; i < MEDS_n; i++)
    for (int j = 0; j < MEDS_n - 1; j++) {
      pmod_mat_w32_t tmp = SUB_w32(
          ADD_w32(sol[MEDS_n * MEDS_n - MEDS_n + i], p),
          GF_mod_w32(MULLO_w32(N[i * (MEDS_n) + MEDS_m - 2 - j],
                           sol[MEDS_n * MEDS_n + MEDS_m * MEDS_m - 2 - j])));
      sol[MEDS_n * MEDS_n - MEDS_n + i] = GF_mod_w32(tmp);
    }

  // complete blocks:

  for (int block = 3; block <= MEDS_n; block++)
    for (int i = 0; i < MEDS_n; i++)
      for (int j = 0; j < MEDS_n; j++) {
        pmod_mat_w32_t tmp = SUB_w32(
            ADD_w32(sol[MEDS_n * MEDS_n + MEDS_m * MEDS_m - block * MEDS_n + i], p),
            GF_mod_w32(
                MULLO_w32(M[i * (MEDS_m + MEDS_m + 2) + MEDS_n + MEDS_n - 1 - j],
                      sol[MEDS_n * MEDS_n + MEDS_m * MEDS_m - 1 -
                          (block - 2) * MEDS_n - j])));
        sol[MEDS_n * MEDS_n + MEDS_m * MEDS_m - block * MEDS_n + i] =
            GF_mod_w32(tmp);
      }

  for (int block = 2; block <= MEDS_n; block++)
    for (int i = 0; i < MEDS_n; i++)
      for (int j = 0; j < MEDS_n; j++) {
        pmod_mat_w32_t tmp =
            SUB_w32(ADD_w32(sol[MEDS_n * MEDS_n - block * MEDS_n + i], p),
                GF_mod_w32(MULLO_w32(N[i * (MEDS_n) + MEDS_m - 1 - j],
                                 sol[MEDS_n * MEDS_n + MEDS_m * MEDS_m - 1 -
                                     (block - 1) * MEDS_n - j])));
        sol[MEDS_n * MEDS_n - block * MEDS_n + i] = GF_mod_w32(tmp);
      }

  for (int i = 0; i < MEDS_m * MEDS_m; i++) A[i] = sol[i + MEDS_n * MEDS_n];

  for (int i = 0; i < MEDS_n * MEDS_n; i++) B_inv[i] = sol[i];

  return valid;
}

void pi_w32(pmod_mat_w32_t *Go, pmod_mat_w32_t *A, pmod_mat_w32_t *B, pmod_mat_w32_t *G0) {
  pmod_mat_w32_t *Go_sub[MEDS_k];
  pmod_mat_w32_t *G0_sub[MEDS_k];

  for (int i = 0; i < MEDS_k; i++) Go_sub[i] = &Go[i * MEDS_m * MEDS_n];
  for (int i = 0; i < MEDS_k; i++) G0_sub[i] = &G0[i * MEDS_m * MEDS_n];

  for (int i = 0; i < MEDS_k; i++) {
    pmod_mat_mul_w32(Go_sub[i], MEDS_m, MEDS_n, A, MEDS_m, MEDS_m, G0_sub[i],
                     MEDS_m, MEDS_n);
    pmod_mat_mul_w32(Go_sub[i], MEDS_m, MEDS_n, Go_sub[i], MEDS_m, MEDS_n, B,
                     MEDS_n, MEDS_n);
  }
}

