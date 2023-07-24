#include "meds.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "api.h"
#include "bitstream.h"
#include "fips202.h"
#include "log.h"
#include "matrixmod.h"
#include "measure.h"
#include "params.h"
#include "randombytes.h"
#include "seed.h"
#include "util.h"
#include "vec.h"

#define CEILING(x, y) (((x) + (y)-1) / (y))

#define START(a) long long a = -cpucycles()
#define END(a) a += cpucycles()
#define MIN(cur, new) \
  if (new < cur) {    \
    cur = new;        \
  }

extern FILE *measure_log;
const char *measure_log_file = "./measure.txt";

int crypto_sign_keypair(unsigned char *pk, unsigned char *sk) {
  clear_measure_log(measure_log_file);
  open_measure_log(measure_log_file);

  double freq = osfreq();
  // LOG_M("Generating keypairs...\n\n");

  /**
    Generate G0
    **/
  START(gen_G0_time);

  uint8_t delta[MEDS_sec_seed_bytes];

  randombytes(delta, MEDS_sec_seed_bytes);

  pmod_mat_t G_data[MEDS_k * MEDS_m * MEDS_n * MEDS_s];
  pmod_mat_t *G[MEDS_s];

  for (int i = 0; i < MEDS_s; i++) G[i] = &G_data[i * MEDS_k * MEDS_m * MEDS_n];

  uint8_t sigma_G0[MEDS_pub_seed_bytes];
  uint8_t sigma[MEDS_sec_seed_bytes];

  XOF((uint8_t *[]){sigma_G0, sigma},
      (size_t[]){MEDS_pub_seed_bytes, MEDS_sec_seed_bytes}, delta,
      MEDS_sec_seed_bytes, 2);

  LOG_VEC(sigma, MEDS_sec_seed_bytes);
  LOG_VEC_FMT(sigma_G0, MEDS_pub_seed_bytes, "sigma_G0");

  rnd_sys_mat(G[0], MEDS_k, MEDS_m * MEDS_n, sigma_G0, MEDS_pub_seed_bytes);

  LOG_MAT(G[0], MEDS_k, MEDS_m * MEDS_n);

  // END(gen_G0_time);
  //  LOG_M("  generate G0:\n");
  //  LOG_M("    %f   (%llu cycles)\n", gen_G0_time / freq, gen_G0_time);

  /**
    Generate public and secret matrices
    **/
  START(gen_all_pm_sm_time);

  pmod_mat_t A_inv_data[MEDS_s * MEDS_m * MEDS_m];
  pmod_mat_t B_inv_data[MEDS_s * MEDS_m * MEDS_m];

  pmod_mat_t *A_inv[MEDS_s];
  pmod_mat_t *B_inv[MEDS_s];

  for (int i = 0; i < MEDS_s; i++) {
    A_inv[i] = &A_inv_data[i * MEDS_m * MEDS_m];
    B_inv[i] = &B_inv_data[i * MEDS_n * MEDS_n];
  }

  long long gen_pm_sm_time_min = (1ll << 62);
  long long gen_G0p_time_min = (1ll << 62);
  long long solve_time_min = (1ll << 62);
  long long check_AB_inv_time_min = (1ll << 62);
  long long pi_time_min = (1ll << 62);
  long long syst_time_min = (1ll << 62);

  for (int i = 1; i < MEDS_s; i++) {
    START(gen_pm_sm_time);

    pmod_mat_t A[MEDS_m * MEDS_m] = {0};
    pmod_mat_t B[MEDS_n * MEDS_n] = {0};

    while (1 == 1)  // redo generation for this index until success
    {
      uint8_t sigma_Ti[MEDS_sec_seed_bytes];
      uint8_t sigma_a[MEDS_sec_seed_bytes];

      XOF((uint8_t *[]){sigma_a, sigma_Ti, sigma},
          (size_t[]){MEDS_sec_seed_bytes, MEDS_sec_seed_bytes,
                     MEDS_sec_seed_bytes},
          sigma, MEDS_sec_seed_bytes, 3);

      pmod_mat_t Ti[MEDS_k * MEDS_k];

      rnd_inv_matrix(Ti, MEDS_k, MEDS_k, sigma_Ti, MEDS_sec_seed_bytes);

      GFq_t Amm;

      {
        keccak_state Amm_shake;
        shake256_absorb_once(&Amm_shake, sigma_a, MEDS_sec_seed_bytes);

        Amm = rnd_GF(&Amm_shake);
      }

      LOG_MAT(Ti, MEDS_k, MEDS_k);
      LOG_VAL(Amm);

      pmod_mat_t G0prime[MEDS_k * MEDS_m * MEDS_n];

      START(gen_G0p_time);

      pmod_mat_mul(G0prime, MEDS_k, MEDS_m * MEDS_n, Ti, MEDS_k, MEDS_k, G[0],
                   MEDS_k, MEDS_m * MEDS_n);

      END(gen_G0p_time);
      MIN(gen_G0p_time_min, gen_G0p_time);

      LOG_MAT(G0prime, MEDS_k, MEDS_m * MEDS_n);

      START(solve_time);

      if (solve(A, B_inv[i], G0prime, Amm) < 0) {
        LOG("no sol");
        continue;
      }

      END(solve_time);
      MIN(solve_time_min, solve_time);

      START(check_AB_inv_time);

      if (pmod_mat_inv(B, B_inv[i], MEDS_n, MEDS_n) < 0) {
        LOG("no inv B");
        continue;
      }

      if (pmod_mat_inv(A_inv[i], A, MEDS_m, MEDS_m) < 0) {
        LOG("no inv A_inv");
        continue;
      }

      END(check_AB_inv_time);
      MIN(check_AB_inv_time_min, check_AB_inv_time);

      LOG_MAT_FMT(A, MEDS_m, MEDS_m, "A[%i]", i);
      LOG_MAT_FMT(A_inv[i], MEDS_m, MEDS_m, "A_inv[%i]", i);
      LOG_MAT_FMT(B, MEDS_n, MEDS_n, "B[%i]", i);
      LOG_MAT_FMT(B_inv[i], MEDS_n, MEDS_n, "B_inv[%i]", i);

      START(pi_time);

      pi(G[i], A, B, G[0]);
      // printf("%d\n", G[i][MEDS_n * MEDS_m]);

      END(pi_time);
      MIN(pi_time_min, pi_time);

      long syst_time = 0;
      syst_time = -cpucycles();

      if (pmod_mat_syst_ct(G[i], MEDS_k, MEDS_m * MEDS_n) != 0) {
        LOG("redo G[%i]", i);
        continue;  // Not systematic; try again for index i.
      }

      END(syst_time);
      MIN(syst_time_min, syst_time);

      LOG_MAT_FMT(G[i], MEDS_k, MEDS_m * MEDS_n, "G[%i]", i);

      // successfull generated G[s]; break out of while loop
      break;
    }

    END(gen_pm_sm_time);
    MIN(gen_pm_sm_time_min, gen_pm_sm_time);
  }

  END(gen_all_pm_sm_time);
  // LOG_M("  generate G, A, B (%d loops):\n", MEDS_s - 1);
  // LOG_M("    %f   (%llu cycles)\n", gen_all_pm_sm_time / freq,
  //       gen_all_pm_sm_time);

  // LOG_M("    max loop time:\n");
  // LOG_M("    loop time:\n");
  // LOG_M("      %f   (%llu cycles)\n", gen_pm_sm_time_min / freq,
  //       gen_pm_sm_time_min);
  //
  // LOG_M("    calculate TG:\n");
  // LOG_M("      %f   (%llu cycles)\n", gen_G0p_time_min / freq,
  //       gen_G0p_time_min);
  //
  // LOG_M("    solve A, B:\n");
  // LOG_M("      %f   (%llu cycles)\n", solve_time_min / freq, solve_time_min);
  //
  // LOG_M("    check A, B invertible:\n");
  // LOG_M("      %f   (%llu cycles)\n", check_AB_inv_time_min / freq,
  //       check_AB_inv_time_min);
  //
  // LOG_M("    calculate pi:\n");
  // LOG_M("      %f   (%llu cycles)\n", pi_time_min / freq, pi_time_min);
  //
  // LOG_M("    calculate syst form:\n");
  // LOG_M("      %f   (%llu cycles)\n", syst_time_min / freq, syst_time_min);

  /**
    Copy public data to bytes
    **/
  START(copy_pk_time);

  // copy pk data
  {
    uint8_t *tmp_pk = pk;

    memcpy(tmp_pk, sigma_G0, MEDS_pub_seed_bytes);
    LOG_VEC(tmp_pk, MEDS_pub_seed_bytes, "sigma_G0 (pk)");
    tmp_pk += MEDS_pub_seed_bytes;

    bitstream_t bs;

    bs_init(&bs, tmp_pk, MEDS_PK_BYTES - MEDS_pub_seed_bytes);

    for (int si = 1; si < MEDS_s; si++) {
      for (int j = (MEDS_m - 1) * MEDS_n; j < MEDS_m * MEDS_n; j++)
        bs_write(&bs, G[si][MEDS_m * MEDS_n + j], GFq_bits);

      for (int r = 2; r < MEDS_k; r++)
        for (int j = MEDS_k; j < MEDS_m * MEDS_n; j++)
          bs_write(&bs, G[si][r * MEDS_m * MEDS_n + j], GFq_bits);

      bs_finalize(&bs);
    }

    LOG_VEC(tmp_pk, MEDS_PK_BYTES - MEDS_pub_seed_bytes, "G[1:] (pk)");
    tmp_pk += MEDS_PK_BYTES - MEDS_pub_seed_bytes;

    LOG_HEX(pk, MEDS_PK_BYTES);

    if (MEDS_PK_BYTES !=
        MEDS_pub_seed_bytes + bs.byte_pos + (bs.bit_pos > 0 ? 1 : 0)) {
      fprintf(
          stderr,
          "ERROR: MEDS_PK_BYTES and actual pk size do not match! %i vs %i\n",
          MEDS_PK_BYTES,
          MEDS_pub_seed_bytes + bs.byte_pos + (bs.bit_pos > 0 ? 1 : 0));
      fprintf(stderr, "%i %i\n", MEDS_pub_seed_bytes + bs.byte_pos,
              MEDS_pub_seed_bytes + bs.byte_pos + (bs.bit_pos > 0 ? 1 : 0));
      return -1;
    }
  }

  // END(copy_pk_time);
  //  LOG_M("  copy public data:\n");
  //  LOG_M("    %f   (%llu cycles)\n", copy_pk_time / freq, copy_pk_time);

  /**
    Copy public data to bytes
    **/
  START(copy_sk_time);

  // copy sk data
  {
    memcpy(sk, delta, MEDS_sec_seed_bytes);
    memcpy(sk + MEDS_sec_seed_bytes, sigma_G0, MEDS_pub_seed_bytes);

    bitstream_t bs;

    bs_init(&bs, sk + MEDS_sec_seed_bytes + MEDS_pub_seed_bytes,
            MEDS_SK_BYTES - MEDS_sec_seed_bytes - MEDS_pub_seed_bytes);

    for (int si = 1; si < MEDS_s; si++) {
      for (int j = 0; j < MEDS_m * MEDS_m; j++)
        bs_write(&bs, A_inv[si][j], GFq_bits);

      bs_finalize(&bs);
    }

    for (int si = 1; si < MEDS_s; si++) {
      for (int j = 0; j < MEDS_n * MEDS_n; j++)
        bs_write(&bs, B_inv[si][j], GFq_bits);

      bs_finalize(&bs);
    }

    LOG_HEX(sk, MEDS_SK_BYTES);
  }

  // END(copy_sk_time);
  //  LOG_M("  copy secret data:\n");
  //  LOG_M("    %f   (%llu cycles)\n", copy_sk_time / freq, copy_sk_time);

  // LOG_M("\nKeypairs generated\n\n");

  close_measure_log();

  return 0;
}

void print_matrix(const char *message, pmod_mat_t *M, int M_r, int M_c) {
  printf(message);

  for (int i = 0; i < M_r; i++) {
    for (int j = 0; j < M_c; j++) {
      printf("%5d", M[i * M_c + j]);
    }
    printf("\n");
  }
  printf("\n");
}

int crypto_sign_vec(unsigned char *sm, unsigned long long *smlen,
                    const unsigned char *m, unsigned long long mlen,
                    const unsigned char *sk) {
  open_measure_log(measure_log_file);

  double freq = osfreq();
  LOG_M("Signing...\n\n");

  /**
    Load secred keypairs
    **/
  START(load_secred_time);

  uint8_t delta[MEDS_sec_seed_bytes];

  randombytes(delta, MEDS_sec_seed_bytes);

  // skip secret seed
  sk += MEDS_sec_seed_bytes;

  pmod_mat_t G_0[MEDS_k * MEDS_m * MEDS_n];

  rnd_sys_mat(G_0, MEDS_k, MEDS_m * MEDS_n, sk, MEDS_pub_seed_bytes);

  sk += MEDS_pub_seed_bytes;

  pmod_mat_t A_inv_data[MEDS_s * MEDS_m * MEDS_m];
  pmod_mat_t B_inv_data[MEDS_s * MEDS_n * MEDS_n];

  pmod_mat_t *A_inv[MEDS_s];
  pmod_mat_t *B_inv[MEDS_s];

  for (int i = 0; i < MEDS_s; i++) {
    A_inv[i] = &A_inv_data[i * MEDS_m * MEDS_m];
    B_inv[i] = &B_inv_data[i * MEDS_n * MEDS_n];
  }

  // Load secret key matrices.
  {
    bitstream_t bs;

    bs_init(&bs, (uint8_t *)sk,
            MEDS_SK_BYTES - MEDS_sec_seed_bytes - MEDS_pub_seed_bytes);

    for (int si = 1; si < MEDS_s; si++) {
      for (int j = 0; j < MEDS_m * MEDS_m; j++)
        A_inv[si][j] = bs_read(&bs, GFq_bits);

      bs_finalize(&bs);
    }

    for (int si = 1; si < MEDS_s; si++) {
      for (int j = 0; j < MEDS_n * MEDS_n; j++)
        B_inv[si][j] = bs_read(&bs, GFq_bits);

      bs_finalize(&bs);
    }

    bs_finalize(&bs);
  }

  for (int i = 1; i < MEDS_s; i++) LOG_MAT(A_inv[i], MEDS_m, MEDS_m);

  for (int i = 1; i < MEDS_s; i++) LOG_MAT(B_inv[i], MEDS_n, MEDS_n);

  LOG_MAT(G_0, MEDS_k, MEDS_m * MEDS_m);

  LOG_VEC(delta, MEDS_sec_seed_bytes);

  /**
    Generate seed tree parameter
    **/
  START(gen_st_param_time);

  uint8_t stree[MEDS_st_seed_bytes * SEED_TREE_size] = {0};
  uint8_t alpha[MEDS_st_salt_bytes];

  uint8_t *rho = &stree[MEDS_st_seed_bytes * SEED_TREE_ADDR(0, 0)];

  XOF((uint8_t *[]){rho, alpha},
      (size_t[]){MEDS_st_seed_bytes, MEDS_st_salt_bytes}, delta,
      MEDS_sec_seed_bytes, 2);

  t_hash(stree, alpha, 0, 0);

  uint8_t *sigma =
      &stree[MEDS_st_seed_bytes * SEED_TREE_ADDR(MEDS_seed_tree_height, 0)];

  /**
    Copy public data to bytes
    **/
  START(gen_rsp_time);

  pmod_mat_t A_tilde_data[MEDS_t][MEDS_m * MEDS_m];
  pmod_mat_t B_tilde_data[MEDS_t][MEDS_n * MEDS_n];
  pmod_mat_t A_tilde_dump[MEDS_m * MEDS_m];
  pmod_mat_t B_tilde_dump[MEDS_m * MEDS_m];

  pmod_mat_t *A_tilde[MEDS_t << 1];
  pmod_mat_t *B_tilde[MEDS_t << 1];

  for (int i = 0; i < MEDS_t; i++) {
    A_tilde[i] = A_tilde_data[i];
    B_tilde[i] = B_tilde_data[i];
  }
  for (int i = MEDS_t; i < MEDS_t << 1; i++) {
    A_tilde[i] = A_tilde_dump;
    B_tilde[i] = B_tilde_dump;
  }

  keccak_state h_shake;
  shake256_init(&h_shake);

  pmod_vec_t G_vec[MEDS_k * MEDS_m * MEDS_n];
  pmod_vec_t G0_vec[MEDS_k * MEDS_m * MEDS_n];

  for (int i = 0; i < MEDS_k * MEDS_m * MEDS_n; i++) G0_vec[i] = SET1(G_0[i]);

  pmod_vec_t A_vec[MEDS_m * MEDS_m];
  pmod_vec_t B_vec[MEDS_n * MEDS_n];

  static pmod_mat_t Gs_data[MEDS_t][MEDS_k * MEDS_m * MEDS_n];
  static pmod_mat_t Gs_dump[MEDS_k * MEDS_m * MEDS_n];
  static pmod_mat_t *Gs[MEDS_t << 1];

  for (int i = 0; i < MEDS_t; i++) Gs[i] = Gs_data[i];
  for (int i = MEDS_t; i < MEDS_t << 1; i++) Gs[i] = Gs_dump;

  int num_valid = 0;
  int num_invalid = 0;

  int num_tried = 0;
  int indexes[MEDS_t << 1];

  for (int i = 0; i < MEDS_t; i++) indexes[i] = i;
  for (int i = MEDS_t; i < (MEDS_t << 1); i++) indexes[i] = 0;

  long long pi_time_min = (1ll << 62);
  long long syst_time_min = (1ll << 62);

  while (num_valid < MEDS_t) {
    for (int t = 0; t < 16; t++) {
      if (t + num_tried >= MEDS_t + num_invalid) break;

      int i = indexes[t + num_tried];

      uint8_t sigma_A_tilde_i[MEDS_st_seed_bytes];
      uint8_t sigma_B_tilde_i[MEDS_st_seed_bytes];

      XOF((uint8_t *[]){sigma_A_tilde_i, sigma_B_tilde_i,
                        &sigma[i * MEDS_st_seed_bytes]},
          (size_t[]){MEDS_st_seed_bytes, MEDS_st_seed_bytes,
                     MEDS_st_seed_bytes},
          &sigma[i * MEDS_st_seed_bytes], MEDS_st_seed_bytes, 3);

      rnd_inv_matrix(A_tilde[i], MEDS_m, MEDS_m, sigma_A_tilde_i,
                     MEDS_st_seed_bytes);
      rnd_inv_matrix(B_tilde[i], MEDS_n, MEDS_n, sigma_B_tilde_i,
                     MEDS_st_seed_bytes);

      LOG_MAT_FMT(A_tilde[i], MEDS_m, MEDS_m, "A_tilde[%i]", i);
      LOG_MAT_FMT(B_tilde[i], MEDS_n, MEDS_n, "B_tilde[%i]", i);
    }

    for (int i = 0; i < MEDS_m * MEDS_m; i++)
      A_vec[i] =
          pmod_mat_entry_vec(A_tilde + num_tried, 1, MEDS_m * MEDS_m, 0, i);
    for (int i = 0; i < MEDS_n * MEDS_n; i++)
      B_vec[i] =
          pmod_mat_entry_vec(B_tilde + num_tried, 1, MEDS_n * MEDS_n, 0, i);

    START(pi_time);
    pi_vec(G_vec, A_vec, B_vec, G0_vec);
    END(pi_time);
    MIN(pi_time_min, pi_time);

    START(syst_simd_time);
    // pmod_vec_mask_t valid =
    //     pmod_mat_syst_ct_vec(G_vec, MEDS_k, MEDS_m * MEDS_n);
    pmod_vec_mask_t valid =
        pmod_mat_syst_ct_vec_inv(G_vec, MEDS_k, MEDS_m * MEDS_n);
    END(syst_simd_time);
    MIN(syst_time_min, syst_simd_time);

    for (int i = 0; i < MEDS_k * MEDS_m * MEDS_n; i++) {
      pmod_mat_set_entry_vec(Gs + num_tried, 1, MEDS_k * MEDS_m * MEDS_n, 0, i,
                             G_vec[i]);
    }

    int invalids = 0;
    int tries = 0;

    for (int t = 0; t < 16; t++) {
      if (t + num_tried >= MEDS_t + num_invalid) break;

      int i = indexes[t + num_tried];
      if (extract_mask(valid, t) == 0) {
        int idx = MEDS_t + num_invalid + invalids;

        indexes[idx] = i;
        Gs[idx] = Gs_data[i];
        A_tilde[idx] = A_tilde_data[i];
        B_tilde[idx] = B_tilde_data[i];

        invalids++;
      }

      tries++;
    }

    num_valid += tries - invalids;
    num_invalid += invalids;
    num_tried += tries;
  }

  // printf("Sqaure\n");
  // for (int i = 0; i < 10; i++) {
  //   for (int j = 0; j < 10; j++)
  //     printf("%6d", Gs_data[0][i * MEDS_m * MEDS_n + j]);
  //   printf("\n");
  // }

  // for (int i = 0; i < MEDS_t; i++) {
  //   printf("G_tilde[%d]\n", i);
  //   print_matrix("", Gs_data[i], MEDS_m, MEDS_m);
  // }

  /*
     Serialize results
     */

  for (int i = 0; i < MEDS_t; i++) {
    LOG_MAT_FMT(G_tilde_ti, MEDS_k, MEDS_m * MEDS_n, "G_tilde[%i]", i);

    bitstream_t bs;
    uint8_t
        bs_buf[CEILING((MEDS_k * (MEDS_m * MEDS_n - MEDS_k)) * GFq_bits, 8)];

    bs_init(&bs, bs_buf,
            CEILING((MEDS_k * (MEDS_m * MEDS_n - MEDS_k)) * GFq_bits, 8));

    for (int r = 0; r < MEDS_k; r++)
      for (int j = MEDS_k; j < MEDS_m * MEDS_n; j++)
        // bs_write(&bs, Gs[i][r * MEDS_m * MEDS_n + j], GFq_bits);
        bs_write(&bs, Gs_data[i][r * MEDS_m * MEDS_n + j], GFq_bits);

    shake256_absorb(
        &h_shake, bs_buf,
        CEILING((MEDS_k * (MEDS_m * MEDS_n - MEDS_k)) * GFq_bits, 8));
  }

  END(gen_rsp_time);
  LOG_M("  generate response:\n");
  LOG_M("    %f   (%llu cycles)\n", gen_rsp_time / freq, gen_rsp_time);

  LOG_M("    calculate pi (simd):\n");
  LOG_M("      %f   (%llu cycles)\n", pi_time_min / freq, pi_time_min);

  LOG_M("    calculate syst form (simd):\n");
  LOG_M("      %f   (%llu cycles)\n", syst_time_min / freq, syst_time_min);

  // LOG_M("\n");

  shake256_absorb(&h_shake, (uint8_t *)m, mlen);

  shake256_finalize(&h_shake);

  uint8_t digest[MEDS_digest_bytes];

  shake256_squeeze(digest, MEDS_digest_bytes, &h_shake);

  LOG_VEC(digest, MEDS_digest_bytes);

  uint8_t h[MEDS_t];

  parse_hash(digest, MEDS_digest_bytes, h, MEDS_t);

  LOG_VEC(h, MEDS_t);

  bitstream_t bs;

  bs_init(&bs, sm,
          MEDS_w * (CEILING(MEDS_m * MEDS_m * GFq_bits, 8) +
                    CEILING(MEDS_n * MEDS_n * GFq_bits, 8)));

  uint8_t *path = sm + MEDS_w * (CEILING(MEDS_m * MEDS_m * GFq_bits, 8) +
                                 CEILING(MEDS_n * MEDS_n * GFq_bits, 8));

  t_hash(stree, alpha, 0, 0);

  stree_to_path(stree, h, path, alpha);

  /**
    Compress mu and nu
     **/

  long long calc_mu_nu_time_min = (1llu << 62);

  START(calc_all_mu_nu_time);

  for (int i = 0; i < MEDS_t; i++) {
    if (h[i] > 0) {
      START(calc_mu_nu_time);

      {
        pmod_mat_t mu[MEDS_m * MEDS_m];

        // pmod_mat_mul(mu, MEDS_m, MEDS_m, A_tilde[i], MEDS_m, MEDS_m,
        //              A_inv[h[i]], MEDS_m, MEDS_m);
        pmod_mat_mul(mu, MEDS_m, MEDS_m, A_tilde_data[i], MEDS_m, MEDS_m,
                     A_inv[h[i]], MEDS_m, MEDS_m);

        LOG_MAT(mu, MEDS_m, MEDS_m);

        for (int j = 0; j < MEDS_m * MEDS_m; j++)
          bs_write(&bs, mu[j], GFq_bits);
      }

      bs_finalize(&bs);

      {
        pmod_mat_t nu[MEDS_n * MEDS_n];

        // pmod_mat_mul(nu, MEDS_n, MEDS_n, B_inv[h[i]], MEDS_n, MEDS_n,
        //              B_tilde[i], MEDS_n, MEDS_n);
        pmod_mat_mul(nu, MEDS_n, MEDS_n, B_inv[h[i]], MEDS_n, MEDS_n,
                     B_tilde_data[i], MEDS_n, MEDS_n);

        LOG_MAT(nu, MEDS_n, MEDS_n);

        for (int j = 0; j < MEDS_n * MEDS_n; j++)
          bs_write(&bs, nu[j], GFq_bits);
      }

      bs_finalize(&bs);

      END(calc_mu_nu_time);
      MIN(calc_mu_nu_time_min, calc_mu_nu_time);
    }
  }

  END(calc_all_mu_nu_time);
  // LOG_M("  calculate mu and nu (%d loops):\n", MEDS_t);
  // LOG_M("    %f   (%llu cycles)\n", calc_all_mu_nu_time / freq,
  //       calc_all_mu_nu_time);

  // LOG_M("    loop time:\n");
  // LOG_M("      %f   (%llu cycles)\n", calc_mu_nu_time_min / freq,
  //       calc_mu_nu_time_min);

  LOG_M("\nMessage signed\n\n");

  memcpy(sm + MEDS_SIG_BYTES - MEDS_digest_bytes - MEDS_st_salt_bytes, digest,
         MEDS_digest_bytes);
  memcpy(sm + MEDS_SIG_BYTES - MEDS_st_salt_bytes, alpha, MEDS_st_salt_bytes);
  memcpy(sm + MEDS_SIG_BYTES, m, mlen);

  *smlen = MEDS_SIG_BYTES + mlen;

  LOG_HEX(sm, MEDS_SIG_BYTES + mlen);

  fclose(measure_log);

  return 0;
}

int crypto_sign(unsigned char *sm, unsigned long long *smlen,
                const unsigned char *m, unsigned long long mlen,
                const unsigned char *sk) {
  open_measure_log(measure_log_file);

  double freq = osfreq();
  LOG_M("Signing...\n\n");

  /**
    Load secred keypairs
    **/
  START(load_secred_time);

  uint8_t delta[MEDS_sec_seed_bytes];

  randombytes(delta, MEDS_sec_seed_bytes);

  // skip secret seed
  sk += MEDS_sec_seed_bytes;

  pmod_mat_t G_0[MEDS_k * MEDS_m * MEDS_n];

  rnd_sys_mat(G_0, MEDS_k, MEDS_m * MEDS_n, sk, MEDS_pub_seed_bytes);

  sk += MEDS_pub_seed_bytes;

  pmod_mat_t A_inv_data[MEDS_s * MEDS_m * MEDS_m];
  pmod_mat_t B_inv_data[MEDS_s * MEDS_n * MEDS_n];

  pmod_mat_t *A_inv[MEDS_s];
  pmod_mat_t *B_inv[MEDS_s];

  for (int i = 0; i < MEDS_s; i++) {
    A_inv[i] = &A_inv_data[i * MEDS_m * MEDS_m];
    B_inv[i] = &B_inv_data[i * MEDS_n * MEDS_n];
  }

  // Load secret key matrices.
  {
    bitstream_t bs;

    bs_init(&bs, (uint8_t *)sk,
            MEDS_SK_BYTES - MEDS_sec_seed_bytes - MEDS_pub_seed_bytes);

    for (int si = 1; si < MEDS_s; si++) {
      for (int j = 0; j < MEDS_m * MEDS_m; j++)
        A_inv[si][j] = bs_read(&bs, GFq_bits);

      bs_finalize(&bs);
    }

    for (int si = 1; si < MEDS_s; si++) {
      for (int j = 0; j < MEDS_n * MEDS_n; j++)
        B_inv[si][j] = bs_read(&bs, GFq_bits);

      bs_finalize(&bs);
    }

    bs_finalize(&bs);
  }

  for (int i = 1; i < MEDS_s; i++) LOG_MAT(A_inv[i], MEDS_m, MEDS_m);

  for (int i = 1; i < MEDS_s; i++) LOG_MAT(B_inv[i], MEDS_n, MEDS_n);

  LOG_MAT(G_0, MEDS_k, MEDS_m * MEDS_m);

  LOG_VEC(delta, MEDS_sec_seed_bytes);

  // END(load_secred_time);
  //  LOG_M("  load secred keypairs:\n");
  //  LOG_M("    %f   (%llu cycles)\n", load_secred_time / freq,
  //  load_secred_time);

  /**
    Generate seed tree parameter
    **/
  START(gen_st_param_time);

  uint8_t stree[MEDS_st_seed_bytes * SEED_TREE_size] = {0};
  uint8_t alpha[MEDS_st_salt_bytes];

  uint8_t *rho = &stree[MEDS_st_seed_bytes * SEED_TREE_ADDR(0, 0)];

  XOF((uint8_t *[]){rho, alpha},
      (size_t[]){MEDS_st_seed_bytes, MEDS_st_salt_bytes}, delta,
      MEDS_sec_seed_bytes, 2);

  t_hash(stree, alpha, 0, 0);

  uint8_t *sigma =
      &stree[MEDS_st_seed_bytes * SEED_TREE_ADDR(MEDS_seed_tree_height, 0)];

  // END(gen_st_param_time);
  //  LOG_M("  generate ST parameters:\n");
  //  LOG_M("    %f   (%llu cycles)\n", gen_st_param_time / freq,
  //        gen_st_param_time);

  /**
    Copy public data to bytes
    **/
  START(gen_all_rsp_time);

  pmod_mat_t A_tilde_data[MEDS_t * MEDS_m * MEDS_m];
  pmod_mat_t B_tilde_data[MEDS_t * MEDS_m * MEDS_m];

  pmod_mat_t *A_tilde[MEDS_t];
  pmod_mat_t *B_tilde[MEDS_t];

  for (int i = 0; i < MEDS_t; i++) {
    A_tilde[i] = &A_tilde_data[i * MEDS_m * MEDS_m];
    B_tilde[i] = &B_tilde_data[i * MEDS_n * MEDS_n];
  }

  keccak_state h_shake;
  shake256_init(&h_shake);

  // print_matrix("G_0\n", G_0, MEDS_m, MEDS_n);
  // print_matrix("G_0\n", G_0, 2, MEDS_n);

  long long gen_rsp_time_min = (1ll << 62);
  long long pi_time_min = (1ll << 62);
  long long syst_time_min = (1ll << 62);

  for (int i = 0; i < MEDS_t; i++) {
    long gen_rsp_time = 0;
    gen_rsp_time = -cpucycles();

    pmod_mat_t G_tilde_ti[MEDS_k * MEDS_m * MEDS_m];

    while (1 == 1) {
      uint8_t sigma_A_tilde_i[MEDS_st_seed_bytes];
      uint8_t sigma_B_tilde_i[MEDS_st_seed_bytes];

      XOF((uint8_t *[]){sigma_A_tilde_i, sigma_B_tilde_i,
                        &sigma[i * MEDS_st_seed_bytes]},
          (size_t[]){MEDS_st_seed_bytes, MEDS_st_seed_bytes,
                     MEDS_st_seed_bytes},
          &sigma[i * MEDS_st_seed_bytes], MEDS_st_seed_bytes, 3);

      rnd_inv_matrix(A_tilde[i], MEDS_m, MEDS_m, sigma_A_tilde_i,
                     MEDS_st_seed_bytes);
      rnd_inv_matrix(B_tilde[i], MEDS_n, MEDS_n, sigma_B_tilde_i,
                     MEDS_st_seed_bytes);

      LOG_MAT_FMT(A_tilde[i], MEDS_m, MEDS_m, "A_tilde[%i]", i);
      LOG_MAT_FMT(B_tilde[i], MEDS_n, MEDS_n, "B_tilde[%i]", i);

      START(pi_time);

      pi(G_tilde_ti, A_tilde[i], B_tilde[i], G_0);

      END(pi_time);
      MIN(pi_time_min, pi_time);

      LOG_MAT_FMT(G_tilde_ti, MEDS_k, MEDS_m * MEDS_n, "G_tilde[%i]", i);

      START(syst_time);

      if (pmod_mat_syst_ct(G_tilde_ti, MEDS_k, MEDS_m * MEDS_n) != 0) continue;
      // if (16 <= i && i < 20) {
      //   printf("%d\n", i);
      //   print_matrix("G_tilde_0\n", G_tilde_ti, 2, MEDS_n);
      // }

      END(syst_time);
      MIN(syst_time_min, syst_time);

      break;
    }

    // memcpy(Gs[i], (const void*)G_tilde_ti, sizeof(pmod_mat_t) * MEDS_k *
    // MEDS_m * MEDS_n);

    LOG_MAT_FMT(G_tilde_ti, MEDS_k, MEDS_m * MEDS_n, "G_tilde[%i]", i);

    bitstream_t bs;
    uint8_t
        bs_buf[CEILING((MEDS_k * (MEDS_m * MEDS_n - MEDS_k)) * GFq_bits, 8)];

    bs_init(&bs, bs_buf,
            CEILING((MEDS_k * (MEDS_m * MEDS_n - MEDS_k)) * GFq_bits, 8));

    for (int r = 0; r < MEDS_k; r++)
      for (int j = MEDS_k; j < MEDS_m * MEDS_n; j++)
        bs_write(&bs, G_tilde_ti[r * MEDS_m * MEDS_n + j], GFq_bits);

    shake256_absorb(
        &h_shake, bs_buf,
        CEILING((MEDS_k * (MEDS_m * MEDS_n - MEDS_k)) * GFq_bits, 8));

    END(gen_rsp_time);
    MIN(gen_rsp_time_min, gen_rsp_time);

    // printf("G_tilde[%d]\n", i);
    // print_matrix("", G_tilde_ti, MEDS_m, MEDS_n);

    // break;
  }

  // for (int i = 0; i < MEDS_t; i++) {
  //   printf("G_tilde[%d]\n", i);
  //   print_matrix("", Gs_data[i], MEDS_m, MEDS_m);
  //   // print_matrix("", Gs_data[i], MEDS_k, MEDS_m * MEDS_n);
  // }

  END(gen_all_rsp_time);
  LOG_M("  generate response (%d loops):\n", MEDS_t);
  LOG_M("    %f   (%llu cycles)\n", gen_all_rsp_time / freq, gen_all_rsp_time);

  // LOG_M("    loop time:\n");
  // LOG_M("      %f   (%llu cycles)\n", gen_rsp_time_min / freq,
  //       gen_rsp_time_min);

  LOG_M("    calculate pi:\n");
  LOG_M("      %f   (%llu cycles)\n", pi_time_min / freq, pi_time_min);

  LOG_M("    calculate syst form:\n");
  LOG_M("      %f   (%llu cycles)\n", syst_time_min / freq, syst_time_min);

  LOG_M("\n");

  // LOG_M("  COMP result 0:\n");
  // for (int i = 0; i < 5; i++) {
  //   LOG_M("    ");
  //   for (int j = 0; j < 12; j++) {
  //     LOG_M("%d ", G_tilde_t0[i * MEDS_m * MEDS_n + j]);
  //   }
  //   LOG_M("\n");
  // }

  shake256_absorb(&h_shake, (uint8_t *)m, mlen);

  shake256_finalize(&h_shake);

  uint8_t digest[MEDS_digest_bytes];

  shake256_squeeze(digest, MEDS_digest_bytes, &h_shake);

  LOG_VEC(digest, MEDS_digest_bytes);

  uint8_t h[MEDS_t];

  parse_hash(digest, MEDS_digest_bytes, h, MEDS_t);

  LOG_VEC(h, MEDS_t);

  bitstream_t bs;

  bs_init(&bs, sm,
          MEDS_w * (CEILING(MEDS_m * MEDS_m * GFq_bits, 8) +
                    CEILING(MEDS_n * MEDS_n * GFq_bits, 8)));

  uint8_t *path = sm + MEDS_w * (CEILING(MEDS_m * MEDS_m * GFq_bits, 8) +
                                 CEILING(MEDS_n * MEDS_n * GFq_bits, 8));

  t_hash(stree, alpha, 0, 0);

  stree_to_path(stree, h, path, alpha);

  /**
    Compress mu and nu
     **/

  long long calc_mu_nu_time_min = (1llu << 62);

  START(calc_all_mu_nu_time);

  for (int i = 0; i < MEDS_t; i++) {
    if (h[i] > 0) {
      START(calc_mu_nu_time);

      {
        pmod_mat_t mu[MEDS_m * MEDS_m];

        pmod_mat_mul(mu, MEDS_m, MEDS_m, A_tilde[i], MEDS_m, MEDS_m,
                     A_inv[h[i]], MEDS_m, MEDS_m);

        LOG_MAT(mu, MEDS_m, MEDS_m);

        for (int j = 0; j < MEDS_m * MEDS_m; j++)
          bs_write(&bs, mu[j], GFq_bits);
      }

      bs_finalize(&bs);

      {
        pmod_mat_t nu[MEDS_n * MEDS_n];

        pmod_mat_mul(nu, MEDS_n, MEDS_n, B_inv[h[i]], MEDS_n, MEDS_n,
                     B_tilde[i], MEDS_n, MEDS_n);

        LOG_MAT(nu, MEDS_n, MEDS_n);

        for (int j = 0; j < MEDS_n * MEDS_n; j++)
          bs_write(&bs, nu[j], GFq_bits);
      }

      bs_finalize(&bs);

      END(calc_mu_nu_time);
      MIN(calc_mu_nu_time_min, calc_mu_nu_time);
    }
  }

  END(calc_all_mu_nu_time);
  // LOG_M("  calculate mu and nu (%d loops):\n", MEDS_t);
  // LOG_M("    %f   (%llu cycles)\n", calc_all_mu_nu_time / freq,
  //       calc_all_mu_nu_time);

  // LOG_M("    loop time:\n");
  // LOG_M("      %f   (%llu cycles)\n", calc_mu_nu_time_min / freq,
  //       calc_mu_nu_time_min);

  LOG_M("\nMessage signed\n\n");

  memcpy(sm + MEDS_SIG_BYTES - MEDS_digest_bytes - MEDS_st_salt_bytes, digest,
         MEDS_digest_bytes);
  memcpy(sm + MEDS_SIG_BYTES - MEDS_st_salt_bytes, alpha, MEDS_st_salt_bytes);
  memcpy(sm + MEDS_SIG_BYTES, m, mlen);

  *smlen = MEDS_SIG_BYTES + mlen;

  LOG_HEX(sm, MEDS_SIG_BYTES + mlen);

  fclose(measure_log);

  return 0;
}

int crypto_sign_open(unsigned char *m, unsigned long long *mlen,
                     const unsigned char *sm, unsigned long long smlen,
                     const unsigned char *pk) {
  LOG_HEX(sm, smlen);

  pmod_mat_t G_data[MEDS_k * MEDS_m * MEDS_n * MEDS_s];
  pmod_mat_t *G[MEDS_s];

  for (int i = 0; i < MEDS_s; i++) G[i] = &G_data[i * MEDS_k * MEDS_m * MEDS_n];

  rnd_sys_mat(G[0], MEDS_k, MEDS_m * MEDS_n, pk, MEDS_pub_seed_bytes);

  {
    bitstream_t bs;

    bs_init(&bs, (uint8_t *)pk + MEDS_pub_seed_bytes,
            MEDS_PK_BYTES - MEDS_pub_seed_bytes);

    for (int i = 1; i < MEDS_s; i++) {
      for (int r = 0; r < MEDS_k; r++)
        for (int c = 0; c < MEDS_k; c++)
          if (r == c)
            pmod_mat_set_entry(G[i], MEDS_k, MEDS_m * MEDS_n, r, c, 1);
          else
            pmod_mat_set_entry(G[i], MEDS_k, MEDS_m * MEDS_n, r, c, 0);

      for (int j = (MEDS_m - 1) * MEDS_n; j < MEDS_m * MEDS_n; j++)
        G[i][MEDS_m * MEDS_n + j] = bs_read(&bs, GFq_bits);

      for (int r = 2; r < MEDS_k; r++)
        for (int j = MEDS_k; j < MEDS_m * MEDS_n; j++)
          G[i][r * MEDS_m * MEDS_n + j] = bs_read(&bs, GFq_bits);

      for (int ii = 0; ii < MEDS_m; ii++)
        for (int j = 0; j < MEDS_n; j++)
          G[i][ii * MEDS_n + j] = ii == j ? 1 : 0;

      for (int ii = 0; ii < MEDS_m - 1; ii++)
        for (int j = 0; j < MEDS_n; j++)
          G[i][MEDS_m * MEDS_n + ii * MEDS_n + j] = (ii + 1) == j ? 1 : 0;

      bs_finalize(&bs);
    }
  }

  for (int i = 0; i < MEDS_s; i++)
    LOG_MAT_FMT(G[i], MEDS_k, MEDS_m * MEDS_n, "G[%i]", i);

  uint8_t *digest =
      (uint8_t *)sm + (MEDS_SIG_BYTES - MEDS_digest_bytes - MEDS_st_salt_bytes);

  uint8_t *alpha = (uint8_t *)sm + (MEDS_SIG_BYTES - MEDS_st_salt_bytes);

  LOG_VEC(digest, MEDS_digest_bytes);

  uint8_t h[MEDS_t];

  parse_hash(digest, MEDS_digest_bytes, h, MEDS_t);

  bitstream_t bs;

  bs_init(&bs, (uint8_t *)sm,
          MEDS_w * (CEILING(MEDS_m * MEDS_m * GFq_bits, 8) +
                    CEILING(MEDS_n * MEDS_n * GFq_bits, 8)));

  uint8_t *path =
      (uint8_t *)sm + MEDS_w * (CEILING(MEDS_m * MEDS_m * GFq_bits, 8) +
                                CEILING(MEDS_n * MEDS_n * GFq_bits, 8));

  uint8_t stree[MEDS_st_seed_bytes * SEED_TREE_size] = {0};

  path_to_stree(stree, h, path, alpha);

  uint8_t *sigma =
      &stree[MEDS_st_seed_bytes * SEED_TREE_ADDR(MEDS_seed_tree_height, 0)];

  pmod_mat_t G_hat_i[MEDS_k * MEDS_m * MEDS_n];

  pmod_mat_t mu[MEDS_m * MEDS_m];
  pmod_mat_t nu[MEDS_n * MEDS_n];

  keccak_state shake;
  shake256_init(&shake);

  for (int i = 0; i < MEDS_t; i++) {
    if (h[i] > 0) {
      for (int j = 0; j < MEDS_m * MEDS_m; j++) mu[j] = bs_read(&bs, GFq_bits);

      bs_finalize(&bs);

      for (int j = 0; j < MEDS_n * MEDS_n; j++) nu[j] = bs_read(&bs, GFq_bits);

      bs_finalize(&bs);

      LOG_MAT_FMT(mu, MEDS_m, MEDS_m, "mu[%i]", i);
      LOG_MAT_FMT(nu, MEDS_n, MEDS_n, "nu[%i]", i);

      pi(G_hat_i, mu, nu, G[h[i]]);

      LOG_MAT_FMT(G_hat_i, MEDS_k, MEDS_m * MEDS_n, "G_hat[%i]", i);

      pmod_mat_syst_ct(G_hat_i, MEDS_k, MEDS_m * MEDS_n);

      LOG_MAT_FMT(G_hat_i, MEDS_k, MEDS_m * MEDS_n, "G_hat[%i]", i);
    } else {
      while (1 == 1) {
        LOG_VEC_FMT(&sigma[i * MEDS_st_seed_bytes], MEDS_st_seed_bytes,
                    "seeds[%i]", i);

        uint8_t sigma_A_tilde_i[MEDS_st_seed_bytes];
        uint8_t sigma_B_tilde_i[MEDS_st_seed_bytes];

        XOF((uint8_t *[]){sigma_A_tilde_i, sigma_B_tilde_i,
                          &sigma[i * MEDS_st_seed_bytes]},
            (size_t[]){MEDS_st_seed_bytes, MEDS_st_seed_bytes,
                       MEDS_st_seed_bytes},
            &sigma[i * MEDS_st_seed_bytes], MEDS_st_seed_bytes, 3);

        pmod_mat_t A_tilde[MEDS_m * MEDS_m];
        pmod_mat_t B_tilde[MEDS_n * MEDS_n];

        LOG_VEC(sigma_A_tilde_i, MEDS_sec_seed_bytes);
        rnd_inv_matrix(A_tilde, MEDS_m, MEDS_m, sigma_A_tilde_i,
                       MEDS_st_seed_bytes);

        LOG_VEC(sigma_B_tilde_i, MEDS_sec_seed_bytes);
        rnd_inv_matrix(B_tilde, MEDS_n, MEDS_n, sigma_B_tilde_i,
                       MEDS_st_seed_bytes);

        LOG_MAT_FMT(A_tilde, MEDS_m, MEDS_m, "A_tilde[%i]", i);
        LOG_MAT_FMT(B_tilde, MEDS_n, MEDS_n, "B_tilde[%i]", i);

        pi(G_hat_i, A_tilde, B_tilde, G[0]);

        LOG_MAT_FMT(G_hat_i, MEDS_k, MEDS_m * MEDS_n, "G_hat[%i]", i);

        if (pmod_mat_syst_ct(G_hat_i, MEDS_k, MEDS_m * MEDS_n) == 0) {
          LOG_MAT_FMT(G_hat_i, MEDS_k, MEDS_m * MEDS_n, "G_hat[%i]", i);
          break;
        }

        LOG_MAT_FMT(G_hat_i, MEDS_k, MEDS_m * MEDS_n, "G_hat[%i]", i);
      }
    }

    bitstream_t bs;
    uint8_t
        bs_buf[CEILING((MEDS_k * (MEDS_m * MEDS_n - MEDS_k)) * GFq_bits, 8)];

    bs_init(&bs, bs_buf,
            CEILING((MEDS_k * (MEDS_m * MEDS_n - MEDS_k)) * GFq_bits, 8));

    for (int r = 0; r < MEDS_k; r++)
      for (int j = MEDS_k; j < MEDS_m * MEDS_n; j++)
        bs_write(&bs, G_hat_i[r * MEDS_m * MEDS_n + j], GFq_bits);

    shake256_absorb(
        &shake, bs_buf,
        CEILING((MEDS_k * (MEDS_m * MEDS_n - MEDS_k)) * GFq_bits, 8));
  }

  shake256_absorb(&shake, (uint8_t *)(sm + MEDS_SIG_BYTES),
                  smlen - MEDS_SIG_BYTES);

  shake256_finalize(&shake);

  uint8_t check[MEDS_digest_bytes];

  shake256_squeeze(check, MEDS_digest_bytes, &shake);

  if (memcmp(digest, check, MEDS_digest_bytes) != 0) {
    fprintf(stderr, "Signature verification failed!\n");

    return -1;
  }

  memcpy(m, (uint8_t *)(sm + MEDS_SIG_BYTES), smlen - MEDS_SIG_BYTES);
  *mlen = smlen - MEDS_SIG_BYTES;

  return 0;
}
