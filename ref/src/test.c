#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "api.h"
#include "measure.h"
#include "meds.h"
#include "params.h"
#include "randombytes.h"

#define COMPARE

#ifdef LOG_MEASURE
const char *measure_log_file = "./measure.txt";
#endif

int main(int argc, char *argv[]) {
#ifdef LOG_MEASURE
  open_measure_log(measure_log_file);
#endif

  printf("paramter set: %s\n\n", MEDS_name);

  long long keygen_time_simd = 0xfffffffffffffff;
  long long sign_time_simd = 0xfffffffffffffff;
  long long verify_time_simd = 0xfffffffffffffff;

  int rounds = 1;

  if (argc > 1)
    rounds = atoi(argv[1]);

  unsigned char entropy_input[48] = {0};

  randombytes_init(entropy_input, NULL, 256);

  char msg[4] = "Test";

  printf("pk:  %i bytes\n", CRYPTO_PUBLICKEYBYTES);
  printf("sk:  %i bytes\n", CRYPTO_SECRETKEYBYTES);
  printf("sig: %i bytes\n", CRYPTO_BYTES);
  printf("\n");

  uint8_t sk[CRYPTO_SECRETKEYBYTES] = {0};
  uint8_t pk[CRYPTO_PUBLICKEYBYTES] = {0};

  uint8_t sig[CRYPTO_BYTES + sizeof(msg)] = {0};
  unsigned long long sig_len = sizeof(sig);

  unsigned char msg_out[4];
  unsigned long long msg_out_len = sizeof(msg_out);

  for (int round = 0; round < rounds; round++) {
    keygen_time_simd = -cpucycles();
    crypto_sign_keypair_vec(pk, sk);
    keygen_time_simd += cpucycles();

    sign_time_simd = -cpucycles();
    crypto_sign_vec(sig, &sig_len, (const unsigned char *)msg, sizeof(msg), sk);
    sign_time_simd += cpucycles();

    int ret;

    verify_time_simd = -cpucycles();
    ret = crypto_sign_open_vec(msg_out, &msg_out_len, sig, sizeof(sig), pk);
    verify_time_simd += cpucycles();

    if (ret == 0)
      printf("success   (SIMD)\n");
    else
      printf("!!! FAILED   (SIMD) !!!\n");
    printf("\n");

#ifdef COMPARE
    long long keygen_time = 0xfffffffffffffff;
    long long sign_time = 0xfffffffffffffff;
    long long verify_time = 0xfffffffffffffff;

    keygen_time = -cpucycles();
    crypto_sign_keypair(pk, sk);
    keygen_time += cpucycles();

    sign_time = -cpucycles();
    crypto_sign(sig, &sig_len, (const unsigned char *)msg, sizeof(msg), sk);
    // crypto_sign_vec(sig, &sig_len, (const unsigned char *)msg, sizeof(msg), sk);
    sign_time += cpucycles();

    verify_time = -cpucycles();
    ret = crypto_sign_open(msg_out, &msg_out_len, sig, sizeof(sig), pk);
    // ret = crypto_sign_open_vec(msg_out, &msg_out_len, sig, sizeof(sig), pk);
    verify_time += cpucycles();

    printf("keypair (normal): %llu\n", keygen_time);
    printf("keypair   (SIMD): %llu\n", keygen_time_simd);
    printf("\n");

    printf("   sign (normal): %llu\n", sign_time);
    printf("   sign   (SIMD): %llu\n", sign_time_simd);
    printf("\n");

    printf(" verify (normal): %llu\n", verify_time);
    printf(" verify   (SIMD): %llu\n", verify_time_simd);
    printf("\n");

    if (ret == 0)
      printf("success (normal)\n");
    else
      printf("!!! FAILED (normal) !!!\n");
    printf("\n");
#endif
  }

  double freq = osfreq();

  printf("Time (min of %i runs):\n", rounds);
  printf("keygen: %f   (%llu cycles)\n", keygen_time_simd / freq, keygen_time_simd);
  printf("sign:   %f   (%llu cycles)\n", sign_time_simd / freq, sign_time_simd);
  printf("verify: %f   (%llu cycles)\n", verify_time_simd / freq, verify_time_simd);

#ifdef LOG_MEASURE
  close_measure_log();
#endif

  return 0;
}
