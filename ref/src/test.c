#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "randombytes.h"

#include "api.h"
#include "meds.h"
#include "params.h"

#include "measure.h"

int main(int argc, char *argv[]) {
  printf("paramter set: %s\n\n", MEDS_name);

  long long time = 0;
  long long keygen_time = 0xfffffffffffffff;
  long long sign_time = 0xfffffffffffffff;
  long long verify_time = 0xfffffffffffffff;

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

  for (int round = 0; round < rounds; round++) {
    uint8_t sk[CRYPTO_SECRETKEYBYTES] = {0};
    uint8_t pk[CRYPTO_PUBLICKEYBYTES] = {0};

    time = -cpucycles();
    crypto_sign_keypair(pk, sk);
    time += cpucycles();

    if (time < keygen_time)
      keygen_time = time;

    uint8_t sig[CRYPTO_BYTES + sizeof(msg)] = {0};
    unsigned long long sig_len = sizeof(sig);

    time = -cpucycles();
    crypto_sign(sig, &sig_len, (const unsigned char *)msg, sizeof(msg), sk);
    time += cpucycles();
    printf("sign (normal): %llu\n", time);

    time = -cpucycles();
    crypto_sign_vec(sig, &sig_len, (const unsigned char *)msg, sizeof(msg), sk);
    time += cpucycles();
    printf("sign   (SIMD): %llu\n", time);

    if (time < sign_time)
      sign_time = time;

    unsigned char msg_out[4];
    unsigned long long msg_out_len = sizeof(msg_out);

    int ret;

    time = -cpucycles();
    ret = crypto_sign_open(msg_out, &msg_out_len, sig, sizeof(sig), pk);
    time += cpucycles();
    printf("verify (normal): %llu\n", time);

    time = -cpucycles();
    ret = crypto_sign_open_vec(msg_out, &msg_out_len, sig, sizeof(sig), pk);
    time += cpucycles();
    printf("verify   (SIMD): %llu\n", time);

    if (time < verify_time)
      verify_time = time;

    if (ret == 0)
      printf("success\n");
    else
      printf("!!! FAILED !!!\n");
  }

  double freq = osfreq();

  printf("\n");
  printf("Time (min of %i runs):\n", rounds);
  printf("keygen: %f   (%llu cycles)\n", keygen_time / freq, keygen_time);
  printf("sign:   %f   (%llu cycles)\n", sign_time / freq, sign_time);
  printf("verify: %f   (%llu cycles)\n", verify_time / freq, verify_time);

  return 0;
}
