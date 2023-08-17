#ifndef FIPS202_W64_H
#define FIPS202_W64_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "vec_w64.h"

#define SHAKE128_RATE 168
#define SHAKE256_RATE 136
#define SHA3_256_RATE 136
#define SHA3_512_RATE 72

#define FIPS202_NAMESPACE(s) meds_fips202_ref_##s

typedef struct {
  pmod_mat_w64_t s[25];
  unsigned int pos;
} keccak_state_w64;

#define shake256_w64_init FIPS202_NAMESPACE(shake256_w64_init)
void shake256_w64_init(keccak_state_w64 *state);

#define shake256_w64_absorb FIPS202_NAMESPACE(shake256_w64_absorb)
void shake256_w64_absorb(keccak_state_w64 *state, const pmod_mat_w64_t *in, size_t inlen);

#define shake256_w64_absorb_raw FIPS202_NAMESPACE(shake256_w64_absorb_test)
void shake256_w64_absorb_raw(keccak_state_w64 *state, const uint8_t **in, size_t inlen, uint32_t offset);

#define shake256_w64_finalize FIPS202_NAMESPACE(shake256_w64_finalize)
void shake256_w64_finalize(keccak_state_w64 *state);

#define shake256_w64_squeeze FIPS202_NAMESPACE(shake256_w64_squeeze)
void shake256_w64_squeeze(pmod_mat_w64_t *out, size_t outlen, keccak_state_w64 *state);

#define shake256_w64_absorb_once FIPS202_NAMESPACE(shake256_w64_absorb_once)
void shake256_w64_absorb_once(keccak_state_w64 *state, const pmod_mat_w64_t *in, size_t inlen);

#define shake256_w64_squeezeblocks FIPS202_NAMESPACE(shake256_w64_squeezeblocks)
void shake256_w64_squeezeblocks(pmod_mat_w64_t *out, size_t nblocks,  keccak_state_w64 *state);

#define shake256_w64 FIPS202_NAMESPACE(shake256_w64)
void shake256_w64(pmod_mat_w64_t *out, size_t outlen, const pmod_mat_w64_t *in, size_t inlen);

#endif

