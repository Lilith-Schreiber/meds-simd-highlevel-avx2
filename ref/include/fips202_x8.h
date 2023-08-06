#ifndef FIPS202_x8_H
#define FIPS202_x8_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

// #include "fips202.h"
#include "vec_x8.h"

#define SHAKE128_RATE 168
#define SHAKE256_RATE 136
#define SHA3_256_RATE 136
#define SHA3_512_RATE 72

#define FIPS202_NAMESPACE(s) meds_fips202_ref_##s

typedef struct {
  pmod_vec_x8_t s[25];
  unsigned int pos;
} keccak_state_x8;

#define shake256_x8_init FIPS202_NAMESPACE(shake256_x8_init)
void shake256_x8_init(keccak_state_x8 *state);

#define shake256_x8_absorb FIPS202_NAMESPACE(shake256_x8_absorb)
void shake256_x8_absorb(keccak_state_x8 *state, const pmod_vec_x8_t *in, size_t inlen);

#define shake256_x8_absorb_raw FIPS202_NAMESPACE(shake256_x8_absorb_test)
void shake256_x8_absorb_raw(keccak_state_x8 *state, const uint8_t **in, size_t inlen, uint32_t offset);

#define shake256_x8_finalize FIPS202_NAMESPACE(shake256_x8_finalize)
void shake256_x8_finalize(keccak_state_x8 *state);

#define shake256_x8_squeeze FIPS202_NAMESPACE(shake256_x8_squeeze)
void shake256_x8_squeeze(pmod_vec_x8_t *out, size_t outlen, keccak_state_x8 *state);

#define shake256_x8_absorb_once FIPS202_NAMESPACE(shake256_x8_absorb_once)
void shake256_x8_absorb_once(keccak_state_x8 *state, const pmod_vec_x8_t *in, size_t inlen);

#define shake256_x8_squeezeblocks FIPS202_NAMESPACE(shake256_x8_squeezeblocks)
void shake256_x8_squeezeblocks(pmod_vec_x8_t *out, size_t nblocks,  keccak_state_x8 *state);

#define shake256_x8 FIPS202_NAMESPACE(shake256_x8)
void shake256_x8(pmod_vec_x8_t *out, size_t outlen, const pmod_vec_x8_t *in, size_t inlen);

#endif

