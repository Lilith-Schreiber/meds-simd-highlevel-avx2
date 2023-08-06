#include "fips202_x8.h"

#include <stddef.h>
#include <stdint.h>

#include "vec_x8.h"

#define NROUNDS 24
#define ROL_x8(a, offset) XOR_x8(SLLI_x8(a, offset), SRLI_x8(a, (64 - offset)))

static pmod_vec_x8_t load64_x8(const pmod_vec_x8_t x[8]) {
  unsigned int i;
  pmod_vec_x8_t r = SET1_x8(0);

  for (i = 0; i < 8; i++) r = OR_x8(r, SLLI_x8(x[i], 8 * i));

  return r;
}

static void store64_x8(pmod_vec_x8_t x[8], pmod_vec_x8_t u) {
  unsigned int i;

  for (i = 0; i < 8; i++) x[i] = SRLI_x8(u, 8 * i);
}

/* Keccak round constants */
static const pmod_vec_x8_t KeccakF_RoundConstants_x8[NROUNDS] = {
    SET1_CT_x8(0x0000000000000001ULL), SET1_CT_x8(0x0000000000008082ULL),
    SET1_CT_x8(0x800000000000808aULL), SET1_CT_x8(0x8000000080008000ULL),
    SET1_CT_x8(0x000000000000808bULL), SET1_CT_x8(0x0000000080000001ULL),
    SET1_CT_x8(0x8000000080008081ULL), SET1_CT_x8(0x8000000000008009ULL),
    SET1_CT_x8(0x000000000000008aULL), SET1_CT_x8(0x0000000000000088ULL),
    SET1_CT_x8(0x0000000080008009ULL), SET1_CT_x8(0x000000008000000aULL),
    SET1_CT_x8(0x000000008000808bULL), SET1_CT_x8(0x800000000000008bULL),
    SET1_CT_x8(0x8000000000008089ULL), SET1_CT_x8(0x8000000000008003ULL),
    SET1_CT_x8(0x8000000000008002ULL), SET1_CT_x8(0x8000000000000080ULL),
    SET1_CT_x8(0x000000000000800aULL), SET1_CT_x8(0x800000008000000aULL),
    SET1_CT_x8(0x8000000080008081ULL), SET1_CT_x8(0x8000000000008080ULL),
    SET1_CT_x8(0x0000000080000001ULL), SET1_CT_x8(0x8000000080008008ULL)};

static void KeccakF1600_StatePermute_x8(pmod_vec_x8_t state[25]) {
  int round;

  static const pmod_vec_x8_t ALL_ONE = SET1_CT_x8(0xFFFFFFFFFFFFFFFFULL);

  pmod_vec_x8_t Aba, Abe, Abi, Abo, Abu;
  pmod_vec_x8_t Aga, Age, Agi, Ago, Agu;
  pmod_vec_x8_t Aka, Ake, Aki, Ako, Aku;
  pmod_vec_x8_t Ama, Ame, Ami, Amo, Amu;
  pmod_vec_x8_t Asa, Ase, Asi, Aso, Asu;
  pmod_vec_x8_t BCa, BCe, BCi, BCo, BCu;
  pmod_vec_x8_t Da, De, Di, Do, Du;
  pmod_vec_x8_t Eba, Ebe, Ebi, Ebo, Ebu;
  pmod_vec_x8_t Ega, Ege, Egi, Ego, Egu;
  pmod_vec_x8_t Eka, Eke, Eki, Eko, Eku;
  pmod_vec_x8_t Ema, Eme, Emi, Emo, Emu;
  pmod_vec_x8_t Esa, Ese, Esi, Eso, Esu;

  // copyFromState(A, state)
  Aba = state[0];
  Abe = state[1];
  Abi = state[2];
  Abo = state[3];
  Abu = state[4];
  Aga = state[5];
  Age = state[6];
  Agi = state[7];
  Ago = state[8];
  Agu = state[9];
  Aka = state[10];
  Ake = state[11];
  Aki = state[12];
  Ako = state[13];
  Aku = state[14];
  Ama = state[15];
  Ame = state[16];
  Ami = state[17];
  Amo = state[18];
  Amu = state[19];
  Asa = state[20];
  Ase = state[21];
  Asi = state[22];
  Aso = state[23];
  Asu = state[24];

  for (round = 0; round < NROUNDS; round += 2) {
    //    prepareTheta
    BCa = XOR_x8(Aba, XOR_x8(Aga, XOR_x8(Aka, XOR_x8(Ama, Asa))));
    BCe = XOR_x8(Abe, XOR_x8(Age, XOR_x8(Ake, XOR_x8(Ame, Ase))));
    BCi = XOR_x8(Abi, XOR_x8(Agi, XOR_x8(Aki, XOR_x8(Ami, Asi))));
    BCo = XOR_x8(Abo, XOR_x8(Ago, XOR_x8(Ako, XOR_x8(Amo, Aso))));
    BCu = XOR_x8(Abu, XOR_x8(Agu, XOR_x8(Aku, XOR_x8(Amu, Asu))));

    // thetaRhoPiChiIotaPrepareTheta(round, A, E)
    Da = XOR_x8(BCu, ROL_x8(BCe, 1));
    De = XOR_x8(BCa, ROL_x8(BCi, 1));
    Di = XOR_x8(BCe, ROL_x8(BCo, 1));
    Do = XOR_x8(BCi, ROL_x8(BCu, 1));
    Du = XOR_x8(BCo, ROL_x8(BCa, 1));

    Aba = XOR_x8(Aba, Da);
    BCa = Aba;
    Age = XOR_x8(Age, De);
    BCe = ROL_x8(Age, 44);
    Aki = XOR_x8(Aki, Di);
    BCi = ROL_x8(Aki, 43);
    Amo = XOR_x8(Amo, Do);
    BCo = ROL_x8(Amo, 21);
    Asu = XOR_x8(Asu, Du);
    BCu = ROL_x8(Asu, 14);
    Eba = XOR_x8(BCa, AND_x8((XOR_x8(BCe, ALL_ONE)), BCi));
    Eba = XOR_x8(Eba, KeccakF_RoundConstants_x8[round]);
    Ebe = XOR_x8(BCe, AND_x8((XOR_x8(BCi, ALL_ONE)), BCo));
    Ebi = XOR_x8(BCi, AND_x8((XOR_x8(BCo, ALL_ONE)), BCu));
    Ebo = XOR_x8(BCo, AND_x8((XOR_x8(BCu, ALL_ONE)), BCa));
    Ebu = XOR_x8(BCu, AND_x8((XOR_x8(BCa, ALL_ONE)), BCe));

    Abo = XOR_x8(Abo, Do);
    BCa = ROL_x8(Abo, 28);
    Agu = XOR_x8(Agu, Du);
    BCe = ROL_x8(Agu, 20);
    Aka = XOR_x8(Aka, Da);
    BCi = ROL_x8(Aka, 3);
    Ame = XOR_x8(Ame, De);
    BCo = ROL_x8(Ame, 45);
    Asi = XOR_x8(Asi, Di);
    BCu = ROL_x8(Asi, 61);
    Ega = XOR_x8(BCa, AND_x8((XOR_x8(BCe, ALL_ONE)), BCi));
    Ege = XOR_x8(BCe, AND_x8((XOR_x8(BCi, ALL_ONE)), BCo));
    Egi = XOR_x8(BCi, AND_x8((XOR_x8(BCo, ALL_ONE)), BCu));
    Ego = XOR_x8(BCo, AND_x8((XOR_x8(BCu, ALL_ONE)), BCa));
    Egu = XOR_x8(BCu, AND_x8((XOR_x8(BCa, ALL_ONE)), BCe));

    Abe = XOR_x8(Abe, De);
    BCa = ROL_x8(Abe, 1);
    Agi = XOR_x8(Agi, Di);
    BCe = ROL_x8(Agi, 6);
    Ako = XOR_x8(Ako, Do);
    BCi = ROL_x8(Ako, 25);
    Amu = XOR_x8(Amu, Du);
    BCo = ROL_x8(Amu, 8);
    Asa = XOR_x8(Asa, Da);
    BCu = ROL_x8(Asa, 18);
    Eka = XOR_x8(BCa, AND_x8((XOR_x8(BCe, ALL_ONE)), BCi));
    Eke = XOR_x8(BCe, AND_x8((XOR_x8(BCi, ALL_ONE)), BCo));
    Eki = XOR_x8(BCi, AND_x8((XOR_x8(BCo, ALL_ONE)), BCu));
    Eko = XOR_x8(BCo, AND_x8((XOR_x8(BCu, ALL_ONE)), BCa));
    Eku = XOR_x8(BCu, AND_x8((XOR_x8(BCa, ALL_ONE)), BCe));

    Abu = XOR_x8(Abu, Du);
    BCa = ROL_x8(Abu, 27);
    Aga = XOR_x8(Aga, Da);
    BCe = ROL_x8(Aga, 36);
    Ake = XOR_x8(Ake, De);
    BCi = ROL_x8(Ake, 10);
    Ami = XOR_x8(Ami, Di);
    BCo = ROL_x8(Ami, 15);
    Aso = XOR_x8(Aso, Do);
    BCu = ROL_x8(Aso, 56);
    Ema = XOR_x8(BCa, AND_x8((XOR_x8(BCe, ALL_ONE)), BCi));
    Eme = XOR_x8(BCe, AND_x8((XOR_x8(BCi, ALL_ONE)), BCo));
    Emi = XOR_x8(BCi, AND_x8((XOR_x8(BCo, ALL_ONE)), BCu));
    Emo = XOR_x8(BCo, AND_x8((XOR_x8(BCu, ALL_ONE)), BCa));
    Emu = XOR_x8(BCu, AND_x8((XOR_x8(BCa, ALL_ONE)), BCe));

    Abi = XOR_x8(Abi, Di);
    BCa = ROL_x8(Abi, 62);
    Ago = XOR_x8(Ago, Do);
    BCe = ROL_x8(Ago, 55);
    Aku = XOR_x8(Aku, Du);
    BCi = ROL_x8(Aku, 39);
    Ama = XOR_x8(Ama, Da);
    BCo = ROL_x8(Ama, 41);
    Ase = XOR_x8(Ase, De);
    BCu = ROL_x8(Ase, 2);
    Esa = XOR_x8(BCa, AND_x8((XOR_x8(BCe, ALL_ONE)), BCi));
    Ese = XOR_x8(BCe, AND_x8((XOR_x8(BCi, ALL_ONE)), BCo));
    Esi = XOR_x8(BCi, AND_x8((XOR_x8(BCo, ALL_ONE)), BCu));
    Eso = XOR_x8(BCo, AND_x8((XOR_x8(BCu, ALL_ONE)), BCa));
    Esu = XOR_x8(BCu, AND_x8((XOR_x8(BCa, ALL_ONE)), BCe));

    //    prepareTheta
    BCa = XOR_x8(Eba, XOR_x8(Ega, XOR_x8(Eka, XOR_x8(Ema, Esa))));
    BCe = XOR_x8(Ebe, XOR_x8(Ege, XOR_x8(Eke, XOR_x8(Eme, Ese))));
    BCi = XOR_x8(Ebi, XOR_x8(Egi, XOR_x8(Eki, XOR_x8(Emi, Esi))));
    BCo = XOR_x8(Ebo, XOR_x8(Ego, XOR_x8(Eko, XOR_x8(Emo, Eso))));
    BCu = XOR_x8(Ebu, XOR_x8(Egu, XOR_x8(Eku, XOR_x8(Emu, Esu))));

    // thetaRhoPiChiIotaPrepareTheta(round+1, E, A)
    Da = XOR_x8(BCu, ROL_x8(BCe, 1));
    De = XOR_x8(BCa, ROL_x8(BCi, 1));
    Di = XOR_x8(BCe, ROL_x8(BCo, 1));
    Do = XOR_x8(BCi, ROL_x8(BCu, 1));
    Du = XOR_x8(BCo, ROL_x8(BCa, 1));

    Eba = XOR_x8(Eba, Da);
    BCa = Eba;
    Ege = XOR_x8(Ege, De);
    BCe = ROL_x8(Ege, 44);
    Eki = XOR_x8(Eki, Di);
    BCi = ROL_x8(Eki, 43);
    Emo = XOR_x8(Emo, Do);
    BCo = ROL_x8(Emo, 21);
    Esu = XOR_x8(Esu, Du);
    BCu = ROL_x8(Esu, 14);
    Aba = XOR_x8(BCa, AND_x8((XOR_x8(BCe, ALL_ONE)), BCi));
    Aba = XOR_x8(Aba, KeccakF_RoundConstants_x8[round + 1]);
    Abe = XOR_x8(BCe, AND_x8((XOR_x8(BCi, ALL_ONE)), BCo));
    Abi = XOR_x8(BCi, AND_x8((XOR_x8(BCo, ALL_ONE)), BCu));
    Abo = XOR_x8(BCo, AND_x8((XOR_x8(BCu, ALL_ONE)), BCa));
    Abu = XOR_x8(BCu, AND_x8((XOR_x8(BCa, ALL_ONE)), BCe));

    Ebo = XOR_x8(Ebo, Do);
    BCa = ROL_x8(Ebo, 28);
    Egu = XOR_x8(Egu, Du);
    BCe = ROL_x8(Egu, 20);
    Eka = XOR_x8(Eka, Da);
    BCi = ROL_x8(Eka, 3);
    Eme = XOR_x8(Eme, De);
    BCo = ROL_x8(Eme, 45);
    Esi = XOR_x8(Esi, Di);
    BCu = ROL_x8(Esi, 61);
    Aga = XOR_x8(BCa, AND_x8((XOR_x8(BCe, ALL_ONE)), BCi));
    Age = XOR_x8(BCe, AND_x8((XOR_x8(BCi, ALL_ONE)), BCo));
    Agi = XOR_x8(BCi, AND_x8((XOR_x8(BCo, ALL_ONE)), BCu));
    Ago = XOR_x8(BCo, AND_x8((XOR_x8(BCu, ALL_ONE)), BCa));
    Agu = XOR_x8(BCu, AND_x8((XOR_x8(BCa, ALL_ONE)), BCe));

    Ebe = XOR_x8(Ebe, De);
    BCa = ROL_x8(Ebe, 1);
    Egi = XOR_x8(Egi, Di);
    BCe = ROL_x8(Egi, 6);
    Eko = XOR_x8(Eko, Do);
    BCi = ROL_x8(Eko, 25);
    Emu = XOR_x8(Emu, Du);
    BCo = ROL_x8(Emu, 8);
    Esa = XOR_x8(Esa, Da);
    BCu = ROL_x8(Esa, 18);
    Aka = XOR_x8(BCa, AND_x8((XOR_x8(BCe, ALL_ONE)), BCi));
    Ake = XOR_x8(BCe, AND_x8((XOR_x8(BCi, ALL_ONE)), BCo));
    Aki = XOR_x8(BCi, AND_x8((XOR_x8(BCo, ALL_ONE)), BCu));
    Ako = XOR_x8(BCo, AND_x8((XOR_x8(BCu, ALL_ONE)), BCa));
    Aku = XOR_x8(BCu, AND_x8((XOR_x8(BCa, ALL_ONE)), BCe));

    Ebu = XOR_x8(Ebu, Du);
    BCa = ROL_x8(Ebu, 27);
    Ega = XOR_x8(Ega, Da);
    BCe = ROL_x8(Ega, 36);
    Eke = XOR_x8(Eke, De);
    BCi = ROL_x8(Eke, 10);
    Emi = XOR_x8(Emi, Di);
    BCo = ROL_x8(Emi, 15);
    Eso = XOR_x8(Eso, Do);
    BCu = ROL_x8(Eso, 56);
    Ama = XOR_x8(BCa, AND_x8((XOR_x8(BCe, ALL_ONE)), BCi));
    Ame = XOR_x8(BCe, AND_x8((XOR_x8(BCi, ALL_ONE)), BCo));
    Ami = XOR_x8(BCi, AND_x8((XOR_x8(BCo, ALL_ONE)), BCu));
    Amo = XOR_x8(BCo, AND_x8((XOR_x8(BCu, ALL_ONE)), BCa));
    Amu = XOR_x8(BCu, AND_x8((XOR_x8(BCa, ALL_ONE)), BCe));

    Ebi = XOR_x8(Ebi, Di);
    BCa = ROL_x8(Ebi, 62);
    Ego = XOR_x8(Ego, Do);
    BCe = ROL_x8(Ego, 55);
    Eku = XOR_x8(Eku, Du);
    BCi = ROL_x8(Eku, 39);
    Ema = XOR_x8(Ema, Da);
    BCo = ROL_x8(Ema, 41);
    Ese = XOR_x8(Ese, De);
    BCu = ROL_x8(Ese, 2);
    Asa = XOR_x8(BCa, AND_x8((XOR_x8(BCe, ALL_ONE)), BCi));
    Ase = XOR_x8(BCe, AND_x8((XOR_x8(BCi, ALL_ONE)), BCo));
    Asi = XOR_x8(BCi, AND_x8((XOR_x8(BCo, ALL_ONE)), BCu));
    Aso = XOR_x8(BCo, AND_x8((XOR_x8(BCu, ALL_ONE)), BCa));
    Asu = XOR_x8(BCu, AND_x8((XOR_x8(BCa, ALL_ONE)), BCe));
  }

  // copyToState(state, A)
  state[0] = Aba;
  state[1] = Abe;
  state[2] = Abi;
  state[3] = Abo;
  state[4] = Abu;
  state[5] = Aga;
  state[6] = Age;
  state[7] = Agi;
  state[8] = Ago;
  state[9] = Agu;
  state[10] = Aka;
  state[11] = Ake;
  state[12] = Aki;
  state[13] = Ako;
  state[14] = Aku;
  state[15] = Ama;
  state[16] = Ame;
  state[17] = Ami;
  state[18] = Amo;
  state[19] = Amu;
  state[20] = Asa;
  state[21] = Ase;
  state[22] = Asi;
  state[23] = Aso;
  state[24] = Asu;
}

static void keccak_init_x8(pmod_vec_x8_t s[25]) {
  static const pmod_vec_x8_t zero = SET1_CT_x8(0);

  unsigned int i;
  for (i = 0; i < 25; i++) s[i] = zero;
}

static unsigned int keccak_absorb_x8_raw(pmod_vec_x8_t s[25], unsigned int pos,
                                          unsigned int r, const uint8_t **in,
                                          size_t inlen, uint32_t offset) {
  unsigned int i;

  while (pos + inlen >= r) {
    for (i = pos; i < r; i++) {
      const uint8_t *tmp = (*in++) + offset;
      static uint64_t buf[8] align64;
      for (int j = 0; j < 8; j++) buf[j] = tmp[j];
      pmod_vec_x8_t tmp_x8 = LOAD_x8(buf);
      s[i / 8] = XOR_x8(s[i / 8], SLLI_x8(tmp_x8, 8 * (i & 0x07)));
    }
    inlen -= r - pos;
    KeccakF1600_StatePermute_x8(s);
    pos = 0;
  }

  for (i = pos; i < pos + inlen; i++) {
    const uint8_t *tmp = (*in++) + offset;
    static uint64_t buf[8] align64;
    for (int j = 0; j < 8; j++) buf[j] = tmp[j];
    pmod_vec_x8_t tmp_x8 = LOAD_x8(buf);
    s[i / 8] = XOR_x8(s[i / 8], SLLI_x8(tmp_x8, 8 * (i & 0x07)));
  }

  return i;
}

static unsigned int keccak_absorb_x8(pmod_vec_x8_t s[25], unsigned int pos,
                                     unsigned int r, const pmod_vec_x8_t *in,
                                     size_t inlen) {
  unsigned int i;

  while (pos + inlen >= r) {
    for (i = pos; i < r; i++) {
      s[i / 8] = XOR_x8(s[i / 8], SLLI_x8(*in++, 8 * (i & 0x07)));
    }
    inlen -= r - pos;
    KeccakF1600_StatePermute_x8(s);
    pos = 0;
  }

  for (i = pos; i < pos + inlen; i++) {
    s[i / 8] = XOR_x8(s[i / 8], SLLI_x8(*in++, 8 * (i & 0x07)));
  }

  return i;
}

static void keccak_finalize_x8(pmod_vec_x8_t s[25], unsigned int pos,
                               unsigned int r, uint8_t p) {
  static const pmod_vec_x8_t pow63 = SET1_CT_x8(1ULL << 63);

  const pmod_vec_x8_t p_x8 = SET1_x8(p);

  s[pos / 8] = XOR_x8(s[pos / 8], SLLI_x8(p_x8, 8 * (pos % 8)));
  s[r / 8 - 1] = XOR_x8(s[r / 8 - 1], pow63);
}

static unsigned int keccak_squeeze_x8(pmod_vec_x8_t *out, size_t outlen,
                                      pmod_vec_x8_t s[25], unsigned int pos,
                                      unsigned int r) {
  static const pmod_vec_x8_t byte_mask = SET1_CT_x8(0xFF);

  unsigned int i;

  while (outlen) {
    if (pos == r) {
      KeccakF1600_StatePermute_x8(s);
      pos = 0;
    }
    for (i = pos; i < r && i < pos + outlen; i++)
      *out++ = AND_x8(SRLI_x8(s[i / 8], 8 * (i % 8)), byte_mask);
    outlen -= i - pos;
    pos = i;
  }

  return pos;
}

static void keccak_absorb_once_x8(pmod_vec_x8_t s[25], unsigned int r,
                                  const pmod_vec_x8_t *in, size_t inlen,
                                  uint8_t p) {
  unsigned int i;

  static const pmod_vec_x8_t zero = SET1_CT_x8(0);
  static const pmod_vec_x8_t pow63 = SET1_CT_x8(1ULL << 63);

  const pmod_vec_x8_t p_x8 = SET1_x8(p);

  for (i = 0; i < 25; i++) s[i] = zero;

  while (inlen >= r) {
    for (i = 0; i < r / 8; i++) s[i] = XOR_x8(s[i], load64_x8(in + 8 * i));
    in += r;
    inlen -= r;
    KeccakF1600_StatePermute_x8(s);
  }

  for (i = 0; i < inlen; i++)
    s[i / 8] = XOR_x8(s[i / 8], SLLI_x8(in[i], 8 * (i % 8)));

  s[i / 8] = XOR_x8(s[i / 8], SLLI_x8(p_x8, 8 * (i % 8)));
  s[(r - 1) / 8] = XOR_x8(s[(r - 1) / 8], pow63);
}

static void keccak_squeezeblocks_x8(pmod_vec_x8_t *out, size_t nblocks,
                                    pmod_vec_x8_t s[25], unsigned int r) {
  unsigned int i;

  while (nblocks) {
    KeccakF1600_StatePermute_x8(s);
    for (i = 0; i < r / 8; i++) store64_x8(out + 8 * i, s[i]);
    out += r;
    nblocks -= 1;
  }
}

void shake256_x8(pmod_vec_x8_t *out, size_t outlen, const pmod_vec_x8_t *in,
                 size_t inlen) {
  size_t nblocks;
  keccak_state_x8 state;

  shake256_x8_absorb_once(&state, in, inlen);
  nblocks = outlen / SHAKE256_RATE;
  shake256_x8_squeezeblocks(out, nblocks, &state);
  outlen -= nblocks * SHAKE256_RATE;
  out += nblocks * SHAKE256_RATE;
  shake256_x8_squeeze(out, outlen, &state);
}

void shake256_x8_init(keccak_state_x8 *state) {
  keccak_init_x8(state->s);
  state->pos = 0;
}

void shake256_x8_absorb_raw(keccak_state_x8 *state, const uint8_t **in,
                             size_t inlen, uint32_t offset) {
  state->pos =
      keccak_absorb_x8_raw(state->s, state->pos, SHAKE256_RATE, in, inlen, offset);
}

void shake256_x8_absorb(keccak_state_x8 *state, const pmod_vec_x8_t *in,
                        size_t inlen) {
  state->pos = keccak_absorb_x8(state->s, state->pos, SHAKE256_RATE, in, inlen);
}

void shake256_x8_finalize(keccak_state_x8 *state) {
  keccak_finalize_x8(state->s, state->pos, SHAKE256_RATE, 0x1F);
  state->pos = SHAKE256_RATE;
}

void shake256_x8_squeeze(pmod_vec_x8_t *out, size_t outlen,
                         keccak_state_x8 *state) {
  state->pos =
      keccak_squeeze_x8(out, outlen, state->s, state->pos, SHAKE256_RATE);
}

void shake256_x8_absorb_once(keccak_state_x8 *state, const pmod_vec_x8_t *in,
                             size_t inlen) {
  keccak_absorb_once_x8(state->s, SHAKE256_RATE, in, inlen, 0x1F);
  state->pos = SHAKE256_RATE;
}

void shake256_x8_squeezeblocks(pmod_vec_x8_t *out, size_t nblocks,
                               keccak_state_x8 *state) {
  keccak_squeezeblocks_x8(out, nblocks, state->s, SHAKE256_RATE);
}
