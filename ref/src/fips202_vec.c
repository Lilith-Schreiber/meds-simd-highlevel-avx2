#include "fips202_vec.h"

#define NROUNDS 24
#define ROL_w64(a, offset) XOR_w64(SLLI_w64(a, offset), SRLI_w64(a, (64 - offset)))

static pmod_mat_w64_t load64_w64(const pmod_mat_w64_t x[8]) {
  unsigned int i;
  pmod_mat_w64_t r = SET1_w64(0);

  for (i = 0; i < 8; i++) r = OR_w64(r, SLLI_w64(x[i], 8 * i));

  return r;
}

static void store64_w64(pmod_mat_w64_t x[8], pmod_mat_w64_t u) {
  unsigned int i;

  for (i = 0; i < 8; i++) x[i] = SRLI_w64(u, 8 * i);
}

/* Keccak round constants */
static const pmod_mat_w64_t KeccakF_RoundConstants_w64[NROUNDS] = {
    SET1_CT_w64(0x0000000000000001ULL), SET1_CT_w64(0x0000000000008082ULL),
    SET1_CT_w64(0x800000000000808aULL), SET1_CT_w64(0x8000000080008000ULL),
    SET1_CT_w64(0x000000000000808bULL), SET1_CT_w64(0x0000000080000001ULL),
    SET1_CT_w64(0x8000000080008081ULL), SET1_CT_w64(0x8000000000008009ULL),
    SET1_CT_w64(0x000000000000008aULL), SET1_CT_w64(0x0000000000000088ULL),
    SET1_CT_w64(0x0000000080008009ULL), SET1_CT_w64(0x000000008000000aULL),
    SET1_CT_w64(0x000000008000808bULL), SET1_CT_w64(0x800000000000008bULL),
    SET1_CT_w64(0x8000000000008089ULL), SET1_CT_w64(0x8000000000008003ULL),
    SET1_CT_w64(0x8000000000008002ULL), SET1_CT_w64(0x8000000000000080ULL),
    SET1_CT_w64(0x000000000000800aULL), SET1_CT_w64(0x800000008000000aULL),
    SET1_CT_w64(0x8000000080008081ULL), SET1_CT_w64(0x8000000000008080ULL),
    SET1_CT_w64(0x0000000080000001ULL), SET1_CT_w64(0x8000000080008008ULL)};

static void KeccakF1600_StatePermute_w64(pmod_mat_w64_t state[25]) {
  int round;

  static const pmod_mat_w64_t ALL_ONE = SET1_CT_w64(0xFFFFFFFFFFFFFFFFULL);

  pmod_mat_w64_t Aba, Abe, Abi, Abo, Abu;
  pmod_mat_w64_t Aga, Age, Agi, Ago, Agu;
  pmod_mat_w64_t Aka, Ake, Aki, Ako, Aku;
  pmod_mat_w64_t Ama, Ame, Ami, Amo, Amu;
  pmod_mat_w64_t Asa, Ase, Asi, Aso, Asu;
  pmod_mat_w64_t BCa, BCe, BCi, BCo, BCu;
  pmod_mat_w64_t Da, De, Di, Do, Du;
  pmod_mat_w64_t Eba, Ebe, Ebi, Ebo, Ebu;
  pmod_mat_w64_t Ega, Ege, Egi, Ego, Egu;
  pmod_mat_w64_t Eka, Eke, Eki, Eko, Eku;
  pmod_mat_w64_t Ema, Eme, Emi, Emo, Emu;
  pmod_mat_w64_t Esa, Ese, Esi, Eso, Esu;

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
    BCa = XOR_w64(Aba, XOR_w64(Aga, XOR_w64(Aka, XOR_w64(Ama, Asa))));
    BCe = XOR_w64(Abe, XOR_w64(Age, XOR_w64(Ake, XOR_w64(Ame, Ase))));
    BCi = XOR_w64(Abi, XOR_w64(Agi, XOR_w64(Aki, XOR_w64(Ami, Asi))));
    BCo = XOR_w64(Abo, XOR_w64(Ago, XOR_w64(Ako, XOR_w64(Amo, Aso))));
    BCu = XOR_w64(Abu, XOR_w64(Agu, XOR_w64(Aku, XOR_w64(Amu, Asu))));

    // thetaRhoPiChiIotaPrepareTheta(round, A, E)
    Da = XOR_w64(BCu, ROL_w64(BCe, 1));
    De = XOR_w64(BCa, ROL_w64(BCi, 1));
    Di = XOR_w64(BCe, ROL_w64(BCo, 1));
    Do = XOR_w64(BCi, ROL_w64(BCu, 1));
    Du = XOR_w64(BCo, ROL_w64(BCa, 1));

    Aba = XOR_w64(Aba, Da);
    BCa = Aba;
    Age = XOR_w64(Age, De);
    BCe = ROL_w64(Age, 44);
    Aki = XOR_w64(Aki, Di);
    BCi = ROL_w64(Aki, 43);
    Amo = XOR_w64(Amo, Do);
    BCo = ROL_w64(Amo, 21);
    Asu = XOR_w64(Asu, Du);
    BCu = ROL_w64(Asu, 14);
    Eba = XOR_w64(BCa, AND_w64((XOR_w64(BCe, ALL_ONE)), BCi));
    Eba = XOR_w64(Eba, KeccakF_RoundConstants_w64[round]);
    Ebe = XOR_w64(BCe, AND_w64((XOR_w64(BCi, ALL_ONE)), BCo));
    Ebi = XOR_w64(BCi, AND_w64((XOR_w64(BCo, ALL_ONE)), BCu));
    Ebo = XOR_w64(BCo, AND_w64((XOR_w64(BCu, ALL_ONE)), BCa));
    Ebu = XOR_w64(BCu, AND_w64((XOR_w64(BCa, ALL_ONE)), BCe));

    Abo = XOR_w64(Abo, Do);
    BCa = ROL_w64(Abo, 28);
    Agu = XOR_w64(Agu, Du);
    BCe = ROL_w64(Agu, 20);
    Aka = XOR_w64(Aka, Da);
    BCi = ROL_w64(Aka, 3);
    Ame = XOR_w64(Ame, De);
    BCo = ROL_w64(Ame, 45);
    Asi = XOR_w64(Asi, Di);
    BCu = ROL_w64(Asi, 61);
    Ega = XOR_w64(BCa, AND_w64((XOR_w64(BCe, ALL_ONE)), BCi));
    Ege = XOR_w64(BCe, AND_w64((XOR_w64(BCi, ALL_ONE)), BCo));
    Egi = XOR_w64(BCi, AND_w64((XOR_w64(BCo, ALL_ONE)), BCu));
    Ego = XOR_w64(BCo, AND_w64((XOR_w64(BCu, ALL_ONE)), BCa));
    Egu = XOR_w64(BCu, AND_w64((XOR_w64(BCa, ALL_ONE)), BCe));

    Abe = XOR_w64(Abe, De);
    BCa = ROL_w64(Abe, 1);
    Agi = XOR_w64(Agi, Di);
    BCe = ROL_w64(Agi, 6);
    Ako = XOR_w64(Ako, Do);
    BCi = ROL_w64(Ako, 25);
    Amu = XOR_w64(Amu, Du);
    BCo = ROL_w64(Amu, 8);
    Asa = XOR_w64(Asa, Da);
    BCu = ROL_w64(Asa, 18);
    Eka = XOR_w64(BCa, AND_w64((XOR_w64(BCe, ALL_ONE)), BCi));
    Eke = XOR_w64(BCe, AND_w64((XOR_w64(BCi, ALL_ONE)), BCo));
    Eki = XOR_w64(BCi, AND_w64((XOR_w64(BCo, ALL_ONE)), BCu));
    Eko = XOR_w64(BCo, AND_w64((XOR_w64(BCu, ALL_ONE)), BCa));
    Eku = XOR_w64(BCu, AND_w64((XOR_w64(BCa, ALL_ONE)), BCe));

    Abu = XOR_w64(Abu, Du);
    BCa = ROL_w64(Abu, 27);
    Aga = XOR_w64(Aga, Da);
    BCe = ROL_w64(Aga, 36);
    Ake = XOR_w64(Ake, De);
    BCi = ROL_w64(Ake, 10);
    Ami = XOR_w64(Ami, Di);
    BCo = ROL_w64(Ami, 15);
    Aso = XOR_w64(Aso, Do);
    BCu = ROL_w64(Aso, 56);
    Ema = XOR_w64(BCa, AND_w64((XOR_w64(BCe, ALL_ONE)), BCi));
    Eme = XOR_w64(BCe, AND_w64((XOR_w64(BCi, ALL_ONE)), BCo));
    Emi = XOR_w64(BCi, AND_w64((XOR_w64(BCo, ALL_ONE)), BCu));
    Emo = XOR_w64(BCo, AND_w64((XOR_w64(BCu, ALL_ONE)), BCa));
    Emu = XOR_w64(BCu, AND_w64((XOR_w64(BCa, ALL_ONE)), BCe));

    Abi = XOR_w64(Abi, Di);
    BCa = ROL_w64(Abi, 62);
    Ago = XOR_w64(Ago, Do);
    BCe = ROL_w64(Ago, 55);
    Aku = XOR_w64(Aku, Du);
    BCi = ROL_w64(Aku, 39);
    Ama = XOR_w64(Ama, Da);
    BCo = ROL_w64(Ama, 41);
    Ase = XOR_w64(Ase, De);
    BCu = ROL_w64(Ase, 2);
    Esa = XOR_w64(BCa, AND_w64((XOR_w64(BCe, ALL_ONE)), BCi));
    Ese = XOR_w64(BCe, AND_w64((XOR_w64(BCi, ALL_ONE)), BCo));
    Esi = XOR_w64(BCi, AND_w64((XOR_w64(BCo, ALL_ONE)), BCu));
    Eso = XOR_w64(BCo, AND_w64((XOR_w64(BCu, ALL_ONE)), BCa));
    Esu = XOR_w64(BCu, AND_w64((XOR_w64(BCa, ALL_ONE)), BCe));

    //    prepareTheta
    BCa = XOR_w64(Eba, XOR_w64(Ega, XOR_w64(Eka, XOR_w64(Ema, Esa))));
    BCe = XOR_w64(Ebe, XOR_w64(Ege, XOR_w64(Eke, XOR_w64(Eme, Ese))));
    BCi = XOR_w64(Ebi, XOR_w64(Egi, XOR_w64(Eki, XOR_w64(Emi, Esi))));
    BCo = XOR_w64(Ebo, XOR_w64(Ego, XOR_w64(Eko, XOR_w64(Emo, Eso))));
    BCu = XOR_w64(Ebu, XOR_w64(Egu, XOR_w64(Eku, XOR_w64(Emu, Esu))));

    // thetaRhoPiChiIotaPrepareTheta(round+1, E, A)
    Da = XOR_w64(BCu, ROL_w64(BCe, 1));
    De = XOR_w64(BCa, ROL_w64(BCi, 1));
    Di = XOR_w64(BCe, ROL_w64(BCo, 1));
    Do = XOR_w64(BCi, ROL_w64(BCu, 1));
    Du = XOR_w64(BCo, ROL_w64(BCa, 1));

    Eba = XOR_w64(Eba, Da);
    BCa = Eba;
    Ege = XOR_w64(Ege, De);
    BCe = ROL_w64(Ege, 44);
    Eki = XOR_w64(Eki, Di);
    BCi = ROL_w64(Eki, 43);
    Emo = XOR_w64(Emo, Do);
    BCo = ROL_w64(Emo, 21);
    Esu = XOR_w64(Esu, Du);
    BCu = ROL_w64(Esu, 14);
    Aba = XOR_w64(BCa, AND_w64((XOR_w64(BCe, ALL_ONE)), BCi));
    Aba = XOR_w64(Aba, KeccakF_RoundConstants_w64[round + 1]);
    Abe = XOR_w64(BCe, AND_w64((XOR_w64(BCi, ALL_ONE)), BCo));
    Abi = XOR_w64(BCi, AND_w64((XOR_w64(BCo, ALL_ONE)), BCu));
    Abo = XOR_w64(BCo, AND_w64((XOR_w64(BCu, ALL_ONE)), BCa));
    Abu = XOR_w64(BCu, AND_w64((XOR_w64(BCa, ALL_ONE)), BCe));

    Ebo = XOR_w64(Ebo, Do);
    BCa = ROL_w64(Ebo, 28);
    Egu = XOR_w64(Egu, Du);
    BCe = ROL_w64(Egu, 20);
    Eka = XOR_w64(Eka, Da);
    BCi = ROL_w64(Eka, 3);
    Eme = XOR_w64(Eme, De);
    BCo = ROL_w64(Eme, 45);
    Esi = XOR_w64(Esi, Di);
    BCu = ROL_w64(Esi, 61);
    Aga = XOR_w64(BCa, AND_w64((XOR_w64(BCe, ALL_ONE)), BCi));
    Age = XOR_w64(BCe, AND_w64((XOR_w64(BCi, ALL_ONE)), BCo));
    Agi = XOR_w64(BCi, AND_w64((XOR_w64(BCo, ALL_ONE)), BCu));
    Ago = XOR_w64(BCo, AND_w64((XOR_w64(BCu, ALL_ONE)), BCa));
    Agu = XOR_w64(BCu, AND_w64((XOR_w64(BCa, ALL_ONE)), BCe));

    Ebe = XOR_w64(Ebe, De);
    BCa = ROL_w64(Ebe, 1);
    Egi = XOR_w64(Egi, Di);
    BCe = ROL_w64(Egi, 6);
    Eko = XOR_w64(Eko, Do);
    BCi = ROL_w64(Eko, 25);
    Emu = XOR_w64(Emu, Du);
    BCo = ROL_w64(Emu, 8);
    Esa = XOR_w64(Esa, Da);
    BCu = ROL_w64(Esa, 18);
    Aka = XOR_w64(BCa, AND_w64((XOR_w64(BCe, ALL_ONE)), BCi));
    Ake = XOR_w64(BCe, AND_w64((XOR_w64(BCi, ALL_ONE)), BCo));
    Aki = XOR_w64(BCi, AND_w64((XOR_w64(BCo, ALL_ONE)), BCu));
    Ako = XOR_w64(BCo, AND_w64((XOR_w64(BCu, ALL_ONE)), BCa));
    Aku = XOR_w64(BCu, AND_w64((XOR_w64(BCa, ALL_ONE)), BCe));

    Ebu = XOR_w64(Ebu, Du);
    BCa = ROL_w64(Ebu, 27);
    Ega = XOR_w64(Ega, Da);
    BCe = ROL_w64(Ega, 36);
    Eke = XOR_w64(Eke, De);
    BCi = ROL_w64(Eke, 10);
    Emi = XOR_w64(Emi, Di);
    BCo = ROL_w64(Emi, 15);
    Eso = XOR_w64(Eso, Do);
    BCu = ROL_w64(Eso, 56);
    Ama = XOR_w64(BCa, AND_w64((XOR_w64(BCe, ALL_ONE)), BCi));
    Ame = XOR_w64(BCe, AND_w64((XOR_w64(BCi, ALL_ONE)), BCo));
    Ami = XOR_w64(BCi, AND_w64((XOR_w64(BCo, ALL_ONE)), BCu));
    Amo = XOR_w64(BCo, AND_w64((XOR_w64(BCu, ALL_ONE)), BCa));
    Amu = XOR_w64(BCu, AND_w64((XOR_w64(BCa, ALL_ONE)), BCe));

    Ebi = XOR_w64(Ebi, Di);
    BCa = ROL_w64(Ebi, 62);
    Ego = XOR_w64(Ego, Do);
    BCe = ROL_w64(Ego, 55);
    Eku = XOR_w64(Eku, Du);
    BCi = ROL_w64(Eku, 39);
    Ema = XOR_w64(Ema, Da);
    BCo = ROL_w64(Ema, 41);
    Ese = XOR_w64(Ese, De);
    BCu = ROL_w64(Ese, 2);
    Asa = XOR_w64(BCa, AND_w64((XOR_w64(BCe, ALL_ONE)), BCi));
    Ase = XOR_w64(BCe, AND_w64((XOR_w64(BCi, ALL_ONE)), BCo));
    Asi = XOR_w64(BCi, AND_w64((XOR_w64(BCo, ALL_ONE)), BCu));
    Aso = XOR_w64(BCo, AND_w64((XOR_w64(BCu, ALL_ONE)), BCa));
    Asu = XOR_w64(BCu, AND_w64((XOR_w64(BCa, ALL_ONE)), BCe));
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

static void keccak_init_w64(pmod_mat_w64_t s[25]) {
  static const pmod_mat_w64_t zero = SET1_CT_w64(0);

  unsigned int i;
  for (i = 0; i < 25; i++) s[i] = zero;
}

static unsigned int keccak_absorb_w64_raw(pmod_mat_w64_t s[25], unsigned int pos,
                                          unsigned int r, const uint8_t **in,
                                          size_t inlen, uint32_t offset) {
  unsigned int i;

  while (pos + inlen >= r) {
    for (i = pos; i < r; i++) {
      const uint8_t *tmp = (*in++) + offset;
      static uint64_t buf[8] aligned;
      for (int j = 0; j < 8; j++) buf[j] = tmp[j];
      pmod_mat_w64_t tmp_w64 = LOAD_w64(buf);
      s[i / 8] = XOR_w64(s[i / 8], SLLI_w64(tmp_w64, 8 * (i & 0x07)));
    }
    inlen -= r - pos;
    KeccakF1600_StatePermute_w64(s);
    pos = 0;
  }

  for (i = pos; i < pos + inlen; i++) {
    const uint8_t *tmp = (*in++) + offset;
    static uint64_t buf[8] aligned;
    for (int j = 0; j < 8; j++) buf[j] = tmp[j];
    pmod_mat_w64_t tmp_w64 = LOAD_w64(buf);
    s[i / 8] = XOR_w64(s[i / 8], SLLI_w64(tmp_w64, 8 * (i & 0x07)));
  }

  return i;
}

static unsigned int keccak_absorb_w64(pmod_mat_w64_t s[25], unsigned int pos,
                                     unsigned int r, const pmod_mat_w64_t *in,
                                     size_t inlen) {
  unsigned int i;

  while (pos + inlen >= r) {
    for (i = pos; i < r; i++) {
      s[i / 8] = XOR_w64(s[i / 8], SLLI_w64(*in++, 8 * (i & 0x07)));
    }
    inlen -= r - pos;
    KeccakF1600_StatePermute_w64(s);
    pos = 0;
  }

  for (i = pos; i < pos + inlen; i++) {
    s[i / 8] = XOR_w64(s[i / 8], SLLI_w64(*in++, 8 * (i & 0x07)));
  }

  return i;
}

static void keccak_finalize_w64(pmod_mat_w64_t s[25], unsigned int pos,
                               unsigned int r, uint8_t p) {
  static const pmod_mat_w64_t pow63 = SET1_CT_w64(1ULL << 63);

  const pmod_mat_w64_t p_w64 = SET1_w64(p);

  s[pos / 8] = XOR_w64(s[pos / 8], SLLI_w64(p_w64, 8 * (pos % 8)));
  s[r / 8 - 1] = XOR_w64(s[r / 8 - 1], pow63);
}

static unsigned int keccak_squeeze_w64(pmod_mat_w64_t *out, size_t outlen,
                                      pmod_mat_w64_t s[25], unsigned int pos,
                                      unsigned int r) {
  static const pmod_mat_w64_t byte_mask = SET1_CT_w64(0xFF);

  unsigned int i;

  while (outlen) {
    if (pos == r) {
      KeccakF1600_StatePermute_w64(s);
      pos = 0;
    }
    for (i = pos; i < r && i < pos + outlen; i++)
      *out++ = AND_w64(SRLI_w64(s[i / 8], 8 * (i % 8)), byte_mask);
    outlen -= i - pos;
    pos = i;
  }

  return pos;
}

static void keccak_absorb_once_w64(pmod_mat_w64_t s[25], unsigned int r,
                                  const pmod_mat_w64_t *in, size_t inlen,
                                  uint8_t p) {
  unsigned int i;

  static const pmod_mat_w64_t zero = SET1_CT_w64(0);
  static const pmod_mat_w64_t pow63 = SET1_CT_w64(1ULL << 63);

  const pmod_mat_w64_t p_w64 = SET1_w64(p);

  for (i = 0; i < 25; i++) s[i] = zero;

  while (inlen >= r) {
    for (i = 0; i < r / 8; i++) s[i] = XOR_w64(s[i], load64_w64(in + 8 * i));
    in += r;
    inlen -= r;
    KeccakF1600_StatePermute_w64(s);
  }

  for (i = 0; i < inlen; i++)
    s[i / 8] = XOR_w64(s[i / 8], SLLI_w64(in[i], 8 * (i % 8)));

  s[i / 8] = XOR_w64(s[i / 8], SLLI_w64(p_w64, 8 * (i % 8)));
  s[(r - 1) / 8] = XOR_w64(s[(r - 1) / 8], pow63);
}

static void keccak_squeezeblocks_w64(pmod_mat_w64_t *out, size_t nblocks,
                                    pmod_mat_w64_t s[25], unsigned int r) {
  unsigned int i;

  while (nblocks) {
    KeccakF1600_StatePermute_w64(s);
    for (i = 0; i < r / 8; i++) store64_w64(out + 8 * i, s[i]);
    out += r;
    nblocks -= 1;
  }
}

void shake256_w64(pmod_mat_w64_t *out, size_t outlen, const pmod_mat_w64_t *in,
                 size_t inlen) {
  size_t nblocks;
  keccak_state_w64 state;

  shake256_w64_absorb_once(&state, in, inlen);
  nblocks = outlen / SHAKE256_RATE;
  shake256_w64_squeezeblocks(out, nblocks, &state);
  outlen -= nblocks * SHAKE256_RATE;
  out += nblocks * SHAKE256_RATE;
  shake256_w64_squeeze(out, outlen, &state);
}

void shake256_w64_init(keccak_state_w64 *state) {
  keccak_init_w64(state->s);
  state->pos = 0;
}

void shake256_w64_absorb_raw(keccak_state_w64 *state, const uint8_t **in,
                             size_t inlen, uint32_t offset) {
  state->pos =
      keccak_absorb_w64_raw(state->s, state->pos, SHAKE256_RATE, in, inlen, offset);
}

void shake256_w64_absorb(keccak_state_w64 *state, const pmod_mat_w64_t *in,
                        size_t inlen) {
  state->pos = keccak_absorb_w64(state->s, state->pos, SHAKE256_RATE, in, inlen);
}

void shake256_w64_finalize(keccak_state_w64 *state) {
  keccak_finalize_w64(state->s, state->pos, SHAKE256_RATE, 0x1F);
  state->pos = SHAKE256_RATE;
}

void shake256_w64_squeeze(pmod_mat_w64_t *out, size_t outlen,
                         keccak_state_w64 *state) {
  state->pos =
      keccak_squeeze_w64(out, outlen, state->s, state->pos, SHAKE256_RATE);
}

void shake256_w64_absorb_once(keccak_state_w64 *state, const pmod_mat_w64_t *in,
                             size_t inlen) {
  keccak_absorb_once_w64(state->s, SHAKE256_RATE, in, inlen, 0x1F);
  state->pos = SHAKE256_RATE;
}

void shake256_w64_squeezeblocks(pmod_mat_w64_t *out, size_t nblocks,
                               keccak_state_w64 *state) {
  keccak_squeezeblocks_w64(out, nblocks, state->s, SHAKE256_RATE);
}
