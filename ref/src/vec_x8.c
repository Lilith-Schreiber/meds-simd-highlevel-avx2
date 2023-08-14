#include "vec_x8.h"

uint64_t extract_vec_x8(__m512i x, int pos) {  
  uint64_t buf[8] align64;
  STORE_x8(buf, x);
  return buf[pos];
}

uint32_t extract_mask_x8(__mmask8 x, int pos) { return ((uint32_t)x) & (1 << pos); }

