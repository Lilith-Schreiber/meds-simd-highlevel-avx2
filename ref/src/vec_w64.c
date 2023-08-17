#include "vec_w64.h"

uint64_t extract_vec_w64(__m512i x, int pos) {  
  uint64_t buf[8] aligned;
  STORE_w64(buf, x);
  return buf[pos];
}

uint32_t extract_mask_w64(__mmask8 x, int pos) { return ((uint32_t)x) & (1 << pos); }

