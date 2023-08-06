#ifndef BITSTREAM_BATCH_H
#define BITSTREAM_BATCH_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

typedef struct {
  uint8_t **data;
  size_t buf_len;
  size_t batch_size;
  uint32_t bit_pos;
  uint32_t byte_pos;
} bitstream_batch_t;

int bs_batch_init(bitstream_batch_t *bs, uint8_t **buf, size_t buf_len,
                  size_t batch_size);
int bs_batch_write(bitstream_batch_t *bs, uint32_t *data, uint32_t data_len, uint32_t offset);
int bs_batch_read(bitstream_batch_t *bs, uint32_t *buf, uint32_t data_len);
int bs_batch_finalize(bitstream_batch_t *bs);

#endif
