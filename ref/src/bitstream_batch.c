#include "bitstream_batch.h"

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

int bs_batch_init(bitstream_batch_t *bs, uint8_t **buf, size_t buf_len,
                  size_t batch_size) {
  bs->data = buf;
  bs->batch_size = batch_size;
  bs->buf_len = buf_len;
  bs->byte_pos = 0;
  bs->bit_pos = 0;

  return 0;
}

int bs_batch_finalize(bitstream_batch_t *bs) {
  bs->byte_pos += (bs->bit_pos > 0 ? 1 : 0);
  bs->bit_pos = 0;
  return bs->byte_pos - 1;
}

int bs_batch_write(bitstream_batch_t *bs, uint32_t *data, uint32_t data_len) {
  if (bs->byte_pos * 8 + bs->bit_pos + data_len > bs->buf_len * 8) {
    fprintf(stderr, "ERROR: bistream - write esceeds buffer!\n");
    return -1;
  }

  if (bs->bit_pos + data_len < 8) {
    for (int i = 0; i < bs->batch_size; i++)
      bs->data[i][bs->byte_pos] |= data[i] << bs->bit_pos;

    bs->bit_pos += data_len;

    if (bs->bit_pos > 7) {
      bs->bit_pos = 0;
      bs->byte_pos += 1;
    }

    return 0;
  }

  if (bs->bit_pos > 0) {
    for (int i = 0; i < bs->batch_size; i++) {
      bs->data[i][bs->byte_pos] |= (data[i] << bs->bit_pos) & 0xFF;
      data[i] >>= 8 - bs->bit_pos;
    }

    data_len -= 8 - bs->bit_pos;

    bs->bit_pos = 0;
    bs->byte_pos += 1;
  }

  while (data_len >= 8) {
    for (int i = 0; i < bs->batch_size; i++) {
      bs->data[i][bs->byte_pos] = data[i] & 0xFF;
      data[i] >>= 8;
    }

    data_len -= 8;

    bs->byte_pos += 1;
  }

  if (data_len > 0) {
    for (int i = 0; i < bs->batch_size; i++)
      bs->data[i][bs->byte_pos] = data[i];

    bs->bit_pos = data_len;
  }

  return 0;
}

int bs_batch_read(bitstream_batch_t *bs, uint32_t *buf, uint32_t data_len) {
  if (bs->byte_pos * 8 + bs->bit_pos + data_len > bs->buf_len * 8) {
    fprintf(stderr, "ERROR: bistream - read esceeds buffer!\n");
    return -1;
  }

  if (bs->bit_pos + data_len < 8) {
    for (int i = 0; i < bs->batch_size; i++)
      buf[i] =
          (bs->data[i][bs->byte_pos] >> bs->bit_pos) & ((1 << data_len) - 1);

    bs->bit_pos += data_len;

    if (bs->bit_pos > 7) {
      bs->bit_pos = 0;
      bs->byte_pos += 1;
    }

    return 0;
  }

  uint32_t off = 0;

  if (bs->bit_pos > 0) {
    for (int i = 0; i < bs->batch_size; i++)
      buf[i] = bs->data[i][bs->byte_pos] >> bs->bit_pos;

    off = 8 - bs->bit_pos;
    data_len -= 8 - bs->bit_pos;

    bs->bit_pos = 0;
    bs->byte_pos += 1;
  }

  while (data_len >= 8) {
    for (int i = 0; i < bs->batch_size; i++)
      buf[i] |= bs->data[i][bs->byte_pos] << off;

    off += 8;
    bs->byte_pos += 1;
    data_len -= 8;
    bs->bit_pos = 0;
  }

  if (data_len > 0) {
    for (int i = 0; i < bs->batch_size; i++)
      buf[i] |= (bs->data[i][bs->byte_pos] & ((1 << data_len) - 1)) << off;

    bs->bit_pos = data_len;
  }

  return 0;
}
