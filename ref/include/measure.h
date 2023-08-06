#ifndef MEASURE_H
#define MEASURE_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>

#define LOG_MEASURE

#ifdef LOG_MEASURE

#define START(a) long long a = -cpucycles()
#define END(s, a)   \
  a += cpucycles(); \
  s += a;

#define LOG_M(...) fprintf(measure_log, __VA_ARGS__)

#else

#define START(a) do { } while (0);
#define END(s, a) do { } while (0);

#define LOG_M(...) do { } while (0);

#endif

double osfreq();
long long cpucycles(void);

void open_measure_log(const char*);
void close_measure_log();

#endif
