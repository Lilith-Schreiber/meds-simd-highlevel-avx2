#ifndef MEASURE_H
#define MEASURE_H

#define LOG_M(...) fprintf(measure_log, __VA_ARGS__)

void open_measure_log(const char*);

void close_measure_log();

void clear_measure_log(const char *);

double osfreq();

long long cpucycles(void);

#endif
