#include <time.h>
#include <stdio.h>

#include "ecsw.h"

void setup(const enum ecsw_mode mode) {
    (void) mode;
}

void teardown(void) {
}

void output(char *message) {
  puts(message);
}

unsigned long long performance_count(void) {
  return clock();
}

char* performance_count_unit(void) {
    static char unit[32];

    snprintf(unit, sizeof(unit),
        "seconds/%llu", (unsigned long long) CLOCKS_PER_SEC);
    return unit;
}
