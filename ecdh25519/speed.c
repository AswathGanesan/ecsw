#include <stdio.h>
#include "smult.h"
#include "../common/ecsw.h"

#define OUTLEN 1024

unsigned char sk[32] = {
  0x57, 0x6c, 0x7c, 0x77, 0x6a, 0xc2, 0x93, 0xc6, 0x78, 0x3a, 0x4a, 0x48, 0xc9, 0x45, 0x20, 0x36, 
  0x7d, 0xb3, 0xd4, 0x8c, 0x66, 0xa0, 0x52, 0xa8, 0xb2, 0xea, 0x09, 0xdc, 0x41, 0x43, 0xc5, 0x61};

unsigned char pk[32];
unsigned char ss[32];

int main(void)
{
  char outstr[128];
  unsigned long long oldcount, newcount;


  setup(ECSW_BENCHMARK);

  output("\n============ IGNORE OUTPUT BEFORE THIS LINE ============\n");

  oldcount = performance_count();
  crypto_scalarmult_base(pk, sk);
  newcount = performance_count()-oldcount;

  sprintf(outstr, "\n%s for scalarmult_base: %u",
    performance_count_unit(), (unsigned) newcount);
  output(outstr);

  oldcount = performance_count();
  crypto_scalarmult(ss, sk, pk);
  newcount = performance_count()-oldcount;

  sprintf(outstr, "\n%s for scalarmult: %u",
    performance_count_unit(), (unsigned) newcount);
  output(outstr);

  teardown();

  return 0;
}
