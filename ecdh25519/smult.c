#include "group.h"
#include "smult.h"


int crypto_scalarmult(unsigned char *ss, const unsigned char *sk, const unsigned char *pk)
{
  group_ge p, r0, r1;
  unsigned char t[32];
  int i;
  unsigned char bit, swap, prevswap = 0;
  uint32_t mask;
  fe25519 tmp;
  int j;

  for(i=0; i<32; i++) {
    t[i] = sk[i];
  }

  t[0] &= 248;
  t[31] &= 127;
  t[31] |= 64;

  if(group_ge_unpack(&p, pk)) return -1;

  // Montgomery ladder - constant time
  r0 = group_ge_neutral;  
  r1 = p;

  for(i=254; i>=0; i--)
  {
    bit = (t[i/8] >> (i & 7)) & 1;
    swap = bit ^ prevswap;
    prevswap = bit;

    // Constant-time conditional swap using XOR trick
    mask = (uint32_t)(-(int32_t)swap);
    for(j=0; j<32; j++) {
      tmp.v[j] = mask & (r0.x.v[j] ^ r1.x.v[j]);
      r0.x.v[j] ^= tmp.v[j];
      r1.x.v[j] ^= tmp.v[j];
      
      tmp.v[j] = mask & (r0.y.v[j] ^ r1.y.v[j]);
      r0.y.v[j] ^= tmp.v[j];
      r1.y.v[j] ^= tmp.v[j];
      
      tmp.v[j] = mask & (r0.z.v[j] ^ r1.z.v[j]);
      r0.z.v[j] ^= tmp.v[j];
      r1.z.v[j] ^= tmp.v[j];
      
      tmp.v[j] = mask & (r0.t.v[j] ^ r1.t.v[j]);
      r0.t.v[j] ^= tmp.v[j];
      r1.t.v[j] ^= tmp.v[j];
    }
    
    group_ge_add(&r0, &r0, &r1);
    group_ge_double(&r1, &r1);
  }
  
  // Final swap
  mask = (uint32_t)(-(int32_t)prevswap);
  for(j=0; j<32; j++) {
    tmp.v[j] = mask & (r0.x.v[j] ^ r1.x.v[j]);
    r0.x.v[j] ^= tmp.v[j];
    r1.x.v[j] ^= tmp.v[j];
    
    tmp.v[j] = mask & (r0.y.v[j] ^ r1.y.v[j]);
    r0.y.v[j] ^= tmp.v[j];
    r1.y.v[j] ^= tmp.v[j];
    
    tmp.v[j] = mask & (r0.z.v[j] ^ r1.z.v[j]);
    r0.z.v[j] ^= tmp.v[j];
    r1.z.v[j] ^= tmp.v[j];
    
    tmp.v[j] = mask & (r0.t.v[j] ^ r1.t.v[j]);
    r0.t.v[j] ^= tmp.v[j];
    r1.t.v[j] ^= tmp.v[j];
  }

  group_ge_pack(ss, &r0);
  return 0;
}

