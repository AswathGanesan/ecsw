#include "group.h"
#include "smult.h"

/* External declarations from group.c */
extern const group_ge group_ge_neutral;
extern const group_ge group_ge_base;

/* Precomputed table for base point multiples - windowed approach */
/* Table structure: base_table[window][multiple] where window goes from 0 to 63 (4 bits per window) */
static group_ge base_table[64][16];
static int base_table_initialized = 0;

/* Constant-time table lookup for side-channel resistance */
static void table_select(group_ge *r, const group_ge *table, int index, int table_size)
{
  int i, b;
  *r = group_ge_neutral;
  
  for(i = 0; i < table_size; i++) {
    /* Constant-time comparison: b = 1 if i == index, 0 otherwise */
    unsigned int cmp = (unsigned int)(i ^ index);
    cmp = (cmp - 1) >> 31;
    b = (int)cmp;
    
    /* Constant-time conditional move for each coordinate */
    group_ge_cmov(r, &table[i], b);
  }
}

/* Initialize precomputed base point table for fast base scalar multiplication */
static void init_base_table(void)
{
  group_ge p;
  int i, j;

  if(base_table_initialized) return;

  p = group_ge_base;

  /* For each 4-bit window position (64 windows for 256 bits) */
  for(i = 0; i < 64; i++) {
    /* Build table for this window: [0]P, [1]P, [2]P, ..., [15]P */
    base_table[i][0] = group_ge_neutral;
    base_table[i][1] = p;
    
    /* Use addition-chain approach for efficiency */
    for(j = 2; j < 16; j++) {
      if(j & 1) {
        /* Odd: add [1]P to previous */
        group_ge_add(&base_table[i][j], &base_table[i][j-1], &p);
      } else {
        /* Even: double the half */
        group_ge_double(&base_table[i][j], &base_table[i][j>>1]);
      }
    }

    /* Move to next window position: multiply P by 2^4 = 16 */
    for(j = 0; j < 4; j++) {
      group_ge_double(&p, &p);
    }
  }

  base_table_initialized = 1;
}

/* Optimized scalar multiplication with 4-bit windowing */
int crypto_scalarmult(unsigned char *ss, const unsigned char *sk, const unsigned char *pk)
{
  group_ge p, k, table[16], temp;
  unsigned char t[32];
  int i, index;

  /* Copy and clamp scalar */
  for(i=0; i<32; i++) {
    t[i] = sk[i];
  }

  t[0] &= 248;
  t[31] &= 127;
  t[31] |= 64;

  /* Unpack the public key point */
  if(group_ge_unpack(&p, pk)) return -1;

  /* Precompute table: [0]P, [1]P, [2]P, ..., [15]P for this specific point */
  table[0] = group_ge_neutral;
  table[1] = p;
  
  for(i = 2; i < 16; i++) {
    if(i & 1) {
      group_ge_add(&table[i], &table[i-1], &p);
    } else {
      group_ge_double(&table[i], &table[i>>1]);
    }
  }

  /* Initialize result to neutral element */
  k = group_ge_neutral;

  /* Process scalar in 4-bit windows from MSB to LSB (64 windows total) */
  for(i = 63; i >= 0; i--) {
    /* Double 4 times for next window */
    group_ge_double(&k, &k);
    group_ge_double(&k, &k);
    group_ge_double(&k, &k);
    group_ge_double(&k, &k);

    /* Extract 4-bit window value */
    int byte_idx = i >> 1;
    int shift = (i & 1) ? 4 : 0;
    index = (t[byte_idx] >> shift) & 0x0F;

    /* Constant-time table lookup and addition */
    table_select(&temp, table, index, 16);
    group_ge_add(&k, &k, &temp);
  }

  /* Pack result */
  group_ge_pack(ss, &k);
  return 0;
}

/* Optimized base point scalar multiplication using precomputed table */
int crypto_scalarmult_base(unsigned char *pk, const unsigned char *sk)
{
  group_ge k, temp;
  unsigned char t[32];
  int i, index;

  /* Initialize precomputed table on first use */
  init_base_table();

  /* Copy and clamp scalar */
  for(i=0; i<32; i++) {
    t[i] = sk[i];
  }

  t[0] &= 248;
  t[31] &= 127;
  t[31] |= 64;

  /* Start with neutral element */
  k = group_ge_neutral;

  /* Process each 4-bit window using precomputed table */
  /* We process from LSB to MSB for the precomputed table */
  for(i = 0; i < 64; i++) {
    int byte_idx = i >> 1;
    int shift = (i & 1) ? 4 : 0;
    index = (t[byte_idx] >> shift) & 0x0F;

    /* Constant-time lookup from precomputed table */
    table_select(&temp, base_table[i], index, 16);
    
    /* Add the precomputed point */
    group_ge_add(&k, &k, &temp);
  }

  /* Pack result */
  group_ge_pack(pk, &k);
  return 0;
}
