#include <string.h>
#include "group.h"
#include "smult.h"

/* Window size for scalar multiplication */
#define WINDOW_SIZE 4
#define WINDOW_MASK ((1 << WINDOW_SIZE) - 1)

/* Identity/neutral element for group operations */
static const group_ge group_ge_neutral = {{0}, {0}, {0}, {0}};

/* Precomputed table for base point multiples (for base point multiplication)
 * Table contains [1]B, [3]B, [5]B, ..., [15]B for each window position
 */

/* Helper function to select from precomputed table in constant time */
// Add this helper function for constant-time table lookup
static void table_select(group_ge *r, const group_ge *table, int index, int table_size)
{
  int i, b;
  *r = group_ge_neutral;
  
  for(i = 0; i < table_size; i++) {
    // Constant-time comparison: b = 1 if i == index, 0 otherwise
    b = ((i ^ index) - 1) >> 31;
    
    // Constant-time conditional move for each coordinate
    for(int j = 0; j < 32; j++) {
      unsigned char mask = -b;
      r->x.v[j] ^= mask & (r->x.v[j] ^ table[i].x.v[j]);
      r->y.v[j] ^= mask & (r->y.v[j] ^ table[i].y.v[j]);
      r->z.v[j] ^= mask & (r->z.v[j] ^ table[i].z.v[j]);
      r->t.v[j] ^= mask & (r->t.v[j] ^ table[i].t.v[j]);
    }
  }
}

// Optimized scalar multiplication with 4-bit windowing
int crypto_scalarmult(unsigned char *ss, const unsigned char *sk, const unsigned char *pk)
{
  group_ge p, k, table[16], temp;
  unsigned char t[32];
  int i, j;

  for(i=0; i<32; i++) {
    t[i] = sk[i];
  }

  t[0] &= 248;
  t[31] &= 127;
  t[31] |= 64;

  if(group_ge_unpack(&p, pk)) return -1;

  // Precompute table: [0]P, [1]P, [2]P, ..., [15]P
  table[0] = group_ge_neutral;
  table[1] = p;
  for(i = 2; i < 16; i++) {
    if(i & 1) {
      group_ge_add(&table[i], &table[i-1], &p);
    } else {
      group_ge_double(&table[i], &table[i>>1]);
    }
  }

  // Start from neutral element
  k = group_ge_neutral;

  // Process 4 bits at a time from MSB to LSB
  for(i = 63; i >= 0; i--) {
    // Double 4 times
    group_ge_double(&k, &k);
    group_ge_double(&k, &k);
    group_ge_double(&k, &k);
    group_ge_double(&k, &k);

    // Extract 4-bit window
    int byte_idx = i >> 1;
    int shift = (i & 1) ? 4 : 0;
    int index = (t[byte_idx] >> shift) & 0x0F;

    // Constant-time table lookup and addition
    table_select(&temp, table, index, 16);
    group_ge_add(&k, &k, &temp);
  }

  group_ge_pack(ss, &k);
  return 0;
}
/* Constant-time conditional negate: if b == 1, negate point */
static void group_ge_cneg(group_ge *r, const group_ge *p, unsigned char b) {
  fe25519 neg_x, neg_t;
  
  fe25519_neg(&neg_x, &p->x);
  fe25519_neg(&neg_t, &p->t);
  
  fe25519_cmov(&r->x, &neg_x, b);
  r->y = p->y;
  r->z = p->z;
  fe25519_cmov(&r->t, &neg_t, b);
}

/* Compute precomputation table for a point P
 * Table will contain: P, 2P, 3P, ..., (2^w-1)P
 */
static void compute_table(group_ge *table, const group_ge *p) {
  int i;
  
  table[0] = group_ge_neutral;
  table[1] = *p;
  
  for (i = 2; i < (1 << WINDOW_SIZE); i++) {
    if (i & 1) {
      /* Odd: add P to previous even entry */
      group_ge_add(&table[i], &table[i-1], p);
    } else {
      /* Even: double the entry at i/2 */
      group_ge_double(&table[i], &table[i/2]);
    }
  }
}

/* Windowed scalar multiplication using constant-time table lookups
 * Processes scalar from MSB to LSB in windows of WINDOW_SIZE bits
 */
static void scalarmult_windowed(group_ge *r, const group_ge *p, const unsigned char *scalar) {
  group_ge table[1 << WINDOW_SIZE];
  group_ge acc, temp;
  int i, j;
  unsigned char window;
  unsigned char is_negative;
  
  /* Compute precomputation table */
  compute_table(table, p);
  
  /* Initialize accumulator to identity */
  acc = group_ge_neutral;
  
  /* Process scalar from MSB to LSB */
  for (i = 255; i >= 0; i -= WINDOW_SIZE) {
    /* Double WINDOW_SIZE times */
    for (j = 0; j < WINDOW_SIZE && i - j >= 0; j++) {
      group_ge_double(&acc, &acc);
    }
    
    /* Extract window */
    if (i >= WINDOW_SIZE - 1) {
      int byte_idx = i / 8;
      int bit_idx = i % 8;
      
      /* Extract WINDOW_SIZE bits starting at position i */
      if (bit_idx >= WINDOW_SIZE - 1) {
        window = (scalar[byte_idx] >> (bit_idx - WINDOW_SIZE + 1)) & WINDOW_MASK;
      } else {
        /* Window spans two bytes */
        window = (scalar[byte_idx] << (WINDOW_SIZE - 1 - bit_idx));
        if (byte_idx + 1 < 32) {
          window |= (scalar[byte_idx + 1] >> (8 - WINDOW_SIZE + 1 + bit_idx));
        }
        window &= WINDOW_MASK;
      }
    } else {
      /* Last partial window */
      window = scalar[0] & ((1 << (i + 1)) - 1);
    }
    
    /* Constant-time table lookup */
    table_select(&temp, table, window, 1 << WINDOW_SIZE);
    
    /* Add selected point to accumulator */
    group_ge_add(&acc, &acc, &temp);
  }
  
  *r = acc;
}

/* Optimized general scalar multiplication using signed windows */
static void scalarmult_signed_windowed(group_ge *r, const group_ge *p, const unsigned char *scalar) {
  group_ge table[1 << (WINDOW_SIZE - 1)];  /* Only store odd multiples */
  group_ge acc, temp, selected;
  int i, j;
  unsigned char window;
  unsigned char is_negative;
  int abs_window;
  
  /* Compute precomputation table: P, 3P, 5P, ..., (2^(w-1)-1)*2*P */
  table[0] = *p;
  group_ge_double(&temp, p);
  
  for (i = 1; i < (1 << (WINDOW_SIZE - 1)); i++) {
    group_ge_add(&table[i], &table[i-1], &temp);
  }
  
  /* Initialize accumulator */
  acc = group_ge_neutral;
  
  /* Process from MSB to LSB */
  for (i = 252; i >= 0; i -= WINDOW_SIZE) {
    /* Double WINDOW_SIZE times */
    for (j = 0; j < WINDOW_SIZE; j++) {
      group_ge_double(&acc, &acc);
    }
    
    /* Extract signed window */
    int byte_idx = i / 8;
    int bit_idx = i % 8;
    
    if (bit_idx >= WINDOW_SIZE - 1) {
      window = (scalar[byte_idx] >> (bit_idx - WINDOW_SIZE + 1)) & WINDOW_MASK;
    } else {
      window = (scalar[byte_idx] << (WINDOW_SIZE - 1 - bit_idx));
      if (byte_idx + 1 < 32) {
        window |= (scalar[byte_idx + 1] >> (8 - WINDOW_SIZE + 1 + bit_idx));
      }
      window &= WINDOW_MASK;
    }
    
    /* Convert to signed representation */
    is_negative = window >> (WINDOW_SIZE - 1);
    abs_window = (window ^ (is_negative * WINDOW_MASK)) + is_negative;
    abs_window = (abs_window >> 1);  /* Map to odd multiples */
    
    /* Constant-time table lookup */
    selected = group_ge_neutral;
    for (j = 0; j < (1 << (WINDOW_SIZE - 1)); j++) {
      unsigned char match = (unsigned char)((j == abs_window) ? 1 : 0);
      group_ge_cmov(&selected, &table[j], match);
    }
    
    /* Conditionally negate */
    temp = selected;
    group_ge_cneg(&temp, &selected, is_negative);
    
    /* Add to accumulator */
    group_ge_add(&acc, &acc, &temp);
  }
  
  *r = acc;
}

/* General scalar multiplication: k * P
 * Uses constant-time windowed method
 */
int crypto_scalarmult(unsigned char *ss, const unsigned char *sk, const unsigned char *pk) {
  group_ge p, result;
  unsigned char clamped_scalar[32];
  int i;
  
  /* Unpack public key point */
  if (group_ge_unpack(&p, pk)) {
    return -1;
  }
  
  /* Clamp scalar (Curve25519 requirement) */
  for (i = 0; i < 32; i++) {
    clamped_scalar[i] = sk[i];
  }
  clamped_scalar[0] &= 248;
  clamped_scalar[31] &= 127;
  clamped_scalar[31] |= 64;
  
  /* Perform constant-time scalar multiplication */
  scalarmult_signed_windowed(&result, &p, clamped_scalar);
  
  /* Pack result */
  group_ge_pack(ss, &result);
  
  /* Clear sensitive data */
  memset(clamped_scalar, 0, 32);
  
  return 0;
}

/* Base point scalar multiplication: k * G
 * Uses larger precomputed table for better performance
 */

/* For base point, we can use a larger table since it's fixed
 * This table would ideally be precomputed offline
 * Here we compute it on first use for simplicity
 */
#define BASE_TABLE_SIZE 64

static int base_table_initialized = 0;
static group_ge base_table[BASE_TABLE_SIZE];

static void init_base_table(void) {
  int i;
  
  if (base_table_initialized) return;
  
  base_table[0] = group_ge_neutral;
  base_table[1] = group_ge_base;
  
  /* Build table: 2G, 3G, 4G, ..., 63G */
  for (i = 2; i < BASE_TABLE_SIZE; i++) {
    if (i & 1) {
      group_ge_add(&base_table[i], &base_table[i-1], &group_ge_base);
    } else {
      group_ge_double(&base_table[i], &base_table[i/2]);
    }
  }
  
  base_table_initialized = 1;
}

int crypto_scalarmult_base(unsigned char *pk, const unsigned char *sk) {
  group_ge result, temp;
  unsigned char clamped_scalar[32];
  int i, j;
  unsigned char window;
  
  /* Initialize base point table if needed */
  init_base_table();
  
  /* Clamp scalar */
  for (i = 0; i < 32; i++) {
    clamped_scalar[i] = sk[i];
  }
  clamped_scalar[0] &= 248;
  clamped_scalar[31] &= 127;
  clamped_scalar[31] |= 64;
  
  /* Scalar multiplication using windowed method with base table */
  result = group_ge_neutral;
  
  /* Process 6 bits at a time (64 = 2^6 table entries) */
  for (i = 252; i >= 0; i -= 6) {
    /* Double 6 times */
    for (j = 0; j < 6; j++) {
      group_ge_double(&result, &result);
    }
    
    /* Extract 6-bit window */
    int byte_idx = i / 8;
    int bit_idx = i % 8;
    
    if (bit_idx >= 5) {
      window = (clamped_scalar[byte_idx] >> (bit_idx - 5)) & 0x3F;
    } else {
      window = (clamped_scalar[byte_idx] << (5 - bit_idx));
      if (byte_idx + 1 < 32) {
        window |= (clamped_scalar[byte_idx + 1] >> (8 - 5 + bit_idx));
      }
      window &= 0x3F;
    }
    
    /* Constant-time table lookup */
    table_select(&temp, base_table, window, BASE_TABLE_SIZE);
    
    /* Add to result */
    group_ge_add(&result, &result, &temp);
  }
  
  /* Pack result */
  group_ge_pack(pk, &result);
  
  /* Clear sensitive data */
  memset(clamped_scalar, 0, 32);
  
  return 0;
}
