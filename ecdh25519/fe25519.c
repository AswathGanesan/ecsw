#include <stdio.h>
#include "fe25519.h"

#if 0
#define WINDOWSIZE 1 /* Should be 1,2, or 4 */
#define WINDOWMASK ((1<<WINDOWSIZE)-1)
#endif

const fe25519 fe25519_zero = {{0,0,0,0,0,0,0,0,0,0}};
const fe25519 fe25519_one  = {{1,0,0,0,0,0,0,0,0,0}};
const fe25519 fe25519_two  = {{2,0,0,0,0,0,0,0,0,0}};

/* sqrt(-1) */
const fe25519 fe25519_sqrtm1 = {{0x20ea0b0, 0x186c9d2, 0x8f189d7, 0x035697f, 0x0bd0c60, 
                                 0x1fbd7a7, 0x2804c9e, 0x1e16569, 0x004fc1d, 0x0ae0c92}};
/* -sqrt(-1) */
const fe25519 fe25519_msqrtm1 = {{ 32595773,
  7943725,
  57730914,
  30054016,
  54719391,
  272472,
  25146209,
  2005654,
  66782178,
  22147949 }};

/* -1 */
const fe25519 fe25519_m1 = {{ 67108844,
  33554431,
  67108863,
  33554431,
  67108863,
  33554431,
  67108863,
  33554431,
  67108863,
  33554431 }}


static uint32_t equal(uint32_t a,uint32_t b) /* 16-bit inputs */
{
  uint32_t x = a ^ b; /* 0: yes; 1..65535: no */
  x -= 1; /* 4294967295: yes; 0..65534: no */
  x >>= 31; /* 1: yes; 0: no */
  return x;
}

static uint32_t ge(uint32_t a,uint32_t b) /* 16-bit inputs */
{
  uint32_t x = a;
  x -= (uint32_t) b; /* 0..65535: yes; 4294901761..4294967295: no */
  x >>= 31; /* 0: yes; 1: no */
  x ^= 1; /* 1: yes; 0: no */
  return x;
}

static uint32_t times19(uint32_t a)
{
  return (a << 4) + (a << 1) + a;
}

static uint32_t times38(uint32_t a)
{
  return (a << 5) + (a << 2) + (a << 1);
}

static void reduce_add_sub(fe25519 *r)
{
  uint32_t t;
  int i,rep;

  for(rep=0;rep<2;rep++)
  {
    t = r->v[31] >> 7;
    r->v[31] &= 127;
    t = times19(t);
    r->v[0] += t;
    for(i=0;i<31;i++)
    {
      t = r->v[i] >> 8;
      r->v[i+1] += t;
      r->v[i] &= 255;
    }
  }
}

static void reduce_mul(fe25519 *r)
{
  uint32_t t;
  int i,rep;

  for(rep=0;rep<3;rep++)
  {
    t = r->v[31] >> 7;
    r->v[31] &= 127;
    t = times19(t);
    r->v[0] += t;
    for(i=0;i<31;i++)
    {
      t = r->v[i] >> 8;
      r->v[i+1] += t;
      r->v[i] &= 255;
    }
  }
}

/* reduction modulo 2^255-19 */
void fe25519_freeze(fe25519 *r) 
{
  int i;
  uint32_t m = equal(r->v[31],127);
  for(i=30;i>0;i--)
    m &= equal(r->v[i],255);
  m &= ge(r->v[0],237);

  m = -m;

  r->v[31] -= m&127;
  for(i=30;i>0;i--)
    r->v[i] -= m&255;
  r->v[0] -= m&237;
}

void fe25519_unpack(fe25519 *r, const unsigned char x[32])
{
  // Unpack 32 bytes into 10 limbs of 25.5 bits
  uint32_t in[8];
  int i;
  
  for(i=0; i<8; i++) {
    in[i] = ((uint32_t)x[4*i]) | ((uint32_t)x[4*i+1] << 8) | 
            ((uint32_t)x[4*i+2] << 16) | ((uint32_t)x[4*i+3] << 24);
  }
  
  r->v[0] = in[0] & 0x3ffffff;
  r->v[1] = (in[0] >> 26) | ((in[1] & 0x7ffff) << 6);
  r->v[2] = (in[1] >> 19) | ((in[2] & 0xfff) << 13);
  r->v[3] = (in[2] >> 12) | ((in[3] & 0x1f) << 20);
  r->v[4] = (in[3] >> 6) & 0x3ffffff;
  r->v[5] = (in[3] >> 32) | ((in[4] & 0x1ffff) << 0);
  r->v[6] = (in[4] >> 25) | ((in[5] & 0xff) << 7);
  r->v[7] = (in[5] >> 18) | ((in[6] & 0x7) << 14);
  r->v[8] = (in[6] >> 11) | ((in[7] & 0x3fff) << 21);
  r->v[9] = (in[7] >> 14) & 0x1ffffff;
  
  r->v[9] &= 0x1ffffff; // Ensure top bit clear
}

/* Assumes input x being reduced below 2^255 */
void fe25519_pack(unsigned char r[32], const fe25519 *x)
{
  fe25519 y = *x;
  uint32_t out[8];
  int i;
  
  // Carry propagation
  uint32_t carry = 0;
  for(i=0; i<10; i++) {
    y.v[i] += carry;
    carry = y.v[i] >> ((i & 1) ? 25 : 26);
    y.v[i] &= (i & 1) ? 0x1ffffff : 0x3ffffff;
  }
  y.v[0] += carry * 19;
  
  // Final reduction
  carry = y.v[0] >> 26;
  y.v[0] &= 0x3ffffff;
  y.v[1] += carry;
  
  // Pack into bytes
  out[0] = y.v[0] | (y.v[1] << 26);
  out[1] = (y.v[1] >> 6) | (y.v[2] << 19);
  out[2] = (y.v[2] >> 13) | (y.v[3] << 12);
  out[3] = (y.v[3] >> 20) | (y.v[4] << 6);
  out[4] = (y.v[4] >> 26) | (y.v[5] << 0) | (y.v[6] << 25);
  out[5] = (y.v[6] >> 7) | (y.v[7] << 18);
  out[6] = (y.v[7] >> 14) | (y.v[8] << 11);
  out[7] = (y.v[8] >> 21) | (y.v[9] << 4);
  
  for(i=0; i<8; i++) {
    r[4*i+0] = out[i] & 0xff;
    r[4*i+1] = (out[i] >> 8) & 0xff;
    r[4*i+2] = (out[i] >> 16) & 0xff;
    r[4*i+3] = (out[i] >> 24) & 0xff;
  }
}

int fe25519_iszero(const fe25519 *x)
{
  return fe25519_iseq(x, &fe25519_zero);
}

int fe25519_isone(const fe25519 *x)
{
  return fe25519_iseq(x, &fe25519_one);
}

// return true if x has LSB set
int fe25519_isnegative(const fe25519 *x)
{
  fe25519 t = *x;
  
  fe25519_freeze(&t);

  return t.v[0] & 1;
}

int fe25519_iseq(const fe25519 *x, const fe25519 *y)
{
  fe25519 t1,t2;
  int i,r=0;

  t1 = *x;
  t2 = *y;
  fe25519_freeze(&t1);
  fe25519_freeze(&t2);
  for(i=0;i<32;i++)
    r |= (1-equal(t1.v[i],t2.v[i]));
  return 1-r;
}

void fe25519_cmov(fe25519 *r, const fe25519 *x, unsigned char b)
{
  uint32_t mask = (uint32_t)(-(int32_t)b);
  int i;
  for(i = 0; i < 32; i++) {
    r->v[i] ^= mask & (r->v[i] ^ x->v[i]);
  }
}

void fe25519_neg(fe25519 *r, const fe25519 *x)
{
  fe25519 t = fe25519_zero;
  fe25519_sub(r, &t, x);
}

void fe25519_add(fe25519 *r, const fe25519 *x, const fe25519 *y)
{
  int i;
  for(i=0;i<10;i++) r->v[i] = x->v[i] + y->v[i];
}
void fe25519_double(fe25519 *r, const fe25519 *x)
{
  int i;
  for(i=0;i<32;i++) r->v[i] = 2*x->v[i];
  reduce_add_sub(r);
}

void fe25519_sub(fe25519 *r, const fe25519 *x, const fe25519 *y)
{
  // Add 2*p before subtracting to avoid underflow
  r->v[0] = x->v[0] + 0x7ffffda - y->v[0];
  r->v[1] = x->v[1] + 0x3fffffe - y->v[1];
  r->v[2] = x->v[2] + 0x7fffffe - y->v[2];
  r->v[3] = x->v[3] + 0x3fffffe - y->v[3];
  r->v[4] = x->v[4] + 0x7fffffe - y->v[4];
  r->v[5] = x->v[5] + 0x3fffffe - y->v[5];
  r->v[6] = x->v[6] + 0x7fffffe - y->v[6];
  r->v[7] = x->v[7] + 0x3fffffe - y->v[7];
  r->v[8] = x->v[8] + 0x7fffffe - y->v[8];
  r->v[9] = x->v[9] + 0x3fffffe - y->v[9];
}
void fe25519_mul(fe25519 *r, const fe25519 *x, const fe25519 *y)
{
  uint32_t x0=x->v[0], x1=x->v[1], x2=x->v[2], x3=x->v[3], x4=x->v[4];
  uint32_t x5=x->v[5], x6=x->v[6], x7=x->v[7], x8=x->v[8], x9=x->v[9];
  uint32_t y0=y->v[0], y1=y->v[1], y2=y->v[2], y3=y->v[3], y4=y->v[4];
  uint32_t y5=y->v[5], y6=y->v[6], y7=y->v[7], y8=y->v[8], y9=y->v[9];
  
  uint64_t t0, t1, t2, t3, t4, t5, t6, t7, t8, t9;
  uint32_t r0, r1, r2, r3, r4, r5, r6, r7, r8, r9;
  
  // Compute products
  t0 = ((uint64_t)x0 * y0);
  t1 = ((uint64_t)x0 * y1) + ((uint64_t)x1 * y0);
  t2 = ((uint64_t)x0 * y2) + ((uint64_t)x1 * y1) + ((uint64_t)x2 * y0);
  t3 = ((uint64_t)x0 * y3) + ((uint64_t)x1 * y2) + ((uint64_t)x2 * y1) + ((uint64_t)x3 * y0);
  t4 = ((uint64_t)x0 * y4) + ((uint64_t)x1 * y3) + ((uint64_t)x2 * y2) + ((uint64_t)x3 * y1) + ((uint64_t)x4 * y0);
  t5 = ((uint64_t)x1 * y4) + ((uint64_t)x2 * y3) + ((uint64_t)x3 * y2) + ((uint64_t)x4 * y1) + ((uint64_t)x5 * y0);
  t6 = ((uint64_t)x2 * y4) + ((uint64_t)x3 * y3) + ((uint64_t)x4 * y2) + ((uint64_t)x5 * y1) + ((uint64_t)x6 * y0);
  t7 = ((uint64_t)x3 * y4) + ((uint64_t)x4 * y3) + ((uint64_t)x5 * y2) + ((uint64_t)x6 * y1) + ((uint64_t)x7 * y0);
  t8 = ((uint64_t)x4 * y4) + ((uint64_t)x5 * y3) + ((uint64_t)x6 * y2) + ((uint64_t)x7 * y1) + ((uint64_t)x8 * y0);
  t9 = ((uint64_t)x5 * y4) + ((uint64_t)x6 * y3) + ((uint64_t)x7 * y2) + ((uint64_t)x8 * y1) + ((uint64_t)x9 * y0);
  
  // Add high terms * 38
  t0 += ((uint64_t)38) * (((uint64_t)x5 * y5) + ((uint64_t)x6 * y4) + ((uint64_t)x7 * y3) + ((uint64_t)x8 * y2) + ((uint64_t)x9 * y1));
  t1 += ((uint64_t)38) * (((uint64_t)x6 * y5) + ((uint64_t)x7 * y4) + ((uint64_t)x8 * y3) + ((uint64_t)x9 * y2));
  t2 += ((uint64_t)38) * (((uint64_t)x7 * y5) + ((uint64_t)x8 * y4) + ((uint64_t)x9 * y3));
  t3 += ((uint64_t)38) * (((uint64_t)x8 * y5) + ((uint64_t)x9 * y4));
  t4 += ((uint64_t)38) * ((uint64_t)x9 * y5);
  
  t0 += ((uint64_t)38) * (((uint64_t)x0 * y5) + ((uint64_t)x1 * y6) + ((uint64_t)x2 * y7) + ((uint64_t)x3 * y8) + ((uint64_t)x4 * y9));
  t1 += ((uint64_t)38) * (((uint64_t)x0 * y6) + ((uint64_t)x1 * y7) + ((uint64_t)x2 * y8) + ((uint64_t)x3 * y9));
  t2 += ((uint64_t)38) * (((uint64_t)x0 * y7) + ((uint64_t)x1 * y8) + ((uint64_t)x2 * y9));
  t3 += ((uint64_t)38) * (((uint64_t)x0 * y8) + ((uint64_t)x1 * y9));
  t4 += ((uint64_t)38) * ((uint64_t)x0 * y9);
  
  // Carry propagation
  r0 = (uint32_t)t0 & 0x3ffffff; t1 += t0 >> 26;
  r1 = (uint32_t)t1 & 0x1ffffff; t2 += t1 >> 25;
  r2 = (uint32_t)t2 & 0x3ffffff; t3 += t2 >> 26;
  r3 = (uint32_t)t3 & 0x1ffffff; t4 += t3 >> 25;
  r4 = (uint32_t)t4 & 0x3ffffff; t5 += t4 >> 26;
  r5 = (uint32_t)t5 & 0x1ffffff; t6 += t5 >> 25;
  r6 = (uint32_t)t6 & 0x3ffffff; t7 += t6 >> 26;
  r7 = (uint32_t)t7 & 0x1ffffff; t8 += t7 >> 25;
  r8 = (uint32_t)t8 & 0x3ffffff; t9 += t8 >> 26;
  r9 = (uint32_t)t9 & 0x1ffffff;
  
  r0 += (uint32_t)((t9 >> 25) * 19);
  r1 += r0 >> 26; r0 &= 0x3ffffff;
  
  r->v[0]=r0; r->v[1]=r1; r->v[2]=r2; r->v[3]=r3; r->v[4]=r4;
  r->v[5]=r5; r->v[6]=r6; r->v[7]=r7; r->v[8]=r8; r->v[9]=r9;
}

void fe25519_square(fe25519 *r, const fe25519 *x)
{
  fe25519_mul(r, x, x);
}

#if 0
void fe25519_invert(fe25519 *r, const fe25519 *x)
{
	fe25519 z2;
	fe25519 z9;
	fe25519 z11;
	fe25519 z2_5_0;
	fe25519 z2_10_0;
	fe25519 z2_20_0;
	fe25519 z2_50_0;
	fe25519 z2_100_0;
	fe25519 t0;
	fe25519 t1;
	int i;
	
	/* 2 */ fe25519_square(&z2,x);
	/* 4 */ fe25519_square(&t1,&z2);
	/* 8 */ fe25519_square(&t0,&t1);
	/* 9 */ fe25519_mul(&z9,&t0,x);
	/* 11 */ fe25519_mul(&z11,&z9,&z2);
	/* 22 */ fe25519_square(&t0,&z11);
	/* 2^5 - 2^0 = 31 */ fe25519_mul(&z2_5_0,&t0,&z9);

	/* 2^6 - 2^1 */ fe25519_square(&t0,&z2_5_0);
	/* 2^7 - 2^2 */ fe25519_square(&t1,&t0);
	/* 2^8 - 2^3 */ fe25519_square(&t0,&t1);
	/* 2^9 - 2^4 */ fe25519_square(&t1,&t0);
	/* 2^10 - 2^5 */ fe25519_square(&t0,&t1);
	/* 2^10 - 2^0 */ fe25519_mul(&z2_10_0,&t0,&z2_5_0);

	/* 2^11 - 2^1 */ fe25519_square(&t0,&z2_10_0);
	/* 2^12 - 2^2 */ fe25519_square(&t1,&t0);
	/* 2^20 - 2^10 */ for (i = 2;i < 10;i += 2) { fe25519_square(&t0,&t1); fe25519_square(&t1,&t0); }
	/* 2^20 - 2^0 */ fe25519_mul(&z2_20_0,&t1,&z2_10_0);

	/* 2^21 - 2^1 */ fe25519_square(&t0,&z2_20_0);
	/* 2^22 - 2^2 */ fe25519_square(&t1,&t0);
	/* 2^40 - 2^20 */ for (i = 2;i < 20;i += 2) { fe25519_square(&t0,&t1); fe25519_square(&t1,&t0); }
	/* 2^40 - 2^0 */ fe25519_mul(&t0,&t1,&z2_20_0);

	/* 2^41 - 2^1 */ fe25519_square(&t1,&t0);
	/* 2^42 - 2^2 */ fe25519_square(&t0,&t1);
	/* 2^50 - 2^10 */ for (i = 2;i < 10;i += 2) { fe25519_square(&t1,&t0); fe25519_square(&t0,&t1); }
	/* 2^50 - 2^0 */ fe25519_mul(&z2_50_0,&t0,&z2_10_0);

	/* 2^51 - 2^1 */ fe25519_square(&t0,&z2_50_0);
	/* 2^52 - 2^2 */ fe25519_square(&t1,&t0);
	/* 2^100 - 2^50 */ for (i = 2;i < 50;i += 2) { fe25519_square(&t0,&t1); fe25519_square(&t1,&t0); }
	/* 2^100 - 2^0 */ fe25519_mul(&z2_100_0,&t1,&z2_50_0);

	/* 2^101 - 2^1 */ fe25519_square(&t1,&z2_100_0);
	/* 2^102 - 2^2 */ fe25519_square(&t0,&t1);
	/* 2^200 - 2^100 */ for (i = 2;i < 100;i += 2) { fe25519_square(&t1,&t0); fe25519_square(&t0,&t1); }
	/* 2^200 - 2^0 */ fe25519_mul(&t1,&t0,&z2_100_0);

	/* 2^201 - 2^1 */ fe25519_square(&t0,&t1);
	/* 2^202 - 2^2 */ fe25519_square(&t1,&t0);
	/* 2^250 - 2^50 */ for (i = 2;i < 50;i += 2) { fe25519_square(&t0,&t1); fe25519_square(&t1,&t0); }
	/* 2^250 - 2^0 */ fe25519_mul(&t0,&t1,&z2_50_0);

	/* 2^251 - 2^1 */ fe25519_square(&t1,&t0);
	/* 2^252 - 2^2 */ fe25519_square(&t0,&t1);
	/* 2^253 - 2^3 */ fe25519_square(&t1,&t0);
	/* 2^254 - 2^4 */ fe25519_square(&t0,&t1);
	/* 2^255 - 2^5 */ fe25519_square(&t1,&t0);
	/* 2^255 - 21 */ fe25519_mul(r,&t1,&z11);
}
#endif

void fe25519_pow2523(fe25519 *r, const fe25519 *x)
{
	fe25519 z2;
	fe25519 z9;
	fe25519 z11;
	fe25519 z2_5_0;
	fe25519 z2_10_0;
	fe25519 z2_20_0;
	fe25519 z2_50_0;
	fe25519 z2_100_0;
	fe25519 t;
	int i;
		
	/* 2 */ fe25519_square(&z2,x);
	/* 4 */ fe25519_square(&t,&z2);
	/* 8 */ fe25519_square(&t,&t);
	/* 9 */ fe25519_mul(&z9,&t,x);
	/* 11 */ fe25519_mul(&z11,&z9,&z2);
	/* 22 */ fe25519_square(&t,&z11);
	/* 2^5 - 2^0 = 31 */ fe25519_mul(&z2_5_0,&t,&z9);

	/* 2^6 - 2^1 */ fe25519_square(&t,&z2_5_0);
	/* 2^10 - 2^5 */ for (i = 1;i < 5;i++) { fe25519_square(&t,&t); }
	/* 2^10 - 2^0 */ fe25519_mul(&z2_10_0,&t,&z2_5_0);

	/* 2^11 - 2^1 */ fe25519_square(&t,&z2_10_0);
	/* 2^20 - 2^10 */ for (i = 1;i < 10;i++) { fe25519_square(&t,&t); }
	/* 2^20 - 2^0 */ fe25519_mul(&z2_20_0,&t,&z2_10_0);

	/* 2^21 - 2^1 */ fe25519_square(&t,&z2_20_0);
	/* 2^40 - 2^20 */ for (i = 1;i < 20;i++) { fe25519_square(&t,&t); }
	/* 2^40 - 2^0 */ fe25519_mul(&t,&t,&z2_20_0);

	/* 2^41 - 2^1 */ fe25519_square(&t,&t);
	/* 2^50 - 2^10 */ for (i = 1;i < 10;i++) { fe25519_square(&t,&t); }
	/* 2^50 - 2^0 */ fe25519_mul(&z2_50_0,&t,&z2_10_0);

	/* 2^51 - 2^1 */ fe25519_square(&t,&z2_50_0);
	/* 2^100 - 2^50 */ for (i = 1;i < 50;i++) { fe25519_square(&t,&t); }
	/* 2^100 - 2^0 */ fe25519_mul(&z2_100_0,&t,&z2_50_0);

	/* 2^101 - 2^1 */ fe25519_square(&t,&z2_100_0);
	/* 2^200 - 2^100 */ for (i = 1;i < 100;i++) { fe25519_square(&t,&t); }
	/* 2^200 - 2^0 */ fe25519_mul(&t,&t,&z2_100_0);

	/* 2^201 - 2^1 */ fe25519_square(&t,&t);
	/* 2^250 - 2^50 */ for (i = 1;i < 50;i++) { fe25519_square(&t,&t); }
	/* 2^250 - 2^0 */ fe25519_mul(&t,&t,&z2_50_0);

	/* 2^251 - 2^1 */ fe25519_square(&t,&t);
	/* 2^252 - 2^2 */ fe25519_square(&t,&t);
	/* 2^252 - 3 */ fe25519_mul(r,&t,x);
}

void fe25519_invsqrt(fe25519 *r, const fe25519 *x)
{
  fe25519 den2, den3, den4, den6, chk, t, t2;
  int b;

  fe25519_square(&den2, x);
  fe25519_mul(&den3, &den2, x);

  fe25519_square(&den4, &den2);
  fe25519_mul(&den6, &den2, &den4);
  fe25519_mul(&t, &den6, x); // r is now x^7

  fe25519_pow2523(&t, &t);
  fe25519_mul(&t, &t, &den3);

  fe25519_square(&chk, &t);
  fe25519_mul(&chk, &chk, x);

  fe25519_mul(&t2, &t, &fe25519_sqrtm1);
  b = 1 - fe25519_isone(&chk);

  fe25519_cmov(&t, &t2, b);

  *r = t;
}
