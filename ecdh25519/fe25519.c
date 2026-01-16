#include <stdio.h>
#include "fe25519.h"

#if 0
#define WINDOWSIZE 1 /* Should be 1,2, or 4 */
#define WINDOWMASK ((1<<WINDOWSIZE)-1)
#endif

typedef uint32_t uint32;
typedef unsigned int uint8;

void fe25519_convert_to_8(const fe25519 *in, uint32_t out[8]) {
  uint32 cur_val = 0;
  unsigned int i;
  for(i=0; i<8; i++) {
    cur_val = ((in->v[i*4+3] << 24) | (in->v[i*4+2] << 16) | (in->v[i*4+1] << 8) | (in->v[i*4]));
    out[i] = cur_val;
  }
}
void fe25519_convert_to_32(const uint32_t in[8], fe25519 *out){
  uint32 cur_val = 0;
  unsigned int i;
  for(i=0; i<8; i++) {
    cur_val = (uint32) in[i];
    out->v[i*4] = (uint32) cur_val & 255;
    out->v[i*4+1] = (uint32) (cur_val >> 8) & 255;
    out->v[i*4+2] = (uint32) (cur_val >> 16) & 255;
    out->v[i*4+3] = (uint32) (cur_val >> 24) & 255;
  }
}

const fe25519 fe25519_zero = {{0}};
const fe25519 fe25519_one  = {{1}};
const fe25519 fe25519_two  = {{2}};

const uint32_t fe25519_zero_8[8]  = {0};
const uint32_t fe25519_one_8[8]  = {1};
const uint32_t fe25519_two_8[8]  = {2};


/* sqrt(-1) */
const fe25519 fe25519_sqrtm1 = {{0xB0, 0xA0, 0x0E, 0x4A, 0x27, 0x1B, 0xEE, 0xC4, 0x78, 0xE4, 0x2F, 0xAD, 0x06, 0x18, 0x43, 0x2F, 
                                 0xA7, 0xD7, 0xFB, 0x3D, 0x99, 0x00, 0x4D, 0x2B, 0x0B, 0xDF, 0xC1, 0x4F, 0x80, 0x24, 0x83, 0x2B}};

/* -sqrt(-1) */
const fe25519 fe25519_msqrtm1 = {{0x3D, 0x5F, 0xF1, 0xB5, 0xD8, 0xE4, 0x11, 0x3B, 0x87, 0x1B, 0xD0, 0x52, 0xF9, 0xE7, 0xBC, 0xD0, 
                                  0x58, 0x28, 0x4, 0xC2, 0x66, 0xFF, 0xB2, 0xD4, 0xF4, 0x20, 0x3E, 0xB0, 0x7F, 0xDB, 0x7C, 0x54}};

uint32 fe25519_msqrtm1_8[8] = {0xB5F15F3D};

/* -1 */
const fe25519 fe25519_m1 = {{236, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 
                             0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7f}};



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

static void reduce_add_sub_8(uint32 r_8[8], uint64_t u) {
  uint32 bitmask = 0xFFFFFFFF;

  for (int i = 0; i< 2; i++) {
    u >>= 31;
    r_8[7] &= 0x7FFFFFFF;
    u = u * 19;
    u += r_8[0];
    r_8[0] = (uint32) u & bitmask;
    for (int j = 0; j < 7; j++)
    {
      u >>= 32;
      u += r_8[j+1];
      r_8[j+1] = (u & bitmask);
    }
  }
}

/* reduction modulo 2^255-19 */
void fe25519_freeze_8(uint32 r_8[8]) {
  fe25519 r;
  fe25519_convert_to_32(r_8, &r);
  fe25519_freeze(&r);
  fe25519_convert_to_8(&r, r_8);
} 

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

void fe25519_unpack_8(uint32 r_8[8], const unsigned char x[32]) {
  fe25519 r;
  fe25519_convert_to_32(r_8, &r);
  fe25519_unpack(&r, x);
  fe25519_convert_to_8(&r, r_8);
} 

void fe25519_unpack(fe25519 *r, const unsigned char x[32])
{
  int i;
  for(i=0;i<32;i++) {
    r->v[i] = x[i];
  }  
  r->v[31] &= 127;
}

void fe25519_pack_8(unsigned char r[32], const uint32 x_8[8]) {
  fe25519 x;
  fe25519_convert_to_32(x_8, &x);
  fe25519_pack(r, &x);
} 

/* Assumes input x being reduced below 2^255 */
void fe25519_pack(unsigned char r[32], const fe25519 *x)
{
  int i;
  fe25519 y = *x;
  fe25519_freeze(&y);
  for(i=0;i<32;i++) 
    r[i] = y.v[i];
}

int fe25519_iszero_8(const uint32 x_8[8])
{
  return fe25519_iseq_8(x_8, fe25519_zero_8);
}

int fe25519_iszero(const fe25519 *x)
{
  return fe25519_iseq(x, &fe25519_zero);
}

int fe25519_isone_8(const uint32 x_8[8])
{
  return fe25519_iseq_8(x_8, fe25519_one_8);
}

int fe25519_isone(const fe25519 *x)
{
  return fe25519_iseq(x, &fe25519_one);
}

int fe25519_isnegative_8(const uint32 x_8[8]) {
  fe25519 x;
  fe25519_convert_to_32(x_8, &x);
  return fe25519_isnegative(&x);
} 

// return true if x has LSB set
int fe25519_isnegative(const fe25519 *x)
{
  fe25519 t = *x;
  
  fe25519_freeze(&t);

  return t.v[0] & 1;
}

int fe25519_iseq_8(const uint32 x_8[8], const uint32 y_8[8])
{
  fe25519 t1,t2;
  int i,r=0;

  fe25519_convert_to_32(x_8, &t1);
  fe25519_convert_to_32(y_8, &t2);

  // t1 = *x;
  // t2 = *y;
  fe25519_freeze(&t1);
  fe25519_freeze(&t2);
  for(i=0;i<32;i++)
    r |= (1-equal(t1.v[i],t2.v[i]));
  return 1-r;
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

  for (int i =0; i< 32; i++) {
    r->v[i] = x->v[i]*b + r->v[i]*(1-b);
  }
  
}

void fe25519_cmov_8(uint32 r_8[8], const uint32 x_8[8], unsigned char b)
{
  // if(b) *r = *x;

  for (int i =0; i< 8; i++) {
    r_8[i] = x_8[i]*b + r_8[i]*(1-b);
  }
}

void fe25519_neg_8(uint32 r_8[8], const uint32 x_8[8])
{
  // uint32 t_8[8] = fe25519_zero_8;
  fe25519_sub_8(r_8, fe25519_zero_8, x_8);
}

void fe25519_neg(fe25519 *r, const fe25519 *x)
{
  fe25519 t = fe25519_zero;
  fe25519_sub(r, &t, x);
}

void fe25519_add(fe25519 *r, const fe25519 *x, const fe25519 *y)
{
  int i;
  for(i=0;i<32;i++) {
    r->v[i] = x->v[i] + y->v[i];
  }
  reduce_add_sub(r);
}

void fe25519_add_8(uint32 r_8[8], const uint32 x_8[8], const uint32 y_8[8])
{
  uint64_t u;
  uint32 bitmask = 0xFFFFFFFF;

  u = 0;
  for (int j = 0; j < 8; j++)
  {
    u >>= 32;
    u += (uint64_t) x_8[j] + y_8[j];
    r_8[j] = (uint32) u & bitmask;
  }

  reduce_add_sub_8(r_8, u);
}

void fe25519_double_8(uint32 r_8[8], const uint32 x_8[8])
{
  fe25519_add_8(r_8, x_8, x_8);
}

void fe25519_double(fe25519 *r, const fe25519 *x)
{
  fe25519_add(r, x, x);
}

void fe25519_sub_8(uint32 r_8[8], const uint32 x_8[8], const uint32 y_8[8])
{
  fe25519 r, x, y;
  
  fe25519_convert_to_32(x_8, &x);
  fe25519_convert_to_32(y_8, &y);
  fe25519_convert_to_32(r_8, &r);

  fe25519_sub(&r, &x, &y);

  fe25519_convert_to_8(&r, r_8);
}

void fe25519_sub(fe25519 *r, const fe25519 *x, const fe25519 *y)
{
  int i;
  
  uint32_t t[32];
  t[0] = x->v[0] + 0x1da;
  // add 000111011010
  t[31] = x->v[31] + 0xfe;
  // add 11111110
  for(i=1;i<31;i++) t[i] = x->v[i] + 0x1fe;
  // add 000111111110
  for(i=0;i<32;i++) r->v[i] = t[i] - y->v[i];
  reduce_add_sub(r);
}

void fe25519_mul(fe25519 *r, const fe25519 *x, const fe25519 *y){
  // fe25519_mul_old(r, x, y);
  
  uint32 r_8[8];
  uint32 x_8[8];
  uint32 y_8[8];

  fe25519_convert_to_8(x, x_8);
  fe25519_convert_to_8(y, y_8);
  fe25519_mul_8(r_8, x_8, y_8);

  fe25519_convert_to_32(r_8, r);

}

void fe25519_mul_8(uint32 r_8[8], const uint32 x_8[8], const uint32 y_8[8]){
  unsigned int i, j;
  uint32_t u_non38[8];
  uint32_t u_38[8];
  uint64_t carry = 0;
  uint64_t u = 0;

  // We want to do the following:
  // hr[0] = x_8[0]*y_8[0] + 38*x_8[1]*y_8[7] + 38*x_8[2]*y_8[6] + 38*x_8[3]*y_8[5] + 38*x_8[4]*y_8[4] + 38*x_8[5]*y_8[3] + 38*x_8[6]*y_8[2] + 38*x_8[7]*y_8[1];   
  // hr[1] = x_8[0]*y_8[1] +    x_8[1]*y_8[0] + 38*x_8[2]*y_8[7] + 38*x_8[3]*y_8[6] + 38*x_8[4]*y_8[5] + 38*x_8[5]*y_8[4] + 38*x_8[6]*y_8[3] + 38*x_8[7]*y_8[2];  
  // hr[2] = x_8[0]*y_8[2] +    x_8[1]*y_8[1] +    x_8[2]*y_8[0] + 38*x_8[3]*y_8[7] + 38*x_8[4]*y_8[6] + 38*x_8[5]*y_8[5] + 38*x_8[6]*y_8[4] + 38*x_8[7]*y_8[3]; 
  // hr[3] = x_8[0]*y_8[3] +    x_8[1]*y_8[2] +    x_8[2]*y_8[1] +    x_8[3]*y_8[0] + 38*x_8[4]*y_8[7] + 38*x_8[5]*y_8[6] + 38*x_8[6]*y_8[5] + 38*x_8[7]*y_8[4];  
  // hr[4] = x_8[0]*y_8[4] +    x_8[1]*y_8[3] +    x_8[2]*y_8[2] +    x_8[3]*y_8[1] +    x_8[4]*y_8[0] + 38*x_8[5]*y_8[7] + 38*x_8[6]*y_8[6] + 38*x_8[7]*y_8[5];  
  // hr[5] = x_8[0]*y_8[5] +    x_8[1]*y_8[4] +    x_8[2]*y_8[3] +    x_8[3]*y_8[2] +    x_8[4]*y_8[1] +    x_8[5]*y_8[0] + 38*x_8[6]*y_8[7] + 38*x_8[7]*y_8[6];  
  // hr[6] = x_8[0]*y_8[6] +    x_8[1]*y_8[5] +    x_8[2]*y_8[4] +    x_8[3]*y_8[3] +    x_8[4]*y_8[2] +    x_8[5]*y_8[1] +    x_8[6]*y_8[0] + 38*x_8[7]*y_8[7];  
  // hr[7] = x_8[0]*y_8[7] +    x_8[1]*y_8[6] +    x_8[2]*y_8[5] +    x_8[3]*y_8[4] +    x_8[4]*y_8[3] +    x_8[5]*y_8[2] +    x_8[6]*y_8[1] +    x_8[7]*y_8[0];  
  // However, there are problems due to x_8[0] * y_8[0] being the multiplication of two 32-bit integers
  // This means that the result of only this part might already fill the entire 64 bits which is why it cannot be done like in the poly version

  // We decided to fix this by doing it in many parts such that the carry can be handled in many different steps as well
  // We first calculate all the parts that are not multiplied with 38, then we calculate it for all the parts that are multiplied by 38
  // Each time, the carry is alraedy being handled
  // Afterwards, we add the two parts together while also handling the carry again.

  // Multiplication for 'non 38-part':
  for (i = 0; i < 8; i++) {
    for (j = 0; j <= i; j++) {
      u += (uint64_t) x_8[j] * y_8[i - j];
      carry += u >> 32;
      u &= 0xFFFFFFFF;
    }
    u_non38[i] = u & 0xFFFFFFFF;
    u = carry & 0xFFFFFFFF;
    carry >>= 32;
  }

  // Multiplication '38-part':
  for (i = 8; i < 16; i++) {
    for (j = i - 7; j < 8; j++) {
      u += (uint64_t) x_8[j] * y_8[i - j];
      carry += u >> 32;
      u &= 0xFFFFFFFF;
    }
    u_38[i-8] = u & 0xFFFFFFFF;
    u = carry & 0xFFFFFFFF;
    carry >>= 32;
  }

  // u = (u_non38[7] >> 31) * 19;
  // u_non38[7] &= 0x7FFFFFFF;
  
  // Now we add the '38 part' to the 'non-38' part
  // This requires actually multiplying the '38 part' by 38
  for(i = 0; i < 8; i++)
  {
      u += u_non38[i];
      u += (uint64_t) u_38[i] * 38;
      u_non38[i] = u & 0xFFFFFFFF;
      u >>= 32;
  }

  // Handle final carry
  u = 38 * u;

  // u += (u_non38[7] >> 31) * 19;
  // u_non38[7] &= 0x7FFFFFFF;
  
  for(i = 0; i < 8; i++)
  {
      u += u_non38[i];
      u_non38[i] = u & 0xFFFFFFFF;
      u >>= 32;
  }

  // Perform the final reduction.

  // Compute r_8 = u_non38 + 19 when MSB of u_non38 is set
  // Idea copied from original reduce_mul function
  // We do not really know why the +19 is needed sometimes
  u = u_non38[7] >> 31;
  u_non38[7] &= 0x7FFFFFFF;
  u *= 19;
  for(i = 0; i < 8; i++)
  {
      u += u_non38[i];
      r_8[i] = u & 0xFFFFFFFF;
      u >>= 32;
  }
}

void fe25519_square(fe25519 *r, const fe25519 *x)
{
  fe25519_mul(r, x, x);
}

void fe25519_square_8(uint32 r_8[8], const uint32 x_8[8])
{
  fe25519_mul_8(r_8, x_8, x_8);
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
  uint32 x_8[8];
  uint32 r_8[8];

  fe25519_convert_to_8(x, x_8);
  fe25519_convert_to_8(r, r_8);
  fe25519_pow2523_8(r_8, x_8);
  fe25519_convert_to_32(r_8, r);
}

void fe25519_pow2523_8(uint32 r_8[8], const uint32 x_8[8])
{

	uint32 z2[8];
	uint32 z9[8];
	uint32 z11[8];
	uint32 z2_5_0[8];
	uint32 z2_10_0[8];
	uint32 z2_20_0[8];
	uint32 z2_50_0[8];
	uint32 z2_100_0[8];
	uint32 t[8];
	int i;
		
	/* 2 */ fe25519_square_8(z2,x_8);
	/* 4 */ fe25519_square_8(t,z2);
	/* 8 */ fe25519_square_8(t,t);
	/* 9 */ fe25519_mul_8(z9,t,x_8);
	/* 11 */ fe25519_mul_8(z11,z9,z2);
	/* 22 */ fe25519_square_8(t,z11);
	/* 2^5 - 2^0 = 31 */ fe25519_mul_8(z2_5_0,t,z9);

	/* 2^6 - 2^1 */ fe25519_square_8(t,z2_5_0);
	/* 2^10 - 2^5 */ for (i = 1;i < 5;i++) { fe25519_square_8(t,t); }
	/* 2^10 - 2^0 */ fe25519_mul_8(z2_10_0,t,z2_5_0);

	/* 2^11 - 2^1 */ fe25519_square_8(t,z2_10_0);
	/* 2^20 - 2^10 */ for (i = 1;i < 10;i++) { fe25519_square_8(t,t); }
	/* 2^20 - 2^0 */ fe25519_mul_8(z2_20_0,t,z2_10_0);

	/* 2^21 - 2^1 */ fe25519_square_8(t,z2_20_0);
	/* 2^40 - 2^20 */ for (i = 1;i < 20;i++) { fe25519_square_8(t,t); }
	/* 2^40 - 2^0 */ fe25519_mul_8(t,t,z2_20_0);

	/* 2^41 - 2^1 */ fe25519_square_8(t,t);
	/* 2^50 - 2^10 */ for (i = 1;i < 10;i++) { fe25519_square_8(t,t); }
	/* 2^50 - 2^0 */ fe25519_mul_8(z2_50_0,t,z2_10_0);

	/* 2^51 - 2^1 */ fe25519_square_8(t,z2_50_0);
	/* 2^100 - 2^50 */ for (i = 1;i < 50;i++) { fe25519_square_8(t,t); }
	/* 2^100 - 2^0 */ fe25519_mul_8(z2_100_0,t,z2_50_0);

	/* 2^101 - 2^1 */ fe25519_square_8(t,z2_100_0);
	/* 2^200 - 2^100 */ for (i = 1;i < 100;i++) { fe25519_square_8(t,t); }
	/* 2^200 - 2^0 */ fe25519_mul_8(t,t,z2_100_0);

	/* 2^201 - 2^1 */ fe25519_square_8(t,t);
	/* 2^250 - 2^50 */ for (i = 1;i < 50;i++) { fe25519_square_8(t,t); }
	/* 2^250 - 2^0 */ fe25519_mul_8(t,t,z2_50_0);

	/* 2^251 - 2^1 */ fe25519_square_8(t,t);
	/* 2^252 - 2^2 */ fe25519_square_8(t,t);
	/* 2^252 - 3 */ fe25519_mul_8(r_8,t,x_8);

}


void fe25519_invsqrt_8(uint32 r_8[8], const uint32 x_8[8])
{
  uint32 den2[8];
  uint32 den3[8];
  uint32 den4[8];
  uint32 den6[8];
  uint32 chk[8];
  uint32 t[8];
  uint32 t2[8];

  int b;

  fe25519_square_8(den2, x_8);
  fe25519_mul_8(den3, den2, x_8);

  fe25519_square_8(den4, den2);
  fe25519_mul_8(den6, den2, den4);
  fe25519_mul_8(t, den6, x_8); // r is now x^7

  fe25519_pow2523_8(t, t);
  fe25519_mul_8(t, t, den3);

  fe25519_square_8(chk, t);
  fe25519_mul_8(chk, chk, x_8);



  const fe25519 fe25519_sqrtm1_temp = {{0xB0, 0xA0, 0x0E, 0x4A, 0x27, 0x1B, 0xEE, 0xC4, 0x78, 0xE4, 0x2F, 0xAD, 0x06, 0x18, 0x43, 0x2F, 0xA7, 0xD7, 0xFB, 0x3D, 0x99, 0x00, 0x4D, 0x2B, 0x0B, 0xDF, 0xC1, 0x4F, 0x80, 0x24, 0x83, 0x2B}};
  uint32 fe25519_sqrtm1_8_temp[8];
  fe25519_convert_to_8(&fe25519_sqrtm1_temp, fe25519_sqrtm1_8_temp);

  fe25519_mul_8(t2, t, fe25519_sqrtm1_8_temp);
  b = 1 - fe25519_isone_8(chk);

  fe25519_cmov_8(t, t2, b);

  for (int i =0; i< 8; i++) {
    r_8[i] = t[i];
  }
}


void fe25519_invsqrt(fe25519 *r, const fe25519 *x)
{
  uint32 den2[8];
  uint32 den3[8];
  uint32 den4[8];
  uint32 den6[8];
  uint32 chk[8];
  uint32 t[8];
  uint32 t2[8];

  uint32 x_8[8];
  uint32 r_8[8];

  fe25519_convert_to_8(x, x_8);
  fe25519_convert_to_8(r, r_8);
  int b;

  fe25519_square_8(den2, x_8);
  fe25519_mul_8(den3, den2, x_8);

  fe25519_square_8(den4, den2);
  fe25519_mul_8(den6, den2, den4);
  fe25519_mul_8(t, den6, x_8); // r is now x^7

  fe25519_pow2523_8(t, t);
  fe25519_mul_8(t, t, den3);

  fe25519_square_8(chk, t);
  fe25519_mul_8(chk, chk, x_8);



  const fe25519 fe25519_sqrtm1_temp = {{0xB0, 0xA0, 0x0E, 0x4A, 0x27, 0x1B, 0xEE, 0xC4, 0x78, 0xE4, 0x2F, 0xAD, 0x06, 0x18, 0x43, 0x2F, 0xA7, 0xD7, 0xFB, 0x3D, 0x99, 0x00, 0x4D, 0x2B, 0x0B, 0xDF, 0xC1, 0x4F, 0x80, 0x24, 0x83, 0x2B}};
  uint32 fe25519_sqrtm1_8_temp[8];
  fe25519_convert_to_8(&fe25519_sqrtm1_temp, fe25519_sqrtm1_8_temp);

  fe25519_mul_8(t2, t, fe25519_sqrtm1_8_temp);
  b = 1 - fe25519_isone_8(chk);

  fe25519_cmov_8(t, t2, b);

  // *r = t;
  fe25519_convert_to_32(t, r);
}
