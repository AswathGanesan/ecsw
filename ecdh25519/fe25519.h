#ifndef FE25519_H
#define FE25519_H

#include <stdint.h>

typedef struct 
{
  uint32_t v[32]; 
}
fe25519;

typedef struct 
{
  uint32_t v[8]; 
}
fe25519_8;

extern const fe25519 fe25519_zero;
extern const uint32_t fe25519_zero_8[8];
extern const fe25519 fe25519_one;
extern const uint32_t fe25519_one_8[8];
extern const fe25519 fe25519_two;

extern const fe25519 fe25519_sqrtm1;
extern const fe25519 fe25519_msqrtm1;
extern const fe25519 fe25519_m1;

void fe25519_convert_to_8(const fe25519 *in, uint32_t out[8]);
void fe25519_convert_to_32(const uint32_t in[8], fe25519 *out);

void fe25519_freeze(fe25519 *r);
void fe25519_freeze_8(uint32_t r_8[8]);

void fe25519_unpack(fe25519 *r, const unsigned char x[32]);
void fe25519_unpack_8(uint32_t r_8[8], const unsigned char x[32]);

void fe25519_pack(unsigned char r[32], const fe25519 *x);
void fe25519_pack_8(unsigned char r[32], const uint32_t x_8[8]);

int fe25519_iszero(const fe25519 *x);
int fe25519_iszero_8(const uint32_t x_8[8]);

int fe25519_isone(const fe25519 *x);
int fe25519_isone_8(const uint32_t x_8[8]);

int fe25519_isnegative(const fe25519 *x);
int fe25519_isnegative_8(const uint32_t x_8[8]);

int fe25519_iseq(const fe25519 *x, const fe25519 *y);
int fe25519_iseq_8(const uint32_t x_8[8], const uint32_t y_8[8]);

int fe25519_iseq_vartime(const fe25519 *x, const fe25519 *y);
int fe25519_iseq_vartime_8(const uint32_t x_8[8], const uint32_t y_8[8]);

void fe25519_cmov(fe25519 *r, const fe25519 *x, unsigned char b);
void fe25519_cmov_8(uint32_t r_8[8], const uint32_t x_8[8], unsigned char b);

void fe25519_setone(fe25519 *r);
void fe25519_setone_8(uint32_t r_8[8]);

void fe25519_setzero(fe25519 *r);
void fe25519_setzero_8(uint32_t r_8[8]);

void fe25519_neg(fe25519 *r, const fe25519 *x);
void fe25519_neg_8(uint32_t r_8[8], const uint32_t x_8[8]);

unsigned char fe25519_getparity(const fe25519 *x);
unsigned char fe25519_getparity_8(const uint32_t x_8[8]);

void fe25519_add(fe25519 *r, const fe25519 *x, const fe25519 *y);
void fe25519_add_8(uint32_t r_8[8], const uint32_t x_8[8], const uint32_t y_8[8]);

void fe25519_double(fe25519 *r, const fe25519 *x);
void fe25519_double_8(uint32_t r_8[8], const uint32_t x_8[8]);

void fe25519_triple(fe25519 *r, const fe25519 *x);
void fe25519_triple_8(uint32_t r_8[8], const uint32_t x_8[8]);

void fe25519_sub(fe25519 *r, const fe25519 *x, const fe25519 *y);
void fe25519_sub_8(uint32_t r_8[8], const uint32_t x_8[8], const uint32_t y_8[8]);

void fe25519_mul(fe25519 *r, const fe25519 *x, const fe25519 *y);
void fe25519_mul_8(uint32_t r_8[8], const uint32_t x_8[8], const uint32_t y_8[8]);

void fe25519_square(fe25519 *r, const fe25519 *x);
void fe25519_square_8(uint32_t r_8[8], const uint32_t x_8[8]);

void fe25519_invert(fe25519 *r, const fe25519 *x);
void fe25519_invert_8(uint32_t r_8[8], const uint32_t x_8[8]);

void fe25519_pow2523(fe25519 *r, const fe25519 *x);
void fe25519_pow2523_8(uint32_t r_8[8], const uint32_t x_8[8]);

void fe25519_invsqrt(fe25519 *r, const fe25519 *x);
void fe25519_invsqrt_8(uint32_t r_8[8], const uint32_t x_8[8]);

void fe25519_print(const fe25519 *x);
void fe25519_print_8(const uint32_t x_8[8]);

#endif
