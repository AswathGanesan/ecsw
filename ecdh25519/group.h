#ifndef GROUP_GE_H 
#define GROUP_GE_H

#include "fe25519.h"

#define GROUP_GE_PACKEDBYTES 32

typedef struct
{	
	fe25519 x;
	fe25519 y;
	fe25519 z;
	fe25519 t;
} group_ge;

typedef struct
{	
	uint32_t x[8];
	uint32_t y[8];
	uint32_t z[8];
	uint32_t t[8];
} group_ge_8;

extern const group_ge group_ge_base;

void group_convert_to_8(const group_ge *in, group_ge_8 *out);
void group_convert_to_32(const group_ge_8 *in, group_ge *out);

int  group_ge_unpack(group_ge *r, const unsigned char x[GROUP_GE_PACKEDBYTES]);
void group_ge_pack(unsigned char r[GROUP_GE_PACKEDBYTES], const group_ge *x);

void group_ge_add(group_ge *r, const group_ge *x, const group_ge *y);
void group_ge_double(group_ge *r, const group_ge *x);

#endif
