all: ecdh25519test.bin ecdh25519speed.bin

include ../common/make.mk

ecdh25519test.elf: fe25519.o fe25519.h group.o group.h smult.o smult.h test.o $(OBJS) $(LDSCRIPT)
		$(LD) -o $@ test.o fe25519.o group.o smult.o $(OBJS) $(LDFLAGS)

ecdh25519speed.elf: fe25519.o fe25519.h group.o group.h smult.o smult.h speed.o $(OBJS) $(LDSCRIPT)
		$(LD) -o $@ speed.o fe25519.o group.o smult.o $(OBJS) $(LDFLAGS)
