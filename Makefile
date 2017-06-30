CC=gcc
#CFLAGS=-O2
CFLAGS=-g -std=gnu89

all: tri-aln-report quad-aln-report find-ancestral pu2fa

tri-aln-report: tri-aln-report.c
	$(CC) $(CFLAGS) -o $@ $<

quad-aln-report: quad-aln-report.c
	$(CC) $(CFLAGS) -o $@ $<

find-ancestral: find-ancestral.c
	$(CC) $(CFLAGS) -o $@ $<

pu2fa: pu2fa.c pileup.c pileup.h
	$(CC) $(CFLAGS) -o $@ $^

.PHONY: tests clean all

tests:
	./test_quad.sh
	./test_pu2fa.sh
clean:
	rm -f tri-aln-report quad-aln-report find-ancestral pu2fa
