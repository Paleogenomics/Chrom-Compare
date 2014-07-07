CC=gcc
#CFLAGS=-O2
CFLAGS=-g

all: tri-aln-report quad-aln-report find-ancestral

tri-aln-report: tri-aln-report.c
	echo "Making tri-aln-report..."
	$(CC) $(CFLAGS) -o tri-aln-report tri-aln-report.c

quad-aln-report: quad-aln-report.c
	echo "Making quad-aln-report..."
	$(CC) $(CFLAGS) -o quad-aln-report quad-aln-report.c

find-ancestral: find-ancestral.c
	echo "Making find-ancestral..."
	$(CC) $(CFLAGS) -o find-ancestral find-ancestral.c

.PHONY:
tests:
	./test_quad.sh

clean:
	echo "Removing targets..."
	rm tri-aln-report quad-aln-report find-ancestral
