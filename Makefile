CC=gcc
#CFLAGS=-O2
CFLAGS=-g

all: tri-aln-report quad-aln-report

tri-aln-report: tri-aln-report.c
	echo "Making tri-aln-report..."
	$(CC) $(CFLAGS) -o tri-aln-report tri-aln-report.c

quad-aln-report: quad-aln-report.c
	echo "Making quad-aln-report..."
	$(CC) $(CFLAGS) -o quad-aln-report quad-aln-report.c

.PHONY:
tests:
	./test_quad.sh

clean:
	echo "Removing targets..."
	rm tri-aln-report quad-aln-report
