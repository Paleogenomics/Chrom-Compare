#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include "pileup.h"

#define MAX_LINE_LEN (10240)
#define BASE_QUAL_CUT (25)
#define MAP_QUAL_CUT (30)
#define COV_CUT (10)
#define LOW_COV_CUT (1)
#define MAX_CHR_LEN (250000000)
#ifndef MAX_ID_LEN
#define MAX_ID_LEN (128)
#endif
#define MAX_FN_LEN (512)
#define MAX_AXT_LINE_LEN (11000000)
#define FASTA_ROW_LEN (60)

typedef struct chr {
  char id[MAX_ID_LEN+1];
  char* seq;
} Chr;
typedef struct chr* ChrP;

void help( void ) {
  printf( "Usage: pu2fa [options] -c <chr name> {pileup file | stdin}\n\n" );
  printf( "Required:\n" );
  printf( "      -c chromosome name to output as fasta.\n" );
  printf( "Options:\n" );
  printf( "      -q <q-score cutoff filename>\n" );
  printf( "      -C <max coverage; inclusive>\n" );
  printf( "      -l <minimum coverage; inclusive>\n" );
  printf( "      -m <map-quality cutoff; default = %d>\n", MAP_QUAL_CUT );
  printf( "      -s <region start> -e <region end>\n" );
  printf( "      -b use most frequently-encountered base rather than random base\n" );
  printf( "      -w use best base, determined by additive quality scores\n\n" );
  printf( "Takes pileup output as produced from samtools mpileup -s and writes a\n" );
  printf( "fasta sequence of the aligned reads from the pileup output for one\n" );
  printf( "chromosome (-c).\n" );
  printf( "Picks randomly the first base to pass the quality score thresholding.\n" );
  printf( "By default, writes the entire fasta sequence, even if only a small\n" );
  printf( "part of it is given as in mpileup input. A small region can\n" );
  printf( "be specified by using the -s and -e options (1-based coordinates)\n" );
  return;
}

void write_fasta( const char* id, const char* desc,
                  const char* seq, const size_t len,
      int reg_start, int reg_end );

int main( int argc, char* argv[] ) {
  extern char* optarg;
  int ich, base_inx, rand_inx;
  size_t chr_len = 0;
  int reg_start = -1;
  int reg_end = -1;

  char qcut_fn[MAX_FN_LEN+1]; // filename of qcut file (-q option)
  // initialize the string so that we can check later whether the -q option was
  // supplied
  qcut_fn[0] = '\0';

  char target_chr[MAX_ID_LEN+1];
  char line[MAX_LINE_LEN+1];
  unsigned int covc = COV_CUT;
  unsigned int lowcovc = LOW_COV_CUT;
  unsigned int mqc  = MAP_QUAL_CUT;
  QcutsP qcp; // struct storing the quality cutoffs
  PulP pp;
  ChrP cp;
  bool best_base = false; // holds option whether to use best base
  bool weighted_best_base = false; // holds option whether to use weighted best

  /* Seed rand() for choosing random bases */
  srand(time(NULL));

  /* Process command line arguments */
  if ( argc == 1 ) {
    help();
    exit( 0 );
  }

  while( (ich=getopt( argc, argv, "q:C:c:l:m:s:e:b:w" )) != -1 ) {
    switch(ich) {
    case 'q' :
      strcpy( qcut_fn, optarg );
      break;
    case 'C' :
      covc = atoi(optarg);
      break;
    case 'l' :
      lowcovc = atoi(optarg);
      break;
    case 'c' :
      strcpy( target_chr, optarg );
      break;
    case 'm' :
      mqc = atoi(optarg);
      break;
    case 'b' :
      best_base = true;
      break;
    case 'w' :
      weighted_best_base = true;
      best_base = true;
      break;
    case 's' :
      reg_start = atoi(optarg);
      break;
    case 'e' :
      reg_end = atoi(optarg);
      break;
    default :
      help();
      exit( 0 );
    }
  }

  /* Init pile-up line structure */
  pp = (PulP)malloc(sizeof(Pul));

  /* Init the chromosome */
  cp = (ChrP)malloc(sizeof(Chr));
  cp->seq = (char*)malloc((MAX_CHR_LEN+1) * sizeof(char));
  memset( cp->seq, 'N', MAX_CHR_LEN );

  /* If qcut file provided, parse it to get the base-specific quality score
   * cutoffs for this library; otherwise, use defaults. */
  qcp = parse_q_score_cut( qcut_fn );

  /* Go through the pileup output, getting bases for
     the comparison guy and human and chimp */
  while( fgets( line, MAX_LINE_LEN, stdin ) != NULL ) {
    if ( (line2pul( line, pp ) == 0) &&
         (strcmp( pp->chr, target_chr ) == 0) ) {
      /* Check the coverage cutoffs first */
      if ( (pp->cov >= lowcovc) &&
           (pp->cov <= covc) ) {
        /* Do I want the BEST base */
        if ( best_base ) {
          base_inx = best_base_from_pul(pp, qcp, mqc, covc, weighted_best_base);
        }
        /* Or a random base? */
        else {
          base_inx = rand_good_base_from_pul( pp, qcp, mqc, covc );
        }
        if ( base_inx >= 0 ) {
          cp->seq[pp->pos - 1] = inx2base( base_inx );
        }
        if ( pp->pos > chr_len ) {
          chr_len = pp->pos;
        }
      }
    }
  }

  /* Make the output */
  write_fasta( target_chr, "", cp->seq, chr_len, reg_start, reg_end );
  exit( 0 );
}


void write_fasta( const char* id, const char* desc,
                  const char* seq, const size_t len,
      int reg_start, int reg_end ) {
  size_t pos;
  size_t line_pos = 0;
  if ( reg_start == -1 ) { // was never set by user, set it to 0
    reg_start = 0;
  }
  else { // it was set by user, but in 1-based coordinates; convert to 0-based
    reg_start -= 1;
  }

  if ( reg_end == -1 ) {// was never set by user, set it to the end
    reg_end = len;
  }
  else { // it was set, but in 1-based coordiantes; convert to 0-based
    reg_end -= 1;
  }

  printf( ">%s %s %d %d", id, desc, reg_start, reg_end );
  for( pos = reg_start; pos < reg_end; pos++ ) {
    if ( line_pos % FASTA_ROW_LEN == 0 ) {
      printf( "\n" );
    }
    printf( "%c", seq[pos] );
    line_pos++;
  }
  printf( "\n" );
  return;
}
