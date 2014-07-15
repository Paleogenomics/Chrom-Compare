#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#define NUM_ALN_CHRS (32)
#define MAX_FN_LEN (512)
#define MAX_ID_LEN (128)
#define MAX_CHR_LEN (250000000)
#define DEFAULT_MIN_DIFF (1)

typedef struct chr {
  char id[MAX_ID_LEN+1];
  char* seq;
  size_t len;
} Chr;
typedef struct chr* ChrP;

typedef struct aln_chrs {
  char id[MAX_ID_LEN];
  ChrP anc_cp; // the ancestral chromosome
  ChrP foc_cp; // the focal chromosome
  ChrP cps[NUM_ALN_CHRS]; // aligned chromosomes
  char* mask;   // which positions to consider, optional
  size_t num_cps;
} Aln_Chrs;
typedef struct aln_chrs* Aln_ChrsP;

FILE * fileOpen(const char *name, char access_mode[]);
void mask_cpg( Aln_ChrsP acsp );
void fasta2chr( const char* fn, ChrP cp );
inline int valid_base( char b );
int cpg_site( const Aln_ChrsP acsp, const size_t pos );
void output_anc_sites( const Aln_ChrsP acsp, 
		       const int min_diff,
		       const int cpg_mask,
		       const char* id_to_assign );
void help( void ) {
  printf( "find-ancestral -a <ancestral fasta> -f <focal fasta> -C [mask CpG sites]\n" );
  printf( "            -i [chromosome ID to put in first column] <aligned chrs...>\n" );
  printf( "            -m [minimum number of sequences that must be different; default = %d]\n", DEFAULT_MIN_DIFF);
  printf( "            [aligned chrs...]\n" );
  printf( "Scans the input chromosomes to find and report positions where the focal\n" );
  printf( "input sequence matches the ancestral input sequence and all other input\n" );
  printf( "sequences match one another but are different. This is meant to find\n" );
  printf( "positions where the focal sequence is an outgroup to all other input\n" );
  printf( "sequences.\n" );
}

int main( int argc, char* argv[] ) {
  extern char* optarg;
  extern int optind;
  char anc_chr_fn[MAX_FN_LEN+1];
  char foc_chr_fn[MAX_FN_LEN+1];
  char id_to_assign[MAX_ID_LEN];
  Aln_ChrsP acsp;
  int ich, input_num;
  int cpg_mask = 0;
  size_t chr_inx = 0;
  int min_diff = DEFAULT_MIN_DIFF;

  if ( argc == 1 ) {
    help();
    exit( 0 );
  }

  while( (ich=getopt( argc, argv, "a:f:i:m:C" )) != -1 ) {
    switch(ich) {
    case 'a' :
      strcpy( anc_chr_fn, optarg );
      break;
    case 'f' :
      strcpy( foc_chr_fn, optarg );
      break;
    case 'i' :
      strcpy( id_to_assign, optarg );
      break;
    case 'm' :
      min_diff = atoi( optarg );
      break;
    case 'C' :
      cpg_mask = 1;
      break;
    default :
      help();
    }
  }

  /* Initialize the Aln_Chrs data structure */
  acsp = (Aln_ChrsP)malloc(sizeof(Aln_Chrs));
  acsp->num_cps = 0;
  strcpy( acsp->id, id_to_assign );

  /* Slurp up the ancestral chromosome */
  acsp->anc_cp = (ChrP)malloc(sizeof(Chr));
  fasta2chr( anc_chr_fn, acsp->anc_cp );

  /* Slurp up the focal chromosome */
  acsp->foc_cp = (ChrP)malloc(sizeof(Chr));
  fasta2chr( foc_chr_fn, acsp->foc_cp );

  /* Init the mask */
  acsp->mask = (char*)malloc((MAX_CHR_LEN+1)*sizeof(char));
  memset( acsp->mask, 1, MAX_CHR_LEN );

  /* Parse each aligned chromsome */
  for ( input_num = optind; input_num < argc; input_num++ ) {
    acsp->cps[chr_inx] = (ChrP)malloc(sizeof(Chr));
    fasta2chr( argv[input_num], acsp->cps[chr_inx] );
    chr_inx++;
  }
  acsp->num_cps = chr_inx;

  /* If CpG mask is requested, do it! */
  if ( cpg_mask ) {
    mask_cpg( acsp );
  }

  /* Output the focal = ancestral != everything else sites */
  output_anc_sites( acsp, min_diff, cpg_mask, id_to_assign );
  exit( 0 );
}

void mask_cpg( Aln_ChrsP acsp ) {
  size_t pos;
  size_t cn;
  for ( pos = 0; pos < acsp->anc_cp->len; pos++ ) {
    if ( (acsp->anc_cp->seq[pos] == 'G') &&
	 (acsp->anc_cp->seq[pos-1] == 'C') ) {
      acsp->mask[pos]   = 0;
      acsp->mask[pos-1] = 0;
    }
    if ( (acsp->foc_cp->seq[pos] == 'G') &&
	 (acsp->foc_cp->seq[pos-1] == 'C') ) {
      acsp->mask[pos]   = 0;
      acsp->mask[pos-1] = 0;
    }
    for( cn = 0; cn < acsp->num_cps; cn++ ) {
      if ( (acsp->cps[cn]->seq[pos] == 'G') &&
	   (acsp->cps[cn]->seq[pos-1] == 'C') ) {
	acsp->mask[pos]   = 0;
	acsp->mask[pos-1] = 0;
      }
    }
  }
}

/** fileOpen **/
FILE * fileOpen(const char *name, char access_mode[]) {
  FILE * f;
  f = fopen(name, access_mode);
  if (f == NULL) {
    fprintf( stderr, "%s\n", name);
    perror("Cannot open file");
    return NULL;
  }
  return f;
}

int cpg_site( const Aln_ChrsP acsp, const size_t pos ) {
  size_t cn;
  for( cn = 0; cn < acsp->num_cps; cn++ ) {
    if ( (acsp->cps[cn]->seq[pos]   == 'C') &&
	 (acsp->cps[cn]->seq[pos+1] == 'G') ) {
      return 1;
    }
    if ( (acsp->cps[cn]->seq[pos-1]   == 'C') &&
	 (acsp->cps[cn]->seq[pos]     == 'G') ) {
      return 1;
    }
  }
  return 0;
}

void fasta2chr( const char* fn, ChrP cp ) {
  FILE * fasta;
  char c;
  size_t i;
  
  fasta = fileOpen( fn, "r" );
  cp->seq = (char*)malloc((MAX_CHR_LEN+1) * sizeof(char));
  c = fgetc( fasta );
  if ( c == EOF ) return;
  if ( c != '>' ) return;

  // get id
  i = 0;
  while( (!isspace( c=fgetc( fasta ) ) &&
          (i < MAX_ID_LEN) ) ) {
    if ( c == EOF ) {
      return;
    }
    cp->id[i] = c;
    i++;
    if ( i == MAX_ID_LEN ) {
      //id is too long - truncate it
      cp->id[i] = '\0';
    }
  }
  cp->id[i] = '\0';

  // everything else on this line is description, if there is anything
  if ( c == '\n' ) {
    ;
  }
  else { // not end of line, so skip past the stupid whitespace...
    while( (c != '\n') &&
           (isspace(c)) ) {
      c = fgetc( fasta );
    }
    ///...everthing else is description
    i = 0;
    ungetc( c, fasta );
    while ( c != '\n' ) {
      c = fgetc( fasta );
    }
  }

  // read sequence
  i = 0;
  c = fgetc( fasta );
  while ( ( c != '>' ) &&
          ( c != EOF ) &&
          (i < MAX_CHR_LEN) ) {
    if ( isspace(c) ) {
      ;
    }
    else {
      c = toupper( c );
      cp->seq[i++] = c;
    }
    c = fgetc( fasta );
  }
  cp->seq[i] = '\0';

  cp->len = i;

  return;
}


void output_anc_sites( const Aln_ChrsP acsp, 
		       const int min_diff,
		       const int cpg_mask,
		       const char* id_to_assign ) {
  int not_anc, diff;
  size_t pos, i;
  char r, anc_b;
  for( pos = 0; pos < acsp->anc_cp->len; pos++ ) {
    not_anc = 0;
    diff    = 0;
    /* Are ancestral and focal bases the same? */
    if (valid_base(acsp->anc_cp->seq[pos]) &&
	(acsp->anc_cp->seq[pos] == acsp->foc_cp->seq[pos]) &&
	(acsp->mask[pos]) ) {
      anc_b = acsp->anc_cp->seq[pos];
      /* Are none of the other bases ancestral? */
      for( i = 0; i < acsp->num_cps; i++ ) {
	if ( valid_base(acsp->cps[i]->seq[pos]) ) {
	  if ( acsp->cps[i]->seq[pos] != anc_b ) {
	    diff++;
	  }
	  else {
	    not_anc = 1; // this base is ancestral - we don't care about this site
	  }
	}
      }
      /* We've looked at all other bases. Is this site ancestral? */
      if ( not_anc ) {
	; // skip it
      }
      else {
	if ( diff >= min_diff ) {
	  printf( "%s %d\n", id_to_assign, (int)pos );
	}
      }
    }
  }
}

inline int valid_base( char b ) {
  char B;
  B = toupper(b);
  switch(B) {
  case 'A' :
    return 1;
  case 'C' :
    return 1;
  case 'G' :
    return 1;
  case 'T' :
    return 1;
  default :
    return 0;
  }
}

