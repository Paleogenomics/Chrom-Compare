#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>

#define MAX_LINE_LEN  (2000000)
#define MAX_ALN_BLOCK (1000000)
#define MAX_CHR_LEN   (250000000)
#define MAX_ID_LEN    (256)
#define MAX_DESC_LEN  (256)
#define MAX_FN_LEN    (2048)
#define FASTA_ROW_LEN (60)

typedef struct tri_aln {
  char* h_seq;  // human sequence
  char* c_seq;  // chimp sequence
  char* n_seq;  // neandertal sequence
  char* mask;   // which positions to consider, optional
  size_t h_len;
  char identifier[MAX_ID_LEN+1];
} Tri_Aln;
typedef struct tri_aln* TRIAP;

/* Prototypes */
FILE * fileOpen(const char *name, char access_mode[]);
void add_fa( char* fn, TRIAP aln, int pos );
char read_fasta_seq( FILE* fp, char* seq );
inline char revcom_base( const char base );
TRIAP init_TRIAP( void );
void output_summary( TRIAP aln );
void output_one_line_summary( const TRIAP aln, 
			      int start, int end );
void output_one_line_windows( const TRIAP aln, int window_size );
inline int valid_base( char b );
void output_bed_windows( const char* mask_fn, const char* chr_mask, const TRIAP aln );
void mask_evidence_code( char* mask, const char evid_code_fn[],
			 const char evid_code_accep[] );
void mask_cpg( TRIAP aln );
int all_diff( const int diverg_counts[] );
int count4kmer( char* kmer, unsigned int kmer_len,
		int* diverg_counts );

void help ( void ) {
  printf( "tri-aln-report -1 [first fa] -2 [second fa] -3 [third fa]\n" );
  printf( "               -o [make one-line simplified output of counts]\n" );
  printf( "               -W [make one-line simplified output over windows of this size\n" );
  printf( "               -b [optional bedfile of regions to consider; must use -c also]\n" );
  printf( "               -N [mask is a negative mask - ignore regions in mask]\n" );
  printf( "               -c [chromosome name to use with bedfile]\n" );
  printf( "               -C [mask CpG sites]\n" );
  printf( "               -e [evidence code filename]\n" );
  printf( "               -E [evidence code string to accept]\n" );
  printf( "               -I [identifier to use as extra, first column in simplified output]\n" );
  printf( "Reports the number of lineage specific differences given three\n" );
  printf( "aligned input fasta files. These can be big, chromosome-sized\n" );
  printf( "input files.\n" );
  exit( 0 );
}

int main( int argc, char* argv[] ) {
  extern char* optarg;
  extern int optind;
  char chr_mask[MAX_ID_LEN+1];
  char mfa_fn[MAX_FN_LEN];
  char mask_fn[MAX_FN_LEN];
  char first_fn[MAX_FN_LEN];
  char second_fn[MAX_FN_LEN];
  char third_fn[MAX_FN_LEN];
  char evid_code_fn[MAX_FN_LEN];
  char evid_code_accept[MAX_ID_LEN];
  int ich;
  int mask = 0;
  int ecm = 0;
  int cpg_mask = 0;
  int neg_mask = 0;
  int mask_seg_width = 0;
  int one_line_out   = 0;
  int window_size    = 0;
  TRIAP aln;

  /* No args - just help */
  if ( argc == 1 ) {
    help();
    exit( 0 );
  }

  /* Init aln; give mask memory and sets all positions to TRUE */
  aln = init_TRIAP();
  if ( aln == NULL ) {
    fprintf( stderr, "Failed to initialize TRIAP\n" );
    exit( 2 );
  }

  /* Process input arguments */
  while( (ich=getopt( argc, argv, "1:2:3:b:W:c:w:e:E:I:CNo" )) != -1 ) {
    switch(ich) {
    case '1' :
      strcpy( first_fn, optarg );
      break;

    case '2' :
      strcpy( second_fn, optarg );
      break;

    case '3' :
      strcpy( third_fn, optarg );
      break;

    case 'w' :
      mask_seg_width = atoi( optarg );
      break;

    case 'W' :
      window_size = atoi( optarg );
      one_line_out = 0;
      break;

    case 'I' :
      strcpy( aln->identifier, optarg );
      break;

    case 'b' :
      strcpy( mask_fn, optarg );
      mask = 1;
      break;
      
    case 'N' :
      neg_mask = 1;
      break;

    case 'c' :
      strcpy( chr_mask, optarg );
      break;

    case 'C' :
      cpg_mask = 1;
      break;

    case 'e' :
      strcpy( evid_code_fn, optarg );
      ecm = 1;
      break;
      
    case 'E' :
      strcpy( evid_code_accept, optarg );
      break;

    case 'o' :
      one_line_out = 1;
      break;

    default :
      help();
      exit( 0 );
    }
  }

  /* Load up the three input files */
  add_fa( first_fn, aln, 0 );
  add_fa( second_fn, aln, 1 );
  add_fa( third_fn, aln, 2 );

  /* Set aln->h_len to the minimum sequence length */
  aln->h_len = strlen( aln->h_seq );
  if ( strlen(aln->c_seq) < aln->h_len ) {
    aln->h_len = strlen( aln->c_seq );
  }
  if ( strlen(aln->n_seq) < aln->h_len ) {
    aln->h_len = strlen( aln->n_seq );
  }

  /* If CpG mask is requested, do it! */
  if ( cpg_mask ) {
    mask_cpg( aln );
  }

    /* Make output */
  if (mask) {
    /* User wants bedfile output over regions */
    output_bed_windows( mask_fn, chr_mask, aln );
    exit( 0 );
  }
  if ( one_line_out ) {
    output_one_line_summary( aln, 0, aln->h_len );
  }
  if ( window_size ) {
    output_one_line_windows( aln, window_size );
  }
  if ( !one_line_out && !window_size ) {
    output_summary( aln );
  }
  exit( 0 );
}

char read_fasta_seq( FILE* fp, char* seq ) {
  char c;
  size_t i = 0;

  c = fgetc(fp);
  while( (c != '>') &&
	 (c != EOF) &&
	 (i < MAX_CHR_LEN) ) {
    if ( isspace(c) ) {
      ;
    }
    else {
      c = toupper(c);
      seq[i++] = c;
    }
    c = fgetc(fp);
  }
  seq[i] = '\0';
  
  return c;
}
  
void mask_cpg( TRIAP aln ) {
  size_t pos;
  int last_ca_C = 0;
  int last_h_C  = 0;
  int last_c_C  = 0;
  /* Remember state, looking for places
     where the last base is C and we see
     a G at current place */
  for( pos = 1; pos < aln->h_len; pos++ ) {

    if ( (aln->h_seq[pos] == 'G') &&
	 (aln->h_seq[pos-1] == 'C') ) {

      /* Here's one! mask it! */
      aln->mask[pos]   = 0;
      aln->mask[pos-1] = 0;
    }
    
    if ( (aln->c_seq[pos]   == 'G') &&
	 (aln->c_seq[pos-1] == 'C') ) {
      /* Here's one! mask it! */
      aln->mask[pos]   = 0;
      aln->mask[pos-1] = 0;
    }

    if ( (aln->n_seq[pos]   == 'G') &&
	 (aln->n_seq[pos-1] == 'C') ) {
      /* Here's one! mask it! */
      aln->mask[pos]   = 0;
      aln->mask[pos-1] = 0;
    }
  }
  return;
}

/* Takes a char* - filename of fasta file to add in
   postions (pos) to the TRIAP aln.
   Dies with error code 1 if the input is bad
*/
void add_fa( char* fn, TRIAP aln, int pos ) {
  FILE* fa;
  char c;
  
  fa = fileOpen( fn, "r" );
  if ( fa == NULL ) {
    fprintf( stderr, "Problem reading input: %s\n", fn );
    exit( 1 );
  }
  c = fgetc(fa);

  if ( c == EOF ) return;
  if ( c != '>' ) return;

  /* Go to end of header line of fasta input sequence */
  while( c != '\n' ) {
    c = fgetc(fa);
  }

  /* Now, fill up the appropriate sequence */
  if ( pos == 0 ) {
    c = read_fasta_seq( fa, aln->c_seq );
    return;
  }
  if ( pos == 1 ) {
    c = read_fasta_seq( fa, aln->h_seq );
    return;
  }
  if ( pos == 2 ) {
    c = read_fasta_seq( fa, aln->n_seq );
    return;
  }
}

/* Mask out regions not mentioned in input file 
   INPUT FILE looks like this:
   chr1 1000 2000
   chr1 2501 34195
   etc.
*/  
void output_bed_windows( const char* mask_fn, const char* chr_mask, const TRIAP aln ) {
  FILE* f;
  char chr[MAX_ID_LEN];
  char line[MAX_LINE_LEN+1];
  int start, end;
  size_t i;

  f = fileOpen( mask_fn, "r" );
  while( fgets( line, MAX_LINE_LEN, f ) != NULL ) {
    if ( sscanf( line, "%s\t%u\t%u", chr, &start, &end ) 
	 == 3 ) {
      /* The right chromosome? */
      if ( strcmp( chr, chr_mask ) == 0 ) {
	/* Range check */
	if ( start < 0 ) {
	  start = 0;
	}
	if ( end > aln->h_len ) {
	  end = aln->h_len;
	}
	output_one_line_summary( aln, start, end );
      }
    }
  }
  fclose(f);
}

/* Takes the TRIAP aln and the user-specified window_size
   Prints the counts of same and divergent positions over
   each window */
void output_one_line_windows( const TRIAP aln, int window_size ) {
  int start = 0;
  int end;

  /* Initialize */
  end = start + window_size;
  
  while( end < aln->h_len ) {
    output_one_line_summary( aln, start, end );
    start += window_size;
    end   = start + window_size;
  }
}

void output_one_line_summary( const TRIAP aln, 
			      int start, int end ) {
  char diverg_kmer[4];
  size_t pos, h, inx;
  int diverg_counts[64];
  char bases[] = { 'A', 'G', 'C', 'T' };

  diverg_kmer[3] = '\0';

  /* Zero the array of diverg_counts */
  for( pos = 0; pos < 64; pos++ ) {
    diverg_counts[pos] = 0;
  }

  /* Go through each position and add a count to the
     appropriate position in the diverg_counts if all
     three have a non-masked valid base */
  /* Check back end */
  if ( end > aln->h_len ) {
    end = aln->h_len;
  }
  for( pos = start; pos < end; pos++ ) {
    if ( aln->mask[pos] && 
	 valid_base( aln->c_seq[pos] ) &&
	 valid_base( aln->h_seq[pos] ) &&
	 valid_base( aln->n_seq[pos] ) ) {
      diverg_kmer[0] = aln->c_seq[pos];
      diverg_kmer[1] = aln->h_seq[pos];
      diverg_kmer[2] = aln->n_seq[pos];
      if ( kmer2inx( diverg_kmer, 3, &inx ) ) {
	diverg_counts[inx]++;
      }
    }
  }

  /* If there's an identifier specified, it means we want 
     it printed here */
  if ( aln->identifier[0] != '\0' ) {
    printf( "%s ", aln->identifier );
  }

  /* Output the region */
  printf( "%d %d ", start, end );

  /* Ouput all the same counts */
  printf( "%d ", count4kmer( "AAA", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "CCC", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "GGG", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "TTT", 3, diverg_counts ) );

  /* Output the "first is different counts" in a sensible
     order */
  printf( "%d ", count4kmer( "CAA", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "GAA", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "TAA", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "ACC", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "GCC", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "TCC", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "AGG", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "CGG", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "TGG", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "ATT", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "CTT", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "GTT", 3, diverg_counts ) );

  /* Output the "second is different counts" in a sensible
     order */
  printf( "%d ", count4kmer( "ACA", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "AGA", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "ATA", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "CAC", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "CGC", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "CTC", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "GAG", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "GCG", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "GTG", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "TAT", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "TCT", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "TGT", 3, diverg_counts ) );

  /* Output the "first is different counts" in a sensible
     order */
  printf( "%d ", count4kmer( "AAC", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "AAG", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "AAT", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "CCA", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "CCG", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "CCT", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "GGA", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "GGC", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "GGT", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "TTA", 3, diverg_counts ) );
  printf( "%d ", count4kmer( "TTC", 3, diverg_counts ) );
  printf( "%d\n", count4kmer( "TTG", 3, diverg_counts ) );

}

/* Takes a kmer, its length, and a populated diverg_counts[]
   Returns the count for the divergence represented by this
   kmer */
int count4kmer( char* kmer, unsigned int kmer_len,
		int* diverg_counts ) {
  size_t inx;
  kmer2inx( kmer, kmer_len, &inx );
  return diverg_counts[inx];
}

void output_summary( TRIAP aln ) {
  char diverg_kmer[4];
  size_t pos, inx, h, h_inx, c_inx, n_inx;
  int diverg_counts[64];
  char bases[] = { 'A', 'G', 'C', 'T' };

  diverg_kmer[3] = '\0';

  for( pos = 0; pos < 64; pos++ ) {
    diverg_counts[pos] = 0;
  }

  for( pos = 0; pos < aln->h_len; pos++ ) {
    if ( aln->mask[pos] && 
	 valid_base( aln->c_seq[pos] ) &&
	 valid_base( aln->h_seq[pos] ) &&
	 valid_base( aln->n_seq[pos] ) ) {
      diverg_kmer[0] = aln->c_seq[pos];
      diverg_kmer[1] = aln->h_seq[pos];
      diverg_kmer[2] = aln->n_seq[pos];
      if ( kmer2inx( diverg_kmer, 3, &inx ) ) {
	diverg_counts[inx]++;
      }
    }
  }
  
  printf( "# All the same, 1st=2nd=3rd\n" );
  for( h = 0; h <= 3; h++ ) {
    diverg_kmer[0] = bases[h];
    diverg_kmer[1] = bases[h];
    diverg_kmer[2] = bases[h];
    kmer2inx( diverg_kmer, 3, &inx );
    printf( "%c %d\n", bases[h], diverg_counts[inx] );
  }
  printf( "# All different, first != second != third: %d\n", 
	  all_diff(diverg_counts) );
  
  printf( "# first second third different\n" );
  kmer2inx( "CAA", 3, &c_inx );
  kmer2inx( "ACA", 3, &h_inx );
  kmer2inx( "AAC", 3, &n_inx );
  printf( "A->C %d %d %d\n", 
	  diverg_counts[c_inx],
	  diverg_counts[h_inx],
	  diverg_counts[n_inx] );

  kmer2inx( "GAA", 3, &c_inx );
  kmer2inx( "AGA", 3, &h_inx );
  kmer2inx( "AAG", 3, &n_inx );
  printf( "A->G %d %d %d\n", 
	  diverg_counts[c_inx],
	  diverg_counts[h_inx],
	  diverg_counts[n_inx] );

  kmer2inx( "TAA", 3, &c_inx );
  kmer2inx( "ATA", 3, &h_inx );
  kmer2inx( "AAT", 3, &n_inx );
  printf( "A->T %d %d %d\n", 
	  diverg_counts[c_inx],
	  diverg_counts[h_inx],
	  diverg_counts[n_inx] );

  kmer2inx( "ACC", 3, &c_inx );
  kmer2inx( "CAC", 3, &h_inx );
  kmer2inx( "CCA", 3, &n_inx );
  printf( "C->A %d %d %d\n", 
	  diverg_counts[c_inx],
	  diverg_counts[h_inx],
	  diverg_counts[n_inx] );

  kmer2inx( "GCC", 3, &c_inx );
  kmer2inx( "CGC", 3, &h_inx );
  kmer2inx( "CCG", 3, &n_inx );
  printf( "C->G %d %d %d\n", 
	  diverg_counts[c_inx],
	  diverg_counts[h_inx],
	  diverg_counts[n_inx] );

  kmer2inx( "TCC", 3, &c_inx );
  kmer2inx( "CTC", 3, &h_inx );
  kmer2inx( "CCT", 3, &n_inx );
  printf( "C->T %d %d %d\n", 
	  diverg_counts[c_inx],
	  diverg_counts[h_inx],
	  diverg_counts[n_inx] );

  kmer2inx( "AGG", 3, &c_inx );
  kmer2inx( "GAG", 3, &h_inx );
  kmer2inx( "GGA", 3, &n_inx );
  printf( "G->A %d %d %d\n", 
	  diverg_counts[c_inx],
	  diverg_counts[h_inx],
	  diverg_counts[n_inx] );

  kmer2inx( "CGG", 3, &c_inx );
  kmer2inx( "GCG", 3, &h_inx );
  kmer2inx( "GGC", 3, &n_inx );
  printf( "G->C %d %d %d\n", 
	  diverg_counts[c_inx],
	  diverg_counts[h_inx],
	  diverg_counts[n_inx] );

  kmer2inx( "TGG", 3, &c_inx );
  kmer2inx( "GTG", 3, &h_inx );
  kmer2inx( "GGT", 3, &n_inx );
  printf( "G->T %d %d %d\n", 
	  diverg_counts[c_inx],
	  diverg_counts[h_inx],
	  diverg_counts[n_inx] );

  kmer2inx( "ATT", 3, &c_inx );
  kmer2inx( "TAT", 3, &h_inx );
  kmer2inx( "TTA", 3, &n_inx );
  printf( "T->A %d %d %d\n", 
	  diverg_counts[c_inx],
	  diverg_counts[h_inx],
	  diverg_counts[n_inx] );

  kmer2inx( "CTT", 3, &c_inx );
  kmer2inx( "TCT", 3, &h_inx );
  kmer2inx( "TTC", 3, &n_inx );
  printf( "T->C %d %d %d\n", 
	  diverg_counts[c_inx],
	  diverg_counts[h_inx],
	  diverg_counts[n_inx] );

  kmer2inx( "GTT", 3, &c_inx );
  kmer2inx( "TGT", 3, &h_inx );
  kmer2inx( "TTG", 3, &n_inx );
  printf( "T->G %d %d %d\n", 
	  diverg_counts[c_inx],
	  diverg_counts[h_inx],
	  diverg_counts[n_inx] );



}

inline char revcom_base( const char base ) {
  switch (base) {
  case 'A' :
    return 'T';
  case 'C' :
    return 'G';
  case 'G' :
    return 'C';
  case 'T' :
    return 'A';
  default :
    return 'N';
  }
}

TRIAP init_TRIAP( void ) {
  TRIAP aln;
  aln = (TRIAP)malloc(sizeof(Tri_Aln));
  aln->h_len  = 0;
  aln->h_seq  = (char*)malloc((MAX_CHR_LEN+1)*sizeof(char));
  aln->c_seq  = (char*)malloc((MAX_CHR_LEN+1)*sizeof(char));
  aln->n_seq  = (char*)malloc((MAX_CHR_LEN+1)*sizeof(char));
  aln->mask   = (char*)malloc((MAX_CHR_LEN+1)*sizeof(char));

  /* initialize to all Ns */
  memset( aln->h_seq, 'N', MAX_CHR_LEN );
  memset( aln->c_seq, 'N', MAX_CHR_LEN );
  memset( aln->n_seq, 'N', MAX_CHR_LEN );

  /* initialize mask to all trues */
  memset( aln->mask, 1, MAX_CHR_LEN );

  /* Set the identifier to NULL. Will be set to appropriate
     value if an identifier is given */
  aln->identifier[0] = '\0';

  return aln;
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
/* kmer2inx
   Args: (1) a pointer to a character string;
             the kmer to find the corresponding index of;
	     might not be null-terminated
	 (2) length of the kmer
	 (3) pointer to size_t to put the index
   Returns: TRUE if the index was set, FALSE if it could not
            be set because of some non A,C,G,T character
   Uses the formula A=>00, C=>01, G=>11, T=>11 to make a
   bit string for the kmer. Any other character is not allowed
   and will cause an error
   The bit string is constructed by reading the kmer from left
   to right. This bit-string is then interpreted as a variable
   of type size_t and is appropriate as an array index
*/
int kmer2inx( const char* kmer,
	      const unsigned int kmer_len,
	      size_t* inx ) {
  size_t l_inx = 0;
  int i = 0;
  char curr_char;

  while( i < kmer_len ) {
    l_inx = l_inx << 2;
    curr_char = toupper(kmer[i]); // Upper case it in case it is not
    switch( curr_char ) {
    case 'A' :
      l_inx += 0;
      break;
    case 'C' :
      l_inx += 1;
      break;
    case 'G' :
      l_inx += 2;
      break;
    case 'T' :
      l_inx += 3;
      break;
    default :
      return 0; // not valid!
    }
    i++;
  }
  *inx = l_inx;
  return 1; // valid!
}

void mask_evidence_code( char* mask, const char evid_code_fn[],
			 const char evid_code_accept[] ) {
  FILE* f;
  char c;
  int accept;
  size_t evl;
  size_t i = 0;
  size_t j = 0;

  f   = fileOpen( evid_code_fn, "r" );
  evl = strlen( evid_code_accept );

  if ( evl == 0 ) {
    fprintf( stderr, "No evidence code accept string given!\n" );
    return;
  }

  c = fgetc( f );
  
  if ( c == EOF ) return;
  if ( c != '>' ) return;
  
  /* Go to end of header line */
  while( c != '\n' ) {
    c = fgetc(f);
  }

  /* Now, read the evidence code and filter out everything 
     that doesn't match one of the allowed values */
  c = fgetc( f );
  while( (c != '>') &&
	 (c != EOF) &&
	 (i < MAX_CHR_LEN) ) {
    if ( isspace(c) ) {
      ;
    }
    else {
      accept = 0;
      for( j = 0; j < evl; j++ ) {
	if ( c == evid_code_accept[j] ) {
	  accept = 1;
	}
      }
      if ( !accept ) {
	mask[i] = 0;
      }
      i++;
    }
    c = fgetc( f );
  }
  return;
}

int all_diff( const int diverg_counts[] ) {
  size_t inx;
  int tot = 0;
  size_t b1, b2, b3;
  char diverg_kmer[4];
  char bases[] = { 'A', 'C', 'G', 'T' };

  for( b1 = 0; b1 <= 3; b1++ ) {
    for ( b2 = 0; b2 <= 3; b2++ ) {
      for ( b3 = 0; b3 <= 3; b3++ ) {
	if (( b1 != b2) &&
	    ( b1 != b3) &&
	    ( b2 != b3) ) {
	  diverg_kmer[0] = bases[b1];
	  diverg_kmer[1] = bases[b2];
	  diverg_kmer[2] = bases[b3];
	  
	  kmer2inx( diverg_kmer, 3, &inx );
	  tot += diverg_counts[inx];
	}
      }
    }
  }
    return tot;
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
