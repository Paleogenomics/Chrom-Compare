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

typedef struct quad_aln {
  char* h1_seq; // human 1 sequence
  char* h2_seq; // human 2 sequence
  char* n_seq;  // neandertal sequence
  char* c_seq;  // chimp sequence
  char* mask;   // which positions to consider, optional
  size_t h_len;
  char identifier[MAX_ID_LEN+1];
} Quad_Aln;
typedef struct quad_aln* QUADAP;

/* Prototypes */
FILE * fileOpen(const char *name, char access_mode[]);
void add_fa( char* fn, QUADAP aln, int pos );
char read_fasta_seq( FILE* fp, char* seq );
inline char revcom_base( const char base );
QUADAP init_QUADAP( void );
void output_summary( const QUADAP aln, const int start, const int end );
void output_one_line_windows( const QUADAP aln, int window_size );
			      
void count_diverg( const QUADAP aln, const int start, const int end, int diverg_counts[] );
inline int valid_base( char b );
void mask_from_fn( char* mask_fn, char* mask, const char* chr_mask,
		   const int mask_seg_width, const int neg_mask );
void mask_evidence_code( char* mask, const char evid_code_fn[],
			 const char evid_code_accep[] );
void mask_cpg( QUADAP aln );
int all_diff( const int diverg_counts[] );
int count4kmer( char* kmer, unsigned int kmer_len,
		int* diverg_counts );

void help ( void ) {
  printf( "quad-aln-report -1 [first fa] -2 [second fa] -3 [third fa] -4 [fourth fa; ancestral]\n" );
  printf( "               -o [make one-line simplified output of counts]\n" );
  printf( "               -W [make one-line simplified output over windows of this size\n" );
  printf( "               -M [optional mask fn]\n" );
  printf( "               -N [mask is a negative mask - ignore regions in mask]\n" );
  printf( "               -c [chromosome name to use with mask]\n" );
  printf( "               -C [mask CpG sites]\n" );
  printf( "               -w [widen mask segments to this width\n" );
  printf( "                   if they are shorter]\n" );
  printf( "               -e [evidence code filename]\n" );
  printf( "               -E [evidence code string to accept]\n" );
  printf( "               -I [identifier to use as extra, first column in simplified output]\n" );
  printf( "Reports the number of positions with various configurations of\n" );
  printf( "matching alleles. The assumption is that the fourth (last) input\n" );
  printf( "sequence is the ancestral sequence. The default output is this:\n" );
  printf( "START END AAAA BAAA ABAA AABA BBBA BBAA BABA ABBA\n" );
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
  char fourth_fn[MAX_FN_LEN];
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
  QUADAP aln;

  /* No args - just help */
  if ( argc == 1 ) {
    help();
    exit( 0 );
  }

  /* Init aln; give mask memory and sets all positions to TRUE */
  aln = init_QUADAP();
  if ( aln == NULL ) {
    fprintf( stderr, "Failed to initialize QUADAP\n" );
    exit( 2 );
  }

  /* Process input arguments */
  while( (ich=getopt( argc, argv, "1:2:3:4:M:W:c:w:e:E:I:CNo" )) != -1 ) {
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

    case '4' :
      strcpy( fourth_fn, optarg );
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

    case 'M' :
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

  /* If mask specified, do it! */
  if ( mask ) {
    mask_from_fn( mask_fn, aln->mask, chr_mask, 
		  mask_seg_width, neg_mask );
  }

  /* Load up the three input files */
  add_fa( first_fn,  aln, 0 );
  add_fa( second_fn, aln, 1 );
  add_fa( third_fn,  aln, 2 );
  add_fa( fourth_fn, aln, 3 );

  /* Set aln->h_len to the minimum sequence length */
  aln->h_len = strlen( aln->h1_seq );
  if ( strlen(aln->c_seq) < aln->h_len ) {
    aln->h_len = strlen( aln->c_seq );
  }
  if ( strlen(aln->n_seq) < aln->h_len ) {
    aln->h_len = strlen( aln->n_seq );
  }
  if ( strlen(aln->h2_seq) < aln->h_len ) {
    aln->h_len = strlen( aln->h2_seq );
  }

  /* If CpG mask is requested, do it! */
  if ( cpg_mask ) {
    mask_cpg( aln );
  }

  /* If evidence code mask is specified, do it */
  if ( ecm ) {
    mask_evidence_code(aln->mask, evid_code_fn, 
		       evid_code_accept );
  }

  /* Make output */
  if ( window_size ) {
    output_one_line_windows( aln, window_size ); 
  }
  else {
    output_summary( aln, 0, aln->h_len );
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
  
void mask_cpg( QUADAP aln ) {
  size_t pos;
  int last_ca_C = 0;
  int last_h_C  = 0;
  int last_c_C  = 0;
  /* Remember state, looking for places
     where the last base is C and we see
     a G at current place */
  for( pos = 1; pos < aln->h_len; pos++ ) {

    if ( (aln->h1_seq[pos] == 'G') &&
	 (aln->h1_seq[pos-1] == 'C') ) {

      /* Here's one! mask it! */
      aln->mask[pos]   = 0;
      aln->mask[pos-1] = 0;
    }

    if ( (aln->h2_seq[pos] == 'G') &&
	 (aln->h2_seq[pos-1] == 'C') ) {

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
   postions (pos) to the QUADAP aln.
   Dies with error code 1 if the input is bad
*/
void add_fa( char* fn, QUADAP aln, int pos ) {
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
    c = read_fasta_seq( fa, aln->h1_seq );
    return;
  }
  if ( pos == 1 ) {
    c = read_fasta_seq( fa, aln->h2_seq );
    return;
  }
  if ( pos == 2 ) {
    c = read_fasta_seq( fa, aln->n_seq );
    return;
  }
  if ( pos == 3 ) {
    c = read_fasta_seq( fa, aln->c_seq );
    return;
  }
}

void count_diverg( const QUADAP aln, int start, int end, int diverg_counts[] ) {
  char diverg_kmer[5]; /* Holds a string of H1, H2, N, C sequence */
  size_t pos, h, inx;
  char bases[] = { 'A', 'G', 'C', 'T' };

  diverg_kmer[4] = '\0';

  /* Zero the array of diverg_counts */
  for( pos = 0; pos < 256; pos++ ) {
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
	 valid_base( aln->h1_seq[pos] ) &&
	 valid_base( aln->h2_seq[pos] ) &&
	 valid_base( aln->n_seq[pos] ) ) {
      diverg_kmer[0] = aln->h1_seq[pos];
      diverg_kmer[1] = aln->h2_seq[pos];
      diverg_kmer[2] = aln->n_seq[pos];
      diverg_kmer[3] = aln->c_seq[pos];
      if ( kmer2inx( diverg_kmer, 4, &inx ) ) {
	diverg_counts[inx]++;
      }
    }
  }
  return;
}

void output_summary( const QUADAP aln, const int start, const int end ) {
  char diverg_kmer[5]; /* Holds a string of H1, H2, N, C sequence */
  size_t pos, h, inx;
  int diverg_counts[256];
  char bases[] = { 'A', 'G', 'C', 'T' };

  count_diverg( aln, start, end, diverg_counts );
  diverg_kmer[4] = '\0';

  /* If there's an identifier specified, it means we want 
     it printed here */
  if ( aln->identifier[0] != '\0' ) {
    printf( "%s ", aln->identifier );
  }

  /* Output the region */
  printf( "%d %d ", start, end );


  /* Print total number of all the same positions */
  printf( "%d ", (count4kmer( "AAAA", 4, diverg_counts ) +
		  count4kmer( "CCCC", 4, diverg_counts ) +
		  count4kmer( "GGGG", 4, diverg_counts ) +
		  count4kmer( "TTTT", 4, diverg_counts )) );

  /* Print total number of first guy is different positions */
  printf( "%d ", (count4kmer( "CAAA", 4, diverg_counts ) +
		  count4kmer( "GAAA", 4, diverg_counts ) +
		  count4kmer( "TAAA", 4, diverg_counts ) +
		  count4kmer( "ACCC", 4, diverg_counts ) +
		  count4kmer( "GCCC", 4, diverg_counts ) +
		  count4kmer( "TCCC", 4, diverg_counts ) +
		  count4kmer( "AGGG", 4, diverg_counts ) +
		  count4kmer( "CGGG", 4, diverg_counts ) +
		  count4kmer( "TGGG", 4, diverg_counts ) +
		  count4kmer( "ATTT", 4, diverg_counts ) +
		  count4kmer( "CTTT", 4, diverg_counts ) +
		  count4kmer( "GTTT", 4, diverg_counts )) );

  /* Print total number of second guy is different positions */
  printf( "%d ", (count4kmer( "ACAA", 4, diverg_counts ) +
		  count4kmer( "AGAA", 4, diverg_counts ) +
		  count4kmer( "ATAA", 4, diverg_counts ) +
		  count4kmer( "CACC", 4, diverg_counts ) +
		  count4kmer( "CGCC", 4, diverg_counts ) +
		  count4kmer( "CTCC", 4, diverg_counts ) +
		  count4kmer( "GAGG", 4, diverg_counts ) +
		  count4kmer( "GCGG", 4, diverg_counts ) +
		  count4kmer( "GTGG", 4, diverg_counts ) +
		  count4kmer( "TATT", 4, diverg_counts ) +
		  count4kmer( "TCTT", 4, diverg_counts ) +
		  count4kmer( "TGTT", 4, diverg_counts )) );

  /* Print total number of third guy is different positions */
  printf( "%d ", (count4kmer( "AACA", 4, diverg_counts ) +
		  count4kmer( "AAGA", 4, diverg_counts ) +
		  count4kmer( "AATA", 4, diverg_counts ) +
		  count4kmer( "CCAC", 4, diverg_counts ) +
		  count4kmer( "CCGC", 4, diverg_counts ) +
		  count4kmer( "CCTC", 4, diverg_counts ) +
		  count4kmer( "GGAA", 4, diverg_counts ) +
		  count4kmer( "GGCG", 4, diverg_counts ) +
		  count4kmer( "GGTG", 4, diverg_counts ) +
		  count4kmer( "TTAT", 4, diverg_counts ) +
		  count4kmer( "TTCT", 4, diverg_counts ) +
		  count4kmer( "TTGT", 4, diverg_counts )) );

  /* Print total number of fourth guy is different positions */
  printf( "%d ", (count4kmer( "AAAC", 4, diverg_counts ) +
		  count4kmer( "AAAG", 4, diverg_counts ) +
		  count4kmer( "AAAT", 4, diverg_counts ) +
		  count4kmer( "CCCA", 4, diverg_counts ) +
		  count4kmer( "CCCG", 4, diverg_counts ) +
		  count4kmer( "CCCT", 4, diverg_counts ) +
		  count4kmer( "GGGA", 4, diverg_counts ) +
		  count4kmer( "GGGC", 4, diverg_counts ) +
		  count4kmer( "GGGT", 4, diverg_counts ) +
		  count4kmer( "TTTA", 4, diverg_counts ) +
		  count4kmer( "TTTC", 4, diverg_counts ) +
		  count4kmer( "TTTG", 4, diverg_counts )) );

  /* Print total number of BBAA positions */
  printf( "%d ", (count4kmer( "AACC", 4, diverg_counts ) +
		  count4kmer( "AAGG", 4, diverg_counts ) +
		  count4kmer( "AATT", 4, diverg_counts ) +
		  count4kmer( "CCAA", 4, diverg_counts ) +
		  count4kmer( "CCGG", 4, diverg_counts ) +
		  count4kmer( "CCTT", 4, diverg_counts ) +
		  count4kmer( "GGAA", 4, diverg_counts ) +
		  count4kmer( "GGCC", 4, diverg_counts ) +
		  count4kmer( "GGTT", 4, diverg_counts ) +
		  count4kmer( "TTAA", 4, diverg_counts ) +
		  count4kmer( "TTCC", 4, diverg_counts ) +
		  count4kmer( "TTGG", 4, diverg_counts )) );

  /* Print total number of BABA positions */
  printf( "%d ", (count4kmer( "ACAC", 4, diverg_counts ) +
		  count4kmer( "AGAG", 4, diverg_counts ) +
		  count4kmer( "ATAT", 4, diverg_counts ) +
		  count4kmer( "CACA", 4, diverg_counts ) +
		  count4kmer( "CGCG", 4, diverg_counts ) +
		  count4kmer( "CTCT", 4, diverg_counts ) +
		  count4kmer( "GAGA", 4, diverg_counts ) +
		  count4kmer( "GCGC", 4, diverg_counts ) +
		  count4kmer( "GTGT", 4, diverg_counts ) +
		  count4kmer( "TATA", 4, diverg_counts ) +
		  count4kmer( "TCTC", 4, diverg_counts ) +
		  count4kmer( "TGTG", 4, diverg_counts )) );

  /* Print total number of ABBA positions */
  printf( "%d\n", (count4kmer( "ACCA", 4, diverg_counts ) +
		  count4kmer( "AGGA", 4, diverg_counts ) +
		  count4kmer( "ATTA", 4, diverg_counts ) +
		  count4kmer( "CAAC", 4, diverg_counts ) +
		  count4kmer( "CGGC", 4, diverg_counts ) +
		  count4kmer( "CTTC", 4, diverg_counts ) +
		  count4kmer( "GAAG", 4, diverg_counts ) +
		  count4kmer( "GCCG", 4, diverg_counts ) +
		  count4kmer( "GTTG", 4, diverg_counts ) +
		  count4kmer( "TAAT", 4, diverg_counts ) +
		  count4kmer( "TCCT", 4, diverg_counts ) +
		  count4kmer( "TGGT", 4, diverg_counts )) );
}



/* Takes the QUADAP aln and the user-specified window_size
   Prints the counts of same and divergent positions over
   each window. */ 
void output_one_line_windows( const QUADAP aln, int window_size ) {
  int start = 0;
  int end;

  /* Initialize */
  end = start + window_size;
  
  while( end <= aln->h_len ) {
    output_summary( aln, start, end );
    start += window_size;
    end   = start + window_size;
  }
  /* write out the remaining fragment of a window */
  if ( start < aln->h_len ) {
    output_summary( aln, start, aln->h_len );
  } 
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

QUADAP init_QUADAP( void ) {
  QUADAP aln;
  aln = (QUADAP)malloc(sizeof(Quad_Aln));
  aln->h_len  = 0;
  aln->h1_seq  = (char*)malloc((MAX_CHR_LEN+1)*sizeof(char));
  aln->h2_seq  = (char*)malloc((MAX_CHR_LEN+1)*sizeof(char));
  aln->c_seq   = (char*)malloc((MAX_CHR_LEN+1)*sizeof(char));
  aln->n_seq   = (char*)malloc((MAX_CHR_LEN+1)*sizeof(char));
  aln->mask    = (char*)malloc((MAX_CHR_LEN+1)*sizeof(char));

  /* initialize to all Ns */
  memset( aln->h1_seq, 'N', MAX_CHR_LEN );
  memset( aln->h2_seq, 'N', MAX_CHR_LEN );
  memset( aln->c_seq, 'N', MAX_CHR_LEN );
  memset( aln->n_seq, 'N', MAX_CHR_LEN );

  /* initialize mask to all trues */
  memset( aln->mask, 1, MAX_CHR_LEN );

  /* Set the identifier to NULL. Will be set to appropriate
     value if an identifier is given */
  aln->identifier[0] = '\0';

  return aln;
}

/* Mask out regions not mentioned in input file 
   INPUT FILE looks like this:
   chr1 1000 2000
   chr1 2501 34195
   etc.
*/  
void mask_from_fn( char* mask_fn, char* mask, const char* chr_mask,
		   const int mask_seg_width, const int neg_mask ) {
  FILE* f;
  char chr[MAX_ID_LEN];
  char line[MAX_LINE_LEN+1];
  unsigned int start, end, midpoint;
  size_t i;

  /* Since there is a mask, figure if it's a positive or negative
     mask. If negative, everything is in except what's mentioned in 
     the mask. If it's not, everything is out by default */
  if ( neg_mask ) {
    memset( mask, 1, MAX_CHR_LEN );
  }
  else  {
    memset( mask, 0, MAX_CHR_LEN );
  }

  f = fileOpen( mask_fn, "r" );
  while( fgets( line, MAX_LINE_LEN, f ) != NULL ) {
    if ( sscanf( line, "%s %u %u", chr, &start, &end ) 
	 == 3 ) {
      /* The right chromosome? */
      if ( strcmp( chr, chr_mask ) == 0 ) {
	if ( (mask_seg_width == 0) ||
	     (mask_seg_width <= (end-start+1)) ) {
	  if ( neg_mask ) {
	    memset( &mask[start], 0, (end-start+1) );
	  }
	  else {
	    memset( &mask[start], 1, (end-start+1) );
	  }
	}
	else {
	  midpoint = (end + start)/2;
	  start = midpoint - mask_seg_width/2;
	  end = start + mask_seg_width - 1;
	  if ( neg_mask ) {
	    memset( &mask[start], 0, (end-start+1) );
	  }
	  else {
	    memset( &mask[start], 1, (end-start+1) );
	  }
	}
      }
    }
  }
  fclose(f);
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
