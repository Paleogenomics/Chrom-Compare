#include "pileup.h"

inline size_t base_inx( char b ) {
  char B;
  B = toupper(b);
  switch(B) {
  case 'A' :
    return 0;
  case 'C' :
    return 1;
  case 'G' :
    return 2;
  case 'T' :
    return 3;
  default :
    return 4;
  }
}

inline char inx2base( const size_t inx ) {
  switch (inx) {
  case 0 :
    return 'A';
  case 1 :
    return 'C';
  case 2 :
    return 'G';
  case 3 :
    return 'T';
  default :
    return 'N';
  }
}

/* Counts the bases that pass the quality cutoffs
   Returns the single most populous base, unless that is
   an N. If there is a two way tie, returns a random one
   of these two. If there is a three-way tie, returns
   -1.

   Arguments:
   PulP pp -- data structure containing data from one position of pileup
   QcutsP qcp -- coverage cutoffs from QCUT file or default
   unsigned int mqc -- minimum mapping quality; that is, don't consider any base
     with a mapping quality score less than this value
   unsigned int covc -- maximum coverage; i.e., return -1 if position coverage
     is greater than this value
   bool weighted -- true if quality scores should be used to weight each base,
     false otherwise

   -1 signifies no good base
*/
int best_base_from_pul( PulP pp, QcutsP qcp, unsigned int mqc,
        unsigned int covc, bool weighted ) {
  size_t base_counts[5];
  size_t cov, i, this_base_index;
  int rand_inx;
  size_t f_counts, s_counts, t_counts, f_inx, s_inx, t_inx;

  cov = pp->cov; // position coverage

  if ( cov > covc ) {
    return -1;
  }

  if ( cov == 0 ) {
    return -1;
  }

  for( i = 0; i <= 4; i++ ) {
    base_counts[i] = 0;
  }

  for ( i = 0; i < pp->cov; i++ ) {
    if ( (pp->map_quals[i] >= mqc) &&
	 qual_check( toupper( pp->bases[i] ),
		     pp->base_quals[i],
		     pp->strands[i],
		     qcp ) ) {
        this_base_index = base_inx(pp->bases[i]);

        // if not using weights, add one to the base count for this base. if
        // using weights, add base and map quality scores to base_counts
        if (weighted) {
          base_counts[this_base_index] += pp->base_quals[i] + pp->map_quals[i];
        }
        else {
          base_counts[this_base_index]++;
        }
    }
  }

  /* Find f_inx & s_inx - the indices of the most populous and
     second most populous bases */
  f_counts = 0;
  s_counts = 0;
  t_counts = 0;
  for( i = 0; i <= 4; i++ ) {
    if ( base_counts[i] >= f_counts ) {
      f_inx = i;
      f_counts = base_counts[f_inx];
    }
  }

  /* No or only N good counts? */
  if ( (f_counts == 0) ||
       (f_inx == 4) ) {
    return -1;
  }

  for( i = 0; i <= 4; i++ ) {
    if ( (i != f_inx) &&
	 (base_counts[i] >= s_counts) ) {
      s_inx = i;
      s_counts = base_counts[s_inx];
    }
  }
  for( i = 0; i <= 4; i++ ) {
    if ( (i != f_inx) &&
	 (i != s_inx) &&
	 (base_counts[i] >= t_counts) ) {
      t_inx = i;
      t_counts = base_counts[t_inx];
    }
  }

  /* Three-way tie? */
  if ( f_counts == t_counts ) {
    return -1;
  }

  /* Single winner? */
  if ( f_counts > s_counts ) {
    return f_inx;
  }

  /* Tie? flip a rand() coin */
  if ( (rand() % 2) == 0 ) {
    return f_inx;
  }
  else {
    return s_inx;
  }
}

/* Makes a random ordering of the bases at this position.
   Goes through the bases in this random ordering and returns the
   index of the first one that passes the map and base/strand
   cutoff criteria. If none do, or if the coverage cutoff is
   exceeded at this position, then -1 is returned
*/
int rand_good_base_from_pul( PulP pp, QcutsP qcp,
			     unsigned int mqc, unsigned int covc ) {
  unsigned int rand_ord[MAX_COV];
  int tmp_inx, rand_inx;
  size_t cov;
  size_t i;
  cov = pp->cov;

  if ( cov > covc ) {
    return -1;
  }

  /* First set up the random order array to be ordered 1..n */
  for( i = 0; i < cov; i++ ) {
    rand_ord[i] = i;
  }

  /* Now, randomize it */
  for( i = 0; i < cov; i++ ) {
    /* Pick a random position */
    rand_inx = (rand() % cov);
    /* Switch current position with the random position */
    tmp_inx = rand_ord[i];
    rand_ord[i] = rand_ord[rand_inx];
    rand_ord[rand_inx] = tmp_inx;
  }

  /* Now, go through the random ordering of bases and return the
     index of the first one that passes the cutoff criteria */
  for( i = 0; i < cov; i++ ) {
    if ( (pp->map_quals[rand_ord[i]] >= mqc) &&
	 qual_check(toupper( pp->bases[rand_ord[i]] ),
		    pp->base_quals[rand_ord[i]],
		    pp->strands[rand_ord[i]], qcp ) ) {
      return base_inx( pp->bases[rand_ord[i]] );
    }
  }
  return -1;
}

/* line2pul
   Args: char* line : a pileup line of output (with -s flag to samtools view)
         PulP pp    : a pointer to a Pul to be populated with info
	              for this line
   Returns: 0 is everything is copacetic ; 1 if there is a problem
   Ignores (skips past) indels
*/
int line2pul( char* line, PulP pp ) {
  char raw_base_field[MAX_FIELD_WIDTH];
  char raw_base_quals[MAX_FIELD_WIDTH];
  char raw_map_quals[MAX_FIELD_WIDTH];
  int indel_len, i, field_len, cur_base_num;
  char base_code;

  /* First try to just get the first parts of the line to 
     check the coverage. It might be too high and smash the
     stack! */
  if ( sscanf( line, "%s\t%u\t%c\t%u\t",
	       pp->chr,
	       &pp->pos,
	       &pp->ref,
	       &pp->cov ) == 4 ) {
    /* Check to make sure there is not more coverage
       than we can handle */
    if ( pp->cov >= MAX_COV ) {
      return 1;
    }
  }
  else { // couldn't even parse the beginning of this line - give up!
    return 1;
  }

  
  if ( sscanf( line, 
	       "%s\t%u\t%c\t%u\t%s\t%s\t%s",
	       pp->chr,
	       &pp->pos,
	       &pp->ref,
	       &pp->cov,
	       raw_base_field,
	       raw_base_quals,
	       raw_map_quals ) == 7 ) {
    
    /* Newer versions of SAMTools now report the last 3 fields as asterisks
       when there is no coverage (old behavior was to just output first
       four fields). Check for this.
    */
    if (pp->cov == 0 && strcmp(raw_base_field, "*") == 0 && 
        strcmp(raw_base_quals, "*") == 0 && strcmp(raw_map_quals, "*") == 0){
        pp->base_quals[0] = 0;
        pp->map_quals[0] = 0;
        return 0;
    }

    /* Parse the raw_base_field */
    field_len = strlen( raw_base_field );
    i = 0;
    cur_base_num = 0;
    while( i < field_len ) {
      base_code = raw_base_field[i];
      switch( base_code ) {
      case '.' : // Reference base on forward strand
	pp->bases[cur_base_num]   = pp->ref;
	pp->strands[cur_base_num] = 1;
	cur_base_num++;
	i++;
	break;

      case ',' : // Reference base on reverse strand
	pp->bases[cur_base_num]   = pp->ref;
	pp->strands[cur_base_num] = -1;
	cur_base_num++;
	i++;
	break;

      case 'A' :
	pp->bases[cur_base_num]   = 'A';
	pp->strands[cur_base_num] = 1;
	cur_base_num++;
	i++;
	break;

      case 'a' :
	pp->bases[cur_base_num]   = 'A';
	pp->strands[cur_base_num] = -1;
	cur_base_num++;
	i++;
	break;

      case 'C' :
	pp->bases[cur_base_num]   = 'C';
	pp->strands[cur_base_num] = 1;
	cur_base_num++;
	i++;
	break;

      case 'c' :
	pp->bases[cur_base_num]   = 'C';
	pp->strands[cur_base_num] = -1;
	cur_base_num++;
	i++;
	break;

      case 'G' :
	pp->bases[cur_base_num]   = 'G';
	pp->strands[cur_base_num] = 1;
	cur_base_num++;
	i++;
	break;

      case 'g' :
	pp->bases[cur_base_num]   = 'G';
	pp->strands[cur_base_num] = -1;
	cur_base_num++;
	i++;
	break;

      case 'T' :
	pp->bases[cur_base_num]   = 'T';
	pp->strands[cur_base_num] = 1;
	cur_base_num++;
	i++;
	break;

      case 't' :
	pp->bases[cur_base_num]   = 'T';
	pp->strands[cur_base_num] = -1;
	cur_base_num++;
	i++;
	break;

      case 'N' :
	pp->bases[cur_base_num]   = 'N';
	pp->strands[cur_base_num] = 1;
	cur_base_num++;
	i++;
	break;

      case 'n' :
	pp->bases[cur_base_num]   = 'N';
	pp->strands[cur_base_num] = -1;
	cur_base_num++;
	i++;
	break;

      case '-' : // deletion
	i++;
	indel_len = 0;
	base_code = raw_base_field[i];
	while( isdigit(base_code) ) {
	  indel_len *= 10;
	  indel_len += base_code - 48;
	  i++;
	  base_code = raw_base_field[i];
	}
	i += indel_len;
	break;

      case '+' : // insertion
	i++;
	indel_len = 0;
	base_code = raw_base_field[i];
	while( isdigit(base_code) ) {
	  indel_len *= 10;
	  indel_len += base_code - 48;
	  i++;
	  base_code = raw_base_field[i];
	}
	i += indel_len;
	break;

      case '$' : // end of read segment
	i++;
	break;
	
      case '^' : // beginning of read segment
	i += 2; // advance past the map-quality character, too
	break;

      case '*' : // deletion marker, treat it like a base
	pp->bases[cur_base_num] = '*';
	pp->strands[cur_base_num] = 0; // no strand for these
	cur_base_num++;
	i++;
	break;

      default :
	fprintf( stderr, "Cannot parse %c in reads field\n",
		 base_code );
	return 1;
      }
    }
    if (pp->cov > 0){
	    /* Check to see if we got all the bases we were
	       expecting */
	    if ( cur_base_num != pp->cov ) {
	      fprintf( stderr, 
		       "Incorrect number of bases read in: %s\n",
		       line );
	      return 1;
	    }
	    
	    /* Check length of raw_base_quals & raw_map_quals
	       fields */
	    if ( (strlen( raw_base_quals ) != cur_base_num) &&
		 (strlen( raw_map_quals )  != cur_base_num) ) {
	      fprintf( stderr, 
		       "Incorrect number of base or map quals in: %s\n",
		       line );
	      return 1;
	    }
    }

    /* parse raw_base_quals & raw_map_quals field */
    for( i = 0; i < pp->cov; i++ ) {
      pp->base_quals[i] = (raw_base_quals[i] - 33);
      pp->map_quals[i]  = (raw_map_quals[i] - 33);
    }
    return 0;
  }
  
  else {
    return 1;
  }
}

/* Returns a random index of a valid base from this pileup line
   or -1 if the randomly picked base doesn't pass the quality
   cutoffs
 */

int base_inx_from_pul( PulP pp, QcutsP qcp,
		       unsigned int mqc, unsigned int covc  ) {
  int base_num;
  base_num = (rand() % pp->cov);

  if ( (pp->cov <= covc) &&
       (pp->map_quals[base_num] >= mqc) &&
       qual_check( toupper( pp->bases[base_num] ),
		   pp->base_quals[base_num], 
		   pp->strands[base_num], qcp ) ) {
    return base_num;
  }
  else {
    return -1;
  }
}

/* Check the quality score cutoff for this base, strand
   combination; return true if it's good enough; bases
   should be upper case!*/
int qual_check( const char base, const unsigned int base_qual,
		const int strand, QcutsP qcp ) {
  char f_base;
  unsigned int qcut;
  double qp;
  int tsd;

  /* Revcom if minus strand base */
  if ( strand < 0 ) {
    f_base = revcom_base(base);
  }
  else {
    f_base = base;
  }

  switch (f_base) {
  case 'A' :
    qcut = qcp->Aqcut;
    qp   = qcp->Ap;
    break;
  case 'C' :
    qcut = qcp->Cqcut;
    qp   = qcp->Cp;
    break;
  case 'G' :
    qcut = qcp->Gqcut;
    qp   = qcp->Gp;
    break;
  case 'T' :
    qcut = qcp->Tqcut;
    qp   = qcp->Tp;
    break;
  default : // If we can't check it, it fails!
    return 0;
  }

  if (base_qual > qcut ) {
    return 1;
  }
  if ( base_qual < qcut ) {
    return 0;
  }

  /* Exactly at the cutoff - probabilisticall allow it to 
     pass */
  tsd = (rand() % 1000); // roll a thousand-sided die
  
  /* If the roll is higher than the qp * 1000, then it fails;
     otherwise, it passes */
  if ( tsd >= (qp * 1000) ) {
    return 0;
  }
  else {
    return 1;
  }
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

QcutsP dummyQcutsP( void ) {
  QcutsP qcp;
  qcp = (QcutsP)malloc(sizeof(Qcuts));
  qcp->Aqcut = 0;
  qcp->Cqcut = 0;
  qcp->Gqcut = 0;
  qcp->Tqcut = 0;
  qcp->Ap = 100.0;
  qcp->Cp = 100.0;
  qcp->Gp = 100.0;
  qcp->Tp = 100.0;
  return qcp;
}

QcutsP parse_q_score_cut( const char* fn ) {
  FILE* f;
  char line[MAX_LINE_LEN];
  char base[4];
  QcutsP qcp;

  qcp = (QcutsP)malloc(sizeof(Qcuts));

  if (fn[0] != '\0') {
    f = fileOpen( fn, "r" );

    /* First line is the A+T- line */
    fgets( line, MAX_LINE_LEN, f );
    if ( sscanf( line, "%s %u %f",
                 base, &qcp->Aqcut, &qcp->Ap ) != 3 ) {
      fprintf( stderr,
               "Problem parsing Q-score cutoff file line: %s\n",
               line );
      exit( 1 );
    }

    /* First line is the C+G- line */
    fgets( line, MAX_LINE_LEN, f );
    if ( sscanf( line, "%s %u %f",
                 base, &qcp->Cqcut, &qcp->Cp ) != 3 ) {
      fprintf( stderr,
               "Problem parsing Q-score cutoff file line: %s\n",
               line );
      exit( 1 );
    }

    /* First line is the G+C- line */
    fgets( line, MAX_LINE_LEN, f );
    if ( sscanf( line, "%s %u %f",
                 base, &qcp->Gqcut, &qcp->Gp ) != 3 ) {
      fprintf( stderr,
               "Problem parsing Q-score cutoff file line: %s\n",
               line );
      exit( 1 );
    }

    /* First line is the T+A- line */
    fgets( line, MAX_LINE_LEN, f );
    if ( sscanf( line, "%s %u %f",
                 base, &qcp->Tqcut, &qcp->Tp ) != 3 ) {
      fprintf( stderr,
               "Problem parsing Q-score cutoff file line: %s\n",
               line );
      exit( 1 );
    }
    fclose( f );
  }
  else {
    qcp->Aqcut = 0; qcp->Ap = 1.0;
    qcp->Cqcut = 0; qcp->Cp = 1.0;
    qcp->Gqcut = 0; qcp->Gp = 1.0;
    qcp->Tqcut = 0; qcp->Tp = 1.0;
  }

  return qcp;
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


/* Mask out regions not mentioned in input file 
   INPUT FILE looks like this:
   chr1 1000 2000
   chr1 2501 34195
   etc.
*/  
void mask_from_fn( const char* mask_fn, char* mask, const char* chr_mask ) {
  FILE* f;
  char chr[MAX_ID_LEN];
  char line[MAX_LINE_LEN+1];
  unsigned int start, end;
  size_t i;

  /* Since there is a mask, reset everything to false to start */
  memset( mask, 0, MAX_CHR_LEN );
  f = fileOpen( mask_fn, "r" );
  while( fgets( line, MAX_LINE_LEN, f ) != NULL ) {
    if ( sscanf( line, "%s %u %u", chr, &start, &end ) 
	 == 3 ) {
      if ( strcmp( chr, chr_mask ) == 0 ) {
	memset( &mask[start], 1, (end-start+1) );
      }
    }
  }
  fclose(f);
}
