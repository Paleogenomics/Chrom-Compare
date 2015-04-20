#ifndef PILEUP_H
#define PILEUP_H

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#define MAX_FIELD_WIDTH (10240)
#define MAX_COV (4096)
#define MAX_LINE_LEN (10240)
#define MAX_ID_LEN (256)
#define MAX_CHR_LEN (250000000)

typedef struct qcuts {
  unsigned int Aqcut;
  unsigned int Cqcut;
  unsigned int Gqcut;
  unsigned int Tqcut;
  float Ap;  
  float Cp;
  float Gp;
  float Tp;
} Qcuts;
typedef struct qcuts* QcutsP;

typedef struct pul {
  char chr[ MAX_FIELD_WIDTH ];
  unsigned int pos;
  char ref;
  unsigned int cov;
  char bases[MAX_COV];
  size_t best_alt_inx;
  unsigned int base_quals[MAX_COV];
  unsigned int map_quals[MAX_COV];
  int strands[MAX_COV];
} Pul;
typedef struct pul* PulP;

typedef struct subcount {
  int refA[4];
  int refC[4];
  int refG[4];
  int refT[4];
} Subcount;
typedef struct subcount* SCP;

inline char revcom_base( const char base );
inline size_t base_inx( char b );
int best_base_from_pul( PulP pp, QcutsP qcp,
			unsigned int mqc, unsigned int covc );
int rand_good_base_from_pul( PulP pp, QcutsP qcp,
			     unsigned int mqc, unsigned int covc );

int line2pul( char* line, PulP pp );
int base_inx_from_pul( PulP pp, QcutsP qcp,
		       unsigned int mqc, unsigned int covc  );
inline int valid_base( char b );
int qual_check( const char base, const unsigned int base_qual,
		const int strand, QcutsP qcp );
QcutsP parse_q_score_cut( const char* fn );
QcutsP dummyQcutsP( void );
FILE * fileOpen(const char *name, char access_mode[]);
void mask_from_fn( const char* mask_fn, char* mask, const char* chr_mask );

#endif /* PILEUP_H */
