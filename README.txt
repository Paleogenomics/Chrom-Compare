30 May 2014
***chrom-compare***
Programs in this set are used to efficiently read in chromosome-sized,
aligned sequences and compare them. Tri- programs read in 3
chromosomes. Quad- programs read in 4 chromosomes.
The chromosomes must already be aligned for these comparisons to make
any sense. Commonly, one will write a genome from one species "on the
coordinates of another genome".

*** Compiling ***
Simply do this:
> make

*** tri-aln-report ***
This program takes three input fasta files. Each file is meant to
represent the aligned sequence of a chromosome (or contig or scaffold,
etc.). The program examines each position and reports the number of
positions where all bases are the same and the number of positions
where there is a base difference.

To see the program options, run it with no options.

The output can be in one of two flavors. By default, it gives a more
verbose count of how many positions all three sequences had an A, C,
G, and T. And, when one of the two differed, what the difference was
and which of the three sequences had it.

By invoking the -o option, one can make a simplified, machine-parsable
one-line output. This output has the following columns:
1. Start of region (0-indexed)
2. End of region (0-indexed, closed-ended, i.e., the first base NOT in
   the region)
3. Counts of all As
4. Counts of all Cs
5. Counts of all Gs
6. Counts of all Ts
7. A->C where 1st sequence is C (other two are A)
8. A->G where 1st sequence is G (other two are A)
9. A->T where 1st sequence it T (other two are A)
10. C->A where 1st sequence is A (other two are C)
11-18. Continuing pattern above 
19-30. Same pattern as above, but 2nd sequence is different.
31-42. Same pattern as above, but 3rd sequence is different.

An extra, first, column is added by specifying the -I (identifier)
option. This is meant to be an identifier in case this program is
being run over many chromosomes to keep track, in the output, of which
chromsome this is (e.g., chr1).

By invoking the -W (windows) option, one gets the one-line simplified
output, but calculated over adjacent, non-overlapping windows of the
specified width in base-pairs. 

The -C (mask CpG sites) option causes the program to disregard sites
that are CpG in any of the three sequences.

Other options are largly self-explanatory. 

*** tri-window-hist.pl ***
This program takes the simplified, one-line output of tri-aln-report
and generates a histogram of divergences. This is a special
divergence. It is the ratio of the 2nd sequence's unique sites to half
the sum of the 1st and 2nd sequence's unique sites. Note, this ignores
the third sequence's unique sites. The idea is that the third sequence
is typically a pseudo-haploid sequence from shotgun data that has many
errors. The first sequence is typically an outgroup sequence. The
divergence reported this way can be stated as how far back the TMRCA
of the second and third sequence is as a percentage of how far back
the TMRCA of the second and first sequence is.

The format of the output is is a space-delimited table with the
following columns:
1. Divergence bin (<= this percent divergence)
2. Number of regions in the bin
3. Fraction of all regions in this bin
4. Cummulative fraction of all regions in the and lower bins
 
*** quad-aln-report ***
This program is similar in concept to tri-aln-report. However, it does
not report the full spectrum of differences between each of the four
input files. Instead, the user is meant to give for option -4 the
ancestral sequence or an outgroup sequence. Defining the base present
in this 4th sequence as "A" (for ancestral) and any derived allele as
"B" (for not "A"), the output is in this format:
START END AAAA BAAA ABAA AABA BBBA BBAA BABA ABBA

The same -W and -I options are available as in tri-aln-report to
specify windows and an identifier to include in the output.

*** find-ancestral ***
find-ancestral -a <ancestral fasta> -f <focal fasta> -C [mask CpG sites]
            -i [chromosome ID to put in first column] <aligned chrs...>
            -m [minimum number of sequences that must be different; default = 1]
            [aligned chrs...]
Scans the input chromosomes to find and report positions where the focal
input sequence matches the ancestral input sequence and all other input
sequences match one another but are different. This is meant to find
positions where the focal sequence is an outgroup to all other input
sequences.
