#!/usr/local/bin/perl

use Getopt::Std;
use vars qw( $opt_w $opt_m );
use strict;

&init();

my ( $l, $first_diff_total, $second_diff_total, $bin, $cu );
my $t = 0; # total bins with any data
my ( @counts, @hist );

open( WI, $opt_w ) or die( "$!: opt_w\n" );
while( chomp( $l = <WI> ) ) {
    ### Skip comments and blank lines ###
    if ( $l =~ /^\#/ ) {
	next;
    }
    if ( $l eq "" ) {
	next;
    }
    @counts = split( " ", $l );

    if ( $#counts == 42 ) { # Extra Identifier we don't care about
	shift( @counts );
    }

    $first_diff_total  = &sum( @counts[6..17] );
    $second_diff_total = &sum( @counts[18..29] );

    if ( (($first_diff_total + $second_diff_total) > 0) &&
	 (($first_diff_total + $second_diff_total) > $opt_m) ) {
	$bin = int(100 * ($second_diff_total/(($first_diff_total + $second_diff_total)/2)));
	$hist[$bin]++;
    }
    $t++;
}
close( WI );

$cu = 0; # cummulative number of bins
for( $bin = 0; $bin <= 100; $bin++ ) {
    $cu += $hist[$bin];
    printf( "%d %d %.5f %.5f\n", $bin, $hist[$bin], $hist[$bin]/$t, $cu/$t );
}

sub sum {
    my $t = 0;
    map { $t += $_ } @_;
    return $t;
}

sub init {
    my $m_DEF = 50;
    getopt( 'w:m:' );
    unless( -f $opt_w ) {
	print( "tri-window-hist.pl -w <window file> -m <minimum counts; DEF=$m_DEF>\n" );
	print( "Makes a histogram of percentage divergence for all windows\n" );
	print( "reported in the input file. This divergence is 2nd/((1st+2nd)/2)\n" );
	print( "Rounded down to the nearest percentile.\n" );
	exit( 0 );
    }
    unless( defined( $m_DEF ) ) {
	$opt_m = $m_DEF;
    }
}
