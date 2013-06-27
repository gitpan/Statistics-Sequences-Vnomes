#!perl -T

use Test::More tests => 1;

BEGIN {
    use_ok( 'Statistics::Sequences::Vnomes', 0.10 ) || print "Bail out!\n";
}

diag( "Testing Statistics::Sequences::Vnomes $Statistics::Sequences::Vnomes::VERSION, Perl $], $^X" );
