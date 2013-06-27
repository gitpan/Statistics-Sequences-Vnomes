use strict;
use warnings;
use Test::More tests => 11;
use constant EPS => 1e-3;

# Tests output of the module by comparing it to published examples of these tests in Tart (197x), Gatlin (1979) and Rushkin et al. (2001) (from NIST) - see the POD for citation:

BEGIN { use_ok('Statistics::Sequences::Vnomes') };

use Statistics::Lite qw(sum);

my $seq = Statistics::Sequences::Vnomes->new();
isa_ok($seq, 'Statistics::Sequences::Vnomes');

my %refdat = (
    tart => {
        observed_sum => 19, # dinomes
        observed_mean => 2.111,
        chisq => 7.053, # for dinomes
        data => [2, 1, 2, 0, 1, 2, 1, 1, 1, 0, 1, 0, 1, 2, 1, 2, 0, 2, 2, 1],
        states => [0, 1, 2],
    },
    nist => {
        psisq => 0.8,
        p_value => 0.67032,
        data => [qw/0 0 1 1 0 1 1 1 0 1/],
        states => [0, 1],
    },
	gatlin => {
        psisq => 36.2909090909091,
        p_value => .45509,
        data => [qw/G A A T A C A G C C T G T C G G T T C T C C G A T G G C A A G T A C T T T A C T G G T T C A G A A T G C G C C C G T A A T C C T T G C A C A G A G G A T C T T A C G C A G T G A A C C G G C T C C G T G G A T T A G C A A C/],
        states => [qw/A C G T/],
    }
);

my $val;

# Tart data test:
eval { $seq->load($refdat{'tart'}->{'data'});};
ok(!$@, $@);

$val = $seq->observed_mean(length => 2, delta => 0, circularize => 0, states => $refdat{'tart'}->{'states'});
ok(equal($val, $refdat{'tart'}->{'observed_mean'}), "observed_mean()  $val = $refdat{'tart'}->{'observed_mean'}");

my @freq = $seq->observed(length => 2, delta => 0, circularize => 0, states => $refdat{'tart'}->{'states'});
my $sum = sum(@freq);
ok(equal($sum, $refdat{'tart'}->{'observed_sum'}), "observed sum of frequences:  $sum = $refdat{'tart'}->{'observed_sum'}");

#ok(equal($val, $refdat{'tart'}->{'chisq'}), "psisq $val = $refdat{'tart'}->{'chisq'}");
   
# Gatlin data:
eval { $seq->load($refdat{'gatlin'}->{'data'});};
ok(!$@, $@);

$val = $seq->psisq(length => 3, delta => 2, circularize => 1, states => $refdat{'gatlin'}->{'states'}, precision_s => 3);
ok(equal($val, $refdat{'gatlin'}->{'psisq'}), "psisq  $val = $refdat{'gatlin'}->{'psisq'}");

$val = $seq->p_value(length => 3, delta => 2, circularize => 1, states => $refdat{'gatlin'}->{'states'}, precision_s => 3);
ok(equal($val, $refdat{'gatlin'}->{'p_value'}), "p_value $val = $refdat{'gatlin'}->{'p_value'}");

# NIST data:
eval { $seq->load($refdat{'nist'}->{'data'});};
ok(!$@, $@);

$val = $seq->psisq(length => 3, delta => 2, states => $refdat{'nist'}->{'states'}, precision_s => 3);
ok(equal($val, $refdat{'nist'}->{'psisq'}), "psisq  $val = $refdat{'nist'}->{'psisq'}");

$val = $seq->p_value(length => 3, delta => 2, states => $refdat{'nist'}->{'states'}, precision_s => 3);
ok(equal($val, $refdat{'nist'}->{'p_value'}), "p_value $val = $refdat{'nist'}->{'p_value'}");

sub equal {
    return 1 if $_[0] + EPS > $_[1] and $_[0] - EPS < $_[1];
    return 0;
}
