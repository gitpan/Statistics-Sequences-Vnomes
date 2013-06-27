package Statistics::Sequences::Vnomes;

use 5.008008;
use strict;
use warnings;
use Carp 'croak';
use vars qw($VERSION @ISA);
use Algorithm::Combinatorics qw(variations_with_repetition);
use Math::Cephes 0.43;
use Statistics::Sequences 0.11;
use Statistics::Lite qw(mean stddev sum);
@ISA = qw(Statistics::Sequences);

$VERSION = '0.10';

=pod

=head1 NAME

Statistics::Sequences::Vnomes - The Serial Test (psi-square) and Generalized Serial Test (delta psi-square) for equiprobability of I<v>-nomes (or I<v>-plets/bits) (Good's and Kendall-Babington Smith's tests)

=head1 SYNOPSIS

 my $vnomes = Statistics::Sequences::Vnomes->new();
 $vnomes->load(qw/1 0 0 0 1 1 0 1 1 0 0 1 0 0 1 1 1 1 0 1/); # data treated categorically - could also be, e.g., 'a' and 'b'
 my $freq_href = $vnomes->observed(length => 3); # returns hashref of frequency distribution for trinomes in the sequence
 my @freq = $vnomes->observed(length => 3); # returns only the observed trinome frequencies (not keyed by the trinomes themselves)
 my $val = $vnomes->observed_mean(length => 3); # mean observed frequencies (2.5); also get their SD by 'wanting' an array
 $val = $vnomes->expected(length => 3); # mean chance expectation for the frequencies (2.5)
 $val = $vnomes->psisq(length => 3); # Good's "second backward differences" psi-square (3.4); option 'delta' gives alternative estimates
 $val = $vnomes->p_value(length => 3, tails => 1); # 1-tailed p-value for psi-square (0.0913)
 $val = $vnomes->z_value(length => 3, tails => 1, ccorr => 1); # inverse-phi of the 1-tailed p-value (1.333)
 my $href = $vnomes->stats_hash(length => 3, values => {observed => 1, p_value => 1}, tails => 1); # include any stat-method (& their options)
 $vnomes->dump(length => 3,
  values => {observed => 1, observed_mean => 1, expected => 1, psisq => 1, p_value => 1}, # what stats to show (or not if => 0)
  format => 'table', flag => 1, precision_s => 3, precision_p => 7, verbose => 1, tails => 1);
  # prints:
 # Vnomes (3) statistics
 #.-----------+-----------+-----------+-----------+-----------.
 #| observed  | expected  | p_value   | psisq     | observed- |
 #|           |           |           |           | _mean     |
 #+-----------+-----------+-----------+-----------+-----------+
 #| '011' = 4 | 2.500     | 0.0913418 | 3.400     | 2.500     |
 #| '010' = 1 |           |           |           |           |
 #| '111' = 2 |           |           |           |           |
 #| '000' = 1 |           |           |           |           |
 #| '101' = 2 |           |           |           |           |
 #| '001' = 3 |           |           |           |           |
 #| '100' = 3 |           |           |           |           |
 #| '110' = 4 |           |           |           |           |
 #'-----------+-----------+-----------+-----------+-----------'
 # these and other methods inherited from Statistics::Sequences and its parent, Statistics::Data:
 $vnomes->dump_data(delim => ','); # a comma-separated single line of the loaded data
 $vnomes->save(path => 'seq.csv'); # serialized for later retrieval by open() method

=head1 DESCRIPTION

Implements tests of the independence of successive elements of a sequence of data: serial tests for I<v>-nomes (a.k.a I<v>-plets or, for binary data, I<v>-bits) - singlets/monobits, dinomes/doublets, trinomes/triplets, etc.. Test are of variations in sub-sequence length I<v> that are equally likely for the sampled sequence. For example, a sequence sampled from a "heads'n'tails" (H and T) distribution can be tested for its equal representation of the trinomes HTH, HTT, TTT, THT, and so on. Counting up these I<v>-nomes at all points in the sequence, permitting overlaps, yields a statistic - psi-square - that is approximately distributed as I<chi>-square; the Kendall-Babington Smith statistic.

However, because these counts are not independent (given the overlaps), Good's Generalized Serial Test is the default test-statistic returned by this module's L<test|test> routine: It computes psi-square by differencing, viz., in relation to not only the specified B<length>, or value of I<v>, but also its value for the first two prior lengths of I<v>, yielding a statistic, delta-square-psi-square (the "second backward difference" measure) that is I<exactly> distributed as I<chi>-square.

The test is suitable for multi-state data, not only the binary, dichotomous sequence suitable for the L<runs|Statistics::Sequences::Runs> or L<joins|Statistics::Sequences::Joins> tests. It can also be used to test that the individual elements in the list are uniformly distributed, that the states are equally represented, i.e., as a I<chi>-square-based frequency test (a.k.a. test of uniformity, equiprobability, equidistribution). This is done by giving a B<length> of 1, i.e., testing for mononomes.

Note that this is I<not> the so-called serial test described by Knuth (1998, Ch. 2), which involves non-overlapping pairs of sequences.

=head1 METHODS

Methods include those described in L<Statistics::Sequences>, which can be used directly from this module, as follows.

=head2 new

 $vnomes = Statistics::Sequences::Vnomes->new();

Returns a new Vnomes object. Expects/accepts no arguments but the classname.

=head2 load

 $vnomes->load(@data); # anonymously
 $vnomes->load(\@data);
 $vnomes->load('sample1' => \@data); # labelled whatever

Loads data anonymously or by name - see L<load|Statistics::Data/load, add, unload> in L<Statistics::Data> for ways data can be loaded and retrieved (more than shown here). Every load unloads all previous loads and any additions to them.

Data for this test of sequences can be categorical or numerical, all being treated categorically. Also, the data do not have to be dichotomous (unlike in tests of L<runs|Statistics::Sequences::Runs> and L<joins|Statistics::Sequences::Joins>.

=head2 add, read, unload

See L<Statistics::Data> for these and other operations on loaded data.

=head2 observed, vnomes_observed, vno

 $href = $vnomes->observed(length => n, circularize => 1|0); # returns keyed distribution; assumes data have already been loaded
 @ari = $vnomes->observed(data => [1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1], length => n, circularize => 1|0); # returns frequencies only

Returns the frequency distribution for the I<observed> number of vnomes of the given B<length> in a sequence. These are counted up as overlapping. Called in array context, returns an array of the counts for each possible pattern; otherwise, where these are keyed by the pattern in a hash reference. A value for B<length> greater than zero is required, and must be no more than the sample-size. So for the sequence

 1 0 0 0 1 0 0 1 0 1 1 0 

there are 12 vnomes of length 1 (mononomes, the number of elements) from the two (0 and 1) that are possible; 11 dinomes from the four (10, 00, 00, 01) that are possible; and 10 trinomes from eight (100, 000, 001, 010, etc.) that are possible.

If the options B<circularize> equals 1 (default), the count continues to loop to the beginning of the sequence until all elements from the end are included. So instead of ending the count for trinomes at (110), the count will also include (101) and (010), increasing the observed value of trinomes to 12.

The data to test can already have been L<load|Statistics::Sequences/load, add, unload>ed, or you send it directly keyed as B<data>.

=cut

sub observed {
    my $self = shift;
    my $args = ref $_[0] ? shift : {@_};
    my $data = ref $args->{'data'} ? $args->{'data'} : $self->read($args);
    my $num = scalar(@$data);
    my $v = $args->{'length'} || croak __PACKAGE__, '::test Must have a v-nome length greater than zero';
    # Sanitize length in proportion to size:
    croak __PACKAGE__, "::test Sequence length v ($v) for testing should be no more than the sample-size ($num) - 2" if $v >= $num - 1;
    #my $lim =  Math::Cephes::floor(Math::Cephes::log2($num)) - 2; # NIST recommmendation
    # Could warn here if $v is too small for $num and $nstates if $v >= $lim ...
    my $circularize = defined $args->{'circularize'} ? $args->{'circularize'} : 1;
    my ($nstates, $states_aref) = _set_states($data, $args->{'states'});# Get list of unique states as given or from sequence itself
    my $v_aref = _get_v_ari($v);
    my ($v_i, @data_i, %stats) = ();
    foreach $v_i(@$v_aref) {# print "v_i = $v_i\n"; # loop through tests of v, v-1 and v-2 lengths ...
        @data_i = @$data;# extend sequence by appending the first $v_i - 1 elements to its end
        push(@data_i, @data_i[0 .. $v_i - 2]) if $circularize; # e.g., if v = 3, append elements [0] & [1]; if v = 2, append only 1st value
        $stats{$v_i} = _frequencies(\@data_i, $v_i, $states_aref);# count up frequencies of each sequence of length $v_i
        last if $v_i == $v;
   }
   return wantarray ? values %{$stats{$v}} : $stats{$v}; 
}
*vnomes_observed = \&observed;
*vno = \&observed;

=head2 observed_mean

 $mean = $vnomes->observed_mean(length => n, data => \@data); # options as for observed()
 ($mean, $stdev) = $vnomes->observed_mean(length => n, data => \@data); # options as for observed()

Returns the mean of the observed frequencies for the v-nomes of the given length. Called in array context, the mean and also the standard deviation of the frequencies are returned. Options as for the L<observed|Statistics::Sequences::Vnomes/observed, vnomes_observed, vno> method. For other descriptives of the observed frequencies, get the array returned from calling L<observed|Statistics::Sequences::Vnomes/observed, vnomes_observed, vno> and use L<Statistics::Lite|Statistics::Lite> methods, which this method depends on.

=cut

sub observed_mean {
    my $self = shift;
    my $args = ref $_[0] ? shift : {@_};
    my @freq = observed($self, $args);
    my $v = $args->{'length'}; # || croak __PACKAGE__, '::test Must have a v-nome length greater than zero';
    return wantarray ? ( mean(@freq), stddev(@freq) ) : mean(@freq);
}

=head2 expected, vnomes_expected, vne

 $count = $vnomes->expected(length => n, circularize => 1|0); # assumes data have already been loaded
 $count = $vnomes->expected(data => [1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1], length => n, circularize => 1|0, states => [0, 1]);

Returns the I<expected> number of observations for each vnome of the given B<length> in a sequence; i.e., the mean chance expectation, assuming that each event is generated by a random uniform process. Options are as for L<observed|Statistics::Sequences::Vnomes/observed>. This expected frequency is given by:

=for html <p>&nbsp;&nbsp;<i>E[V]</i> = <i>N</i><i>k</i><sup>&ndash;<i>v</i></sup>

where I<k> is the number of possible states (alternatives; as read from the data or as explicitly given as an array in B<states>), I<v> is the vnome length, and I<N> is the length of the sequence, less I<v> + 1 if the count is not to be circularized.

Another way to think of this is as the number of mononome observations in the sequence divided by the number of possible permutations of its states for the given length. So, for a sequence made up of 0s and 1s, there are four possible variations of length 2 (00, 10, 01 and 11), so that the expected frequency for each of these variations in a sequence of 20 values is  20 / 4, i.e., 5.

=cut

sub expected {
    my $self = shift;
    my $args = ref $_[0] ? shift : {@_};
    my $data = ref $args->{'data'} ? $args->{'data'} : $self->read($args);
    my $num = scalar(@$data);
    my $v = $args->{'length'} || croak __PACKAGE__, '::test Must have a v-nome length greater than zero';
    # Sanitize length in proportion to size:
    croak __PACKAGE__, "::test Sequence length v ($v) for testing should be no more than the sample-size ($num) - 2" if $v >= $num - 1;
    #my $lim =  Math::Cephes::floor(Math::Cephes::log2($num)) - 2; # NIST recommmendation
    # Could warn here if $v is too small for $num amd $nstates if $v >= $lim ... At present, only use $lim for "verbose" dumps.
    my $circularize = defined $args->{'circularize'} ? $args->{'circularize'} : 1;
    my ($nstates, $states_aref) = _set_states($data, $args->{'states'});# Get list of unique states as given or from sequence itself
    my $v_aref = _get_v_ari($v);
    # While looping through tests of v, v-1 and v-2 (= $v_i) sequence lengths ...
    my ($v_i, $n_i, $x_i, $nx_i, $freq, @data_i, %stats) = ();
    foreach $v_i(@$v_aref) { # print "v_i = $v_i\n";
        @data_i = @$data;# extend sequence by appending the first $v_i - 1 elements to its end
        push(@data_i, @data_i[0 .. $v_i - 2]) if $circularize; # e.g., if v = 3, append elements [0] & [1]; if v = 2, append only 1st value
        $n_i = scalar(@$data);
        # Compute expected number of any form of sequences of length $v_i:
        $nx_i = $circularize ? $n_i : ($num - $v_i + 1);
        $x_i = $nx_i * ( $nstates**(-1 * $v_i) );
        $stats{$v_i} = $x_i;
    }
    return $stats{$v};
}
*vnomes_expected = \&expected;
*vne = \&expected;

=head2 variance

(Not defined for this package.)

=cut

sub variance {
   croak 'variance() method is not defined for ' . __PACKAGE__;
}

=head2 obsdev, observed_deviation

(Not defined for this package.)

=cut

sub obsdev {
   croak 'obsdev() method is not defined for ' . __PACKAGE__;
}
*observed_deviation = \&obsdev;

=head2 stdev, standard_deviation

(Not defined for this package.)

=cut

sub stdev {
   croak 'stdev() method is not defined for ' . __PACKAGE__;
}
*standard_deviation = \&stdev;

=head2 psisq

 $vnomes->psisq(length => n, delta => '2|1|0', circularize => '1|0', states => [qw/A C G T/]);

Performs Good's Generalized Serial Test (by default), of I<v>-nomes on the given or named distribution, yielding a I<psi-square> statistic. Actually, this is I<not> the raw psi-square value for sub-sequences of length I<v> (Kendall-Babington Smith statistic), because, unless I<length> (I<v>) = 1, this value is not asymptotically distributed as I<chi>-square. Instead, this method returns, by default (without specifiying a value for B<delta>, or if B<delta> is not 1 or 0) the "second backward differences" psi-square (I<delta^2-psi^2>). This incorporates psi-square values for backwardly adjacent values of I<length> (I<v>), i.e., for sub-sequences of length I<v>, I<v> - 1, and I<v> - 2. It is not only asymptotically I<chi>-square distributed, but uses statistically independent counts of all the possible variations of sequences of the given B<length> (Good, 1953).

To get the "I<first> backward differences" of psi-square, which is the difference between the psi-square values for sub-sequences of length I<v> and length I<v> - 1, specify B<delta> => 1. While it is I<chi>-square distributed, counts of first-differences are not statistically independent (Good, 1953; Good & Gover, 1967), and "the sequence of second differences forms a much better set of statistics for testing the hypothesis of flat-randomness" (Good & Gover, 1967, p. 104). To incoroporate no backward differences in the calculation, specifiy B<delta> => 0. Leave B<delta> undefined to use the second-backward differences.

The implemented algorithm is that given by Good (1953, Eq. 1); benchmarking shows no reliable speed difference to alternative forms of the equation, as given by Good (1957, Eq. 2) and in the NIST test suite (Rukhin et al., 2001). Good's original algorithm can also be found in individual papers describing the application of the Serial Test (e.g., Davis & Akers, 1974).

By default, the I<p>-value associated with the test-statistic is 2-tailed. See the L<Statistics::Sequences|Statistics::Sequences/test> manpage for generic options other than the following Vnome test-specific ones.

=cut

sub psisq {
    my $self = shift;
    my $args = ref $_[0] ? shift : {@_};
    my $data = ref $args->{'data'} ? $args->{'data'} : $self->read($args);
    my $num = scalar(@$data);
    my $v = $args->{'length'} || croak __PACKAGE__, '::test Must have a v-nome length greater than zero';
    # Sanitize length in proportion to size:
    croak __PACKAGE__, "::test Sequence 'length' ($v) for testing should be no more than the sample-size ($num) - 2" if $v >= $num - 1;
    #my $lim =  Math::Cephes::floor(Math::Cephes::log2($num)) - 2; # NIST recommmendation
    # Could warn here if $v is too small for $num amd $nstates if $v >= $lim ... At present, only use $lim for "verbose" dumps.
    my $circularize = defined $args->{'circularize'} ? $args->{'circularize'} : 1;
    my ($nstates, $states_aref) = _set_states($data, $args->{'states'});# Get list of unique states as given or from sequence itself
    my $v_aref = _get_v_ari($v);
    # While looping through tests of v, v-1 and v-2 (= $v_i) sequence lengths ...
    my ($v_i, $n_i, $x_i, $nx_i, $freq, @data_i, %stats) = ();
    foreach $v_i(@$v_aref) {
        @data_i = @$data;# extend sequence by appending the first $v_i - 1 elements to its end
        push(@data_i, @data_i[0 .. $v_i - 2]) if $circularize; # e.g., if v = 3, append elements [0] & [1]; if v = 2, append only 1st value
        $n_i = scalar(@$data);
        # Count up frequencies of each sequence of length $v_i - used to get observed()
        $freq = _frequencies(\@data_i, $v_i, $states_aref);
        # Compute expected number of any form of sequences of length $v_i:
        $nx_i = $circularize ? $n_i : ($num - $v_i + 1);
        $x_i = $nx_i * ( $nstates**(-1 * $v_i) );
        #print "$v_i: nstates**vi/vi = ", (($nstates**$v_i) / $n_i) , "\n"; "$v_i: sum = ", sum( map{ ($_ - $x_i)**2 } values %{$freq}), "\n";
        # Compute psi^2 itself for this v-nome length ($v_i) (Good, 1953, Eq. 1):
        $stats{$v_i} = ( ($nstates**$v_i) / $n_i) * sum( map{ ($_ - $x_i)**2 } values %{$freq});
    }
    my $delta = defined $args->{'delta'} ? $args->{'delta'} : 2;
    my ($psisq, $df) = _get_psisq_and_df($v, $nstates, \%stats, $delta);
    return wantarray ? ($psisq, $df) : $psisq;
}

=head2 z_value, vnomes_zscore, vzs, zscore

 $val = $vnomes->z_value(); # data already loaded, use default windows and prob
 $val = $vnomes->z_value(data => $aref);
 ($zvalue, $pvalue) =  $vnomes->z_value(data => $aref, tails => 2); # same but wanting an array, get the p-value too

Returns the zscore not from a direct test of deviation but from the (inverse-phi) of the L<p_value|Statistics::Sequences::Vnomes/p_value>, using L<Math::Cephes|Math::Cephes> C<ndtri>. Called in array context, also returns the I<p>-value itself. Same options and conditions as above. Other options are B<precision_s> (for the z_value) and B<precision_p> (for the p_value).

=cut

sub z_value {
   my $self = shift;
   my $args = ref $_[0] ? shift : {@_};   
   my $pval = $self->p_value($args);
   require Statistics::Zed;
   my $zed = Statistics::Zed->new();
   $args->{'tails'} ||= 1;
   my $zval = $zed->p2z(value => $pval, tails => $args->{'tails'});
   return wantarray ? ($zval, $pval) : $zval;
}
*vzs = \&z_value;
*vnomes_zscore = \&z_value;
*zscore = \&z_value;

=head2 p_value, test, vnomes_test, vnt

 $p = $vnomes->p_value(); # using loaded data and default args
 $p = $vnomes->p_value(data => [1, 0, 1, 1, 0], exact => 1); #  using given data (by-passing load and read)
 $p = $vnomes->p_value(trials => 20, observed => 10); # without using data

Returns probability of obtaining the psisq value for data already L<load|Statistics::Sequences/load, add, unload>ed, or directly keyed as B<data>. The I<p>-value is read off the complemented chi-square distribution (incomplete gamma integral) using L<Math::Cephes|Math::Cephes> C<igamc>.

=cut

sub p_value {
   my $self = shift;
   my $args = ref $_[0] ? shift : {@_};
   my ($psisq, $df, $p_value) = $self->psisq($args);
   $p_value = Math::Cephes::igamc($df/2, $psisq/2);
   $p_value /= 2 if defined $args->{'tails'} and $args->{'tails'} == 1;
   return $p_value;
}
*test = \&p_value;
*vnomes_test = \&p_value;
*vnt = \&p_value;

=head2 dump

 $vnomes->dump(length => 3, values => {psisq => 1, p_value => 1}, format => 'table|labline|csv', flag => 1, precision_s => 3, precision_p => 7, verbose => 1, tails => 1);

Print Vnome-test results to STDOUT. See L<dump|Statistics::Sequences/dump> in the Statistics::Sequences manpage for details. For comparability with other modules in the Statistics::Sequences package, the I<Z>-value associated (by Math::Cephes::ndtri) with the obtained I<p>-value is reported. If B<text> => 2, then you get a verbose dump, including (1) the actual test-statistic depending on the value of B<delta> tested (I<delta2_psi2> for the second difference measure (default), I<delta_psi2> for the first difference measure, and I<psi2> for the raw measure), followed by degrees-of-freedom in parentheses; and (2) a warning, if relevant, that your B<length> value might be too large with respect to the sample size (see NIST reference, above, in discussing B<length>). If B<text> => 1, you just get the average observed and expected frequencies for each I<v>-nome, the I<Z>-value, and its associated I<p>-value.

=cut

sub dump {
    my $self = shift;
    my $args = ref $_[0] ? $_[0] : {@_};
    $args->{'stat'} = 'vnomes';
    $args->{'stat'} .= " ($args->{'length'})" if defined $args->{'length'};
    $self->SUPER::dump($args);
    return $self;
}

=head1 OPTIONS

Options common to the above stats methods.

=head2 length

The length of the I<v>-nome, i.e., the value of I<v>. Must be an integer greater than or equal to 1, and smaller than than the sample-size.

What is a meaningful maximal value of B<length>? As a I<chi>-square test, it could be held that there should be an expected frequency of at least 5 for each I<v>-nome. This is "conventional wisdom" recommended by Knuth (1988) but can be judged to be too conservative (Delucchi, 1993). The NIST documentation on the serial test (Rukhin et al., 2001) recommends that length should be less than the floored value of log2 of the sample-size, minus 2. No tests are here made of these recommendations.

=head2 circularize

By default, circularizes the data sequence; i.e., the datum after the last element is the first element. This affects (and slightly simplifies) the calculation of the expected frequency of each I<v>-nome, and so the value of each psi-square. Circularizing ensures that the expected frequencies are accurate; otherwise, they might only be approximate. As Good and Gover (1967) offer, "It is convenient to circularize in order to get exact checks of the arithmetic and also in order to simplify some of the theoretical formulae" (p. 103).

=head2 states

A referenced array listing the unique states (or 'events' or 'letters') in the population from which the sequence was sampled. This is useful to specify if the sequence itself is likely not to include all the possible states; it might even include only one of them. If this array is not specified, the unique states are identified from the sequence itself. If giving a list of states, a check in each test is made to ensure that the sequence contains I<only> those elements in the list.

=cut

# PRIVATMETHODEN

sub _set_states {
    my ($data, $states) = @_;
    if (! ref $states) { # Get states from the data themselves:
        my %hash = map { $_, 1 } @{$data};
        $states = [keys %hash];
    }
    else { # Ensure that the data only contain states in the given list:
        my ($g, $h) = ();
        DATA:
        foreach $g(@{$data}) {
            foreach $h(@{$states}) {
                if ($h eq $g) {
                    next DATA;
                }
            }
            croak __PACKAGE__, "::test The element $g in the data is not represented in the given states"; 
        }
    }
    my $nstates = scalar(@{$states});
    #croak __PACKAGE__, '::test At least two different values must be in the sequence to test its sub-sequences' if $nstates <= 1;
    return ($nstates, $states);
}

sub _frequencies {
    my ($data_i, $v_i, $states_aref) = @_;

    # Get a list of all possible combinations of states at the current length ($v_i):
    my @variations = variations_with_repetition($states_aref, $v_i);

    # Count up the frequency of each variation in the data:
    my ($i, $probe_str, $test_str, %freq) = ();
    foreach (@variations) {
        $probe_str = join'', @{$_};
        $freq{$probe_str} = 0;
        for ($i = 0; $i < scalar(@{$data_i}) - $v_i + 1; $i++) {
            $test_str = join'', @{$data_i}[$i .. ($i + $v_i - 1)];
            $freq{$probe_str}++ if $probe_str eq $test_str;
        }
    }#    print "FREQ:\n"; while (my($key, $val) = each %freq) { print "\t$key = $val\n"; } print "\n";
    return \%freq;
}
 
sub _get_psisq_and_df {# get psisq and its df
    my ($v, $nstates, $stats, $delta) = @_;
    my ($psisq, $df) = ();
    
    if ($v == 1) { # psisq is asymptotically distributed chisq, can use psisq for chisq distribution:
        $psisq = $stats->{1};
        $df = $nstates - 1;
    }
    else {
        if ($delta == 0) { # Raw psisq:
            $psisq = $stats->{$v};
            $df = $nstates**$v - 1;
        }
        elsif ($delta == 1) { # First backward difference:
            $psisq = $stats->{$v} - ($v - 1 <= 0 ? 0 : $stats->{$v - 1});
            $df = $nstates**$v - $nstates**($v - 1);
        }
        else { # $delta == 2 # Second backward difference (default):
            $psisq = $stats->{$v} - ( 2 * ($v - 1 <= 0 ? 0 : $stats->{$v - 1}) ) + ($v - 2 <= 0 ? 0 : $stats->{$v - 2});
            $df = ( $nstates**($v - 2) ) * ( $nstates - 1)**2;
        }
    }
    return ($psisq, $df);
}

sub _get_v_ari {
    my $v = shift;    # Init a hash to keep the psi-square values for the v, v-1, and v-2 ( = $v_i) sequence lengths, where relevant:
    my @ari = ();
    foreach (0 .. 2) {
        $v > $_ ? push @ari, $v - $_ : last;
    } # print "v ari = ", join(' ', @ari), "\n";  #push @ari, $v - 1 if $v >= 2; #push @ari, $v - 2 if $v >= 3;
    return \@ari;    
}

__END__

=head1 EXAMPLE

=head2 Seating at the diner

This is the data from Swed and Eisenhart (1943) also given as an example for the L<Runs test|Statistics::Sequences::Runs/EXAMPLE> and L<Turns test|Statistics::Sequences::Turns/EXAMPLE>. It lists the occupied (O) and empty (E) seats in a row at a lunch counter.
Have people taken up their seats on a random basis - or do they show some social phobia, or are they trying to pick up? What does the test of Vnomes reveal?

 use Statistics::Sequences::Vnomes;
 my $vnomes = Statistics::Sequences::Vnomes->new();
 my @seating = (qw/E O E E O E E E O E E E O E O E/);
 $vnomes->load(\@seating);
 $vnomes->dump(length => 3, values => {z_value => 1, p_value => 1}, format => 'labline', flag => 1, precision_s => 3, precision_p => 3, tails => 1);

This prints: 

 z_value = 2.015, p_value = 0.022*

That is, the observed frequency of each possible trio of seating arrangements (the dinomes OOO, OOE, OEE, EEE, etc.) differed significantly from that expected. Look up the observed frequencies for each possible trinome to see if this is because there are more empty or occupied neighbouring seats (phobia or philia):

 $vnomes->dump(length => 3, values => {observed => 1}, format => 'labline');

This prints:

 observed = ('OEE' = 4,'EEO' = 4,'EEE' = 2,'OEO' = 1,'EOE' = 5,'OOO' = 0,'OOE' = 0,'EOO' = 0)

As the chance-expected frequency is 2.5 (using the L<expected|Statistics::Sequences::Vnomes/expected, vnomes_expected, vne> method), there are clearly more than expected trinomes involving empty seats than occupied seats - suggesting a non-random factor like social phobia or body odour is at work in sequencing people's seating here. Noting that the sequencing isn't significant for dinomes (with B<length> => 2) might also tell us something about what's going on.

=head1 REFERENCES

Davis, J. W., & Akers, C. (1974). Randomization and tests for randomness. I<Journal of Parapsychology>, I<38>, 393-407.

Delucchi, K. L. (1993). The use and misuse of chi-square: Lewis and Burke revisited. I<Psychological Bulletin>, I<94>, 166-176.

Gatlin, L. L. (1979). A new measure of bias in finite sequences with applications to ESP data. I<Journal of the American Society for Psychical Research>, I<73>, 29-43. (Used for one of the reference tests in the CPAN distribution.)

Good, I. J. (1953). The serial test for sampling numbers and other tests for randomness. I<Proceedings of the Cambridge Philosophical Society>, I<49>, 276-284.

Good, I. J. (1957). On the serial test for random sequences. I<Annals of Mathematical Statistics>, I<28>, 262-264.

Good, I. J., & Gover, T. N. (1967). The generalized serial test and the binary expansion of [square-root]2. I<Journal of the Royal Statistical Society A>, I<130>, 102-107.

Kendall, M. G., & Babington Smith, B. (1938). Randomness and random sampling numbers. I<Journal of the Royal Statistical Society>, I<101>, 147-166.

Knuth, D. E. (1998). I<The art of computer programming> (3rd ed., Vol. 2 Seminumerical algorithms). Reading, MA, US: Addison-Wesley.

Rukhin, A., Soto, J., Nechvatal, J., Smid, M., Barker, E., Leigh, S., et al. (2001). A statistical test suite for random and pseudorandom number generators for cryptographic applications. Retrieved September 4 2010, from L<http://csrc.nist.gov/groups/ST/toolkit/rng/documents/SP800-22b.pdf>.

=head1 SEE ALSO

L<Statistics::Sequences|Statistics::Sequences> sub-modules for other tests of sequences, and for sharing data between these tests.

=head1 TO DO/BUGS

Handle non-overlapping I<v>-nomes.

=head1 AUTHOR/LICENSE

=over 4

=item Copyright (c) 2006-2013 Roderick Garton

rgarton AT cpan DOT org

This program is free software. It may be used, redistributed and/or modified under the same terms as Perl-5.6.1 (or later) (see L<http://www.perl.com/perl/misc/Artistic.html>).

=back

=head1 DISCLAIMER

To the maximum extent permitted by applicable law, the author of this module disclaims all warranties, either express or implied, including but not limited to implied warranties of merchantability and fitness for a particular purpose, with regard to the software and the accompanying documentation.

=head1 END

This ends documentation of the Perl implementation of the I<psi>-square statistic, Kendall-Babington Smith test, and Good's Generalized Serial Test, for randomness in a sequence.

=cut
