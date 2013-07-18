package Statistics::Sequences::Vnomes;

use 5.008008;
use strict;
use warnings;
use Carp 'croak';
use vars qw($VERSION @ISA);
use Algorithm::Combinatorics qw(variations_with_repetition);
use Math::Cephes 0.43;
use Statistics::Sequences 0.11;
use Statistics::Lite qw(sum);
@ISA = qw(Statistics::Sequences);

$VERSION = '0.11';

=pod

=head1 NAME

Statistics::Sequences::Vnomes - The Serial Test (psi-square) and Generalized Serial Test (delta psi-square) for equiprobability of v-nomes (or I<v>-plets/bits) (Good's and Kendall-Babington Smith's tests)

=head1 SYNOPSIS

 use Statistics::Sequences::Vnomes 0.11;
 my $vnomes = Statistics::Sequences::Vnomes->new();
 $vnomes->load(qw/1 0 0 0 1 1 0 1 1 0 0 1 0 0 1 1 1 1 0 1/); # data treated categorically - could also be, e.g., 'a' and 'b'
 my $r_freq_href = $vnomes->observed(length => 3); # returns hashref of frequency distribution for trinomes in the sequence
 my @freq = $vnomes->observed(length => 3); # returns only the observed trinome frequencies (not keyed by the trinomes themselves)
 $val = $vnomes->expected(length => 3); # mean chance expectation for the frequencies (2.5)
 $val = $vnomes->psisq(length => 3); # Good's "second backward differences" psi-square (3.4); option 'delta' gives alternative estimates
 $val = $vnomes->p_value(length => 3, tails => 1); # 1-tailed p-value for psi-square (0.0913)
 $val = $vnomes->z_value(length => 3, tails => 1, ccorr => 1); # inverse-phi of the 1-tailed p-value (1.333)
 my $href = $vnomes->stats_hash(length => 3, values => {psisq => 1, p_value => 1}, tails => 1); # include any stat-method (& their options)
 $vnomes->dump(length => 3,
  values => {observed => 1, expected => 1, psisq => 1, p_value => 1}, # what stats to show (or not if => 0)
  format => 'table', flag => 1, precision_s => 3, precision_p => 7, verbose => 1, tails => 1);
  # prints:
   # Vnomes (3) statistics
   #.-----------+-----------+-----------+-----------.
   #| observed  | expected  | p_value   | psisq     |
   #|           |           |           |           |
   #+-----------+-----------+-----------+-----------+
   #| '011' = 4 | 2.500     | 0.0913418 | 3.400     |
   #| '010' = 1 |           |           |           |
   #| '111' = 2 |           |           |           |
   #| '000' = 1 |           |           |           |
   #| '101' = 2 |           |           |           |
   #| '001' = 3 |           |           |           |
   #| '100' = 3 |           |           |           |
   #| '110' = 4 |           |           |           |
   #'-----------+-----------+-----------+-----------'
   #psisq is Good's delta^2-psi^2, calculated with second backward differences, and has 2 degrees of freedom.
   #p-value is 1-tailed.
 # these and other methods inherited from Statistics::Sequences and its parent, Statistics::Data:
 $vnomes->dump_data(delim => ','); # comma-separated single line of the loaded data
 $vnomes->save(path => 'seq.csv'); # serialized for retrieval by open() method

=head1 DESCRIPTION

Implements tests of the independence of successive elements of a sequence of data: serial tests for v-nomes (a.k.a I<v>-plets or, for binary data, I<v>-bits) - singlets/monobits, dinomes/doublets, trinomes/triplets, etc.. Test are of variations in sub-sequence length I<v> that are equally likely for the sampled sequence. For example, a sequence sampled from a "heads'n'tails" (H and T) distribution can be tested for its equal representation of the trinomes HTH, HTT, TTT, THT, and so on. Counting up these v-nomes at all points in the sequence, permitting overlaps, yields a statistic - psi-square - that is approximately distributed as I<chi>-square; the Kendall-Babington Smith statistic.

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

 $href = $vnomes->observed(length => int, circularize => 1|0); # returns keyed distribution; assumes data have already been loaded
 @ari = $vnomes->observed(data => [1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1], length => int, circularize => 1|0); # returns frequencies only

Returns the frequency distribution for the I<observed> number of v-nomes of the given B<length> in a sequence. These are counted up as overlapping. Called in array context, returns an array of the counts for each possible pattern; otherwise, where these are keyed by the pattern in a hash reference. For descriptives of the observed frequencies, try calling L<Statistics::Lite|Statistics::Lite> methods with this array.

A value for B<length> greater than zero is required, and must be no more than the sample-size. So for the sequence

 1 0 0 0 1 0 0 1 0 1 1 0 

there are 12 v-nomes of length 1 (mononomes, the number of elements) from the two (0 and 1) that are possible; 11 dinomes from the four (10, 00, 00, 01) that are possible; and 10 trinomes from eight (100, 000, 001, 010, etc.) that are possible.

By default, the sequence is counted up for v-nomes as a cyclic sequence, treating the first element of the sequence as following the last one. So the count loops to the beginning of the sequence until all elements from the end are included, and instead of ending the count for trinomes in this sequence at '110', the count includes '101' and '010', increasing the observed sum of trinomes to 12. Set L<circularize|circularize> => 0 if the count is not to be made cyclically.

The data to test can already have been L<load|Statistics::Sequences/load, add, unload>ed, or you send it directly keyed as B<data>.

In the code for this and/or other methods, following the work of Good (1953, 1957; Good & Gover, 1967), I<v> is used to define variables referring to the length of the subsequence (I<mono>nome, I<di>nome, I<tri>nome, etc.), I<t> defines the number of states (events, letters, etc.) in the sequence, and I<r> identifies each possible subsequence of length I<v>. 

=cut

sub observed {
    my $self = shift;
    my $args = ref $_[0] ? shift : {@_};
    my $data = ref $args->{'data'} ? $args->{'data'} : $self->read($args);
    my $num = _get_n($self, $args);
    my $v = $args->{'length'} || croak __PACKAGE__, '::test Must have a v-nome length greater than zero';
    # Sanitize length in proportion to size:
    croak __PACKAGE__, "::test Sequence length v ($v) for testing should be no more than the sample-size ($num) - 2" if $v >= $num - 1;
    #my $lim =  Math::Cephes::floor(Math::Cephes::log2($num)) - 2; # re NIST recommmendation that $v is too small for $num and $t if $v >= $lim
    my $circularize = defined $args->{'circularize'} ? $args->{'circularize'} : 1;
    my $states_aref = _get_stateslist($data, $args->{'states'});# Get list of unique states as given or from sequence itself
    my $v_aref = _get_v_ari($v);
    my ($v_i, @data_i, %psisq_v) = ();
    foreach $v_i(@$v_aref) {# print "v_i = $v_i\n"; # loop through tests of v, v-1 and v-2 lengths ...
        @data_i = @$data;# extend sequence by appending the first $v_i - 1 elements to its end
        push(@data_i, @data_i[0 .. $v_i - 2]) if $circularize; # e.g., if v = 3, append elements [0] & [1]; if v = 2, append only 1st value
        $psisq_v{$v_i} = _frequencies(\@data_i, $v_i, $states_aref);# count up frequencies of each sequence of length $v_i
        last if $v_i == $v;
   }
   return wantarray ? values %{$psisq_v{$v}} : $psisq_v{$v}; 
}
*vnomes_observed = \&observed;
*vno = \&observed;

=head2 expected, vnomes_expected, vne

 $count = $vnomes->expected(length => int, circularize => 1|0); # assumes data have already been loaded
 $count = $vnomes->expected(data => [1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1], length => int, circularize => 1|0, states => [0, 1]);

Returns the I<expected> number of observations for each v-nome of the given B<length> in a sequence; i.e., the mean chance expectation, assuming that each event is generated by a random uniform process. Options are as for L<observed|Statistics::Sequences::Vnomes/observed>. This expected frequency is given by:

=for html <p>&nbsp;&nbsp;<i>E[V]</i> = <i>N</i><i>t</i><sup>&ndash;<i>v</i></sup>

where I<t> is the number of possible states (alternatives; as read from the data or as explicitly given as an array in B<states>), I<v> is the v-nome length, and I<N> is the length of the sequence, less I<v> + 1 if the count is not to be circularized (Good, 1953, Eq. 12).

Another way to think of this is as the number of mononome observations in the sequence divided by the number of possible permutations of its states for the given length. So, for a sequence made up of 0s and 1s, there are four possible variations of length 2 (00, 10, 01 and 11), so that the expected frequency for each of these variations in a sequence of 20 values is  20 / 4, i.e., 5.

=cut

sub expected {
    my $self = shift;
    my $args = ref $_[0] ? shift : {@_};
    my $num = _get_n($self, $args);
    my $v = $args->{'length'} || croak __PACKAGE__, '::test Must have a v-nome length greater than zero';
    # Sanitize length in proportion to size:
    croak __PACKAGE__, "::test Sequence length v ($v) for testing should be no more than the sample-size ($num) - 2" if $v >= $num - 1; # or permit $v to extend to limit of circularized data?
    #my $lim =  Math::Cephes::floor(Math::Cephes::log2($num)) - 2; # re NIST recommmendation that $v is too small for $num and $t if $v >= $lim
    my $circularize = defined $args->{'circularize'} ? $args->{'circularize'} : 1; # circularize by default
    # Compute expected freq of variations of $t of length $v:
    $num = _count_lim($num, $v) if ! $circularize;
    my $p = $self->prob_r(length => $v);
    return $num * $p;
}
*vnomes_expected = \&expected;
*vne = \&expected;

=head2 variance

 $var = $vnomes->variance(length => int, circularize => 1|0); # assumes data have already been loaded
 $var = $vnomes->variance(data => [1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1], length => int, circularize => 1|0, states => [0, 1]);

Returns the variance in the expected frequency of each v-nome, as per Good (1953, Eqs. 11 and 14).

=cut

sub variance {
    #croak 'variance() method is not defined for ' . __PACKAGE__;
    my $self = shift;
    my $args = ref $_[0] ? shift : {@_};
    my $num = _get_n($self, $args);
    my $v = $args->{'length'} || croak __PACKAGE__, '::test Must have a v-nome length greater than zero';
    # Sanitize length in proportion to size:
    croak __PACKAGE__, "::test Sequence length v ($v) for testing should be no more than the sample-size ($num) - 2" if $v >= $num - 1; # or permit $v to extend to limit of circularized data?
    #my $lim =  Math::Cephes::floor(Math::Cephes::log2($num)) - 2; # re NIST recommmendation that $v is too small for $num and $t if $v >= $lim
    # Compute expected freq of variations of $t of length $v:
    my $pr = $self->prob_r(length => $v);
    my $var = $pr;
    my ($sum_1, $sum_2, $m, $pm) = (0, 0);
    foreach $m(1 .. $v - 1) {
        $pm = $self->prob_r(length => $m);
        $sum_1 += $pm;
        $sum_2 += ($m + $v - 1) * $pm;
    }
    $var += 2 * $pr * $sum_1;
    $var -= (2 * $v - 1) * $pr**2;
    my $circularize = defined $args->{'circularize'} ? $args->{'circularize'} : 1;
    $var += $num**-1 * $pr * ( ( ($v - 1) * (3 * $v - 1) * $pr ) - ($v - 1) - 2 * $sum_2 ) if !$circularize;
    return $var;
}

sub obsdev {
   croak 'obsdev() method is not defined for ' . __PACKAGE__;
}
*observed_deviation = \&obsdev;

=head2 stdev, standard_deviation

 $var = $vnomes->stdev(length => int, circularize => 1|0); # assumes data have already been loaded
 $var = $vnomes->stdev(data => [1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1], length => int, circularize => 1|0, states => [0, 1]);

Returns the standard deviation in the expected frequency of each v-nome, as per Good (1953, Eqs. 11 and 14).

=cut

sub stdev {
   return sqrt(variance(@_));
}
*standard_deviation = \&stdev;

=head2 psisq

 $vnomes->psisq(length => int, delta => '2|1|0', circularize => '1|0', states => [qw/A C G T/]);

Performs Good's Generalized Serial Test (by default), of v-nomes on the given or named distribution, yielding a I<psi-square> statistic. The option B<delta> specifies if backward differencing.

=over 4

=item delta =E<gt> 0

The raw psi-square value for sub-sequences of length I<v> (Kendall-Babington Smith statistic), i.e., without backward differencing.

=for html <p>&nbsp;&nbsp;&nbsp; &Psi;<sup>2</sup><sub>v</sub> = &sum;{<i>r</i>} ( <i>n<sub>r</sub></i>&ndash; (<i>N</i> &ndash; <i>v</i> + 1) / <i>t<sup>v</sup></i> )<sup>2</sup> ) / ( (<i>N</i> &ndash; <i>v</i> + 1) / <i>t<sup>v</sup></i> )</p>

for uncircularized sequences (Good, 1953, Eq. 1), and

=for html <p>&nbsp;&nbsp;&nbsp; <span style="text-decoration:overline;">&Psi;</span><sup>2</sup><sub>v</sub> = &sum;{<i>r</i>} ( <i><span style="text-decoration:overline;">n</style><sub>r</sub></i>&ndash; <i>Nt<sup>&ndash;v</sup></i> ) / <i>Nt<sup>&ndash;v</sup></i></p>

for circularized sequences (Good, 1953, Eq. 2), where an overline indicates circularization, and where I<N> is the length of the sequence, I<v> is the v-nome length, I<t> is the number of unique states that each element of the sequence can take, and I<r> is the number of variations of length I<v> for the given number of unique states.

This statistic is only asymptotically distributed as I<chi>-square if I<length> (I<v>) = 1, and is therefore not used as the default.

=item delta =E<gt> 1

The "I<first> backward differences" of psi-square, which is the difference between the psi-square values for sub-sequences of length I<v> and length I<v> - 1.

=for html <p>&nbsp;&nbsp;&nbsp; &Delta;&Psi;<sup>2</sup><sub>v</sub> = &Psi;<sup>2</sup><sub>v</sub> &ndash; &Psi;<sup>2</sup><sub>v&ndash;1</sub> (<i>v</i> &ge; 1)</p>

While it is I<chi>-square distributed, counts of first-differences are not statistically independent (Good, 1953; Good & Gover, 1967), and "the sequence of second differences forms a much better set of statistics for testing the hypothesis of flat-randomness" (Good & Gover, 1967, p. 104).

=item delta =E<gt> 2 (or undef)

This method returns by default (without specifying a value for B<delta>), or if B<delta> is not 1 or 0, the "second backward differences" psi-square (I<delta^2-psi^2>). This incorporates psi-square values for backwardly adjacent values of I<length> (I<v>), i.e., for sub-sequences of length I<v>, I<v> - 1, and I<v> - 2.

=for html <p>&nbsp;&nbsp;&nbsp; &Delta;<sup>2</sup>&Psi;<sup>2</sup><sub>v</sub> = &Psi;<sup>2</sup><sub>v</sub> &ndash; 2&Psi;<sup>2</sup><sub>v&ndash;1</sub> &ndash; &Psi;<sup>2</sup><sub>v&ndash;2</sub> (<i>v</i> &ge; 2)</p>

It is not only asymptotically I<chi>-square distributed, but uses statistically independent counts of all the possible variations of sequences of the given B<length> (Good, 1953).

=back

See also Good's algorithm in individual papers describing application of the Serial Test (e.g., Davis & Akers, 1974).

=cut

sub psisq {
    my $self = shift;
    my $args = ref $_[0] ? shift : {@_};
    my $data = ref $args->{'data'} ? $args->{'data'} : $self->read($args);
    my $num = _get_n($self, $args);
    my $v = $args->{'length'} || croak __PACKAGE__, '::test Must have a v-nome length greater than zero';
    # Sanitize length in proportion to size:
    croak __PACKAGE__, "::test Sequence 'length' ($v) for testing should be no more than the sample-size ($num) - 2" if $v >= $num - 1;
    #my $lim =  Math::Cephes::floor(Math::Cephes::log2($num)) - 2; # re NIST recommmendation that $v is too small for $num and $t if $v >= $lim
    my $circularize = defined $args->{'circularize'} ? $args->{'circularize'} : 1;
    my $delta = defined $args->{'delta'} ? $args->{'delta'} : 2;
    my $t = _get_t($self, $args);
    my $states_aref = _get_stateslist($data, $args->{'states'});# Get list of unique states as given or from sequence itself
    my $v_aref = _get_v_ari($v);
    # While looping through tests of v, v-1 and v-2 (= $v_i) sequence lengths ...
    my ($v_i, $r_freq, @data_i, %psisq_v) = ();
    foreach $v_i(@$v_aref) {
        # counting vnome frequencies: handle circularization by firstly appending elements from the start to the sequence's end:
        if ($circularize) {
            @data_i = @$data; # copy the original sequence
            push(@data_i, @data_i[0 .. $v_i - 2]); # e.g., if v = 3, append elements [0] & [1]; if v = 2, append only 1st value
            $r_freq = _frequencies(\@data_i, $v_i, $states_aref); # frequencies of each sequence of length $v_i
            $psisq_v{$v_i} = _psisq_uncirc(scalar @data_i, $t, $v_i, $r_freq);
        }
        else {
            $r_freq = _frequencies($data, $v_i, $states_aref);
            $psisq_v{$v_i} = _psisq_uncirc(scalar @$data, $t, $v_i, $r_freq);
        }
    }
    my ($psisq, $df) = _sel_psisq_and_df($v, $t, \%psisq_v, $delta);
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
   $args->{'tails'} ||= 1;
   my $pval = $self->p_value($args);
   require Statistics::Zed;
   my $zed = Statistics::Zed->new();
   my $zval = $zed->p2z(value => $pval, %$args);
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
   $args->{'tails'} ||= 1;
   $p_value /= 2 if $args->{'tails'} == 1;
   return $p_value;
}
*test = \&p_value;
*vnomes_test = \&p_value;
*vnt = \&p_value;

=head2 dump

 $vnomes->dump(length => 3, values => {psisq => 1, p_value => 1}, format => 'table|labline|csv', flag => 1, precision_s => 3, precision_p => 7, verbose => 1, tails => 1);

Print Vnome-test results to STDOUT. See L<dump|Statistics::Sequences/dump> in the Statistics::Sequences manpage for details. If B<verbose> => 1, then you get (1) the actual test-statistic depending on the value of B<delta> tested (I<delta^2-psi^2> for the second difference measure (default), I<delta-psi^2> for the first difference measure, and I<psi2> for the raw measure), followed by degrees-of-freedom in parentheses; and (2) a warning, if relevant, that your B<length> value might be too large with respect to the sample size (see NIST reference, above, in discussing B<length>). If B<text> => 1, you just get the average observed and expected frequencies for each v-nome, the I<Z>-value, and its associated I<p>-value.

=cut

sub dump {
    my $self = shift;
    my $args = ref $_[0] ? $_[0] : {@_};
    $args->{'stat'} = 'vnomes';
    $args->{'stat'} .= " ($args->{'length'})" if defined $args->{'length'};
    $self->SUPER::dump($args);
    if ($args->{'verbose'} && $args->{'format'} eq 'table') {
        if ($args->{'values'}->{'psisq'}) {
            my $delta = defined $args->{'delta'} ? $args->{'delta'} : 2;
            my ($psisq, $df) = $self->psisq($args);
            my $df_str = 'degree';
            $df_str .= 's' if $df != 1; 
            if ($delta == 0) { # Raw psisq:
                print "psisq is the Kendall-Babington Smith statistic, calculated without backward differencing, and has $df $df_str of freedom.\n";
            }
            elsif ($delta == 1) { # First backward difference:
                print "psisq is Good's delta-psi^2, calculated with first backward differences, and has $df $df_str of freedom.\n";
            }
            else { 
                print "psisq is Good's delta^2-psi^2, calculated with second backward differences, and has $df $df_str of freedom.\n";
            }
        }
        if ($args->{'values'}->{'p_value'}) {
            print "p-value is $args->{'tails'}-tailed.\n";
        }
    }
    return $self;
}

=head2 nnomes

 $r = $vnomes->nnomes(length => int, states => \@ari); # minimal option required; the "v" value itself, with number of states
 $r = $vnomes->nnomes(length => int); # minimal option required, assuming data are loaded so can read its number of states
 $r = $vnomes->nnomes(length => int, data=> \@data); # minimal option required; the "v" value itself, with these data

Returns the number of possible subsequences of the given length (I<v>) for the given number of states (I<t>). This is quantity denoted as I<r> in Good's (1953, 1957) papers; i.e.,

=for html <p>&nbsp;&nbsp;<i>r[V]</i> = <i>t</i><sup><i>v</i></sup>

The routine needs to know two things: the "v" value itself, i.e., the length of the possible subsequences to test (I<mono>nomes, I<di>nomes, I<tri>nomes, etc.), and the number of states (events, letters, etc.) that the process generating the data could take (from 1 to whatever). The former is always required to be specified by the named argument B<length>. The latter (the number of states) can be directly taken from the named argument B<nstates>, indirectly from the size of the array referenced to the named argument B<states>, from whatever states exist in the array referenced to the named argument B<data>, or from whatever states exist in data already L<load|Statistics::Sequences::Vnomes/load>ed. 

=cut

sub nnomes {
    my $self = shift;
    my $args = ref $_[0] ? $_[0] : {@_};
    my $t = _get_t($self, $args);
    my $v = $args->{'length'} || croak __PACKAGE__, '::nnomes Must define argument \'length\' with a value greater than zero';
    return $t**$v;
}

=head2 prob_r

 $Pr = $vnomes->prob_r(length => $v); # length is 1 (mononomes) by default

Returns the probability of the occurrence of any of the individual elements ("digits") in the sequence (I<v> = 1), or of the given B<length>, assuming they are equally likely and independent.

=for html <p>&nbsp;&nbsp;<i>P<sub>r</sub></i> = <i>t</i><sup><i>&ndash;v</i></sup>

=cut

sub prob_r {
    my $self = shift;
    my $args = ref $_[0] ? $_[0] : {@_};
    my $t = _get_t($self, $args);
    my $v = $args->{'length'} || 1;
    return $t**(-1 * $v);
}

=head1 OPTIONS

Options common to the above stats methods.

=head2 length

This is currently a I<required> "option", giving the length of the v-nome of interest, i.e., the value of I<v> - an integer greater than or equal to 1, and smaller than than the sample-size.

What is a meaningful maximal value of B<length>? As a I<chi>-square test, it is conventionally required that the expected frequency is at least 5 for each v-nome (Knuth, 1988). This can be judged to be too conservative (Delucchi, 1993). The NIST documentation on the serial test (Rukhin et al., 2010) recommends that length should be less than the floored value of log2 of the sample-size, minus 2. No tests are here made of these recommendations.

=head2 circularize

By default, L<observed|Statistics::Sequences::Vnomes/observed, vnomes_observed, vno> and L<expected|Statistics::Sequences::Vnomes/expected, vnomes_expected, vne> counts, and the value of L<psisq|Statistics::Sequences::Vnomes/psisq>, are made by treating the sequence as a cyclic one, where the first element of the sequence follows the last one. This affects (and simplifies) the calculation of the expected frequency of each v-nome, and so the value of each psi-square. Also, circularizing ensures that the expected frequencies are accurate; otherwise, they might only be approximate. As Good and Gover (1967) state, "It is convenient to circularize in order to get exact checks of the arithmetic and also in order to simplify some of the theoretical formulae" (p. 103). These methods, however, can also treat the sequence non-cyclically by calling them with B<circularize> => 0.

=head2 states

Optionally send a referenced array listing the unique states (or 'events', 'letters') in the population from which the sequence was sampled, e.g., B<states> => [qw/A C G T/]. This is useful if the sequence itself might not include all the possible states. If this is not specified, the states are identified from the sequence itself. If giving a list of states, a check in each test is made to ensure that the sequence contains I<only> those elements in the list.

=cut

# PRIVATMETHODEN
sub _count_lim { # N - v + 1; if N = 5 and v = 3, starting the count-up for v-nomes must stop from the first 3 elements of the data
    return $_[0] - $_[1] + 1;
}

sub _get_n { # the size of the sequence
    my ($self, $args) =  @_;
    my $data = ref $args->{'data'} ? $args->{'data'} : $self->read($args);
    return scalar @$data;
}

sub _get_t { # the number of unique states; e.g., 2 for a binary sequences; 10 for a typical random digit sequence (0 .. 9)
    my ($self, $args) = @_;
    my $t;
    if ($args->{'nstates'}) {
        $t = $args->{'nstates'};
    }
    elsif (ref $args->{'states'} and ref $args->{'states'} eq 'ARRAY') {
        $t = scalar(@{$args->{'states'}});
    }
    else {
        my $data = ref $args->{'data'} ? $args->{'data'} : $self->read($args);
        my %hash = map { $_, 1 } @{$data};
        my $states = [keys %hash];
        $t = scalar(@{$states});
    }
    return $t;
}

sub _get_stateslist {
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
                next DATA if $h eq $g;
            }
            croak __PACKAGE__, "::test The element $g in the data is not represented in the given states"; 
        }
    }
    #croak __PACKAGE__, '::test At least two different values must be in the sequence to test its sub-sequences' if $t <= 1;
    return $states;
}

sub _frequencies {
    my ($data_i, $v_i, $states_aref) = @_;
    # Get a list of all possible combinations of states at the current length ($v_i):
    my @variations = variations_with_repetition($states_aref, $v_i);
    # Count up the frequency of each variation in the data:
    my $num = scalar(@{$data_i});
    my ($i, $probe_str, $test_str, %r_freq) = ();
    foreach (@variations) {
        $probe_str = join'', @{$_};
        $r_freq{$probe_str} = 0;
        for ($i = 0; $i < _count_lim($num, $v_i); $i++) {
            $test_str = join'', @{$data_i}[$i .. ($i + $v_i - 1)];
            $r_freq{$probe_str}++ if $probe_str eq $test_str;
        }
    }#    print "FREQ:\n"; while (my($key, $val) = each %freq) { print "\t$key = $val\n"; } print "\n";
    return \%r_freq;
}
 
sub _sel_psisq_and_df {# get psisq and its df
    my ($v, $t, $psisq_v, $delta) = @_;
    my ($psisq, $df) = ();
    
    if ($v == 1) { # psisq is asymptotically distributed chisq, can use psisq for chisq distribution:
        $psisq = $psisq_v->{1};
        $df = $t - 1;
    }
    else {
        if ($delta == 0) { # Raw psisq:
            $psisq = $psisq_v->{$v};
            $df = $t**$v - 1; # Good (1957, Eq. 6) - if circularized or not
        }
        elsif ($delta == 1) { # First backward difference:
            $psisq = $psisq_v->{$v} - ($v - 1 <= 0 ? 0 : $psisq_v->{$v - 1});
            $df = $t**$v - $t**($v - 1);
        }
        else { # $delta == 2 # Second backward difference (default):
            $psisq = $psisq_v->{$v} - ( 2 * ($v - 1 <= 0 ? 0 : $psisq_v->{$v - 1}) ) + ($v - 2 <= 0 ? 0 : $psisq_v->{$v - 2});
            $df = ( $t**($v - 2) ) * ( $t - 1)**2;
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

sub _psisq_uncirc {# Compute psi^2 for uncircularized sequence from freq of each variation of length $v for given states(Good, 1953, Eq. 1):
    my ($n, $t, $v, $r_freq) = @_;
    my $k = ($n - $v + 1) / $t**$v;
    my $sum = sum( map{ ($_ - $k)**2 } values %{$r_freq});#foreach (keys %{$r_freq}) {$psisq += ($r_freq->{$_} - $k)**2;}
    return $sum / $k; 
}

sub _psisq_circ {
    my ($n, $t, $v, $r_freq) = @_;# Compute psi^2 for circularized sequence from freq of each variation of length $v for given states(Good, 1953, Eq. 2):
    my $k = $n * $t**(-1 * $v);
    my $sum = sum( map{ ($_ - $k) } values %{$r_freq});#foreach (keys %{$r_freq}) { $sum += ($r_freq->{$_} - $k)**2;}
    return $sum / $k;
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

That is, the observed frequency of each possible trio of seating arrangements (the trinomes OOO, OOE, OEE, EEE, etc.) differed significantly from that expected. Look up the observed frequencies for each possible trinome to see if this is because there are more empty or occupied neighbouring seats ("phobia" or "philia"):

 $vnomes->dump(length => 3, values => {observed => 1}, format => 'labline');

This prints:

 observed = ('OEE' = 4,'EEO' = 4,'EEE' = 2,'OEO' = 1,'EOE' = 5,'OOO' = 0,'OOE' = 0,'EOO' = 0)

As the chance-expected frequency is 2.5 (from the L<expected|Statistics::Sequences::Vnomes/expected, vnomes_expected, vne> method), there are clearly more than expected trinomes involving empty seats than occupied seats - suggesting a non-random factor like social phobia (or body odour?) is at work in sequencing people's seating here. Noting that the sequencing isn't significant for dinomes (with B<length> => 2) might also tell us something about what's going on. What happens for v-nomes of 4 or more in length? Maybe the L<runs|Statistics::Sequences::Runs> or L<pot|Statistics::Sequences::Pot> test might be a better summary of what's going on.

=head1 REFERENCES

Davis, J. W., & Akers, C. (1974). Randomization and tests for randomness. I<Journal of Parapsychology>, I<38>, 393-407.

Delucchi, K. L. (1993). The use and misuse of chi-square: Lewis and Burke revisited. I<Psychological Bulletin>, I<94>, 166-176.

Gatlin, L. L. (1979). A new measure of bias in finite sequences with applications to ESP data. I<Journal of the American Society for Psychical Research>, I<73>, 29-43. (Used for one of the reference tests in the CPAN distribution.)

Good, I. J. (1953). The serial test for sampling numbers and other tests for randomness. I<Mathematical Proceedings of the Cambridge Philosophical Society>, I<49>, 276-284.

Good, I. J. (1957). On the serial test for random sequences. I<Annals of Mathematical Statistics>, I<28>, 262-264.

Good, I. J., & Gover, T. N. (1967). The generalized serial test and the binary expansion of [square-root]2. I<Journal of the Royal Statistical Society A>, I<130>, 102-107.

Kendall, M. G., & Babington Smith, B. (1938). Randomness and random sampling numbers. I<Journal of the Royal Statistical Society>, I<101>, 147-166.

Knuth, D. E. (1998). I<The art of computer programming> (3rd ed., Vol. 2 Seminumerical algorithms). Reading, MA, US: Addison-Wesley.

Rukhin, A., Soto, J., Nechvatal, J., Smid, M., Barker, E., Leigh, S., et al. (2010). A statistical test suite for random and pseudorandom number generators for cryptographic applications. Retrieved September 4 2010, from L<http://csrc.nist.gov/groups/ST/toolkit/rng/documents/SP800-22b.pdf>, and July 17, 2013, from L<http://csrc.nist.gov/publications/nistpubs/800-22-rev1a/SP800-22rev1a.pdf|SP800-22rev1a.pdf> (revised).

=head1 SEE ALSO

L<Statistics::Sequences|Statistics::Sequences> sub-modules for other tests of sequences, and for sharing data between these tests.

=head1 TO DO/BUGS

Handle non-overlapping v-nomes.

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
