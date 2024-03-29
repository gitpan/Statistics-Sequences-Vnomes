use 5.006;
use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'Statistics::Sequences::Vnomes',
    license             => 'perl',
    dist_author         => q{Roderick Garton <rgarton@cpan.org>},
    dist_version_from   => 'lib/Statistics/Sequences/Vnomes.pm',
    configure_requires => {
        'Module::Build' => '0.36',
    },
    build_requires => {
        'Test::More' => 0,
        'Algorithm::Combinatorics' => '0.27',
        'Math::Cephes' => '0.43',
        'Statistics::Lite' => '3.2',
        'Statistics::Zed' => '0.072',
        'Statistics::Sequences' => '0.11',
    },
    requires => {
        'perl' => 5.006,
        'Algorithm::Combinatorics' => '0.27',
        'Math::Cephes' => '0.43',
        'Statistics::Lite' => '3.2',
        'Statistics::Zed' => '0.072',
        'Statistics::Sequences' => '0.11',
    },
    add_to_cleanup      => [ 'Statistics-Sequences-Vnomes-0.11-*' ],
    create_makefile_pl => 'traditional',
);

$builder->create_build_script();
