#!perl

# LICENSE: GNU GENERAL PUBLIC LICENSE GPL3.0
# Joseph Ryan <joseph.ryan@whitney.ufl.edu>

use strict;
use warnings;
use Data::Dumper;

my $dir = $ARGV[0] or die "usage: $0 INDIR OUTDIR\n"; # asks user to type indir name. should be OrthoFinder/Sequences dir
my $outdir = $ARGV[1] or die "usage: $0 INDIR OUTDIR\n"; #asks user to type outdir name
die "$outdir exists" if (-d $outdir);
mkdir $outdir or die "cannot mkdir $outdir:$!";

opendir DIR, $dir or die "cannot opendir $dir$!"; #open the indir

my @seq_files = grep { /pruned3$/ } readdir DIR; #put all the files that end in .fa into array @seq_file

foreach my $sf (@seq_files) {
    die "unexpected" unless ($sf =~ m/\.pruned3/);
    open OUT, ">$outdir/$sf" or die "cannot open >$outdir/$sf:$!";
    open IN, "$dir/$sf" or die "cannot open $dir/$sf:$!";
    while (my $line = <IN>) {
        chomp $line;
        #Ceri_amer_fl|38426
        #Leptogorgia_sp|82766
        #Stro_purp|03831
        if ($line =~ m/^>[A-Z][a-z]+_[a-z]+\|[0-9]+/) {
            $line =~ s/^(>[A-Z][a-z]+_[a-z]+)\|[0-9]+/$1/;
            print OUT "$line\n";
        } elsif ($line =~ m/^>[A-Z][a-z]+_[a-z]+_[a-z]+\|[0-9]+/) {
            $line =~ s/^(>[A-Z][a-z]+_[a-z]+_[a-z]+)\|[0-9]+/$1/;
            print OUT "$line\n";
        } else {
            $line =~ s/\*\s*$//;
            print OUT "$line\n";
        }
    }
}

