#!/usr/bin/perl

use strict;
use warnings;
use JFR::Fasta;
use JFR::GFF3;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

# this version fixes transcript bug (0.02 did not expect transcripts )
our $VERSION = 0.05;

MAIN: {
    my $rh_o = process_options();
    my $rh_map = get_map($rh_o->{'gff'});
    make_new_gff($rh_o->{'gff'},$rh_map,$rh_o->{'out'});
    replace_deflines($rh_map,$rh_o->{'out'},$rh_o,'aa');
    replace_deflines($rh_map,$rh_o->{'out'},$rh_o,'codingseq');
    replace_deflines($rh_map,$rh_o->{'out'},$rh_o,'mrna') if (-f "$rh_o->{'pre'}.mrna");
    
}

sub make_new_gff {
    my $file = shift;
    my $rh_m = shift;
    my $out  = shift;
    open IN, $file or die "cannot open $file:$!";
    open OUT, ">$out.gff" or die "cannot open $out.gff:$!";
    while (my $line = <IN>) {
        if ($line =~ m/((g\d+)\.t\d+)/) {
            my $id = $1;
            my $gid = $2;
            die "unexpected" unless ($rh_m->{$id});
            $rh_m->{$id} =~ m/(.*)\.t\d+$/ or die "unexpected";
            my $new_gid = $1;
            $line =~ s/$id/$rh_m->{$id}/;
            $line =~ s/Parent=$gid/Parent=$new_gid/;
            die "not expect multiple ids per line" if ($line =~ m/(g\d+t\d+)/);
        } elsif ($line =~ m/start gene (g\d+)/ || $line =~ m/end gene (g\d+)/) {
            my $id = $1;
            $rh_m->{"${id}\.t1"} =~ m/(.*)\.t\d+$/ or die "unexpected: $id";
            my $new_gid = $1;
            $line =~ s/$id/$new_gid/;
        }
        print OUT $line;
    }
}

sub replace_deflines {
    my $rh_map = shift;
    my $out = shift;
    my $rh_o = shift;
    my $suf = shift;
    my $fa = "$rh_o->{'pre'}.$suf";

    open OUT, ">$out.$suf" or die "cannot open $out.$suf:$!";
    my $fp = JFR::Fasta->new($fa);
    while (my $rec = $fp->get_record()) {
        my $id = JFR::Fasta->get_def_w_o_gt($rec->{'def'});
        $id =~ s/\.cds\d+$//;
        die "unexpected: $id" unless ($rh_map->{$id});
        print OUT ">$rh_map->{$id}\n$rec->{'seq'}\n";
    }
}

sub get_map {
    my $gff = shift;
    my %map = ();
    my %ids = ();
    open IN, $gff or die "cannot open $gff:$!";
    while (my $line = <IN>) {
        next unless ($line =~ m/AUGUSTUS\ttranscript/);
        chomp $line;
        my @f = split /\t/, $line;
        $f[8] =~ m/(g\d+)\.t(\d+)/ or die "unexpected: $line";
        my $g_id = $1;
        my $t_num = $2;
        $ids{$f[0]}++ if ($t_num == 1);
        my $g_num = $ids{$f[0]};
        my $id = $g_id . '.t' . $t_num;
#        my $id = 'g' . $g . '.t' . $t_num;
        $map{$id} = $f[0] . 'g' . $g_num . '.t' . $t_num;
    }
    return \%map;
}

sub usage {
    print "usage: $0 [--version] [--help] --gff=GFF_FILE --out=OUTPREFIX --pre=INPUTPREFIX\n";
    exit;
}

sub process_options {
    my $rh_o = {};
    my $opt_results = Getopt::Long::GetOptions(
                              "version" => \$rh_o->{'version'},
                                "gff=s" => \$rh_o->{'gff'},
                                "out=s" => \$rh_o->{'out'},
                                "pre=s" => \$rh_o->{'pre'},
                                 "help" => \$rh_o->{'help'});
    die "$VERSION\n" if ($rh_o->{'version'});
    pod2usage({-exitval => 0, -verbose => 2}) if $rh_o->{'help'};
    unless ($rh_o->{'out'} && $rh_o->{'gff'} && $rh_o->{'pre'} &&
            -f "$rh_o->{'pre'}.aa" &&
            -f "$rh_o->{'pre'}.codingseq") {
        if ($rh_o->{'pre'}) {
            warn "missing $rh_o->{'pre'}.aa" unless -f "$rh_o->{'pre'}.aa";
            warn "missing $rh_o->{'pre'}.codingseq" unless -f "$rh_o->{'pre'}.codingseq";
        }
        warn "--pre is required" unless $rh_o->{'pre'};
        warn "--gff is required" unless $rh_o->{'gff'};
        warn "--out is required" unless $rh_o->{'out'};
        usage();
    }
    return $rh_o;
}

