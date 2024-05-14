#!/usr/bin/perl
use strict;
use FindBin;

open IN,"<$FindBin::RealBin/../data/dfToName.tsv" or die;
my %dfToName = ();
while (<IN>) {
  chomp;
  my ($df, $name) = split(/\t/);
  $dfToName{$df} = $name;
}
close IN;


open IN,"<sequences.fasta" or die;
my @seqs = ();
my $seq = "";
my $id = "";
while (<IN>) {
  if (/^>(\S+)/) {
    my $tmp = $1;
    if ($seq) {
      push @seqs, [$id, $seq];
    }
    $id = $tmp;
    $seq = "";
    next;
  } 
  s/[\n\r\s]//g;
  $seq .= $_;
}
if ( $seq ) {
  push @seqs, [$id, $seq];
}
close IN;

open OUT,">tmpSeqs.fa" or die;
while ( <> ) {
  # Key: 1126, Neighbor: Neighbor { min_distance: 0.18799205, label_count: 2 }
  if ( /^Key: (\d+), Neighbor: Neighbor { min_distance: (\d+\.\d+), label_count: (\d+) }/ ) {
    my $subj_idx = $1;
    my $distance = $2;
    my $label_count = $3;
    print OUT ">$dfToName{$seqs[$subj_idx]->[0]}  $seqs[$subj_idx]->[0]\n$seqs[$subj_idx]->[1]\n";
  }
}
close OUT;

system("/home/rhubley/projects/RepeatModeler/util/align.pl tmpSeqs.fa query.fasta");
 

exit;


