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

#my $maxFamilies = 5000;
my $maxFamilies = 105000;
my $famIdx = 0;
open OUT,">tmpSeqs.fa" or die;
while ( <> ) {
  # Key: 1126, Neighbor: Neighbor { min_distance: 0.18799205, label_count: 2 }
  if ( /^Key: (\d+), Neighbor: Neighbor { min_distance: (\d+\.\d+), label_count: (\d+) }/ ) {
    my $subj_idx = $1;
    my $distance = $2;
    my $label_count = $3;
    my $label = "";
    if ( $seqs[$subj_idx]->[0] =~ /(\S+)\#(\S+)/ ) {
      my $tID = $1;
      my $tCL = $2;
      if ( $tID =~ /(D[FR])\d\d(\d\d\d\d\d\d\d)/ ) {
        $tID = $1 . $2;
      }
      if ( exists $dfToName{$tID} ) {
        $label = $dfToName{$tID} . "#" . $tCL;
      }else {
        $label = $dfToName{$seqs[$subj_idx]->[0]};
      }
    }else {
      $label = $dfToName{$seqs[$subj_idx]->[0]};
    }
    if ( $label eq "" ) {
      $label = $seqs[$subj_idx]->[0];
    }
    if ( $label_count > 4 ) {
    print OUT ">$label  $seqs[$subj_idx]->[0]\n$seqs[$subj_idx]->[1]\n";
    $famIdx++;
    }

    last if ( $famIdx == $maxFamilies );
    last if ( $distance > 0.15 );

  }
}
close OUT;

system("/home/rhubley/projects/RepeatModeler/util/align.pl -force tmpSeqs.fa -gap_init 25 -extension 5 -minmatch 7 query.fasta");
 

exit;


