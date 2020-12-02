#!/usr/bin/perl -w
###################################################################################################
#                                            xml2msa.pl
#
# This is a perl script used to extract the multiple sequence alignment (MSA) from the PSI-BLAST
# result. MSA file can be used to generate profile for protein design used in EvoDesign.
####################################################################################################
use strict;
use warnings;

if(@ARGV != 4){
  print "This script is used to convert PSI-BLAST xml format to MSA file\n";
  print "Usage: ./xml2msa.pl -in blast.xml -out msa_psiblast.txt\n";

  # exit with error code -1
  exit(-1);
}

my $xmlfile = $ARGV[1];
my $msafile = $ARGV[3];

open my $fh, "$xmlfile" or die "Cannot open the blast file\n";
local $/ = "<\/Iteration>\n";

open FILE2, ">$msafile";
while(<$fh>){
  my $query = $1 if ($_ =~ /<Iteration_query-def>(\S+)<\/Iteration_query-def>/);
  my $query_len = $1 if ($_ =~ /<Iteration_query-len>(\S+)<\/Iteration_query-len>/);
  while ($_ =~ /<Hit>(.*?)<\/Hit>/sg){
    chomp(my $hit = $1);
    my $id = $1 if ($hit =~ /<Hit_id>(\S+)<\/Hit_id>/);
    my $def = $1 if ($hit =~ /<Hit_def>(.*?)<\/Hit_def>/);
    my $num = 0;
    while ($hit =~ /<Hsp>(.*?)<\/Hsp>/sg){
      $num++;
      my $hsp = $1;
      my $query_from = $1 if ($hsp =~ /<Hsp_query-from>(\S+)<\/Hsp_query-from>/);
      my $query_to = $1 if ($hsp =~ /<Hsp_query-to>(\S+)<\/Hsp_query-to>/);
      my $hit_frame = $1 if ($hsp =~ /<Hsp_query-frame>(\S+)<\/Hsp_query-frame>/);
      my $qseq = $1 if ($hsp =~ /<Hsp_qseq>(\S+)<\/Hsp_qseq>/);
      my $hseq = $1 if ($hsp =~ /<Hsp_hseq>(\S+)<\/Hsp_hseq>/);
      my $align_len = $1 if ($hsp =~ /<Hsp_align-len>(\S+)<\/Hsp_align-len>/);
      #print ">$query;orf$num len=$align_len frame:$hit_frame start:$query_from end:$query_to $id $def\n";
      
      my @qseq_con = split(//, $qseq);
      my @hseq_con = split(//, $hseq);
      
      ##1. check if the hit sequence have >98% sequence identity to the query sequence
      my $nonspacechar = 0;
      my $sameaa = 0;
      for(my $j = 0; $j < @qseq_con; $j++){
        if($qseq_con[$j] ne '-'){
          if($hseq_con[$j] eq 'B' or $hseq_con[$j] eq 'J' or $hseq_con[$j] eq 'O' or $hseq_con[$j] eq 'U' or $hseq_con[$j] eq 'X' or $hseq_con[$j] eq 'Z'){
            ;
          }
          else{ # normal amino acid
            $nonspacechar++;
            if($hseq_con[$j] eq $qseq_con[$j]){
              $sameaa++;
            }
          }
        }
      }
      ##determine if the sequence will be accepted or not.
      ##1) reject the hit sequence if there are too many gaps
      ##too many gaps mean the coverage is low => the hit may not be reliable
      if($nonspacechar/$query_len < 0.7){
        next;
      }
      ##2) reject the hit sequence if it is quite similar to the query
      #elsif($nonspacechar/$query_len>0.99 && $sameaa/$nonspacechar>0.99){
      #  next;
      #}
      
      
      ##2. output the hit sequence if it is accepted based on the above filters
      ##print '-' at the begining if any
      for(my $j = 1; $j < $query_from; $j++){
        printf FILE2 '-';
      }
      
      for(my $j = 0; $j < @qseq_con; $j++){
        if($qseq_con[$j] ne '-'){
          if($hseq_con[$j] eq 'B' or $hseq_con[$j] eq 'J' or $hseq_con[$j] eq 'O' or $hseq_con[$j] eq 'U' or $hseq_con[$j] eq 'X' or $hseq_con[$j] eq 'Z'){
            printf FILE2 '-';
          }
          else{ #normal amino acid
            printf FILE2 $hseq_con[$j];
          }
        }
      }
      ##print '-' in the end if any
      for(my $j = $query_to; $j < $query_len; $j++){
        printf FILE2 '-';
      }
      printf FILE2 "\n";
    }
  }
}
close $fh;
close FILE2;

