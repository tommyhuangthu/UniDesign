#!/usr/bin/perl

###############################################################################
# GenerateProfile.pl
#
# Prepares profile file of a given target protein for subsequent design process
################################################################################

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
    Getopt::Long::Configure qw(gnu_getopt);
use Data::Dumper;
use File::Basename;
use Cwd 'abs_path';
use Benchmark qw(:hireswallclock);


my $starttime = Benchmark->new;
my $finishtime;
my $timespent;
my $hostName = `hostname`;
my $timeStart = localtime();

# initiate variables
my $cutoff=0.5;
my $pdbset="/nfs/amino-library/PDB";
if($hostName =~ /comet/){
  $pdbset="/oasis/projects/nsf/mia181/zhanglab/library/PDB";
}
my $bin_path = dirname(abs_path(__FILE__));

# parse options
GetOptions(
    'pdblib|l=s' => \$pdbset,
    'tmcutoff|t=f' => \$cutoff,
) or die "Illegal arguments are used for GenerateProfile.pl";

# parse input pdb
my $querypdb = $ARGV[0];

if(-e "prf.txt"){goto POS4; }
if(-e "msa.txt"){goto POS3; }
if(-e "msa_tmalgin.txt"){ goto POS1; }

print "#starting GenerateProfile.pl at $timeStart\n";
print "host: $hostName";

# start
print "structure alignment will be done against PDB structures available in: $pdbset\n";
my @decoys = glob("$pdbset/*.pdb");
my $ntmp = scalar @decoys;
printf "number of redundant pdb templates:  %d\n", $ntmp;
$ntmp = 0; # reset to zero

print "starting the structure alignment process, with TM-score cutoff: $cutoff \n";
print "TM-Scores with each template are given in first column...\n";

open TMFILE,  ">msa_tmalign.txt" or die "Error in creating msa_tmalign.txt\n";
open TMFASTA, ">msa_tmalign.fas" or die "Error in creating msa_tmalign.fas\n";

for(my $i = 0; $i<@decoys; ++$i){
  chomp(my $dec = $decoys[$i]);
  my @pdbpath=split(/\//, $dec);
  my $PDBID = $pdbpath[$#pdbpath];
  $PDBID =~s/\.pdb//g;
  my $idLength = length $PDBID;
  
  # run TM-align for full chains only; skip fragment domains
  my @tmo;
  if($idLength<6){
    $ntmp++; # counter
    @tmo = `$bin_path/TMalign $dec $querypdb `;
  }

  my $tmScore;
  my $homology = 0;
  my $seqid;
  for(my $j=0; $j<@tmo; $j++) {
    if($tmo[$j] =~ /TM-score=/) {
      # extract TM-score from TM-align output
      my @dat1=split(/,/, $tmo[$j]);
      my @dat2=split(/=/, $dat1[2]);
      # Printing TM-Scores for each decoy
      $tmScore = $dat2[1];
      print "$dat2[1]\t$dec\n";
      `echo "$dat2[1]\t$PDBID" >> tmscore.txt`;

      # print the sequence identity
      my @dat3=split(/=/, $dat1[3]);
      chomp($seqid = $dat3[1]);
      `echo "$seqid\t$PDBID" >> seqid.txt`;

      if($dat2[1]>$cutoff){ $homology=1; }
      else{ $homology=0; }
      last;
    }
  }
 
  # skip nonhomologues
  if($homology==0) {next;}
  `echo "$seqid\t$PDBID" >> msa.seqid`;

  my ($dali,$tali);
  for(my $j=0; $j<@tmo; $j++) {
    if($tmo[$j] =~ /denotes residue/) {
      chomp($dali = $tmo[$j+1]);
      chomp($tali = $tmo[$j+3]);
      last;
    }
  }

  my ($lena,$ti);
  my @ali;
  $lena = length $tali;
  my $sameaa = 0;
  my $notgap = 0;
  for(my $j=0; $j<$lena; $j++){
    my $char_t = substr($tali, $j, 1);
    if($char_t ne '-') {
      my $char_d = substr($dali, $j, 1);
      $ali[$i][$ti++] = $char_d;
      if($char_d ne '-'){ $notgap++; }
      if($char_t eq $char_d){ $sameaa++; }
    }
  }

  # set filters to get high-quality template structural homologs
  # print "Sequence identity (comparted to query): %f\n", $sameaa/$ti;
  # 1. reject sequence with too many gaps
  if($notgap/$ti < 0.7){ next; }
  # 2. reject sequences nearly identical to the query sequence
  #elsif($notgap/$ti>0.99 && $sameaa/$notgap>0.99){
  #  next;
  #}

  #Print-Decoy-Sequences 
  print TMFASTA "> $PDBID; TM-score=$tmScore \n";
  for(my $j=0; $j<$ti; $j++) {
    print TMFILE $ali[$i][$j];
    print TMFASTA $ali[$i][$j];
  }
  print TMFILE "\n";
  print TMFASTA "\n";
}
    
close TMFILE;
close TMFASTA;

# total non-redundant PDBs in database 
print "\n";
printf "number of non-redundant PDB templates searched: %d\n", $ntmp;

# sort TM-score and sequence identity from high to low 
`sort -k1 -n -r tmscore.txt > tmscore_sorted.txt`;
`sort -k1 -n -r seqid.txt > seqid_sorted.txt`;

POS1:

`cat msa_tmalign.txt > msa.txt`;

if(-e "msa_psiblast.txt"){ goto POS2; }
# if the number of TM-align hits is <10, include psi-blast hits
my $MIN_ANALOG_NUM = 10;
my $analog_num = `cat msa_tmalign.txt | wc -l`;

if($analog_num<$MIN_ANALOG_NUM){
  my $database = "/nfs/amino-library/nr/nr";
  if($hostName =~ /comet/){
    $database = "/oasis/projects/nsf/mia181/zhanglab/library/nr/nr";
  }
  # make fasta
  `$bin_path/PDB2FAS $querypdb > query.fasta`;
  my $arguments = "-query query.fasta -out blast.xml -db $database -outfmt 5 -evalue 1e-4";
  # do psi-blast search and write in .xml format
  `$bin_path/ncbi-blast-2.7.1+/bin/psiblast $arguments`;
  `$bin_path/xml2msa.pl -in blast.xml -out msa_psiblast.txt`;

POS2:

  # combine the msa obtained from tmalign & psi-blast
  my (@msa1, @msa2);
  open MSA1, "<msa_tmalign.txt";
  @msa1 = <MSA1>; chomp(@msa1);
  close MSA1;
  open MSA2, "<msa_psiblast.txt";
  @msa2 = <MSA2>; chomp(@msa2);
  close MSA2;
  for(my $i=0; $i<@msa2; $i++){
    my $same = 0;
    for(my $j=0; $j<@msa1; $j++){
      if($msa1[$j] eq $msa2[$i]){ $same=1; last;}
    }
    if($same == 0){ `echo $msa2[$i] >> msa.txt`; }
  }

  `rm query.fasta`;
}

POS3:
print "generating position specific scoring matrix... \n";
system "$bin_path/mkprf msa.txt > prf.txt";

POS4:
$finishtime = Benchmark->new;
$timespent = timediff($finishtime,$starttime);
print "generating Profile, computational time spent: ". timestr($timespent);
my $timeEnd = localtime();
print "#normal ending of GenerateProfile.pl at $timeEnd\n";
exit;
