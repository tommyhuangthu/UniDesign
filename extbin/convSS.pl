#!/usr/bin/perl
################################################################################
#                                                                              #
#                                   convSS.pl                                  #
#                                                                              #
# Converts the secondary structure of a protein into the I-TASSER format.      #
#                                                                              #
# ARGUMENT:                                                                    #
# ---> path to the file containing the protein's secondary structure.          #
#      The secondary structure is a horizontal sequence of consecutive         #
#      characters, the i-th of which refers to the i-th residue (i=0,...,N-1,  #
#      where N is the number of residues.)                                     #
#                                                                              #
# OUTPUT:                                                                      #
#  The I-TASSER representation of the secondary structure, printed to screen.  #
#  The i-th character refers to the i-th residue (i=0,...,N-1).                #
#  The translation rule is the following:                                      #
#  H ---> 1                                                                    #
#  E ---> 2                                                                    #
#  X ---> 3                                                                    #
#  where X stands for any character different from H and E.                    #
#                                                                              #
################################################################################

@data=`cat dssp.txt`;
$i=0;
while(1) {
	if($data[$i] =~ m/#/) { last; }
        $i++;
}
$i++;
$k=0;
open(FW,">sa.txt") or die;
open(FW1,">phi-psi.txt") or die;
for(;$i<=$#data;$i++) {
	$aa=substr $data[$i],13,1;
	$asa=substr $data[$i],34,5;
        $tot=sander($aa);
        if($tot==0) {
        	print "\n$aa\n";
                next;
        }
        $rsa=$asa/$tot;
        if($rsa<0.09)                   {       print FW "1"; } #'B'    
        elsif($rsa>=0.09 and $rsa<0.64) {       print FW "2"; } #'I'    
        if($rsa>=0.64)                  {       print FW "3"; } #'E'
	$k++;

	$phi = substr $data[$i],103,6;
	$psi = substr $data[$i],109,6;
	print FW1 "$phi $psi\n";
}
print FW "\n";
print FW1 "\n";
close(FW);
close(FW1);


open(INSS, "$ARGV[0]");
chomp($inss = <INSS>);
$len = length $inss;
for($i=0; $i<$len; $i++) {
	$ci = substr($inss, $i, 1);
	if($ci eq 'H') {
		print "1";
	} elsif($ci eq 'E') {
		print "2";
	} else {
		print "3";
	}
}
print "\n";
close(INSS);

if($len != $k) {
	print "Length error in ss.txt and sa.txt $length $k\n";
	exit(0);
}



#ABCDEFGHIKLMNPQRSTVWXYZ
#float SANDER[23]={106,160,135,163,194,197,84,184,169,205,164,188,157,136,198,248,130,142,142,227,180,222,196}

sub sander
{
        $aa=$_[0];
        if($aa eq 'A')  { return 106; }
        if($aa eq 'B')  { return 160; }
        if($aa eq 'C')  { return 135; }
        if($aa eq 'D')  { return 163; }
        if($aa eq 'E')  { return 194; }
        if($aa eq 'F')  { return 197; }
        if($aa eq 'G')  { return 84; }
        if($aa eq 'H')  { return 184; }
        if($aa eq 'I')  { return 169; }
        if($aa eq 'K')  { return 205; }
        if($aa eq 'L')  { return 164; }
        if($aa eq 'M')  { return 188; }
        if($aa eq 'N')  { return 157; }
        if($aa eq 'P')  { return 136; }
        if($aa eq 'Q')  { return 198; }
        if($aa eq 'R')  { return 248; }
        if($aa eq 'S')  { return 130; }
        if($aa eq 'T')  { return 142; }
        if($aa eq 'V')  { return 142; }
        if($aa eq 'W')  { return 227; }
        if($aa eq 'X')  { return 180; }
        if($aa eq 'Y')  { return 222; }
        if($aa eq 'Z')  { return 196; }

        return 0;
}

