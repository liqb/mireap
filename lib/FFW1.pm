package FFW1;
use strict;
use Exporter;
our @ISA=qw(Exporter);
our @EXPORT=qw(ffw1 define_precursor);

use RNA; # Vienna RNA package perl interface
use Struct;

=head1 NAME

  FFW1 - Fold and Filter With 1 mature
  A method used to check whether a candidate sequence like a real pre-miRNA

=head1 SYNOPSIS

  use FFW1;
  $seq="ATGCTTCCGGCCTGTTCCCTGAGACCTCAAGTGTGAGTGTACTATTGATGCTTCACACCTGGGCTCTCCGGGTACCAGGACGGTTTGAGCAGAT";
  $tag="TCCCTGAGACCTCAAGTGTGA";
  $fold={};
  $pass=ffw1($seq,$fold,$tag,'max_energy' => -18, 'min_pair' => 14, 'size_diff' => 4,
    'max_bulge' => 4, 'asymmetry' => 5, 'min_space'=>5, 'max_space' => 55, 'flank' => 10);

=head1 AUTHOR
liqb <liqb@genomics.org.cn> 2008-10-20

=cut


# fold and filter with one small RNA
sub ffw1 {
	my($seq,$fold,$tag)=@_;
	my $pass=0;
	$fold->{pass}=0;
	$fold->{method}="ffw1";

	# parameters
	my $MAX_ENERGY=-18;
	my $MAX_UNPAIR=6; #
	my $MIN_PAIR=14;
	my $MAX_SIZEDIFF=4;
	my $MAX_BULGE=4;
	my $ASYMMETRY=5;
	my $MIN_UNPAIR=0;
	my $MIN_SPACE=5;
	my $MAX_SPACE=55;
	my $FLANK=10;
	while (@_) {
		my $argument = shift @_;
		if ($argument=~/max_energy/i) {$MAX_ENERGY=shift @_}
		if ($argument=~/max_unpair/i) {$MAX_UNPAIR=shift @_}
		if ($argument=~/min_pair/i) {$MIN_PAIR=shift @_}
		if ($argument=~/size_diff/i) {$MAX_SIZEDIFF=shift @_}
		if ($argument=~/max_bulge/i) {$MAX_BULGE=shift @_}
		if ($argument=~/asymmetry/i) {$ASYMMETRY=shift @_}
		if ($argument=~/min_unpair/i) {$MIN_UNPAIR=shift @_}
		if ($argument=~/min_space/i) {$MIN_SPACE=shift @_}
		if ($argument=~/max_space/i) {$MAX_SPACE=shift @_}
		if ($argument=~/flank/i) {$FLANK=shift @_}
	}
	
	# skip if the quality of reference if low
	my $N_count=$seq=~tr/N//;
	if ($N_count > 5) {
		$fold->{reason}="$N_count Ns (5)";
		return 0;
	}

	my $seq_length=length $seq;
	my $tag_length=length $tag;

	# position
	my $tag_beg=index($seq,$tag,0)+1;
	if ($tag_beg < 1) {
		 warn "[ffw1] coordinate error.\n";
		 $fold->{reason}="coordinate error";
		 return $pass;
	}
	my $tag_end=$tag_beg+length($tag)-1;
	
	
	# define candidate precursor by hybrid short arm to long arm, not solid enough
	my($beg,$end)=define_precursor($seq,$tag);
	if (not defined $beg) {
		$fold->{reason}="no precursor candidate";
		return $pass;
	}
	if (not defined $end) {
		$fold->{reason}="no precursor candidate";
		return $pass;
	}
	$fold->{beg}=$beg;
	$fold->{end}=$end;
	$seq=substr($seq,$beg-1,$end-$beg+1);
	$seq_length=length $seq;
	$fold->{seq}=$seq;
	$fold->{length}=$seq_length;
	
	
	# fold
	my ($struct,$mfe)=RNA::fold($seq);
	$mfe=sprintf "%.2f",$mfe;
	$fold->{struct}=$struct;
	$fold->{mfe}=$mfe;
	if ($mfe > $MAX_ENERGY) {
		$pass=0;
		$fold->{reason}="mfe=$mfe ($MAX_ENERGY)";
		return $pass;
	}

	# reposition
	$tag_beg=index($seq,$tag,0)+1;
	if ($tag_beg < 1) {
		 warn "[ffw1] coordinate error.\n";
		 $fold->{reason}="coordinate error";
		 return 0;
	}
	$tag_end=$tag_beg+length($tag)-1;

	my $tag_struct=substr($struct,$tag_beg-1,$tag_length);
	my $tag_arm=which_arm($tag_struct);
	my $tag_unpair=$tag_struct=~tr/.//;
	my $tag_pair=$tag_length-$tag_unpair;
	my $tag_max_bulge=biggest_bulge($tag_struct);
	if ($tag_arm eq "-") {$fold->{reason}="mature not in arm"; return $pass}
#	if ($tag_unpair > $MAX_UNPAIR) {$fold->{reason}="unpair=$tag_unpair ($MAX_UNPAIR)"; return $pass}
	if ($tag_pair < $MIN_PAIR) {$fold->{reason}="pair=$tag_pair ($MIN_PAIR)"; return $pass}
	if ($tag_max_bulge > $MAX_BULGE) {$fold->{reason}="maxbulge=$tag_max_bulge ($MAX_BULGE)"; return $pass}
	if ($tag_arm eq "5p") {
		$fold->{m5}=$tag;
		$fold->{m5pair}=$tag_pair;
		$fold->{m5status}="cloned";
	}
	else {
		$fold->{m3}=$tag;
		$fold->{m3pair}=$tag_pair;
		$fold->{m3status}="cloned";
	}

	# build base pairing table
	my %pairtable;
	&parse_struct($struct,\%pairtable); # coords count from 1
	
	# get star
	my ($star_beg,$star_end)=get_star(\%pairtable,$tag_beg,$tag_end);
	my $star=substr($seq,$star_beg-1,$star_end-$star_beg+1);
	my $star_length=$star_end-$star_beg+1;
	my $star_struct=substr($struct,$star_beg-1,$star_end-$star_beg+1);
	my $star_arm=which_arm($star_struct);
	my $star_unpair=$star_struct=~tr/.//;
	my $star_pair=$star_length-$star_unpair;
	my $star_max_bulge=biggest_bulge($star_struct);
	if ($star_arm eq "-") {$fold->{reason}="star not in arm"; return $pass}
#	if ($star_unpair > $MAX_UNPAIR) {$fold->{reason}="unpair=$star_unpair ($MAX_UNPAIR)"; return $pass}
	if ($star_pair < $MIN_PAIR) {$fold->{reason}="pair=$star_pair ($MIN_PAIR)"; return $pass}
	if ($star_max_bulge > $MAX_BULGE) {$fold->{reason}="maxbulge=$star_max_bulge ($MAX_BULGE)"; return $pass}
	if ($star_arm eq "5p") {
		$fold->{m5}=$star;
		$fold->{m5pair}=$star_pair;
		$fold->{m5status}="putative";
	}
	else {
		$fold->{m3}=$star;
		$fold->{m3pair}=$star_pair;
		$fold->{m3status}="putative";
	}

	if ($tag_arm eq $star_arm) {$fold->{reason}="two mature in same arm";return $pass}
	
	# space size between miR and miR*
	my $space;
	if ($tag_beg < $star_beg) {
		$space=$star_beg-$tag_end-1;
	}
	else {
		$space=$tag_beg-$star_end-1;
	}
	$fold->{space}=$space;
	if ($space < $MIN_SPACE) {$fold->{reason}="space=$space ($MIN_SPACE)"; return $pass}
	if ($space > $MAX_SPACE) {$fold->{reason}="space=$space {$MAX_SPACE}"; return $pass}

	# size diff
	my $size_diff=abs($tag_length-$star_length);
	$fold->{sizediff}=$size_diff;
	if ($size_diff > $MAX_SIZEDIFF) {$fold->{reason}="size_diff=$size_diff ($MAX_SIZEDIFF)"; return $pass}
	
	# asymmetry
	my $asy=get_asy(\%pairtable,$tag_beg,$tag_end);
	$fold->{asymmetry}=$asy;
	if ($asy > $ASYMMETRY) {$fold->{reason}="asymmetry=$asy ($ASYMMETRY)"; return $pass}
	
	#
	$fold->{DroshaOverhang}="-";
	$fold->{DicerOverhang}="-";

#	# GC content
#	my $hairpin_gc_n=$seq=~tr/GC//;
#	my $hairpin_gc_p=sprintf "%.2f", 100*$hairpin_gc_n/$seq_length;
#	my $duplex=$tag.$star;
#	my $duplex_gc_n=$duplex=~tr/GC//;
#	my $duplex_gc_p=sprintf "%.2f", 100*$duplex_gc_n/length($duplex);
	
	$pass=1;
	$fold->{pass}=1;	
	return $pass;
}



sub define_precursor {
	my $seq=shift;
	my $tag=shift;
	my $flank=shift || 10;

	my $seq_length=length $seq;
	my $tag_length=length $tag;
	my $tag_beg=index($seq,$tag,0)+1;
	my $tag_end=$tag_beg+$tag_length-1;

	my $MIN_PAIR=15; # minimal pair number required between short arm and long arm
	
	# split the candidate region into short arm and long arm
	my $tag_arm;
        my ($larm,$larm_beg,$larm_end);
	my ($sarm,$sarm_beg,$sarm_end);
        if ($tag_beg-1 < $seq_length-$tag_end) { # on 5' arm
                $sarm=substr($seq,0,$tag_end);
		$larm=substr($seq,$tag_end);
		$sarm_beg=1;
		$sarm_end=$tag_end;
		$larm_beg=$tag_end+1;
		$larm_end=$seq_length;
		$tag_arm="5p";
        }
        else {
                $larm=substr($seq,0,$tag_beg-1); # on 3' arm
		$sarm=substr($seq,$tag_beg-1);
		$larm_beg=1;
		$larm_end=$tag_beg-1;
		$sarm_beg=$tag_beg;
		$sarm_end=$seq_length;
		$tag_arm="3p";
        }
	
#	print "$sarm_beg,$sarm_end $sarm\n";
#	print "$larm_beg,$larm_end $larm\n";

	# clipping short arm
	if ($tag_arm eq "5p") {
		$sarm_beg=$tag_beg-$flank; $sarm_beg=1 if $sarm_beg < 1;
		$sarm=substr($seq,$sarm_beg-1,$sarm_end-$sarm_beg+1);
	}
	else {
		$sarm_end=$tag_end+$flank; $sarm_end=$seq_length if $sarm_end > $seq_length;
		$sarm=substr($seq,$sarm_beg-1,$sarm_end-$sarm_beg+1);
	}
#	print "$sarm_beg,$sarm_end $sarm\n";
#	print "$larm_beg,$larm_end $larm\n";

	# define the precursor by hybriding short arm to long arm
	my $duplex=RNA::duplexfold($sarm,$larm);
	my $struct=$duplex->{structure};
	my $energy=sprintf "%.2f", $duplex->{energy};
	my ($str1,$str2)=split(/&/,$struct);
	my $pair=$str1=~tr/(//;
#	print "pair=$pair\n";
	my $beg1=$duplex->{i}+1-length($str1);
	my $end1=$duplex->{i};
	my $beg2=$duplex->{j};
	my $end2=$duplex->{j}+length($str2)-1;
#	print "$beg1:$end1 $beg2:$end2\n";
	# transform coordinates
	$beg1=$beg1+$sarm_beg-1;
	$end1=$end1+$sarm_beg-1;
	$beg2=$beg2+$larm_beg-1;
	$end2=$end2+$larm_beg-1;
#	print "$beg1:$end1 $beg2:$end2\n";
	
	my $off5p=$beg1-$sarm_beg;
	my $off3p=$sarm_end-$end1;
	$beg2-=$off3p; $beg2=1 if $beg2 < 1;
	$end2+=$off5p; $end2=$seq_length if $end2 > $seq_length;

#	print "$beg1:$end1 $beg2:$end2\n";
	
	my $beg=$sarm_beg < $beg2 ? $sarm_beg : $beg2;
	my $end=$sarm_end > $end2 ? $sarm_end : $end2;
	
	return if $pair < $MIN_PAIR;
#	print "$beg,$end\n";
	return ($beg,$end);
}

1;

