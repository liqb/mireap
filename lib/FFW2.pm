package FFW2;
use strict;
use Exporter;
our @ISA=qw(Exporter);
our @EXPORT=qw(ffw2 likeMirDuplex1);

use RNA; # Vienna RNA package perl Interface
use Struct;


=head1 NAME

  FFW2 - Fold and Filter With 2 short small RNAs
  A method used to evaluate whether the input sequence is a real pre-miRNA

=head1 SYNOPSIS

  use FFW2;
  $seq = "ATGCTTCCGGCCTGTTCCCTGAGACCTCAAGTGTGAGTGTACTATTGATGCTTCACACCTGGGCTCTCCGGGTACCAGGACGGTTTGAGCAGAT";
  $tag1 = "TCCCTGAGACCTCAAGTGTGA";
  $tag2 = "ACACCTGGGCTCTCCGGGTACC";
  $fold = {};
  $pass = ffw2($seq, $fold, $tag1, $tag2, 'max_bulge' => 4,'asymmetry' => 3);
  $pass = ffw2($seq, $fold, $tag1, $tag2, 'max_energy' => -18, 'min_pair' => 14, 'size_diff' => 4,
    'max_bulge' => 4, 'asymmetry' => 5, 'min_space'=>5, 'max_space' => 55, 'flank' => 10);


=head1 AUTHOR
liqb <liqb@genomics.cn> 

=head1 BUG

2014-10-10
Kira Neller <nellerk@yorku.ca> reports a bug at line 109:  my declaration masks earlier declaration in same scope

=cut


sub ffw2{
	my ($seq,$fold,$tag1,$tag2)=@_;
	my $pass=0;
	$fold->{pass}=0;
	$fold->{method}="ffw2";

	# parameters
	my $MAX_ENERGY=-18;
	my $MAX_UNPAIR=6;
	my $MIN_PAIR=14;
	my $MAX_SIZEDIFF=4;
	my $MAX_BULGE=4;
	my $ASYMMETRY=5;
	my $MIN_UNPAIR=0;
	my $MIN_SPACE=5;
	my $MAX_SPACE=55;
	my $FLANK=15;
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

	my $N_count=$seq=~tr/N//;
	if ($N_count > 5) {
		$fold->{reason}="$N_count Ns (5)";
		return 0;
	}

	my $seq_length=length $seq;
	# position tag1 and tag2
	my $tag1_beg=index($seq,$tag1,0)+1;
	if ($tag1_beg < 1) {
		warn "[ffw2] coordinate error.\n";
		$fold->{reason}="coordinate error";
		return 0;
	}
	my $tag2_beg=index($seq,$tag2,0)+1;
	if ($tag2_beg < 1) {
		warn "[ffw2] coordinate error.\n";
		$fold->{reason}="coordinate error";
		return 0;
	}
	if ($tag2_beg < $tag1_beg) {
		# swap tag1 and tag2
		($tag1,$tag2)=($tag2,$tag1);
		($tag1_beg,$tag2_beg)=($tag2_beg,$tag1_beg);
	}
	my $tag1_end=$tag1_beg+length($tag1)-1;
	my $tag2_end=$tag2_beg+length($tag2)-1;

	# clipping
	my $beg=$tag1_beg-$FLANK; $beg=1 if $beg < 1;
	my $end=$tag2_end+$FLANK; $end=$seq_length if $end > $seq_length;
	$seq=substr($seq,$beg-1,$end-$beg+1);
	$seq_length=length $seq;
	$fold->{beg}=$beg;
	$fold->{end}=$end;
	$fold->{seq}=$seq;

	# reposition
	$tag1_beg=index($seq,$tag1,0)+1;
	if ($tag1_beg < 1) {
		warn "[ffw2] coordinate error.\n";
		$fold->{reason}="coordinate error";
		return 0;
	}
	$tag2_beg=index($seq,$tag2,0)+1; # a bug fixed here
	if ($tag2_beg < 1) {
		warn "[ffw2] coordinate error.\n";
		$fold->{reason}="coordinate error";
		return 0;
	}
	$tag1_end=$tag1_beg+length($tag1)-1;
	$tag2_end=$tag2_beg+length($tag2)-1;
	
	# fold
	my($struct,$mfe)=RNA::fold($seq);
	$mfe=sprintf "%.2f", $mfe;
	$fold->{struct}=$struct;
	$fold->{mfe}=$mfe;
	if ($mfe > $MAX_ENERGY) {$fold->{reason}="mfe=$mfe ($MAX_ENERGY)"; return $pass}

	# tag1
	my $tag1_length=$tag1_end-$tag1_beg+1;
	my $tag1_struct=substr($struct,$tag1_beg-1,$tag1_length);
	my $tag1_arm=which_arm($tag1_struct);
	my $tag1_unpair=$tag1_struct=~tr/.//;
	my $tag1_pair=$tag1_length-$tag1_unpair;
	my $tag1_max_bulge=biggest_bulge($tag1_struct);
	if ($tag1_arm ne "5p") {$fold->{reason}="mature5 not in arm"; return $pass} # tag not in stem
#	if ($tag1_unpair > $MAX_UNPAIR) {$fold->{reason}="unpair=$tag1_unpair ($MAX_UNPAIR)"; return $pass}
	if ($tag1_pair < $MIN_PAIR) {$fold->{reason}="pair=$tag1_pair ($MIN_PAIR)"; return $pass}
	if ($tag1_max_bulge > $MAX_BULGE) {$fold->{reason}="maxbulge=$tag1_max_bulge ($MAX_BULGE)"; return $pass}
	$fold->{m5}=$tag1;
	$fold->{m5pair}=$tag1_pair;
	$fold->{m5status}="cloned";

	# tag2
	my $tag2_length=$tag2_end-$tag2_beg+1;
	my $tag2_struct=substr($struct,$tag2_beg-1,$tag2_length);
	my $tag2_arm=which_arm($tag2_struct);
	my $tag2_unpair=$tag2_struct=~tr/.//;
	my $tag2_pair=$tag2_length-$tag2_unpair;
	my $tag2_max_bulge=biggest_bulge($tag2_struct);
	if ($tag2_arm ne "3p") {$fold->{reason}="mature not in arm"; return $pass} # star not in stem
#	if ($tag2_unpair > $MAX_UNPAIR) {$fold->{reason}="unpair=$tag2_unpair ($MAX_UNPAIR)"; return $pass}
	if ($tag2_pair < $MIN_PAIR) {$fold->{reason}="pair=$tag2_pair ($MIN_PAIR)"; return $pass}
	if ($tag2_max_bulge > $MAX_BULGE) {$fold->{reason}="maxbulge=$tag2_max_bulge ($MAX_BULGE)"; return $pass}
	$fold->{m3}=$tag2;
	$fold->{m3pair}=$tag2_pair;
	$fold->{m3status}="cloned";

	# space size between miR and miR*
	my $space=$tag2_beg-$tag1_end-1;
	$fold->{space}=$space;
	if ($space < $MIN_SPACE) {$fold->{reason}="space=$space ($MIN_SPACE)"; return $pass}
	if ($space > $MAX_SPACE) {$fold->{reason}="space=$space ($MAX_SPACE)"; return $pass}
	
	# size diff of miR and miR*
	my $size_diff=abs($tag1_length-$tag2_length);
	$fold->{sizediff}=$size_diff;
	if ($size_diff > $MAX_SIZEDIFF) {$fold->{reason}="size_diff=$size_diff ($MAX_SIZEDIFF)"; return $pass}

	# build base pairing table
	my %pairtable;
	&parse_struct($struct,\%pairtable); # coords count from 1

	my $asy1=get_asy(\%pairtable,$tag1_beg,$tag1_end);
	my $asy2=get_asy(\%pairtable,$tag2_beg,$tag2_end);
	my $asy=($asy1 < $asy2) ? $asy1 : $asy2;
	$fold->{asymmetry}=$asy;
	if ($asy > $ASYMMETRY) {$fold->{reason}="asymmetry=$asy ($ASYMMETRY)"; return $pass}
	

	# duplex fold, determine whether two matures like a miR/miR* ike duplex
	my ($like_mir_duplex1,$duplex_pair,$overhang1,$overhang2)=likeMirDuplex1($tag1,$tag2);
	# parse hairpin, determine whether two matures form miR/miR* duplex in hairpin context
	my ($like_mir_duplex2,$duplex_pair2,$overhang_b,$overhang_t)=likeMirDuplex2(\%pairtable,$tag1_beg,$tag1_end,$tag2_beg,$tag2_end);
	if ($like_mir_duplex1==0 && $like_mir_duplex2==0) {
		$fold->{reason}="likeMirDuplex=0";
		return $pass;
	}
	$fold->{DroshaOverhang}=$overhang_b;
	$fold->{DicerOverhang}=$overhang_t;

	$pass=1;
	$fold->{pass}=1;
	return $pass;	
#	# GC content
#	my $hairpin_gc_n=$seq=~tr/GC//;
#	my $hairpin_gc_p=sprintf "%.2f", 100*$hairpin_gc_n/$seq_length;
#	my $duplex=$tag1.$tag2;
#	my $duplex_gc_n=$duplex=~tr/GC//;
#	my $duplex_gc_p=sprintf "%.2f", 100*$duplex_gc_n/length($duplex);

}

# duplex fold, judge whether two short seqs like a miRNA/miRNA* duplex
sub likeMirDuplex1 {
	my $seq1=shift;
	my $seq2=shift;
	my $like_mir_duplex=1;
	
	my $length1=length $seq1;
	my $length2=length $seq2;
	my $duplex=RNA::duplexfold($seq1, $seq2);
	my $duplex_struct=$duplex->{structure};
	my $duplex_energy=sprintf "%.2f", $duplex->{energy};
	my ($str1,$str2)=split(/&/,$duplex_struct);
	my $beg1=$duplex->{i}+1-length($str1);
	my $end1=$duplex->{i};
	my $beg2=$duplex->{j};
	my $end2=$duplex->{j}+length($str2)-1;
	
	# revise beg1, end1, beg2, end2
	$str1=~/^(\.*)/;
	$beg1+=length($1);
	$str1=~/(\.*)$/;
	$end1-=length($1);
	$str2=~/^(\.*)/;
	$beg2+=length($1);
	$str2=~/(\.*)$/;
	$end2-=length($1);

	my $pair_num=$str1=~tr/(//;
	my $overhang1=($length2-$end2)-($beg1-1); # 3' overhang at hairpin bottom
	my $overhang2=($length1-$end1)-($beg2-1); # 3' overhang at hairpin neck
#	print $pair_num,"\n";
#	print $overhang1,"\n";
#	print $overhang2,"\n";
	if ($pair_num < 13) {
		$like_mir_duplex=0;
	}
	if ($overhang1 < 0 || $overhang2 < 0 ) {
		$like_mir_duplex=0;
	}
	if ($overhang1 > 4 || $overhang2 > 4) {
		$like_mir_duplex=0;
	}
	return ($like_mir_duplex,$pair_num,$overhang1,$overhang2);
}

# judge whether two matures form miR/miR* duplex, in hairpin context
sub likeMirDuplex2 {
	my ($table,$beg1,$end1,$beg2,$end2)=@_;
	my $like_mir_duplex=1;

#	   s1         e1
#   5 ----------------------------3	
#      | | |||| |||               |
#3 -------------------------------5
#      e2         s2

	my $pair_num;
	my $overhang1;
	my $overhang2;
	my ($s1,$e1,$s2,$e2);
	foreach my $i ($beg1..$end1) {
		if (defined $table->{$i}) {
			my $j=$table->{$i};
			if ($j <= $end2 && $j >= $beg2) {
				$s1=$i;
				$e2=$j;
				last;
			}
		}
	}
	foreach my $i (reverse ($beg1..$end1)) {
		if (defined $table->{$i}) {
			my $j=$table->{$i};
			if ($j <= $end2 && $j >= $beg2) {
				$e1=$i;
				$s2=$j;
				last;
			}
		}
	}

#	print "$beg1,$end1 $s1,$e1\n";
#	print "$beg2,$end2 $s2,$e2\n";

	foreach my $i ($beg1..$end1) {
		if (defined $table->{$i}) {
			my $j=$table->{$i};
			if ($j <= $end2 && $j >= $beg2) {
				++$pair_num;
			}
		}
	}
	if (defined $s1 && defined $e2) {
		$overhang1=($end2-$e2)-($s1-$beg1);
	}
	if (defined $e1 && defined $s2) {
		$overhang2=($end1-$e1)-($s2-$beg2);
	}
	
	if ($pair_num < 13) {
		$like_mir_duplex=0;
	}
	if ($overhang1 < 0 && $overhang2 < 0) {
		$like_mir_duplex=0;
	}
	return ($like_mir_duplex,$pair_num,$overhang1,$overhang2);
}

1;

