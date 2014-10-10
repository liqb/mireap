#!/usr/bin/perl
# Reap microRNAs from deeply sequenced small RNA library

use strict;
use FFW1;
use FFW2;
use Struct;
use Getopt::Std;

# ------------------------------------------------------------------
# Options And Usage
# ------------------------------------------------------------------
use vars qw($opt_i $opt_m $opt_r $opt_o $opt_t $opt_A $opt_B $opt_a
    $opt_b $opt_u $opt_e $opt_d $opt_p $opt_v $opt_s $opt_f $opt_h);
getopts('i:m:r:o:t:A:B:a:b:u:e:d:p:v:s:f:h');

my $prog_name="mireap";
my $version="0.2";
my $seq_file=$opt_i;
my $map_file=$opt_m;
my $ref_file=$opt_r;
my $out_dir=$opt_o ? $opt_o : "./";
my $label=$opt_t ? $opt_t : "xxx";

my $MIN_LEN=$opt_A ? $opt_A : 18;  # minimal miRNA sequence length
my $MAX_LEN=$opt_B ? $opt_B : 26;  # maximal miRNA sequence length
my $MIN_ref_LEN=$opt_a ? $opt_a : 20; # minimal miRNA reference sequence length
my $MAX_ref_LEN=$opt_b ? $opt_b : 24; # minimal miRNA reference sequence length
my $UNIQNESS=$opt_u ? $opt_u : 20; # uniqueness

my $MAX_ENERGY=$opt_e ? $opt_e : -18;
my $MIN_SPACE=5; # minimal space size between miRNA and miRNA*
my $MAX_SPACE=$opt_d ? $opt_d : 35; # maximal space size between miRNA and miRNA*
my $MIN_PAIR=$opt_p ? $opt_p : 14; # minimal pair of miRNA and miRNA*
my $MAX_BULGE=$opt_v ? $opt_v : 4; # maximal bulge of miRNA and miRNA*
my $MAX_SIZEDIFF=4; # maximal size different of miRNA and miRNA*
my $MAX_ASY=$opt_s ? $opt_s : 5; # maximal asymmetry of miRNA/miRNA* duplex
my $FLANK=$opt_f ? $opt_f : 10;
my $help=$opt_h ? 1 : 0;
my $test=0;

# constant
my $DOMINANCE=0.9;
my $MATURE_SIZE=22; # canonical size of an miRNA mature sequence

# usage
unless (-e $seq_file) {
	usage(); exit;
}
unless (-e $map_file) {
	usage(); exit;
}
unless (-e $ref_file) {
	usage(); exit;
}
if ($help) {
	usage(); exit;
}
unless (-e $out_dir) {
	mkdir $out_dir;
}
$out_dir .="/" unless $out_dir =~/\/$/;


my $log_file=$out_dir."$prog_name-$label.log";
open LOG, ">$log_file" || die $!;
my $start_time=localtime;
print LOG "$prog_name start at $start_time\n\n";
print LOG "Parameters:\n";
print LOG "  minimal miRNA length: $MIN_LEN\n";
print LOG "  maximal miRNA length: $MAX_LEN\n";
print LOG "  minimal miRNA(reference) length: $MIN_ref_LEN\n";
print LOG "  maximal miRNA(reference) length: $MAX_ref_LEN\n";
print LOG "  uniqueness of miRNA: $UNIQNESS\n";
print LOG "  maximal energy: $MAX_ENERGY\n";
print LOG "  minimal space: $MIN_SPACE\n";
print LOG "  maximal space: $MAX_SPACE\n";
print LOG "  minimal mature pair: $MIN_PAIR\n";
print LOG "  maximal mature bulge: $MAX_BULGE\n";
print LOG "  maximal duplex asymmetry: $MAX_ASY\n";
print LOG "  flank sequence length: $FLANK\n\n";

# ------------------------------------------------------------------
# load smrna file (fasta format, with count)
# ------------------------------------------------------------------
print LOG "Load short tag file: $seq_file\n";
my %tagId_seq;
my %tagId_count;
my $uniq_tag_num;
my $total_tag_count;
open IN, $seq_file || die $!;
while (<IN>) {
	if (/^>(\S+)/) {
		my $tagId=$1;
		my $count=0;
		if (/>\S+\s+(\d+)/) {
			$count=$1;
		}
		my $seq=<IN>; chomp $seq;
		
		my $length=length $seq;
		if ($length >=$MIN_LEN && $length <=$MAX_LEN) {
			$tagId_seq{$tagId}=$seq;
			$tagId_count{$tagId}=$count;
			++$uniq_tag_num;
			$total_tag_count+=$count;
		}
	}
}
close IN;
print LOG "  small RNAs [$MIN_LEN - $MAX_LEN nt]\n";
print LOG "  unique: $uniq_tag_num\n";
print LOG "  total:  $total_tag_count\n\n";


# ------------------------------------------------------------------
# load mapping file (text file)
# ------------------------------------------------------------------
print LOG "Load mapping file: $map_file\n";
my %tagId_hitn;
my %refId_tags;
my $mapped_c;
my $mapped_n;
open IN, $map_file || die $!;
while (<IN>) {
	next if /^#/;
	chomp;
	my @d=split;
	my $length=$d[3]-$d[2]+1;
	next if ($length < $MIN_LEN || $length > $MAX_LEN);
	my $tag=[$d[0],$d[2],$d[3],$d[4]]; # [tagId, beg, end, +/-]
	++$tagId_hitn{$d[0]};
	my $refId=join(":",$d[1],$d[4]); # "refId:+/-"
	push @{$refId_tags{$refId}},$tag;
}
close IN;
foreach my $tagId (keys %tagId_hitn) {
	$mapped_c++;
	$mapped_n+=$tagId_count{$tagId};
}
print LOG "  small RNAs mapped [$MIN_LEN - $MAX_LEN nt]\n";
print LOG "  unique: $mapped_c\n";
print LOG "  total:  $mapped_n\n\n";


# ------------------------------------------------------------------
# load reference file (fasta format)
# ------------------------------------------------------------------
print LOG "Load reference file: $ref_file\n";
my $refId_seq={};
my $ref_entry_n=readfasta($ref_file,$refId_seq);
print LOG "  total $ref_entry_n entries\n\n";


# ------------------------------------------------------------------
# Identify miRNA candidates
# ------------------------------------------------------------------
print LOG "miRNA discovery\n";
my $gff_file=$out_dir."$prog_name-$label.gff";
my $aln_file=$out_dir."$prog_name-$label.aln";
my $tmp_file=$out_dir."$prog_name-$label.tmp";
open GFF, ">$gff_file" || die;
open ALN, ">$aln_file" || die;
if ($test) {
    open TMP, ">$tmp_file" || die;
}

my $c=0;
foreach (sort keys %refId_tags) {
	next unless ($refId_tags{$_});
	my @tags=@{$refId_tags{$_}};
	my ($refId,$strand)=split(/:/,$_);
	# sort tags by start position
	@tags=sort {$a->[1] <=> $b->[1]} @tags;
	
	# -----------------------------------------------------------------------
	# Step 1:  Screen for candidates miRNA sites
	# -----------------------------------------------------------------------
	# 5' breakpoints defined by small RNA mapping
	my %break_count; # total count of tags that define this breakpoint
	my %break_indexes; # indexes of tags that define this breakpoint
	foreach my $k (0..$#tags) {
		my $tagId=$tags[$k][0];
		my $break;
		if ($strand eq "+") {
			$break=$tags[$k][1];
		}
		else {
			$break=$tags[$k][2];
		}
		$break_count{$break}+=$tagId_count{$tagId};
		push @{$break_indexes{$break}},$k;
	}
	my $break_n=scalar keys %break_count;
#	print LOG "\t$refId:$strand breakpoints: $total_break_n\n";

	# screen for Drosha/Dicer cutting sites
	my @csites;
	my %csite_indexes;  # indexes of @tags
	my %csite_matureix;  # mature index
	my %csite_count;
	foreach my $break (sort {$a <=> $b} keys %break_count) {
		my $count=$break_count{$break};
		
		my $likeCS=1; # like Drosha/Dicer cutting site
		# only consider stable breakpoints
		$likeCS=0 if $break_count{$break} < 3;
		next if $likeCS==0;
		# Is the dominant breakpoint from 20 nt upstream to 20 nt downstream
		foreach my $i (1..20) {
			if ($break_count{$break+$i} > $count) {
				$likeCS=0;
			}
			if ($break-$i > 0 && $break_count{$break-$i} >= $count) {
				$likeCS=0;
			}
		}
		next if $likeCS==0;
		# Drosha and Dicer processing are not absolutely accurate, so we also take
		# breakpoint 2 nt upstream or 2nt downstream as valid Dicer/Drosha cutting site,
		# while sites beyond this range as random breakpoints.
		my $count_csites=$count;
		my $count_region=$count;
		foreach my $i (1..20) {
			$count_region+=$break_count{$break+$i};
			if ($break-$i > 0) {
				$count_region+=$break_count{$break-$i};
			}
		}
		foreach my $i (1..2) {
			$count_csites+=$break_count{$break+$i};
			if ($break-$i > 0) {
				$count_csites+=$break_count{$break-$i};
			}
		}
		if ($count_csites/$count_region < $DOMINANCE) {
			$likeCS=0;
		}
		if ($count/$count_csites < 0.5) { # 080912
			$likeCS=0;
		}
		next if $likeCS==0;
		
		# tags belong to this candidate site
		my @indexes; # indexes of tags belong to this candidate site
		my $sum_count;
		foreach my $i (-2..2) {
			next if $break+$i <= 0;
			if (defined $break_indexes{$break+$i}) {
				foreach my $i (@{$break_indexes{$break+$i}}) {
					push @indexes,$i;
					my $tagId=$tags[$i][0];
					$sum_count+=$tagId_count{$tagId};
				}
			}
		}
		# dominant tag
		my $max_index; # index of the dominant tag
		my $max_count=0;
		foreach my $i (@indexes) {
			my $tagId=$tags[$i][0];
			my $count=$tagId_count{$tagId};
			if ($count > $max_count) {
				$max_count=$count;
				$max_index=$i;
			}
		}
		my $matureix=$max_index;
		my $mature=$tags[$matureix];
		my $mature_length=$mature->[2]-$mature->[1]+1;
		if ($mature_length < $MIN_ref_LEN || $mature_length > $MAX_ref_LEN) {
			$likeCS=0;
		}
		next if $likeCS==0;

		push @csites,$break if $likeCS;
		$csite_indexes{$break}=\@indexes;
		$csite_matureix{$break}=$matureix;
		$csite_count{$break}=$sum_count;
	}
	my $csites_n = scalar @csites;
#	print LOG "\t$refId:$strand candidate csites: $csites_n\n";

	# remove cutting site defined by highly repeated tags
	my %inrepeat;
	my @csites_repeat; # cutting sites defined by highly repeated tags
	my @csites_remain;
	foreach my $csite (@csites) {
		my $matureix=$csite_matureix{$csite};
		my $mature=$tags[$matureix]; # mature is the representative tag of this csite
		my $mature_Id=$mature->[0];
		my $hitn=$tagId_hitn{$mature_Id};
		my $mature_seq=$tagId_seq{$mature_Id};
		my $mature_count=$tagId_count{$mature_Id};
		if ($hitn > $UNIQNESS) {
			$inrepeat{$csite}=1;
			push @csites_repeat,$csite;
		}
		else {
			$inrepeat{$csite}=0;
			push @csites_remain,$csite;
		}
	}
	my $csites_repeat_n=scalar @csites_repeat;
	my $csites_remain_n=scalar @csites_remain;
#	print LOG "\t$refId:$strand csites_repeat: $csites_repeat_n csites_remain: $csites_remain_n\n";
	
	
	# ----------------------------------------------------------------------------
	# Step 2: Structure validation and Deep testing
	# ----------------------------------------------------------------------------
	# fold and fiter
	my %csite_pass; # csites pass filter
	foreach my $i (0..$#csites) {
		my $csite=$csites[$i];
		next if $inrepeat{$csite}; # skip csites locate at repeat region
		next if $csite_pass{$csite}; # avoid repeated prediction
		my $matureix=$csite_matureix{$csite};
		my $mature=$tags[$matureix];
		my $mature_Id=$mature->[0];
		my $mature_beg=$mature->[1];
		my $mature_end=$mature->[2];
		my $mature_seq=$tagId_seq{$mature_Id};
		my $mature_seqref=subseq($refId_seq,$refId,$mature_beg,$mature_end,$strand);
		my $mature_count=$tagId_count{$mature_Id};
		
		### downstream extension
		my $dsext=0; # flag, 0 means no good precursor, has to try again
		my $dsfold; #
		# case 1: has a closely located downstream cutting site, test by ffw2 only
		if ($dsext==0 && defined $csites[$i+1]) {
			my $csite_next=$csites[$i+1];
			my $matureix_next=$csite_matureix{$csite_next};
			my $mature_next=$tags[$matureix_next];
			my $mature_next_Id=$mature_next->[0];
			my $mature_next_beg=$mature_next->[1];
			my $mature_next_end=$mature_next->[2];
			my $mature_next_seq=$tagId_seq{$mature_next_Id};
			my $mature_next_seqref=subseq($refId_seq,$refId,$mature_next_beg,$mature_next_end,$strand);
			my $dist=$mature_next_beg-$mature_end-1;
			if ($dist <= $MAX_SPACE && $dist >= $MIN_SPACE) {
				my $p1=$mature_beg-$FLANK; $p1=1 if $p1 < 1;
				my $p2=$mature_next_end+$FLANK;
				my $seq=subseq($refId_seq,$refId,$p1,$p2,$strand);
				$p2=$p1+length($seq)-1;
				$dsext=1;
				my $fold={};
				ffw2($seq,$fold,$mature_seqref,$mature_next_seqref,'max_energy'=>$MAX_ENERGY,'min_pair'=>$MIN_PAIR,'size_diff'=>$MAX_SIZEDIFF,'max_bulge'=>$MAX_BULGE,'asymmetry'=>$MAX_ASY,'min_space'=>$MIN_SPACE,'max_space'=>$MAX_SPACE,'flank'=>$FLANK);
				if ($test) {print TMP "$refId:$p1:$p2:$strand\t","ffw2($seq,$mature_seqref,$mature_next_seqref)\t",$fold->{pass},"\t",$fold->{reason},"\n"}
				$fold->{location}="$refId:$p1:$p2:$strand";
				if ($fold->{pass}) {
					$dsfold=$fold;
				}
			}
		}
		# case 2: has downstream tags, may be mature parter, test by ffw2
		if ($dsext==0) {
			# collect downstream tags
			my @indexes_after;
			for (my $j=$matureix+1; $j<$#tags; $j++) {
				my $tagj=$tags[$j];
				my $tagj_Id=$tagj->[0];
				my $tagj_beg=$tagj->[1];
				my $tagj_end=$tagj->[2];
				my $tagj_seq=$tagId_seq{$tagj_Id};
				my $tagj_count=$tagId_count{$tagj_Id};
				my $dist=$tagj_beg-$mature_end-1;
				if ($dist < $MIN_SPACE) {
					next;
				}
				elsif ($dist > $MAX_SPACE) {
					last;
				}
				else {
					push @indexes_after,$j;
				}
			}
			if (@indexes_after) {
				# sort downstream tags by count
				@indexes_after=sort {$tagId_count{$tags[$b][0]} <=> $tagId_count{$tags[$a][0]}} @indexes_after;
				# only consider the dominant one
				my $j=$indexes_after[0];
				my $tagj=$tags[$j];
				my $tagj_Id=$tagj->[0];
				my $tagj_beg=$tagj->[1];
				my $tagj_end=$tagj->[2];
				my $tagj_seq=$tagId_seq{$tagj_Id};
				my $tagj_seqref=subseq($refId_seq,$refId,$tagj_beg,$tagj_end,$strand);
				my $tagj_count=$tagId_count{$tagj_Id};
				my $p1=$mature_beg-$FLANK; $p1=1 if $p1 < 1;
				my $p2=$tagj_end+$FLANK;
				my $seq=subseq($refId_seq,$refId,$p1,$p2,$strand);
				$p2=$p1+length($seq)-1;
				my $fold={};
				ffw2($seq,$fold,$mature_seqref,$tagj_seqref,'max_energy'=>$MAX_ENERGY,'min_pair'=>$MIN_PAIR,'size_diff'=>$MAX_SIZEDIFF,'max_bulge'=>$MAX_BULGE,'asymmetry'=>$MAX_ASY,'min_space'=>$MIN_SPACE,'max_space'=>$MAX_SPACE,'flank'=>$FLANK);
				if ($test) {print TMP "$refId:$p1:$p2:$strand\t","ffw2($seq,$mature_seqref,$tagj_seqref)\t",$fold->{pass},"\t",$fold->{reason},"\n"}
				$fold->{location}="$refId:$p1:$p2:$strand";
				if ($fold->{pass}) {
					$dsfold=$fold;
					$dsext=1;
				}
			}
		}
		# case 3: case 2 failed or do not have any downstream tags, test by ffw1
		if ($dsext==0) {
			my $p1=$mature_beg-$FLANK; $p1=1 if $p1 < 1;
			my $p2=$mature_end+$MAX_SPACE+$MATURE_SIZE+$FLANK;
			my $seq=subseq($refId_seq,$refId,$p1,$p2,$strand);
			$p2=$p1+length($seq)-1;
			$dsext=1;
			my $fold={};
#			print "ffw1($seq,$fold,$mature_seqref)\n";
			ffw1($seq,$fold,$mature_seqref,'max_energy'=>$MAX_ENERGY,'min_pair'=>$MIN_PAIR,'size_diff'=>$MAX_SIZEDIFF,'max_bulge'=>$MAX_BULGE,'asymmetry'=>$MAX_ASY,'min_space'=>$MIN_SPACE,'max_space'=>$MAX_SPACE,'flank'=>$FLANK);
			if ($test) {print TMP "$refId:$p1:$p2:$strand\t","ffw1($seq,$mature_seqref)\t",$fold->{pass},"\t",$fold->{reason},"\n"}
			$fold->{location}="$refId:$p1:$p2:$strand";
			if ($fold->{pass}) {
				$dsfold=$fold;
			}
		}
		
		### upstream extension
		my $usext=0; # flag, 0 means no good precursor, has to try again
		my $usfold;
		# case 1: has a closely located upstream cutting site, test by ffw2 only
		if ($usext==0 && $i-1 >= 0 && defined $csites[$i-1]) {
			my $csite_front=$csites[$i-1];
			my $matureix_front=$csite_matureix{$csite_front};
			my $mature_front=$tags[$matureix_front];
			my $mature_front_Id=$mature_front->[0];
			my $mature_front_beg=$mature_front->[1];
			my $mature_front_end=$mature_front->[2];
			my $mature_front_seq=$tagId_seq{$mature_front_Id};
			my $mature_front_seqref=subseq($refId_seq,$refId,$mature_front_beg,$mature_front_end,$strand);
			my $dist=$mature_beg-$mature_front_end-1;
			if ($dist <= $MAX_SPACE && $dist >= $MIN_SPACE) {
				my $p1=$mature_front_beg-$FLANK; $p1=1 if $p1 < 1;
				my $p2=$mature_end+$FLANK;
				my $seq=subseq($refId_seq,$refId,$p1,$p2,$strand);
				$p2=$p1+length($seq)-1;
				$usext=1;
				my $fold={};
#				print "ffw2($seq,$fold,$mature_front_seqref,$mature_seqref)\n";
				ffw2($seq,$fold,$mature_seqref,$mature_front_seqref,'max_energy'=>$MAX_ENERGY,'min_pair'=>$MIN_PAIR,'size_diff'=>$MAX_SIZEDIFF,'max_bulge'=>$MAX_BULGE,'asymmetry'=>$MAX_ASY,'min_space'=>$MIN_SPACE,'max_space'=>$MAX_SPACE,'flank'=>$FLANK);
				if ($test) {print TMP "$refId:$p1:$p2:$strand\t","ffw2($seq,$mature_seqref,$mature_front_seqref)\t",$fold->{pass},"\t",$fold->{reason},"\n"}
				$fold->{location}="$refId:$p1:$p2:$strand";
				if ($fold->{pass}) {
					$usfold=$fold;
				}
			}
		}
		# case 2:  has upstream tags, may be mature parters, test by ffw2
		if ($usext==0) {
			# collect upstream tags
			my @indexes_before;
			for (my $j=$matureix-1; $j >=0; $j--) {
				my $tagj=$tags[$j];
				my $tagj_Id=$tagj->[0];
				my $tagj_beg=$tagj->[1];
				my $tagj_end=$tagj->[2];
				my $tagj_seq=$tagId_seq{$tagj_Id};
				my $tagj_count=$tagId_count{$tagj_Id};
				my $dist=$mature_beg-$tagj_end-1;
				if ($dist < $MIN_SPACE) {
					next;
				}
				elsif ($dist > $MAX_SPACE) {
					last;
				}
				else {
					push @indexes_before,$j;
				}
			}
			# sort tags by count
			if (@indexes_before) {
				@indexes_before=sort {$tagId_count{$tags[$b][0]} <=> $tagId_count{$tags[$a][0]}} @indexes_before;			
				my $j=$indexes_before[0];
				my $tagj=$tags[$j];
				my $tagj_Id=$tagj->[0];
				my $tagj_beg=$tagj->[1];
				my $tagj_end=$tagj->[2];
				my $tagj_seq=$tagId_seq{$tagj_Id};
				my $tagj_seqref=subseq($refId_seq,$refId,$tagj_beg,$tagj_end,$strand);
				my $tagj_count=$tagId_count{$tagj_Id};
				my $p1=$tagj_beg-$FLANK; $p1=1 if $p1 < 1;
				my $p2=$mature_end+$FLANK;
				my $seq=subseq($refId_seq,$refId,$p1,$p2,$strand);
				$p2=$p1+length($seq)-1;
				my $fold={};
#				print "ffw2($seq,$fold,$mature_seqref,$tagj_seqref)\n";
				ffw2($seq,$fold,$mature_seqref,$tagj_seqref,'max_energy'=>$MAX_ENERGY,'min_pair'=>$MIN_PAIR,'size_diff'=>$MAX_SIZEDIFF,'max_bulge'=>$MAX_BULGE,'asymmetry'=>$MAX_ASY,'min_space'=>$MIN_SPACE,'max_space'=>$MAX_SPACE,'flank'=>$FLANK);
				if ($test) {print TMP "$refId:$p1:$p2:$strand\t","ffw2($seq,$mature_seqref,$tagj_seqref)\t",$fold->{pass},"\t",$fold->{reason},"\n"}
				$fold->{location}="$refId:$p1:$p2:$strand";
				if ($fold->{pass}) {
					$usfold=$fold;
					$usext=1;
				}
			}
		}
		# case 3: case 2 failed or do not have upstream tags, test by ffw1
		if ($usext==0) {
			my $p1=$mature_beg-$MAX_SPACE-$MATURE_SIZE-$FLANK; $p1=1 if $p1 < 1;
			my $p2=$mature_end+$FLANK;
			my $seq=subseq($refId_seq,$refId,$p1,$p2,$strand);
			$p2=$p1+length($seq)-1;
			$usext=1;
			my $fold={};
#			print "ffw1($seq,$fold,$mature_seqref)\n";
			ffw1($seq,$fold,$mature_seqref,'max_energy'=>$MAX_ENERGY,'min_pair'=>$MIN_PAIR,'size_diff'=>$MAX_SIZEDIFF,'max_bulge'=>$MAX_BULGE,'asymmetry'=>$MAX_ASY,'min_space'=>$MIN_SPACE,'max_space'=>$MAX_SPACE,'flank'=>$FLANK);
			if ($test) {print TMP "$refId:$p1:$p2:$strand\t","ffw2($seq,$mature_seqref)\t",$fold->{pass},"\t",$fold->{reason},"\n"}
			$fold->{location}="$refId:$p1:$p2:$strand";
			if ($fold->{pass}) {
				$usfold=$fold;
			}
		}
		
		### transform coors
		if ($dsfold) {
			my @locates=split(/:/,$dsfold->{location});
			my $beg_local=$dsfold->{beg};
			my $end_local=$dsfold->{end};
			my ($beg,$end);
			if ($locates[3] eq "+") {
				$beg=$locates[1]+$beg_local-1;
				$end=$locates[1]+$end_local-1;
			}
			else {
				$beg=$locates[2]-$end_local+1;
				$end=$locates[2]-$beg_local+1;
			}
			$dsfold->{refId}=$locates[0];
			$dsfold->{strand}=$locates[3];
			$dsfold->{beg}=$beg;
			$dsfold->{end}=$end;
			$dsfold->{location}=join(":",$locates[0],$beg,$end,$locates[3]);
		}
		if ($usfold) {
			my @locates=split(/:/,$usfold->{location});
			my $beg_local=$usfold->{beg};
			my $end_local=$usfold->{end};
			my ($beg,$end);
			if ($locates[3] eq "+") {
				$beg=$locates[1]+$beg_local-1;
				$end=$locates[1]+$end_local-1;
			}
			else {
				$beg=$locates[2]-$end_local+1;
				$end=$locates[2]-$beg_local+1;
			}
			$usfold->{refId}=$locates[0];
			$usfold->{strand}=$locates[3];
			$usfold->{beg}=$beg;
			$usfold->{end}=$end;
			$usfold->{location}=join(":",$locates[0],$beg,$end,$locates[3]);
		}
		
		next unless ($dsfold->{pass} || $usfold->{pass});
		
		### collect tags
		foreach my $fold ($dsfold,$usfold) {
			next unless $fold->{pass};
			my $beg=$fold->{beg};
			my $end=$fold->{end};
			my $strand=$fold->{strand};
			my @indexes;
			for (my $i=$matureix; $i<=$#tags; $i++) {
				my $tag=$tags[$i];
				if ($tag->[1] >= $beg && $tag->[2] <= $end ) {
					push @indexes,$i;
				}
				else {
					last;
				}
			}
			for (my $i=$matureix-1; $i>=0; $i--) {
				my $tag=$tags[$i];
				if ($tag->[1] >= $beg && $tag->[2] <= $end ) {
					push @indexes,$i;
				}
				else {
					last;
				}
			}
			@indexes=sort {$a <=> $b} @indexes;
			foreach my $i (@indexes) {
				my $tag=$tags[$i];
				my $tagId=$tag->[0];
				my $tag_seq=$tagId_seq{$tagId};
				my $tag_count=$tagId_count{$tagId};
				my ($beg_local,$end_local);
				# transform coors
				if ($tag->[3] eq "+") {
					$beg_local=$tag->[1]-$beg+1;
					$end_local=$tag->[2]-$beg+1;
				}
				else {
					$beg_local=$end-$tag->[2]+1;
					$end_local=$end-$tag->[1]+1;
				}
				my $tmp=[$tagId,$beg_local,$end_local,$tag->[3],$tag_count,$tag_seq];
				push @{$fold->{tags}},$tmp;
			}
		}
		
		### Deep testing
		# Majority of Small RNAs mapped to miRNA precursor should be products of 
		# sequential processing of Dicer and Drosha with rare ones may be random
		# degradation species. The pattern small RNAs distribute along the 
		# precursor of known miRNA genes agree this obervation well, indicating
		# we can use is it as a rule to exploit miRNAs from depp sequencing
		# datasets.
		if ($dsfold->{pass}) {
			dptest($dsfold);
		}
		if ($usfold->{pass}) {
			dptest($usfold);
		}

		### make decision
		my $fold_pf;
		# dsfold and usfold both pass filter, choose the best one, need further test
		if ($dsfold->{pass} && $usfold->{pass}) {
			if ($dsfold->{method} eq "ffw2" || $usfold->{method} eq "ffw2") {
				if ($dsfold->{method} eq "ffw2" && $usfold->{method} eq "ffw2") {
					if ($dsfold->{mfe} < $usfold->{mfe}) {
						$fold_pf=$dsfold;
					}
					else {
						$fold_pf=$usfold;
					}
				}
				elsif ($dsfold->{method} eq "ffw2") {
					$fold_pf=$dsfold;
				}
				elsif ($usfold->{method} eq "ffw2") {
					$fold_pf=$usfold;
				}
			}
			else {
				if ($dsfold->{mfe} < $usfold->{mfe}) {
					$fold_pf=$dsfold;
				}
				else {
					$fold_pf=$usfold;
				}
			}
		}
		# only dsfold pass filter
		elsif ($dsfold->{pass}) {
			$fold_pf=$dsfold;
		}
		# only usfold pass filter
		elsif ($usfold->{pass}) {
			$fold_pf=$usfold;
		}

		### write result
		if ($fold_pf->{pass}) {
			++$c;
#			print Dumper($fold_pf),"\n";
			my $seqId=sprintf "%.4d",$c;
			$seqId=$label."-m".$seqId;
			$fold_pf->{seqId}=$seqId;
			my $aln=mkaln($fold_pf);
			my $gff=mkgff($fold_pf);
			print ALN $aln;
			print GFF $gff;
		}
		### change site status 081029
		if ($fold_pf->{pass}) {
			$csite_pass{$csite}=1;
			my $beg=$fold_pf->{beg};
			my $end=$fold_pf->{end};
			for (my $j=$i+1; $j<=$#csites; $j++) {
				my $csitej=$csites[$j];
				if ($csitej < $end && $csitej > $beg) {
					$csite_pass{$csitej}=1;
				}
				elsif ($csitej > $end) {
					last;
				}
			}
			for (my $j=$i-1; $j >=0; $j--) {
				my $csitej=$csites[$j];
				if ($csitej < $end && $csitej > $beg) {
					$csite_pass{$csitej}=1;
				}
				elsif ($csitej < $beg) {
				    last;
				}
			}
		}
	}
	my $csites_pass_n=scalar keys %csite_pass;
	print LOG "\t$refId:$strand $break_n $csites_n($csites_repeat_n) $csites_pass_n\n";
}

my $end_time=localtime;
print LOG "\n$prog_name finish at $end_time\n";
print LOG "  $c miRNA genes are identified\n";
close GFF;
close ALN;
close LOG;

# test the autheticity of a miRNA gene by deep sequencing data
sub dptest {
	my $fold=shift;
	my $seq=$fold->{seq};
	my $m5=$fold->{m5};
	my $m3=$fold->{m3};
	my $m5_beg=index($seq,$m5,0)+1;
	my $m5_end=$m5_beg+length($m5)-1;
	my $m3_beg=index($seq,$m3,0)+1;
	my $m3_end=$m3_beg+length($m3)-1;
	my @tags=@{$fold->{tags}};

	my $count_m5=0;
	my $count_5p=0;
	my $count_m3=0;
	my $count_3p=0;
	my $count_total=0;
	my $count_loop=0;
	my $count_other=0;
	foreach my $tag (@tags) {
		my $tagId=$tag->[0];
		my $tag_beg=$tag->[1];
		my $tag_end=$tag->[2];
		my $tag_count=$tag->[4];
		$count_total+=$tag_count;
		# classify srnas into four categories
		# 1. 5' matures
		# 2. 3' matures
		# 3. srnas on loop
		# 4. chimeric srnas
		# the first three categories are products of Dicer and Drosha processing, while the fourth not
		if (abs($tag_beg-$m5_beg) <= 3) {
			$count_m5+=$tag_count;
		}
		elsif (abs($tag_beg-$m3_beg) <= 3) {
			$count_m3+=$tag_count;
		}
		elsif ($tag_beg >= $m5_end-3 && $tag_end <= $m3_beg+3 ) {
			$count_loop+=$tag_count;
		}
		else {
			$count_other+=$tag_count;
		}

		if ($tag_end <= $m5_end+5) {
			$count_5p+=$tag_count;
		}
		elsif ($tag_beg >= $m3_beg-5) {
			$count_3p+=$tag_count;
		}
	}
	$fold->{count}=$count_total;
	$fold->{m5count}=$count_m5;
	$fold->{m3count}=$count_m3;
	# >95% srnas on 5p arm should be mature
	if ($count_5p-$count_m5 > 1) {
		$fold->{pass}=0 if ($count_5p-$count_m5 >= int($count_5p*0.05));
	}
	# >95% srnas on 3p arm should be mature
	if ($count_3p-$count_m3 > 1) {
		$fold->{pass}=0 if ($count_3p-$count_m3 >= int($count_3p*0.05));
	}
	# chimeric srnas should less than 1%
	if ($count_other > 1) {
		$fold->{pass}=0 if ($count_other >= int($count_total*0.01));
	}
}

# produce gff format output
sub mkgff {
	my $fold=shift;
	my $gff;
	
	my $seq=$fold->{seq};
	my $seqId=$fold->{seqId};
	my $seq_len=length $fold->{seq};
	my $m5=$fold->{m5};
	my $m3=$fold->{m3};
	my $m5Id=$seqId."-5p";
	my $m3Id=$seqId."-3p";
	my $m5_beg=index($seq,$m5,0)+1;
	my $m5_end=$m5_beg+length($m5)-1;
	my $m3_beg=index($seq,$m3,0)+1;
	my $m3_end=$m3_beg+length($m3)-1;
	
	# transform coors
	my ($refId,$beg,$end,$strand)=split(/:/,$fold->{location});
	my ($m5_beg_abs,$m5_end_abs,$m3_beg_abs,$m3_end_abs);
	if ($strand eq "+") {
		$m5_beg_abs=$beg+$m5_beg-1;
		$m5_end_abs=$beg+$m5_end-1;
		$m3_beg_abs=$beg+$m3_beg-1;
		$m3_end_abs=$beg+$m3_end-1;
	}
	else {
		$m5_beg_abs=$end-$m5_end+1;
		$m5_end_abs=$end-$m5_beg+1;
		$m3_beg_abs=$end-$m3_end+1;
		$m3_end_abs=$end-$m3_beg+1;
	}
	
	$gff .= "$refId\t$prog_name\tprecursor\t$beg\t$end\t.\t$strand\t.\tID=$seqId;Count=$fold->{count};mfe=$fold->{mfe}\n";
	if ($fold->{m5count} > 0) {
		$gff .= "$refId\t$prog_name\tmature-5p\t$m5_beg_abs\t$m5_end_abs\t.\t$strand\t.\tID=$m5Id;Parent=$seqId;Count=$fold->{m5count};Seq=$m5\n";
	}
	if ($fold->{m3count} > 0) {
		$gff .= "$refId\t$prog_name\tmature-3p\t$m3_beg_abs\t$m3_end_abs\t.\t$strand\t.\tID=$m3Id;Parent=$seqId;Count=$fold->{m3count};Seq=$m3\n";
	}

	return $gff;
}

sub mkaln {
	my $fold=shift;
	my $aln;

	my $seq=$fold->{seq};
	my $seqId=$fold->{seqId};
	my $seq_len=length $fold->{seq};
	my $m5=$fold->{m5};
	my $m3=$fold->{m3};
	my $m5_str=newstr($m5,$seq);
	my $m3_str=newstr($m3,$seq);
	my $m5Id=$seqId."-5p";
	my $m3Id=$seqId."-3p";

	$aln="$prog_name\n";
	$aln.="$fold->{seqId} $fold->{location} $seq_len(nt) $fold->{mfe}(kcal/mol)\n";
	$aln.="$fold->{seq} $fold->{seqId} $fold->{count}\n";
	$aln.="$fold->{struct}\n";
#	$aln.="$m5_str $m5Id $fold->{m5count}\n";
#	$aln.="$m3_str $m3Id $fold->{m3count}\n";
	$aln.="$m5_str $m5Id $fold->{m5count}\n" if ($fold->{m5count} > 0);
	$aln.="$m3_str $m3Id $fold->{m3count}\n" if ($fold->{m3count} > 0);
	foreach my $tag (sort {$a->[1] <=> $b->[1] or $a->[2] <=> $b->[2]} @{$fold->{tags}}) {
		my $tagId=$tag->[0];
		my $tag_beg=$tag->[1];
		my $tag_end=$tag->[2];
		my $tag_count=$tag->[4];
		my $tag_seq=$tag->[5];
		my $str="-" x $seq_len;
		substr($str,$tag_beg-1,$tag_end-$tag_beg+1)=$tag_seq;
		$aln.="$str $tagId $tag_count\n";
	}
	$aln.="//\n";
	return $aln;
}

sub readfasta {
	my $infile=shift;
	my $seqId_seq=shift;
	
	my $c=0;
	open IN, $infile || die $!;
	my $seqId;
	while (<IN>) {
		if (/^>(\S+)/) {
			$seqId=$1;
			$c++;
		}
		else {
			$_=~s/\s//g;
			$seqId_seq->{$seqId}.=$_;
		}
	}
	close IN;
	return $c;
}

sub subseq{
	my $seqId_seq=shift;
	my $seqId=shift;
	my $beg=shift;
	my $end=shift;
	my $strand=shift;
	
	my $subseq=substr($seqId_seq->{$seqId},$beg-1,$end-$beg+1);
	if ($strand eq "-") {
		$subseq=revcom($subseq);
	}
	return uc $subseq;
}

sub revcom{
	my $seq=shift;
	$seq=~tr/ATCGatcg/TAGCtagc/;
	$seq=reverse $seq;
	return uc $seq;
}


sub newstr{
	my $substr=shift;
	my $str=shift;
	my $start=index($str, $substr, 0);
	warn "[error: newstr] no substr in str\n" if ($start < 0);

	my $strlen=length $str;
	my $substrlen=length $substr;
	my $newstr= '*' x $strlen;
	substr($newstr,$start,$substrlen)=$substr;
	return $newstr;
}

sub overlen {
	my ($p1,$p2,$t1,$t2)=@_;
		if ($p1>$p2) {
			my $tmp=$p1;
				$p1=$p2;
				$p2=$tmp;
		}
		if ($t1>$t2) {
				my $tmp=$t1;
				$t1=$t2;
				$t2=$tmp;
		}
		my $m = $p1<$t1 ? $t1 : $p1;
		my $n = $p2<$t2 ? $p2 : $t2;
		if ($m > $n) {
			return 0;
		}
		else {
			return $n-$m+1;
		}
}

sub usage {
	my $usage = << "USAGE";
Program: MIREAP (Discover miRNAs from deeply sequenced smRNA library)
Version: $version
Contact: Li Qibin <liqb\@genomics.org.cn>

Usage: mireap.pl -i <smrna.fa> -m <map.txt> -r <reference.fa> -o <outdir>
  Options:
  -i <file>  Small RNA library, fasta format, forced
  -m <file>  Mapping file, tabular format, forced
  -r <file>  Reference file, fasta format, forced
  -o <dir>   Directory where results produce (current directory)
  -t <str>   Sample label (xxx)
  -A <int>   Minimal miRNA sequence length (18)
  -B <int>   Maximal miRNA sequence length (26)
  -a <int>   Minimal miRNA reference sequence length (20)
  -b <int>   Maximal miRNA reference sequence length (24)
  -u <int>   Maximal copy number of miRNAs on reference (20)
  -e <folat> Maximal free energy allowed for a miRNA precursor (-18 kcal/mol)
  -d <int>   Maximal space between miRNA and miRNA* (35)
  -p <int>   Minimal base pairs of miRNA and miRNA*
  -v <int>   Maximal bulge of miRNA and miRNA* (4)
  -s <int>   Maximal asymmetry of miRNA/miRNA* duplex
  -f <int>   Flank sequence length of miRNA precursor (10)
  -h         Help

USAGE
	print $usage;
}

