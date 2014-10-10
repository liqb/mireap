package Struct;
use strict;

use Exporter;
our @ISA=qw(Exporter);
our @EXPORT=qw(parse_struct which_arm get_star biggest_bulge get_asy);

=head1 NAME

Struct - methods used to parse RNA structure

=head1 SYNOPSIS

  use Struct;
  $struct="";
  %pair;
  parse_struct($struct,\%pair);
  
  $substruct="";
  $arm=which_arm($substruct); # 5p/3p/-

  ($star_beg,$star_end)=get_star(\%pair,$beg,$end);

  ($biggest_bulge)=biggest_bulge($substruct);
  ($asymmetry)=get_asy(\%pair,$beg,$end);

=head1 AUTHOR
liqb <liqb@genomics.org.cn>

=cut


# build base pair table, coors count from 1
sub parse_struct {
	my $struct=shift;
	my $table=shift;

	my @t=split('',$struct);
	my @lbs; # left brackets
	foreach my $k (0..$#t) {
		if ($t[$k] eq "(") {
			push @lbs, $k+1;
		}
		elsif ($t[$k] eq ")") {
			my $lb=pop @lbs;
			my $rb=$k+1;
			$table->{$lb}=$rb;
			$table->{$rb}=$lb;
		}
	}
	if (@lbs) {
		warn "unbalanced RNA struct.\n";
	}
}

# define which arm tag reside, return (5p/3p/-)
sub which_arm {
	my $substruct=shift;
	my $arm;
	if ($substruct=~/\(/ && $substruct=~/\)/) {
		$arm="-";
	}
	elsif ($substruct=~/\(/) {
		$arm="5p";
	}
	else {
		$arm="3p";
	}
	return $arm;
}

# given a sub-region, get its star region, 2 nt 3' overhang
sub get_star {
	my($table,$beg,$end)=@_;
	
	my ($s1,$e1,$s2,$e2); # s1 pair to s2, e1 pair to e2
	foreach my $i ($beg..$end) {
		if (defined $table->{$i}) {
			my $j=$table->{$i};
			$s1=$i;
			$s2=$j;
			last;
		}
	}
	foreach my $i (reverse ($beg..$end)) {
		if (defined $table->{$i}) {
			my $j=$table->{$i};
			$e1=$i;
			$e2=$j;
			last;
		}
	}
#	print "$s1,$e1 $s2,$e2\n";
	
	# correct terminus
	my $off1=$s1-$beg;
	my $off2=$end-$e1;
	$s2+=$off1;
	$s2+=2; # 081009
	$e2-=$off2; $e2=1 if $e2 < 1;
	$e2+=2; $e2=1 if $e2 < 1; # 081009
	($s2,$e2)=($e2,$s2) if ($s2 > $e2);
	return ($s2,$e2);
}

# compute biggest bulge size
sub biggest_bulge {
	my $struct=shift;
	my $bulge_size=0;
	my $max_bulge=0;
	while ($struct=~/(\.+)/g) {
		$bulge_size=length $1;
		if ($bulge_size > $max_bulge) {
			$max_bulge=$bulge_size;
		}
	}
	return $max_bulge;
}

# compute asymmetry
sub get_asy {
	my($table,$a1,$a2)=@_;
	my ($pre_i,$pre_j);
	my $asymmetry=0;
	foreach my $i ($a1..$a2) {
		if (defined $table->{$i}) {
			my $j=$table->{$i};
			if (defined $pre_i && defined $pre_j) {
				my $diff=($i-$pre_i)+($j-$pre_j);
				$asymmetry += abs($diff);
			}
			$pre_i=$i;
			$pre_j=$j;
		}
	}
	return $asymmetry;
}

1;


