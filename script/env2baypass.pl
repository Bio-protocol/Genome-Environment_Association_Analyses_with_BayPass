## This script transforms environmental dataset to BayPass covariate data file
## run as: cat {env_data} | perl env2baypass.pl {pop} {prefix}
## Kaichi Huang 2022 May

use warnings;
use strict;

my $out_prefix = $ARGV[1];
my $first = 1;
my @all_pop;
my @all_cov;
my %all_info; # values of all covariates of all populations

# Read in population order from genetic data
open POP, $ARGV[0];
while (<POP>){
	chomp;
	push @all_pop, $_;
}
close POP;

# Store all info
while (<STDIN>) {
	chomp;
	my @a = split(/\t/, $_);
	if ($first) {
		$first = 0;
		foreach my $i (1..$#a) {
			push @all_cov, $a[$i];
			$all_info{$a[$i]} = {};
		}
		next;
	}
	my $pop = $a[0];
	foreach my $i (1..$#a) {
		$all_info{$all_cov[$i-1]}{$pop} = $a[$i];
	}
}

# Output
open OUT, ">$out_prefix.txt";
open OUT_COV, ">$out_prefix.cov";
foreach my $cov (@all_cov) {
	my @values = values $all_info{$cov};
	print OUT_COV "$cov\n";
	my @line;
	foreach my $pp (@all_pop) {
		push @line, "$all_info{$cov}{$pp}";
	}
	print OUT join(" ",@line),"\n";
}
close OUT;
close OUT_COV;
