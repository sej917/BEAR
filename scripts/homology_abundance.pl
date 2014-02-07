#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

=head1 NAME

homology_abundance.pl

=head1 AUTHOR(S)

Brett Trost

=head1 SYNOPSIS

homology_abundance.pl [-b bit_score_threshold] [-r] search_file mapping_file fasta_file

=head1 DESCRIPTION

Takes as input a file containing search results, as well as a file mapping accession numbers to organism lineages,
and outputs two files - a detailed file and a summary file - giving the taxonomic assignment given to each query sequence.
The lineage for a given query is calculated by taking all the lineages for all of the search results hat are different by a certain bit score threshold
from the result with the highest bit score for that query, and finding the lowest common ancestor (LCA) among all of those lineages.  This is similar
to what is done in the MEGAN software package.

=head1 INPUT

search_file: a tab-delimited table with the following format:

query_identifier	hit_identifier	bit_score
		
The database used for searching must have been refseq_protein

mapping_file: a table in tab-delimited text format where the first column is a refseq accession number, and the second column is the lineage of the organism associated with that sequence.
fasta_file: the fasta file that was used as input to the search program.  This is needed because some search programs (i.e. rapsearch) don't put any mention of non-matching queries in its results file, so we need to add those to the results manually.

=head1 OUTPUT

<name_of_search_file>.lca_map: a file in tab-delimited text format where the first column is the identifier of a read in the RAPSearch2 results file, and the second column is the lineage assigned to that query
<name_of_search_file>.lca_summary: a file in tab-delimited text format where the second column is a lineage, the first column is the number queries assigned to that lineage, and the third column is a list of the queries assigned to that lineage
<name_of_search_file>.lca_summary.normalized: a file in tab-delimited text format where the second column is a classification, and the first column is the frequency of that classification

=head1 ASSUMPTIONS
Assumes that refseq_protein was used as the database when searching.

=head1 DATA

=head1 EXAMPLE
If using a file from rapsearch, you need to delete the first five lines and then extract the first, second, and 12th columns; e.g.

tail -n +6  alpha_core1_02cm_OT1.clipped.trimmedQ17.fasta.rapsearch.m8 | cut -f 1,2,12 >  alpha_core1_02cm_OT1.clipped.trimmedQ17.fasta.rapsearch.m8.cut

If using a file from summarize_blast.pl, you need to delete the first line and then extract the first, fifth, and 13th columns.  You also need to delete all the records
that didn't have significant hits.  E.g.:


=head1 SEE ALSO
make_organism_to_accession_map.pl: script used to create the organism to accession map.

=cut

my $bit_score_parameter = 10;
my $print_reads = 0;
GetOptions("-b=i" => \$bit_score_parameter, "-r" => \$print_reads);

sub get_map($) {
	my $mapping_file = shift;

	my %map;
	open MAPPING_FILE, $ARGV[1] or die "Could not open mapping file $ARGV[1]: $!";
	while (defined(my $line = <MAPPING_FILE>)) {
		chomp($line);
		my ($accession, $lineage) = split "\t", $line;
		$map{$accession} = $lineage;
	}
	return %map;
}

my %lca_summary;
my %accession_to_lineage_map = get_map($ARGV[1]);
my $highest_bit_score;
my ($current_query_sequence, $hit_sequence, $bit_score);
my @this_set;
my $previous_query_sequence;
my %queries_found;


sub lowest_common_ancestor(@) {
	my @lineage_strings = @_;
	
	my @lineages;
	
	foreach my $lineage_string (@lineage_strings) {
		push @lineages, [ split /;\s*/, $lineage_string ];
	}
	
	my $lowest_common_ancestor = "Root";
	
	for my $i (0..$#{$lineages[0]}) {
		for my $j (0..($#lineages - 1)) {
			
			
			#if (!$lineages[$j][$i] or !$lineages[$j+1][$i]) {
			#	print "Current query sequence is $current_query_sequence and First is ($lineages[$j][$i]) and second is ($lineages[$j+1][$i])\n";
			#	print "Something's amiss!  Lineage strings are:\n", join("\n", @lineage_strings), "\n";
			#	die;
			#}
			if (!$lineages[$j][$i] or !$lineages[$j+1][$i] or $lineages[$j][$i] ne $lineages[$j+1][$i]) {
				return $lowest_common_ancestor;
			}
		}
		
		$lowest_common_ancestor .= "; $lineages[0][$i]";
	}
	
	return $lowest_common_ancestor;
}


sub process_query() {
	my @lineage_strings;
	my @eligible_sets = grep { $_->{BIT_SCORE} >= $highest_bit_score - $bit_score_parameter } @this_set;
		
	foreach my $set (@eligible_sets) {
		warn "Missing mapping entry for hit accession $set->{HIT_ACCESSION}\n" if !defined($accession_to_lineage_map{$set->{HIT_ACCESSION}});
		push @lineage_strings, $accession_to_lineage_map{$set->{HIT_ACCESSION}};
	}
	
	my $lca = lowest_common_ancestor(@lineage_strings);
	push @{$lca_summary{$lca}}, $previous_query_sequence;
	print LCA_MAP "$previous_query_sequence\t$lca\n";
			
	$highest_bit_score = $bit_score;
	
}

###########################################################################################
###########################################################################################



my %lowest_common_ancestors;

open LCA_MAP, ">$ARGV[0].lca_map";
open LCA_SUMMARY, ">$ARGV[0].lca_summary";
open LCA_SUMMARY_NORMALIZED, ">$ARGV[0].lca_summary.normalized";

open SEARCH_FILE, $ARGV[0] or die "Could not open search results file $ARGV[0]: $!";

# Process the first line
chomp(my $line = <SEARCH_FILE>);
($current_query_sequence, $hit_sequence, $bit_score) = split "\t", $line;
$queries_found{$current_query_sequence} = 1;
my $hit_accession;
if ($hit_sequence =~ /ref\|(\w+)/) {	# If the whole identifier is given (not just the accession number), then extract it...
	$hit_accession = $1;
} else {
	$hit_accession = $hit_sequence;	# Otherwise, assume that the identifier consists entirely of the accession number
}

@this_set = ({ HIT_ACCESSION => $hit_accession, BIT_SCORE => $bit_score });
$previous_query_sequence = $current_query_sequence;
$highest_bit_score = $bit_score;

# Process the remaining lines
while (defined($line = <SEARCH_FILE>)) {
	chomp($line);
	($current_query_sequence, $hit_sequence, $bit_score) = split "\t", $line;
	$queries_found{$current_query_sequence} = 1;
	if ($hit_sequence =~ /ref\|(\w+)/) {	# If the whole identifier is given (not just the accession number), then extract it...
		$hit_accession = $1;
	} else {
		$hit_accession = $hit_sequence;	# Otherwise, assume that the identifier consists entirely of the accession number
	}

	if ($current_query_sequence ne $previous_query_sequence) { # Starting with a new set of query sequences, so process the old ones
		process_query();
		@this_set = ({ HIT_ACCESSION => $hit_accession, BIT_SCORE => $bit_score });
	} else {
		push @this_set, { HIT_ACCESSION => $hit_accession, BIT_SCORE => $bit_score };
	}
	$previous_query_sequence = $current_query_sequence;	
}

process_query();

my $total_length = `grep ">" $ARGV[2] | wc -l | tr -d " \\n"`;

open FASTA_FILE, $ARGV[2] or die "Could not open FASTA file: $!";
while (defined(my $line = <FASTA_FILE>)) {
	next if $line !~ /^>/;
	(my $accession) = $line =~ />(\S+)/;
	if (!$queries_found{$accession}) {
		print LCA_MAP "$accession\tNomatch\n";
		push @{$lca_summary{Nomatch}}, $previous_query_sequence;
	}
}
close FASTA_FILE;

# Now print out the summary
foreach my $lca (sort { @{$lca_summary{$b}} <=> @{$lca_summary{$a}} } keys %lca_summary) {
	my $num = @{$lca_summary{$lca}};
	my $frequency = $num / $total_length;
	if ($print_reads) {	
		print LCA_SUMMARY "$num\t$lca\t", join(", ", @{$lca_summary{$lca}}), "\n";
	} else {
		print LCA_SUMMARY "$num\t$lca\n";
	}
	print LCA_SUMMARY_NORMALIZED "$frequency\t$lca\n";
}

close SEARCH_FILE;
close LCA_MAP;
close LCA_SUMMARY;
close LCA_SUMMARY_NORMALIZED;

