#!/usr/bin/perl

=head1 name

normalize_abundances.pl

=head1 AUTHOR(s)

Stephen Johnson

=head1 SYNOPSIS

./normalize_abundances.pl <lca_summary_file> <genomes_fasta_file>

=head1 DESCRIPTION

Takes as inputa file in tab-delimited text format where the second column is a lineage and the first column is the number queries assigned to that lineage, as well as a multi-FASTA file of genome sequences

=head1 INPUT

lca_summary_file: a tab-delimited file with the following format:

number_of_hits	organism_linege

genomes_fasta_file: A fasta file composed of at least one genomic sequence

=head1 OUTPUT

To stdout, tab-delimited data of the following format:

fasta_header	relative_abundance_percentage


=cut

open(MYFILE, $ARGV[0]);
open(FASTA, $ARGV[1]);
my $sum = 0;
my @counts;
my @headers;
while(<MYFILE>){
	if($_ =~  /(\d+)([^;]+;){7,}[^;]*/){ 
		push(@counts, $1);	   
	}	
}
close(MYFILE);
my $i=0;
my $abund=0;
while(<FASTA>){
	if($_ =~ /^>.*/){
		chomp;
		push(@headers, $_);
	}
}
my $i=0;
for (@headers){
	$sum += $counts[$i];
	$i++;
}
$i=0;
for (@headers){
	$abund = $counts[$i]/$sum;
	print "$_\t$abund\n";
	$i++;
}
