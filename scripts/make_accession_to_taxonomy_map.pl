#!/usr/bin/perl

use strict;
use warnings;

=head1 NAME

make_accession_to_taxonomy_map.pl

=head1 AUTHOR(S)

Brett Trost

=head1 SYNOPSIS

make_accession_to_taxonomy_map.pl genpept_format_file

=head1 DESCRIPTION

Takes as input a file of sequence information in genpept format, and parses each record
so to extract the accession number of that record, as well as the taxonomic lineage
of the organism to which that record corresponds.  

=head1 INPUT

genpept_format_file: file in genpept format

=head1 OUTPUT

Outputs to standard output a table in tab-delimited text format, where the first column
is the accession number of a record and the second column is the lineage of the organism to which
that record corresponds.

=head1 EXAMPLE

gunzip completeNZ_ABFL.protein.gpff.gz
make_accession_to_taxonomy_map.pl completeNZ_ABFL.protein.gpff

=head1 SEE ALSO

homology_abundance.pl

=cut

my $genpept_format_file = $ARGV[0];
open FILE, $genpept_format_file;
undef $/;
my $file_contents = <FILE>;
close FILE;

my @records = split "\n//\n", $file_contents;

foreach my $record (@records) {
	(my $accession) = $record =~ /ACCESSION\s+(\S+)/;
	
	my ($species, $lineage_raw) = $record =~ /ORGANISM\s+(.+?)\n\s*((?:Bacteria|Viruses|Archaea|Eukaryota).*?)\./si;
	
	$species =~ s/\n/ /s;
	
	if (!$accession or !$species or !$lineage_raw) {
		warn "ERROR FOR RECORD: $record\n";
		next;
	}
		
	my $completed_lineage = "$lineage_raw; $species";
	
	$completed_lineage =~ s/\n//g;
	$completed_lineage =~ s/\s{2,}/ /g;
	
	print "$accession\t$completed_lineage\n";
} 
