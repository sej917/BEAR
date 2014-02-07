#!/usr/bin/perl

=head1 NAME

parametric_abundance.pl

=head1 AUTHOR(S)

Stephen Johnson

=head1 SYNOPSIS

./parametric_abundance.pl <genomes_file> <low|med|high>

=head1 DESCRIPTION

Takes as input a multi-FASTA file of genome sequences and a qualifier determining community complexity ("low", "medium", or "high" species complexity)
and generates a tab-delimited abundance profile based on the user-specific complexity and FASTA file

=head1 INPUT

genomes_file: A multi-fasta file of multiple genomes

low|med|high: A measure of species complexity, can be either "low", "med", or "high". E.g., low species complexity would mean there are few dominant organisms in the dataset

=head1 OUTPUT

To stdout, tab-delimited data of the following format:

fasta_header    relative_abundance_percentage

=cut

my $denom;
open(MYFILE, $ARGV[0]);
while(<MYFILE>){$denom++ if $_ =~ /^>.*/ ;}
close(MYFILE);
my $constant=0;
my $exponent=0;
my $complexity = $ARGV[1];

if( $complexity eq "low" ){
	$constant = 31.4034355;
	$exponent = -1.2868576;
}elsif ( $complexity eq "med" ){
	$constant = 21.2301963;
	$exponent = -1.0582662;
}elsif ( $complexity eq "high" ){
	$constant = 2.08348893;
	$exponent = -0.233125;
}else{
	print "Invalid complexity\n"; die;
}
my $abundance;
my $abundance_total=0;
my $species_num = 1;
open (MYFILE, $ARGV[0]);

for(my $i=1; $i <= $denom; $i++){
	$abundance_total += ($constant * (($i)**$exponent));
}

while( <MYFILE>){
	chomp;
	if ($_ =~ m/^>/){
		$abundance = ($constant * (($species_num)**$exponent));
		$abundance = $abundance / $abundance_total;
		print "$_\t$abundance\n";
		$species_num = $species_num + 1;
	}
}
