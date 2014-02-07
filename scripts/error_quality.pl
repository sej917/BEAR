#!/usr/bin/perl


#This script takes the UCLUST summary output from DRISEE and determines the quality scores of erroneous nucleotides
#for each position within the read

use strict;
use warnings;
use Bio::SeqIO;


sub analyze_cluster{
	my %err_hash;
	my %sub_matr;
	my %ins_matr;

	my $length = $_[3];

	my @bases=('A', 'T', 'G', 'C');
	for my $base1 (@bases){
		for my $base2 (@bases){
			$sub_matr{$base1}{$base2} = 0;
			$ins_matr{$base1}{$base2} = 0;
		}
	}

	#for (my $i = 0; $i < $length; $i++){
	for my $i(0..$length){
		$err_hash{'D'}{$i} = 0;
		$err_hash{'I'}{$i} = 0;
	#	print "$i  $length\n"; 
		@{$err_hash{'A'}{$i}} =  @{$err_hash{'T'}{$i}} =  @{$err_hash{'G'}{$i}} =  @{$err_hash{'C'}{$i}} =  @{$err_hash{'X'}{$i}} = (0);
	}	
	my @cluster = split /\n/, $_[0];
	my $seqRef = $_[1];
	my $qualRef = $_[2];
	my $s_id = (split /\t/, $cluster[0])[8];
	shift @cluster;
	foreach(@cluster){
		my( $h_ident, $h_align, $h_id ) = (split /\t/, $_)[3, 7, 8];
		if($h_ident < 100.0){ #only concerned with erroneous nucletotides
			my $h_pos = my $s_pos = 0;
			my $h_tmp = my $s_tmp = "";
			my @h_tmpqual = my @s_tmpqual = (0);
			my @h_qual = split(/ /, $qualRef->{$h_id});
			my @s_qual = split(/ /, $qualRef->{$s_id});
			$h_align =~ s/([MID])/$1 /g;
			my @align_array = split(/ /, $h_align);
			foreach my $sec (@align_array){
				if( $sec =~ /^(\d*)([MID])/){
					my $num = ($1 ne "") ? int($1) : 1;
        				my $type = $2;
					if($type eq "M"){
						#print "$seqRef->{$h_id}, $h_pos, $num";
						my $h_subseq = substr($seqRef->{$h_id}, $h_pos, $num);
						my $s_subseq = substr($seqRef->{$s_id}, $s_pos, $num);
						$h_tmp .= $h_subseq;
						$s_tmp .= $s_subseq;
						
						my @h_subqual = @h_qual[$h_pos..$h_pos+$num];
						my @s_subqual = @s_qual[$s_pos..$s_pos+$num];
						push(@h_tmpqual, @h_subqual);
						push(@s_tmpqual, @s_subqual);
						
						for(my $i = 0; $i < $num; $i++){ 
							my $h_char = substr($h_subseq, $i, 1);
							my $s_char = substr($s_subseq, $i, 1);
							if( $h_char ne $s_char ){
								my $subqual = $h_subqual[$i] <= $s_subqual[$i] ? $h_subqual[$i] : $s_subqual[$i];
								my $cor_nuc = $h_subqual[$i] <= $s_subqual[$i] ? $s_char : $h_char;
								my $bad_nuc = $h_subqual[$i] <= $s_subqual[$i] ? $h_char : $s_char;
								my $subpos  = $i + $s_pos;	
								push(@{$err_hash{$bad_nuc}{$subpos}}, $subqual);
								$sub_matr{$cor_nuc}{$bad_nuc} += 1; 
							}
						}
						
						$h_pos += $num;
						$s_pos += $num;
					}elsif($type eq "D"){
						my $s_subseq = substr($seqRef->{$s_id}, $s_pos, $num);
						for( my $i = 0; $i < $num; $i++){
							my $d_char = substr($s_subseq, $i, 1);
							$err_hash{'D'}{$i+$s_pos} += 1;
						}
						$s_pos += $num;
					}elsif($type eq "I"){
						#if($num > 1){
						#	print "HELLOOO! $num\n";
						#}
						$s_tmp .= "X" x $num;
						my @h_subqual = @h_qual[$h_pos..$h_pos+$num];
						my @s_subqual = @s_qual[$s_pos..$s_pos+$num];
						#my @h_subseq = $h_pos > 0 ? substr($seqRef->{$h_id}, $h_pos-1, $num) : substr($seqRef->{$h_id}, $h_pos-1, $num);
						my @h_subseq = substr($seqRef->{$h_id}, $h_pos, $num);
						my $sublen = length(@h_subseq);
						for( my $i = 0; $i < $num; $i++){
							#if($h_pos > 0 && $i +1 < $sublen){
							#	$ins_matr{$h_subseq[$i]}{$h_subseq[$i+1]}+=1;
							#}
							my $ins_qual = $h_subqual[$i];
							push(@{$err_hash{'X'}{$i+$h_pos}}, $ins_qual);
							$err_hash{'I'}{$i+$h_pos} += 1;
						}
						my $prev_char = $h_pos > 0 ? substr($seqRef->{$h_id}, $h_pos-1, 1) : $bases[rand @bases]; #value for "prev char" should go between ? and :
						for my $char (@h_subseq){
							$ins_matr{$prev_char}{$char} +=1;
							$prev_char = $char;
						}
						$h_pos += $num;
					}
					
				}
			}
		}
	}
	return (\%err_hash, \%sub_matr,\%ins_matr);
}

my ($fastq_file, $summary_file) = @ARGV;
my $cur_line = 1;

#my $fastq_file = 'ecoli_test100k.fastq';
#my $summary_file = 'ecoli_test100k.fastq.uc';
my $quality_output = $fastq_file . ".err.qual";
my $matrix_output = $fastq_file . ".err.matr";
open(SUMMARY, $summary_file);

my @clusters;
my $lines = "";
my @seq_ids;
my $test_count = 0;
while(<SUMMARY>){
	if($_ !~ m/^[C|#]/i){
		$lines .= $_;
		my $seq_id = (split /\t/)[8];
		push(@seq_ids, $seq_id);
	}
}
my %id_hash = map { $_ => 1 } @seq_ids;
close(SUMMARY);

my $seq_in= Bio::SeqIO->new(
				-format => 'fastq',
				-file   => $fastq_file,
			);

my %seqs;
my %quals;
my $hash_counter = 0;
my $seq_counter = 0;
print "Loading sequences..\n";

my $max_length = 0;
while(my $seq = $seq_in->next_seq() ){
	$seq_counter += 1;
	if(exists($id_hash{$seq->id()})){
		$hash_counter += 1;
		$seqs{$seq->id()} = $seq->seq();
		$quals{$seq->id()} = $seq->qual_text();
		if (length($seq->seq()) > $max_length){
			$max_length = length($seq->seq());
		}
	}
	if($hash_counter % 4 == 0){
		print "$hash_counter sequences hashed, $seq_counter sequences processed.\r";
	}
}
print "\nLoading complete.\n";
system("rm -f $fastq_file.tmp");
@clusters = split(/\n(?=S)/, $lines);


my $numelements = @clusters;
print "Num: $numelements\n";
my $counter = 0;
my %big_hash;

my @bases = ('A', 'T', 'C', 'G', 'X','D','I');
foreach my $base(@bases){
	foreach my $i (0..$max_length){
		if($base eq 'D' || $base eq 'I'){
			$big_hash{$base}{$i} = 1;
		}else{
			@{$big_hash{$base}{$i}} = ();
		}
	}
}


my %sub_matr;
my %ins_matr;
my @nucl = ('A', 'T', 'G', 'C');
foreach my $base1 (@nucl){
	foreach my $base2 (@nucl){
		$sub_matr{$base1}{$base2} = 1;
        	$ins_matr{$base1}{$base2} = 1;
        }
}

my $small_hash_ref;
my $sub_matr_ref;
my $ins_matr_ref;

foreach my $cluster (@clusters){
	$counter += 1;
	($small_hash_ref, $sub_matr_ref, $ins_matr_ref) = &analyze_cluster($cluster, \%seqs, \%quals, $max_length);
	my %small_hash = %$small_hash_ref;
	my %sub_matr_small = %$sub_matr_ref;
	my %ins_matr_small = %$ins_matr_ref;
	foreach my $base (@bases){
		if (!($base eq 'D' || $base eq 'I')){
			foreach my $i (0..$max_length){
				foreach my $qual (@{$small_hash{$base}{$i}}){
					if($qual ne 0){
						push(@{$big_hash{$base}{$i}}, $qual);
					}
				}
			}
			if ($base ne 'X'){
				foreach my $base2(@nucl){
					#print "$base $base2 $sub_matr_small{$base}{$base2}\n";
					$sub_matr{$base}{$base2} += $sub_matr_small{$base}{$base2};
					$ins_matr{$base}{$base2} += $ins_matr_small{$base}{$base2};	
				}
			}
		}else{
			foreach my $i(0..$max_length){
				#print "$base $i $max_length $small_hash{$base}{$i}\n";
				$big_hash{$base}{$i} += $small_hash{$base}{$i};
			}
		}
	}
}

my %avg_quals;

my @nucl_x = ('A', 'T', 'G', 'C', 'X');

foreach my $base(@nucl_x){
	foreach my $i(0..$max_length){
		if(@{$big_hash{$base}{$i}}){
			my $total = 0;
			my $count = 0;
			foreach my $qual (@{$big_hash{$base}{$i}}){
				if ($qual <= 40 && $base ne 'X'){
					$total += $qual;
					$count += 1.0;
				}
				if($base eq 'X'){
					$total += $qual;
					$count += 1.0;
				}
			}
			if($total eq 0 || $count eq 0){
				$avg_quals{$base}{$i} = 0;
			}else{
                                $avg_quals{$base}{$i} = $total/$count;
			}
		}else{
			$avg_quals{$base}{$i} = 0;
		}
	}
}

##now average insertion and deletion rates x_i/(x_i+y_i)
foreach my $i(0..$max_length){
	my $sum = ($big_hash{'D'}{$i} + $big_hash{'I'}{$i} );
	$big_hash{'D'}{$i} /= $sum;
	$big_hash{'I'}{$i} /= $sum;
#	print "$i\t$big_hash{'I'}{$i}\t$big_hash{'D'}{$i}\n";
}


##now fix averages 
open(QUAL_FILE, ">", $quality_output) || die "Could not open $quality_output"; 
foreach my $i (0..$max_length){
	print QUAL_FILE "$i\t$avg_quals{'A'}{$i}\t$avg_quals{'T'}{$i}\t$avg_quals{'G'}{$i}\t$avg_quals{'C'}{$i}\t$avg_quals{'X'}{$i}\t$big_hash{'D'}{$i}\t$big_hash{'I'}{$i}\n";
}
close(QUAL_FILE);
my %ins_total;
my %sub_total;
$ins_total{'A'} = $ins_total{'T'} = $ins_total{'G'} = $ins_total{'C'} = 0;
$sub_total{'A'} = $sub_total{'T'} = $sub_total{'G'} = $sub_total{'C'} = 0;
foreach my $base(@nucl){
	foreach my $base2 (@nucl){
		$ins_total{$base} += $ins_matr{$base}{$base2};
		$sub_total{$base} += $sub_matr{$base}{$base2};
	}
}


open(MATR_FILE, ">", $matrix_output) || die "Could not open $matrix_output";
#ATGC x ATGC
print MATR_FILE "#Insertion matrix\n#format {A,T,G,C}x{A,T,G,C}\n";
foreach my $base1(@nucl){
	foreach my $base2(@nucl){
		my $perc = $ins_matr{$base1}{$base2} / $ins_total{$base1};
		print MATR_FILE "$perc\t";
	}
	print MATR_FILE "\n";
}
print MATR_FILE "\n";
#ATGC x ATGC
print MATR_FILE "#Substitution matrix\n#format {A,T,G,C}x{A,T,G,C}\n";
foreach my $base1(@nucl){
        foreach my $base2(@nucl){
                my $perc = $sub_matr{$base1}{$base2} / $sub_total{$base1};
                print MATR_FILE "$perc\t";
        }
        print MATR_FILE "\n";
}
close(MATR_FILE);
