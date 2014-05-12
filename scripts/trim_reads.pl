#!/usr/bin/perl 

use warnings;
use strict;
use Bio::SeqIO;
use Getopt::Long;
use List::Util qw(max);
my $input="";
my $fasta="";
my $output="";
my $err_rate="";
my $err_qual="";
my $matrix_file="";
my $h_ins=0;
my $h_del=0;

my @lengths;
my @freq;

my @seeds; #first quality value of each sequence

GetOptions ("i=s" => \$input,
	    "f=s" => \$fasta,
	    "o=s" => \$output,
	    "r=s" => \$err_rate,
	    "q=s" => \$err_qual,
	    "m=s" => \$matrix_file,
	    "hi=f" => \$h_ins,
	    "hd=f" => \$h_del);
open(DRISEE_ERR_FILE, $err_rate) unless ($err_rate eq "0");
open(DRISEE_QUAL_FILE, $err_qual) unless ($err_qual eq "0");

my %rates;
my %quals;
if($err_rate ne "0"){
	open(DRISEE_ERR_FILE, $err_rate);
	system("error_models.py $err_rate > $err_rate.r.out");
	close(DRISEE_ERR_FILE);
	open(ERR_RATE_FILE, "$err_rate.r.out"); 
	while(<ERR_RATE_FILE>){
		my ($m, $b);
		my $next;
		chomp;
		if($_ =~ /(.*)([ATGCX])(\.)(nonzero)(.*)/){
			my $base = $2;
			$next = <ERR_RATE_FILE>;
			chomp $next;
			if ($next =~ /(\s*)(-?\d*\.*\d*)(\s*)(-?\d*\.*\d*)/){
				($b, $m) = ($2, $4);
				$rates{$base}{'slope'} = $m;
				$rates{$base}{'intercept'} = $b;
				print $rates{$base}{'slope'};
				print $rates{$base}{'intercept'};
			}
		}
	}
	close(ERR_RATE_FILE);
}

#assume insertion/deletion rates are equal unless otherwise specified
my %indel_rates;
for my $init (0..1200){
	$indel_rates{'D'}{$init} = 0.5;
	$indel_rates{'I'}{$init} = 1.0;
}
if($err_qual ne "0"){
	open(DRISEE_QUAL_FILE, $err_qual);
	my $i = 0;
	while(<DRISEE_QUAL_FILE>){
		($indel_rates{'D'}{$i}, $indel_rates{'I'}{$i}) = (split("\t", $_))[6, 7];
		$indel_rates{'I'}{$i} = $indel_rates{'D'}{$i} + $indel_rates{'I'}{$i}; #cumulative probability
		$i++;
	}
	system("error_quality.py $err_qual > $err_qual.r.out");
	close(DRISEE_QUAL_FILE);
	open(ERR_QUAL_FILE, "$err_qual.r.out");
	while(<ERR_QUAL_FILE>){
	        my ($inter, $coef1, $coef2);
	        my $next;
	        chomp;
	        if( $_ =~ /(.*)([ATGCX])(\.)(nonzero)(.*)/ ){
	        	my $base = $2;
	        	$next = <ERR_QUAL_FILE>;
	        	chomp $next;
	        	if($next =~ /(\s*)(-?\d*\.*\d*)(\s*)(-?\d*\.*\d*)(\s*)(-?\d*\.*\d*)/){
					($inter, $coef1, $coef2) = ($2, $4, $6);
					$quals{$base}{'intercept'} = $inter;
					$quals{$base}{'coef1'} = $coef1;
					$quals{$base}{'coef2'} = $coef2;	        		
	        	}
	        }
	}
	close(ERR_QUAL_FILE);
}

my @ins_values;
my @sub_values;
if($matrix_file ne "0"){
	open(MATRIX_FILE, $matrix_file) or die "Could not open $matrix_file\n";
	my $mat_count = 0; #should only read in 8 lines of input (4 lines insertion matrix, 4 lines substitution matrix
	while(<MATRIX_FILE>){
		chomp;
		if(($_ =~ /(\d).*/) && ($mat_count < 8)){
			my @line = split("\t", $_);
			foreach my $val (@line){
				$mat_count < 4 ? push(@ins_values, $val) : push(@sub_values, $val);
			}
			$mat_count++;	
		}
	}
	close(MATRIX_FILE);
}


my %ins_matrix;
my %sub_matrix;
my @bases=('A', 'T', 'G', 'C');
if($matrix_file ne "0"){
	foreach my $base1 (@bases){
		foreach my $base2(@bases){
			$ins_matrix{$base1}{$base2} = shift @ins_values;
			$sub_matrix{$base1}{$base2} = shift @sub_values;
		}
	}

	#convert to cumulative probabilities
	foreach my $base1 (@bases){
		my $prev_base = ' ';
		foreach my $base2 (@bases){
			if($prev_base ne ' '){
				$ins_matrix{$base1}{$base2} += $ins_matrix{$base1}{$prev_base};
				$sub_matrix{$base1}{$base2} += $ins_matrix{$base1}{$prev_base};
			}
			$prev_base = $base2;
		}
	}
}

my %hash;

my $min_q = 1;
my $max_q = 40;

#hash{pos}{pos-1}{qual}
for( my $pos = 0; $pos <= 1200; $pos++){
	for (my $row = $min_q; $row <= $max_q; $row++){
		for (my $col = $min_q; $col <= $max_q; $col++){
			$hash{$pos}{$row}{$col} = 1; #uninformed prior
		}
	}
}
open(MYFILE, ">$output");
my $seq_in= Bio::SeqIO->new(
                                -format => 'fastq',
                                -file   => $input,
                        );

my $prev_q1 = my $prev_q2 = $min_q;

#populate the markov chain based on previous value + current state
my $counter=0;
while( my $cur_seq = $seq_in->next_seq() ){
	my $prev_states = 0;
	my $prev_avg = 0;
	push( @lengths, length($cur_seq->seq()) ); 
	$counter += 1;
	if ($counter % 1000 == 0){
		print "Training: $counter\n";
	}
	my $position = 0;
	foreach my $cur_qual (split / /, $cur_seq->qual_text()){
		if($prev_states == 0){ 
			push(@seeds, $cur_qual);
			$prev_states += 1;
		}else{
			$prev_q2 = $prev_q1;
			$prev_states += 1;	
		}
		$prev_q1 = $cur_qual;
		$hash{$position}{$prev_q2}{$cur_qual} += 10;
		$position += 1;
	}
}			

my %rowsum;
my $length = max(@lengths);

for (my $pos = 0; $pos <= $length; $pos++){
	for (my $row = $min_q; $row <= $max_q; $row++){
		$rowsum{$pos}{$row} = 0;
	}

}

for (my $pos = 0; $pos <= $length; $pos++){ 
	for (my $row = $min_q; $row <= $max_q; $row++){
		for ( my $col = $min_q; $col <= $max_q; $col++){
			$rowsum{$pos}{$row} += $hash{$pos}{$row}{$col};
		}
	}
}

#correct for divide by zero error
for (my $pos = 0; $pos <= $length; $pos++){
	for (my $row = $min_q; $row <= $max_q; $row++){
		if ($rowsum{$pos}{$row} == 0){
			$rowsum{$pos}{$row} = 1;
		}
	}
}

#convert counts to probabilities
for (my $pos = 0; $pos <= $length; $pos++){
	for(my $row = $min_q; $row <= $max_q; $row++){
		for (my $col = $min_q; $col <= $max_q; $col++){

			$hash{$pos}{$row}{$col} = $hash{$pos}{$row}{$col} / $rowsum{$pos}{$row};
		}
	}
}

#convert counts to cumulative probabilities
for (my $pos = 0; $pos <= $length; $pos++){
	for (my $row = $min_q; $row <= $max_q; $row++){
		for (my $col = $min_q+1; $col <= $max_q; $col++){
			$hash{$pos}{$row}{$col} += $hash{$pos}{$row}{$col-1};
		}
	}
}


my $seqs = 0;
$counter=0;

open(MYFILE3, $fasta);
my $header = ' ';
my $curseq = ' ';
#generate sequences

my @nucleotides = ('A', 'T', 'G', 'C');
while(<MYFILE3>){
	my $line = $_;
        chomp($line);
        if($line =~ s/>/@/){
                $header = $line;
		next;
        }else{
		$curseq = $line;
	}
	$counter += 1;
	if ($counter % 1000 == 0){
		print "Generating $counter\n";
	}
	my $rand_per = rand();
	my $markov_qual = 0;

	$markov_qual = $seeds[rand @seeds];
	my @qual_string=();	
	push(@qual_string, chr($markov_qual+33));
	my $state_num = 0;
	my $prev1 = my $prev2 = my $prev3 = my $prev4 = my $prev5 = my $prev6 = my $avg = $min_q;

	my $rand_length = $lengths[rand @lengths];
	my $rand2 = rand();
	my $new_seq = substr($curseq, 0, $rand_length);
	for (my $i = 1; $i < $rand_length; $i++){ #$i=0 is already on the string
		my $avg2 = $avg;
		$rand_per = rand();
		my $cur_nuc = substr($new_seq, $i, 1);
		my $sub_check = rand();
		my $indel_check = rand();
		if ($state_num == 0){
			$avg = $min_q;
			$state_num += 1;
		}elsif ($state_num == 1){
			$prev2 = $prev1;
			$avg = $prev1;
			$state_num += 1;
		}elsif ($state_num == 2){
			$prev3 = $prev2;
			$prev2 = $prev1;
			$avg = int( ($prev2 + $prev3 ) / 2.0);
			$state_num += 1;
		}elsif ($state_num == 3){
			$prev4 = $prev3;
			$prev3 = $prev2;
			$prev2 = $prev1;
			$avg = int ( ($prev2 + $prev3 + $prev4) / 3.0);
			$state_num += 1;
		}elsif($state_num == 4){
			$prev5 = $prev4;
			$prev4 = $prev3;
			$prev3 = $prev2;
			$prev2 = $prev1;
			$avg = int ( ($prev2 + $prev3 + $prev4 + $prev5) / 4.0);
			$state_num += 1;
		}else{
			$prev6 = $prev5;
			$prev5 = $prev4;
			$prev4 = $prev3;
			$prev3 = $prev2;
			$prev2 = $prev1;
			$avg = int( ($prev2*.38) + ($prev3*.29) + ($prev4*.20) + ($prev5*.11) + ($prev6*.02));
		}
		my $nuc_check;
		if( $cur_nuc ~~ @nucleotides){ #handle ambiguous nucleotides, smartmatch requires perl 5.10
			$nuc_check = $cur_nuc;
		}else{
			$nuc_check = $nucleotides[rand @nucleotides];
		}
		$prev1 = $markov_qual;
		my $sub_rate_check = $err_rate ne "0" ? exp( ($rates{$nuc_check}{'slope'} * $i) + $rates{$nuc_check}{'intercept'}) : 10**((-$markov_qual)/10);
		#my $sub_rate_check = $err_rate ne "0" ? $rates{$nuc_check}{'slope'} * exp( exp($rates{$nuc_check}{'intercept'}) * $i) : 10**((-$markov_qual)/10);
		my $sub_rate=0;
		if($sub_rate_check > 1.0){
			$sub_rate = 1.0;
		}else{
			$sub_rate = $sub_rate_check;
		}
		#my $sub_rate = $err_rate ne "0" ? exp( ($rates{$nuc_check}{'slope'} * $i) + $rates{$nuc_check}{'intercept'})/100 : 10**((-$markov_qual)/10);
		#my $indel_rate = $err_rate ne "0" ? exp( ($rates{'X'}{'slope'} * $i) + $rates{'X'}{'intercept'})/100 : 0;
		#my $indel_rate_check = $err_rate ne "0" ? $rates{'X'}{'slope'} * exp( exp($rates{'X'}{'intercept'}) * $i) : 0;
		my $indel_rate_check = $err_rate ne "0" ? exp( ($rates{'X'}{'slope'} * $i) + $rates{'X'}{'intercept'}) : 0;
		my $indel_rate = 0;
		if($indel_rate_check > 1.0){
			$indel_rate = 1.0; 
		}else{
			$indel_rate = $indel_rate_check;
		}
		my $found = 0;
		my $test_col = 0;
		$rand_per = rand();
		my $seq_chunk = substr($new_seq, 0, $i);
		if($i != scalar(@qual_string)){
			my $els = scalar(@qual_string);
			print "POS $i\tELS $els\n";
			print "$new_seq\n",@qual_string,"\n";
			die;
		}
		$avg = $prev1;	
		my $type=rand();
		if( ($sub_check < $sub_rate)){
			if($sub_check < $sub_rate){
				if($err_rate ne "0"){
					my $cur_base = substr($new_seq, $i, 1);
					my $nuc_prob = rand();
					for my $base_check (@bases){
						if($nuc_prob < $sub_matrix{$cur_base}{$base_check}){
							$cur_base = $base_check;
							last;
						}
					}
					substr($new_seq, $i, 1) = $cur_base;
					#my $sub_qual = int( ($quals{$cur_nuc}{'slope'} * $i) + $quals{$cur_nuc}{'intercept'});
					my $new_sub_qual = int($quals{$cur_nuc}{'intercept'} + ($quals{$cur_nuc}{'coef1'} * $i) + (($quals{$cur_nuc}{'coef2'}) * ($i**2)));
					my $sub_qual;
					if ($new_sub_qual < $min_q){
						$sub_qual = $min_q;
					}elsif($new_sub_qual > $max_q){
						$sub_qual = $max_q;
					}else{
						$sub_qual = $new_sub_qual;
					}
					push(@qual_string, chr($sub_qual+33));
					$markov_qual = $sub_qual;
				}else{
					substr($new_seq, $i, 1) = $nucleotides[rand @nucleotides];
					my $sub_qual = $qual_string[$#qual_string];
					push(@qual_string, $sub_qual);
					$markov_qual = ord($sub_qual)-33;
				}
			}
		}else{
			for(my $col = $max_q; $col >= ($min_q); $col--){
				$avg2 = $avg;	
				if( ($rand_per == $hash{$i}{$avg}{$col}) ){
					push(@qual_string, chr($col+33));
					$markov_qual = $col;
					last;
				}elsif ( ($col != $min_q) && ($rand_per < $hash{$i}{$avg}{$col}) && ($hash{$i}{$avg}{$col-1} <= $rand_per)){
					push(@qual_string, chr($col+33));
					$markov_qual = $col;
					last;
				}elsif( ($col == $min_q+1) && ($hash{$i}{$avg}{$col-1} > $rand_per)) {
					push(@qual_string, chr($col-1+33));
					$markov_qual = $col-1;
					last;
				}elsif( ($col == $min_q)){ #&& ($hash{$i}{$avg}{$col} > $rand_per)){
					push(@qual_string, chr($col+33));
					$markov_qual = $col;
					last;

				}else{
					#print "RANDPER: $rand_per\n POS: $i\nAVG: $avg \n QUAL: $col $blah \n COL $hash{$i}{$avg}{$col}\n COL-1 $hash{$i}{$avg}{$col-1}\n";

				}
			}
			my $temp_els = scalar(@qual_string);
			my $cur_base = substr($new_seq, $i, 1);
			my $prev_base = substr($new_seq, $i-1, 1);
			my $h_total = $h_ins + $h_del;
			my $h_i_per = $h_ins / ($h_ins + $h_del);
			my $h_d_per = $h_del / ($h_ins + $h_del);
			#if($cur_base eq $prev_base){
				my $h_check = rand();
				my $h_max = 12;
				my $h_cur = 0;
				while($h_check < $h_total && $h_cur < $h_max){
					my $h_err_type = rand();
					if($h_err_type < $h_i_per){
					#homopolymer insertion
						#print "OLD SEQ: $new_seq\n";
						substr($new_seq, $i+1, 0, $cur_base);
						#print "NEW SEQ: $new_seq\n";
						my $new_ins_qual = int($quals{'X'}{'intercept'} + ($quals{'X'}{'coef1'} * $i) + (($quals{'X'}{'coef2'}) * ($i**2)));
						my $ins_qual;
						if($new_ins_qual > $max_q){
							$ins_qual = $max_q;
						}elsif($new_ins_qual < $min_q){
							$ins_qual = $min_q;
						}else{
							$ins_qual = $new_ins_qual;
						}
						push(@qual_string, chr($ins_qual+33));
						$i++;
						$rand_length++;
						$markov_qual = $ins_qual;
					}else{
						substr($new_seq, $i, 1) = "";
						$i--;
						$rand_length--;
						$markov_qual = ord(pop(@qual_string))-33;
					}
					$h_check = rand();
					$h_cur++;
				}
			#}


			if($indel_check < $indel_rate){
				if($type>$indel_rates{'D'}{$i}){ #insertion
					my $nuc_prob = rand();
					for my $base_check (@bases){
						if($nuc_prob < $ins_matrix{$cur_base}{$base_check}){
							$cur_base = $base_check;
							last;
						}
					}
					substr($new_seq, $i+1, 0, $cur_base);
					my $new_ins_qual = int($quals{'X'}{'intercept'} + ($quals{'X'}{'coef1'} * $i) + (($quals{'X'}{'coef2'}) * ($i**2)));
					my $ins_qual;
					if($new_ins_qual > $max_q){
						$ins_qual = $max_q;
					}elsif($new_ins_qual < $min_q){
						$ins_qual = $min_q;
					}else{
						$ins_qual = $new_ins_qual;
					}
					#my $ins_qual = int( ($quals{'X'}{'slope'} * $i) + $quals{'X'}{'intercept'});
                    			push(@qual_string, chr($ins_qual+33));
                    			$i++;
                    			$rand_length++;
                    			$markov_qual = $ins_qual;
				}else{ #deletion
					substr($new_seq, $i, 1) = "";
					$i--;
					$rand_length--;
					$markov_qual = ord(pop(@qual_string))-33;
				}
			}
		}
	}
	$prev1 = 0;
	$prev2 = 0;
	$prev3 = 0;
	$prev4 = 0;
	$prev5 = 0;
	$avg = 0;
	$seqs += 1;

	my $fasta_length = length($new_seq);
	my $qual_length = scalar(@qual_string);

	my $trimseq='';
	my $trimqual='';
	my $temp_qual= join "", @qual_string;	
        
	$trimseq = $new_seq;
    	$trimqual = $temp_qual;
	print MYFILE "$header\n";
	print MYFILE "$trimseq";
	print MYFILE  "\n+\n";
	print MYFILE "$trimqual";
	print MYFILE "\n";
}
close(MYFILE);
