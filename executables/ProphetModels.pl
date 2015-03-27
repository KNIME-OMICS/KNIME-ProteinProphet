#!/usr/bin/perl -w
#############################################################################
# Program       : ProphetModels.pl                                          #
# Author        : Henry Lam and David Shteynberg                            #
# Date          : 7.28.08                                                   #
#                                                                           #
#                                                                           #
#                                                                           #
#                                                                           #
# Copyright (C) 2008 Henry Lam , David Shteynberg                           # 
#                                                                           #
# This library is free software; you can redistribute it and/or             #
# modify it under the terms of the GNU Lesser General Public                #
# License as published by the Free Software Foundation; either              #
# version 2.1 of the License, or (at your option) any later version.        #
#                                                                           #
# This library is distributed in the hope that it will be useful,           #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU         #
# General Public License for more details.                                  #
#                                                                           #
# You should have received a copy of the GNU Lesser General Public          #
# License along with this library; if not, write to the Free Software       #
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA #
#                                                                           #
# Insitute for Systems Biology                                              #
# 1441 North 34th St.                                                       #
# Seattle, WA  98103  USA                                                   #
# akeller@systemsbiology.org                                                #
#                                                                           #
#############################################################################

use strict;
use Getopt::Long;

$, = ' ';
$\ = "\n";

# grab our tpplib exports from the same directory as this script
use File::Basename;
use Cwd qw(realpath cwd);
use lib realpath(dirname($0));
use tpplib_perl; # exported TPP lib function points

my $posobs_color = 4;
my $negobs_color = 3;
my $pos_color = 3;
my $neg_color = 1;
my $gnuplot_cmd = tpplib_perl::getGnuplotBinary();

#my $gnuplot_cmd = "gnuplot ";
my %ip_dec;

my %pp_dec;

my %options;


if (scalar(@ARGV) <= 0) {
    print "ProphetModels.pl by Henry Lam (maintained by David Shteynberg)";
    print "  Utility to plot PeptideProphet and iProphet models using gnuplot.";
    print '';
    print "Usage: ProphetModels.pl <options>";
    print "Options: -i <FILE> -- Specify pepXML file to processed.";
    print "         -d <STR>  -- Specify protein prefix indicating decoy to be counted.";
    print "         -x <STR>  -- Specify protein prefix indicating decoy to be excluded.";
    print "         -r <NUM>  -- Specify decoy ratio. Will guess from P<0.001 hits if not specified.";
    print "         -w <NUM>  -- Specify window for Prob/Prob plot. Default 100";
    print "         -k        -- Don't delete intermediate files.";
    print "         -M        -- Make Prob/Prob and Corr/Corr plots.";
    print "         -u        -- Consider only top probabilities of unique peptide sequences for iProphet.";
    print "         -q        -- Consider only top probabilities of unique peptide sequences for PeptideProphet.";
    print "         -n        -- Consider only top probabilities of unique peptide ions for PeptideProphet.";
    print "         -T        -- Consider only top probabilities for each PSMs";
    print "         -P        -- Provide simplistic protein ROC plot. Only uniquely-mapped proteins counted with prob = max prob of mapping peptides.";
    print "";

    exit;
}


GetOptions( \%options, 'decoy_string|d=s', 'exclude_string|x=s', 'input_file|i=s', 'decoy_ratio|r=s', 'window_prob|w=s', 'Make_comp_plots|M', 'keep_intermediate_files|k', 'uniq_iproph_peps|u', 'uniq_pproph_peps|n', 'uniq_pproph_pepseqs|q', 'uniq_psm|T', 'Protein_ROC|P');


my $pwd = cwd();
chomp $pwd;

my $infile = '';
if ($options{'input_file'}) {
    $infile = $options{'input_file'};
} elsif ($ARGV[-1] =~ /\.xml(\.gz)?/) {
    $infile = $ARGV[-1];
} else {
    die "No input xml file specified.";
}

if ($options{"uniq_pproph_pepseqs"} || $options{"uniq_iproph_peps"}) {
    $options{"uniq_pproph_peps"} = 1;
}

if ($options{"uniq_psm"}) {
    $options{"uniq_psm"} = 1;
}

#my $tmpinfile = tpplib_perl::uncompress_to_tmpfile($infile); # decompress .gz if needed



print "Analyzing $infile ...";
open(INFILE, $infile) or die "Cannot open file $infile";


(my $infile_pfx = $infile) =~ s/(.*)\.xml(\.gz)?/$1/g;

open(FVOUT, ">$infile_pfx\_FVAL\.tsv") or die "Cannot open file $infile_pfx\_FVAL\.tsv";


my %roc;
my %model;

my $spec = "";
my $decoy_string = $options{decoy_string} || 'DECOY';
my $decoy_ratio = $options{decoy_ratio} || -1;
my $pp_decoy_ratio = $options{decoy_ratio} || -1;
my $ip_decoy_ratio = $options{decoy_ratio} || -1;
my $prob_window =  $options{window_prob} || 100;
my $exclude_string = $options{exclude_string} || 'nEbRaSkA';

if ($decoy_ratio >= 0) {
    print "Assuming decoy ratio of $decoy_ratio";
}

# initializing arrays
my @pp_prob_array;
my @pp_prob_array_decoy;
my @pp_cum_neg_array;
my @ip_prob_array;
my @ip_prob_array_decoy;
my @ip_cum_neg_array;
my %ip_peps;
my %pp_peps;
my %psms;
my %ipsms;
my %pp_prots;
my %ip_prots;

for (my $ppp = 0; $ppp <= 10000; $ppp++) {
    $pp_prob_array[$ppp] = 0;
    $pp_prob_array_decoy[$ppp] = 0;
    $pp_cum_neg_array[$ppp] = 0.0;
    $ip_prob_array[$ppp] = 0;
    $ip_prob_array_decoy[$ppp] = 0;
    $ip_cum_neg_array[$ppp] = 0.0;
}

my $base_name = '';
my $line = '';
my $cur_model = '';
my $is_decoy = -1;
my $is_nondecoy = -1;
my $has_points = 0;
my $num_decoy = 0; # total number of decoys read so far
my $num_decoy_in_file = 0; # number of decoys read in one run (between <msms_summary> tags)
my $num_exclude = 0;
my $num_exclude_in_file = 0;
my $num_hits = 0;
my $num_hits_in_file = 0;
my $num_hits_set = 0;
my $pp_total_neg = 0.0;
my $ip_total_neg = 0.0;
my $num_ip_hits = 0;
my $num_ip_decoy = 0;

my $exclude = 0;
my $reject = 0;
my $charge = 0;
my $has_iprophet = 0;
my $has_fval = 0;
my $secondary_hit = 0;

my $pep = "";
my $pepion = "";
my $prot = "";
my $num_prot = 0;
my $pepchg = "";
my $pepcalcmass = "";
my $got_top_hit = 0;
my $has_obs = 0;
my $topcat = 1;

my $pp_prob = -100;

my $ip_prob = -100;

while ($line = <INFILE>) {

    if ($line =~ /^<roc_data_point min_prob=\"(.*)\" sensitivity=\"(.*)\" error=\"(.*)\" num_corr=\"(.*)\" num_incorr=\"(.*)\"/) {

	$roc{$1}->{"error"} = $3;
	$roc{$1}->{"ncorr"} = $4;
	$roc{$1}->{"nincorr"} = $5;

    }

    if ($line =~ /^<error_point error=\"(.*)\" min_prob=\"(.*)\" num_corr=\"(.*)\" num_incorr=\"(.*)\"/) {

	$roc{$2}->{"error"} = $1;
	$roc{$2}->{"ncorr"} = $3;
	$roc{$2}->{"nincorr"} = $4;
    }

    if ($line =~ /^<distribution_point fvalue=\"(.*)\" obs_1_distr=\"(.*)\" model_1_pos_distr=\"(.*)\" model_1_neg_distr=\"(.*)\" obs_2_distr=\"(.*)\" model_2_pos_distr=\"(.*)\" model_2_neg_distr=\"(.*)\" obs_3_distr=\"(.*)\" model_3_pos_distr=\"(.*)\" model_3_neg_distr=\"(.*)\" obs_4_distr=\"(.*)\" model_4_pos_distr=\"(.*)\" model_4_neg_distr=\"(.*)\" obs_5_distr=\"(.*)\" model_5_pos_distr=\"(.*)\" model_5_neg_distr=\"(.*)\" obs_6_distr=\"(.*)\" model_6_pos_distr=\"(.*)\" model_6_neg_distr=\"(.*)\" obs_7_distr=\"(.*)\" model_7_pos_distr=\"(.*)\" model_7_neg_distr=\"(.*)\"\/>/) {

	$has_fval = 1;
	print FVOUT $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22;

    }

    if ($line =~ /^<mixture_model precursor_ion_charge=\"([\d]+)\"/) {
	$topcat = 0;
	$charge = $1;
    }

    if (!$topcat && $line =~ /^<\/mixture_model>/) {
	$charge = 0;
    }
    if ($line =~ /^<mixturemodel name=\"TopCat\">/) {
	$topcat = 1;
    }
    if ($line =~ /^<mixturemodel name=\"(.*)\" pos_bandwidth=\"(.*)\" neg_bandwidth=\"(.*)\">/) {
	$topcat = 0;
	my $ch_str = '';
	$ch_str = "+" . $charge if ($charge > 0); 
	print "Reading $1 model $ch_str ...";
	$cur_model = $1;
	$cur_model = "AccMass" if ($cur_model eq 'Accurate Mass Model');
	$cur_model = "RT" if ($cur_model eq 'kernel density SSRCalc RT [RT]');
	$cur_model = "PI" if ($cur_model eq 'kernel density calc pI [pI]');

	if ($charge != 0) {
	    $cur_model = $cur_model . '_' . $charge;
	}

	open(PLOTOUT, ">$infile_pfx" . "_$cur_model" . ".tsv") or die "Cannot open $infile_pfx" . "_$cur_model" . ".tsv";
	$model{$cur_model}->{"posbw"} = $2;
	$model{$cur_model}->{"negbw"} = $3;
	$has_points = 0;
    }

    
    if ($line =~ /^<point value=\"(.*)\" pos_dens=\"(.*)\" neg_dens=\"(.*)\" neg_obs_dens=\"(.*)\" pos_obs_dens=\"(.*)\"\/>/) {

	my $pos = $2;
	$pos = 0.0000000001 if ($pos < 0.0000000001);
	my $neg = $3;
	$neg = 0.0000000001 if ($neg < 0.0000000001);
	my $negobs = $4;
	my $posobs = $5;
	print PLOTOUT $1, $2, $3, log($pos / $neg), $4, $5;
	$has_points = 1;
	$has_obs = 1;
    }
    elsif ($line =~ /^<point value=\"(.*)\" pos_dens=\"(.*)\" neg_dens=\"(.*)\"\/>/) {

	my $pos = $2;
	$pos = 0.0000000001 if ($pos < 0.0000000001);
	my $neg = $3;
	$neg = 0.0000000001 if ($neg < 0.0000000001);
	print PLOTOUT $1, $2, $3, log($pos / $neg);
	$has_points = 1;
	$has_obs = 0;
    }

    if (!$topcat && $line =~ /^<\/mixturemodel>/) {
	close(PLOTOUT);
	
	if ($has_points == 1) {

	    open(GPOUT, ">$infile_pfx" . "_$cur_model" . ".gp") or die "Cannot open $infile_pfx" . "_$cur_model" . ".pg";
	    print GPOUT "set terminal png truecolor ";
	    print GPOUT "set output \"$infile_pfx" . "_$cur_model" . ".png\"";
	    print GPOUT "set size 1.0, 1.0";
	    print GPOUT "set border 1";
	    print GPOUT "set xtics border nomirror out";
	    print GPOUT "set ytics border nomirror out";
	    print GPOUT "set y2tics border nomirror out";

#	    print GPOUT "set y2zeroaxis linetype -1 linewidth 1.0";
	    print GPOUT "set xzeroaxis linetype -1 linewidth 1.0";
	    print GPOUT "set x2zeroaxis linetype 0 linewidth 1.0";
	    print GPOUT "set xlabel \"$cur_model\"";
	    print GPOUT "set ylabel \"Density\"";
	    print GPOUT "set y2label \"ln(P/N)\"";
	    print GPOUT "set mxtics";
	    print GPOUT "set mytics";
	    print GPOUT "set my2tics";
	    print GPOUT "set key bottom right";
	    
	    # print GPOUT "set yrange [0.0:1.0]";
	    # print GPOUT, "set label \"$label\" at $x,$y front";
	    
	    if ($has_obs == 1) {
		print GPOUT "plot \"$infile_pfx" . "_$cur_model" . ".tsv\" using 1:2 title \"Positive Model (P)\" with line lw 2 lc $pos_color , \"$infile_pfx" . "_$cur_model" . ".tsv\" using 1:3 title \"Negative Model (N)\" with line lw 2 lc $neg_color,  \"$infile_pfx" . "_$cur_model" . ".tsv\" using 1:6 title \"Positive Observed\" with boxes lc $pos_color fs transparent solid 0.2 noborder, \"$infile_pfx" . "_$cur_model" . ".tsv\" using 1:5 title \"Negative Observed\" with boxes lc $neg_color fs transparent solid 0.2 noborder,  \"$infile_pfx" . "_$cur_model" . ".tsv\" using 1:4 title \"ln(P/N)\" axes x1y2 with line lc -1";
	    }
	    else {
		print GPOUT "plot \"$infile_pfx" . "_$cur_model" . ".tsv\" using 1:2 title \"Positive Model (P)\" with line lw 2 lc $pos_color , \"$infile_pfx" . "_$cur_model" . ".tsv\" using 1:3 title \"Negative Model (N)\" with line lw 2 lc $neg_color, \"$infile_pfx" . "_$cur_model" . ".tsv\" using 1:4 title \"ln(P/N)\" axes x1y2 with line lc -1";
	    }
	    print GPOUT "unset output";
	    close(GPOUT);
	    
	    `$gnuplot_cmd $infile_pfx\_$cur_model\.gp`;

#	    `sed -i 's|_$cur_model\.png|_$cur_model\_ZOOM.png|g' $infile_pfx\_$cur_model\.gp`;
#	    `sed -i 's|plot|plot [] [0.0 to 0.01]|g' $infile_pfx\_$cur_model\.gp`;
#           `$gnuplot_cmd $infile_pfx\_$cur_model\.gp`;

	    unless ($options{"keep_intermediate_files"}) {
		`rm -f $infile_pfx\_$cur_model\.gp $infile_pfx\_$cur_model\.tsv`;
	    }
	} else {
	    `rm -f $infile_pfx\_$cur_model\.tsv`;
	}
    }
    
    if ($line =~ /^<msms_run_summary.* base_name=\"(.*?)\"/) {
	$base_name = $1;
    }

    if ($line =~ /^<search_summary.* search_engine=\"(.*?)\"/) {
	print "Parsing search results \"$base_name ($1)\"...";
	$num_hits_in_file = 0;
	$num_decoy_in_file = 0;
	$num_exclude_in_file = 0;
    }

    if ($line =~ /^<spectrum_query (.*)/) {
	my $fields = $1;
	$spec = $1 if ($fields =~ /spectrum=\"(.*?)\"/);
	$pepchg = $1 if ($fields =~ /assumed_charge=\"([^\"]*)\"/);
	$got_top_hit = 0;
	$ip_prob = -100;
	$pp_prob = -100;
    }
    
    if ($line =~ /^<search_hit (.*)hit_rank=\"1\"(.*)/) {
	my $fields = $1 . $2;
	
	if ($got_top_hit == 0) {

	    $got_top_hit = 1;
	    $secondary_hit = 0;
	    $has_iprophet = 0;
	    $pep = $1 if ($fields =~ /peptide=\"([^\"]*)\"/);
	    $prot = $1 if ($fields =~ /protein=\"([^\"]*)\"/);
	    $num_prot = 1;

	    $pepion = $pepchg . $pep;
	    if ($prot !~ /^$exclude_string/) {
		$exclude = 0;
		if ($prot =~ /^$decoy_string/) {
		    $is_decoy = 1;
		} else {
		    $is_decoy = 0;
		}	
	    } else {
		$exclude = 1;
	    }	    
	    

	} else {
	    # already have a top hit
	    $got_top_hit++;
	}

	$reject = 0; $reject = 1 if ($fields =~ /is_rejected=\"1\"/);
	    
    }

    # check for secondary hit -- these will be ignored
    if ($line =~ /^<search_hit / && $line !~ /hit_rank=\"1\"/) {
	$secondary_hit = 1;
	$reject = 1;
    }

    # replace pep sequence with modified pep sequence if found
    if ($got_top_hit == 1 && $secondary_hit == 0 && $line =~ /^<modification_info.* modified_peptide=\"([^\"]*)\"/) {
	$pepion = $pepchg . $1;
    }

    if ($got_top_hit == 1 && $secondary_hit == 0 && $line =~ /^\s*<alternative_protein.* protein=\"(.*?)\".*>/) {
    	my $pt = $1;
	$num_prot++;

	#A decoy considered decoy when all alternatives are also decoys
	if ($pt !~ /^$exclude_string/) {
	    if ($pt !~ /^$decoy_string/) {
		$is_decoy = 0;
	    }
	    else {
		if ($is_decoy < 0) {
		    $is_decoy = 1 ;		    
		}
		else {
		    $is_decoy = 1 && $is_decoy; # no change, but make logic easier to grasp, to be decoy all must be decoy or excluded proteins
		}
	    }
	    
	    $exclude = 0;
	}
	else {
	    $exclude = 1 && $exclude; #To exclude all proteins matched must be excluded
	}	
    }

    # SpectraST case -- use lib_remark field to see if it is a decoy. Since this line occurs after all protein fields, this means
    # whether the lib_remark field is decoy will override what is indicated by protein names
    if ($got_top_hit == 1 && $secondary_hit == 0 && $line =~ /^<search_score name=\"lib_remark\" value=\"(.*)\"/) {
	my $pt = $1;
	
	if ($pt !~ /^$exclude_string/) {
	    if ($pt !~ /^$decoy_string/) {
		#### HENRY: also check alternative proteins to see if there is any non-decoy entries
		$is_decoy = 0;
	    } else {
		if ($is_decoy < 0) {
		    $is_decoy = 1 ;
		} else {
		    $is_decoy = 1 && $is_decoy;
		}
	    }	    
	    $exclude = 0;
	}
	else {
	    $exclude = 1 && $exclude; #To exclude all proteins matched must be excluded
	}
	
    }
 
    # end of search_hit, now count
    if (!$reject && $got_top_hit == 1 && $secondary_hit == 0 && $line =~ /^<\/search_hit>/) {

	if ($options{"uniq_pproph_pepseqs"}) {
	    $pepion = $pep;
	}

	$pp_prob = 1.0 if ($pp_prob < 0.0); # this means no peptideprophet probability is found, somehow set it to 1
	$ip_prob = 1.0 if ($has_iprophet && $ip_prob < 0.0);

	if ($exclude == 0) {

	    my $ignore_repeat_pp = 0;
	    my $ignore_repeat_ip = 0;

	    if (exists($psms{$spec})) {
		# this psm is already seen

		if ($psms{$spec} < $pp_prob) { # the current prob is higher than those observed for this psm thus far

		    if ($options{"uniq_psm"}) {
			# back out the probability entered before for this psm
			$pp_prob_array[int($psms{$spec} * 10000)]--;
			$pp_total_neg -= (1 - $psms{$spec});
			$pp_cum_neg_array[int( $psms{$spec} * 10000)] -= (1 -  $psms{$spec});
			$pp_prob_array_decoy[int($psms{$spec} * 10000)]-- if ($is_decoy == 1);
			$num_hits--;
			$num_hits_in_file--;
			if ($is_decoy == 1) {
			    $num_decoy--;
			    $num_decoy_in_file--;
			}
		    }

                    # update the max prob for this peptide
		    $psms{$spec} = $pp_prob; 

		} else { # the current prob is lower than those observed for this psm thus far
		    $ignore_repeat_pp = 1 if ($options{"uniq_psm"});
		}

	    } else {
		$psms{$spec} = $pp_prob;
	    }

	    if ($has_iprophet) {
		if (exists($ipsms{$spec})) {

		    if ($ipsms{$spec} < $ip_prob) { # the current prob is higher than those observed for this psm thus far

			if ($options{"uniq_psm"}) {
			    # back out the probability entered before for this psm
			    $ip_prob_array[int($ipsms{$spec} * 10000)]--;
			    $ip_total_neg -= (1 - $ipsms{$spec});
			    $ip_cum_neg_array[int( $ipsms{$spec} * 10000)] -= (1 -  $ipsms{$spec});
			    $ip_prob_array_decoy[int($ipsms{$spec} * 10000)]-- if ($is_decoy == 1);
			    $num_ip_hits--;
			    $num_ip_decoy-- if ($is_decoy == 1);
			}
			
			# update the max prob for this peptide
			$ipsms{$spec} = $ip_prob; 

		    } else { # the current prob is lower than those observed for this psm thus far
			$ignore_repeat_ip = 1 if ($options{"uniq_psm"});
		    }

		} else {
		    $ipsms{$spec} = $ip_prob;
		}
	    }

	    if (exists($pp_peps{$pepion})) {
		# this peptide is already seen

		if ($pp_peps{$pepion} < $pp_prob) { # the current prob is higher than those observed for this peptide thus far

		    if ($options{"uniq_pproph_peps"} or $options{"uniq_pproph_pepseqs"}) {
			# back out the probability entered before for this peptide
			$pp_prob_array[int($pp_peps{$pepion} * 10000)]--;
			$pp_total_neg -= (1 - $pp_peps{$pepion});
			$pp_cum_neg_array[int( $pp_peps{$pepion} * 10000)] -= (1 -  $pp_peps{$pepion});
			$pp_prob_array_decoy[int($pp_peps{$pepion} * 10000)]-- if ($is_decoy == 1);
			$num_hits--;
			$num_hits_in_file--;
			if ($is_decoy == 1) {
			    $num_decoy--;
			    $num_decoy_in_file--;
			}
		    }			
		    
                    # update the max prob for this peptide
		    $pp_peps{$pepion} = $pp_prob; 

		} else { # the current prob is lower than those observed for this peptide thus far
		    $ignore_repeat_pp = 1 if ($options{"uniq_pproph_peps"} or $options{"uniq_pproph_pepseqs"});
		}

	    } else {
		$pp_peps{$pepion} = $pp_prob;
		$pp_dec{$pepion} = $is_decoy;
	    }
	    
	    if ($has_iprophet) {
	
		if (exists($ip_peps{$pepion})) {

		    if ($ip_peps{$pepion} < $ip_prob) { # the current prob is higher than those observed for this peptide thus far

			if ($options{"uniq_iproph_peps"} or $options{"uniq_iproph_pepseqs"}) {
			    # back out the probability entered before for this peptide
			    $ip_prob_array[int($ip_peps{$pepion} * 10000)]--;
			    $ip_total_neg -= (1 - $ip_peps{$pepion});
			    $ip_cum_neg_array[int( $ip_peps{$pepion} * 10000)] -= (1 -  $ip_peps{$pepion});
			    $ip_prob_array_decoy[int($ip_peps{$pepion} * 10000)]-- if ($is_decoy == 1);
			    $num_ip_hits--;
			    $num_ip_decoy-- if ($is_decoy == 1);
			}			
			
			# update the max prob for this peptide
			$ip_peps{$pepion} = $ip_prob; 
			
		    } else { # the current prob is lower than those observed for this peptide thus far
			$ignore_repeat_ip = 1 if ($options{"uniq_iproph_peps"} or $options{"uniq_iproph_pepseqs"});
		    }
		    
		} else {
		    $ip_peps{$pepion} = $ip_prob;
		    $ip_dec{$pepion} = $is_decoy;
		}
	    }

	    if ($options{"Protein_ROC"} and $num_prot == 1) {

		$pp_prots{$prot} = $pp_prob unless (exists($pp_prots{$prot}) and $pp_prots{$prot} > $pp_prob);
	    
		if ($has_iprophet) {

		    $ip_prots{$prot} = $ip_prob unless (exists($ip_prots{$prot}) and $ip_prots{$prot} > $ip_prob);
		
		}

	    }

	    unless ($ignore_repeat_pp) {
		$num_hits++;
		$num_hits_in_file++;
		$pp_prob_array[int($pp_prob * 10000)]++;
		$pp_total_neg += (1 - $pp_prob);
		$pp_cum_neg_array[int($pp_prob * 10000)] += (1 - $pp_prob);
		if ($is_decoy == 1) {
		    $pp_prob_array_decoy[int($pp_prob * 10000)]++;
		    $num_decoy++;
		    $num_decoy_in_file++;
		}	
	    }

	       
	    if ($has_iprophet) {
		unless ($ignore_repeat_ip) {
		    $num_ip_hits++ if ($has_iprophet);
		    $ip_prob = 1.0 if ($ip_prob < 0.0);
		    $ip_prob_array[int($ip_prob * 10000)]++;
		    $ip_total_neg += (1 - $ip_prob);
		    $ip_cum_neg_array[int($ip_prob * 10000)] += (1 - $ip_prob);
		    if ($is_decoy == 1) {
			$ip_prob_array_decoy[int($ip_prob * 10000)]++;
			$num_ip_decoy++;
		    }	
		}	    
	    }

	} else {
	    # exclude = 1
	    $num_exclude++;
	    $num_exclude_in_file++;
	}

	$exclude = -1;
	$is_decoy = -1;
	$pep = "";
	$pepion = "";
	$prot = "";
	$num_prot = 0;
    }

    if ($line =~ /^<peptideprophet_result probability=\"([e\d\.\-]*)\"/) {
	$pp_prob = $1;		  
    }

    if ($line =~ /^<interprophet_result probability=\"([e\d\.\-]*)\"/) {
	$has_iprophet = 1;
	$ip_prob = $1;
    }

    if ($line =~ /^<\/spectrum_query>/) {
	$exclude = -1;
	$is_decoy = -1;
	$spec = "";
	$pep = "";
	$pepion = "";
	$prot = "";
	$num_prot = 0;
	$reject = 0;
    }
    
    if ($line =~ /^<\/msms_run_summary>/) {
	if ($num_hits_set != -1) {
	    my $thing = "hits";
	    $thing = "unique PSMs" if ($options{"uniq_psm"});
	    $thing = "unique peptide ions" if ($options{"uniq_pproph_peps"});
	    $thing = "unique peptide sequences" if ($options{"uniq_pproph_pepseqs"} || $options{"uniq_iproph_peps"});

	    print "  => Found $num_hits_in_file $thing. ($num_decoy_in_file decoys, $num_exclude_in_file excluded)";
	    print "  => Total so far: $num_hits $thing. ($num_decoy decoys, $num_exclude excluded)";
        }
        $base_name = '';
    }

}
     
close(INFILE);
#unlink($tmpinfile) if ($tmpinfile ne $infile); # did we decompress xml.gz?

close (FVOUT);



if ($has_fval == 1) {

    open(GPOUT, ">$infile_pfx\_FVAL.gp") or die "Cannot open $infile_pfx" . "_FVAL.gp";
    
    print GPOUT "set terminal png";
    print GPOUT "set size 1.0, 1.0";
    print GPOUT "set border 1";
    print GPOUT "set xtics border nomirror out";
    print GPOUT "set ytics border nomirror out";
    print GPOUT "set xzeroaxis linetype -1 linewidth 1.0";
    print GPOUT "set xlabel \"F-VALUE\"";
    print GPOUT "set ylabel \"Density\"";
    print GPOUT "set mxtics";
    print GPOUT "set mytics";
    
    print GPOUT "set yrange [0.0:]";
    
    for ($charge = 1; $charge <= 7; $charge++) {
	
	my $obs_col = $charge * 3 - 1;
	my $pos_col = $charge * 3;
	my $neg_col = $charge * 3 + 1;
	
	print GPOUT "set output \"$infile_pfx" . "_FVAL_$charge\.png\"";
	print GPOUT "plot \"$infile_pfx" . "_FVAL.tsv\" using 1:$obs_col title \"Observed\" with line lc -1 , \"$infile_pfx" . "_FVAL.tsv\" using 1:$pos_col title \"Model Pos\" with line lc $pos_color , \"$infile_pfx" . "_FVAL.tsv\" using 1:$neg_col title \"Model Neg\" with line lc $neg_color";
	
    }

    close(GPOUT);
    
    `$gnuplot_cmd $infile_pfx\_FVAL\.gp 2> /dev/null`;

    unless ($options{"keep_intermediate_files"}) {
	`rm -f $infile_pfx\_FVAL\.gp $infile_pfx\_FVAL\.tsv`;
    }

} else {
    `rm -f $infile_pfx\_FVAL\.tsv`;
}

if ($num_decoy > 0  and ($decoy_ratio > 0 or $pp_prob_array[0] > 0) ) {
    
    open(PLOTOUT, ">$infile_pfx" . "_FDR.tsv") or die "Cannot open $infile_pfx" . "_FDR.tsv";
    
    # foreach my $prob_cutoff (sort(keys(%roc))) {
    
    $pp_decoy_ratio = $decoy_ratio;
    $ip_decoy_ratio = $decoy_ratio;
    
    if ($decoy_ratio < 0) {
	$pp_decoy_ratio = ($pp_prob_array_decoy[0] / $pp_prob_array[0]);
	print "Using decoy ratio of $pp_decoy_ratio for PeptideProphet.";
	if ($has_iprophet == 1) {
	    if ($ip_prob_array[0] > 0) {
		$ip_decoy_ratio = ($ip_prob_array_decoy[0] / $ip_prob_array[0]);
		print "Using decoy ratio of $ip_decoy_ratio for iProphet.";
	    } else {
		# HENRY - somehow there is no decoys at P<=0.02, cannot use these hits to determine decoy ratio
		# just use that of PeptideProphet
		$ip_decoy_ratio = $pp_decoy_ratio;
		print "Cannot determine decoy ratio from iProphet. Use that of PeptideProphet: $ip_decoy_ratio .";
	    }

	}
    }

    my $fdr_ip_t = 1.0;
    my $fdr_ip_decoy_t = 1.0;
    my $fdr_pp_t = 1.0;
    my $fdr_pp_decoy_t = 1.0;
    my $num_corr_pp_decoy_t = -1;
    my $num_corr_ip_decoy_t = -1;
    my $print = 0;
    for (my $prob_cutoff = 0; $prob_cutoff <= 1; $prob_cutoff += 0.0001) {
	my $pp_num_decoy_below_prob = 0;
	my $pp_num_below_prob = 0;
	my $pp_cum_neg_below_prob = 0;
	my $ip_num_decoy_below_prob = 0;
	my $ip_num_below_prob = 0;
	my $ip_cum_neg_below_prob = 0;

	$print = 0;
	
	for (my $pp = 0; $pp < int($prob_cutoff * 10000); $pp++) {
	    $pp_num_decoy_below_prob += $pp_prob_array_decoy[$pp];
	    $pp_num_below_prob += $pp_prob_array[$pp];
	    $pp_cum_neg_below_prob += $pp_cum_neg_array[$pp];
	    
	    $ip_num_decoy_below_prob += $ip_prob_array_decoy[$pp];
	    $ip_num_below_prob += $ip_prob_array[$pp];
	    $ip_cum_neg_below_prob += $ip_cum_neg_array[$pp];
	}
	
	my $fdr_ip_roc = $roc{$prob_cutoff}->{"error"} || -1.0; 
	
	my $fdr_ip = 0.0;
	my $num_pos_ip =  $num_ip_hits - $ip_num_below_prob;
	my $num_corr_ip = $num_pos_ip;

	if ($has_iprophet == 1 and $num_pos_ip > 0) {
	    $fdr_ip = ($ip_total_neg - $ip_cum_neg_below_prob) / $num_pos_ip;
	    $num_corr_ip = $num_pos_ip * (1 - $fdr_ip);
	}
	
	my $fdr_pp = 0.0;
	my $num_pos_pp = $num_hits - $pp_num_below_prob;
	my $num_corr_pp = $num_pos_pp;

	if ($num_pos_pp > 0) {
	    $fdr_pp = ($pp_total_neg - $pp_cum_neg_below_prob) / $num_pos_pp;
	    $num_corr_pp = $num_pos_pp * (1 - $fdr_pp);
	}
	
	my $fdr_ip_decoy = 0.0;
	my $num_corr_ip_decoy = $num_pos_ip;
	if ($has_iprophet == 1 and $num_pos_ip > 0) {
	    my $num_incorr_ip = ($num_ip_decoy - $ip_num_decoy_below_prob) / $ip_decoy_ratio;
	    $fdr_ip_decoy = $num_incorr_ip / $num_pos_ip;
	    $num_corr_ip_decoy = $num_pos_ip - $num_incorr_ip;
	}

	my $fdr_pp_decoy = 0.0;
	my $num_corr_pp_decoy = $num_pos_pp;
	if ($num_hits - $pp_num_below_prob > 0) {
	    my $num_incorr_pp = ($num_decoy - $pp_num_decoy_below_prob) / $pp_decoy_ratio;
	    $fdr_pp_decoy = $num_incorr_pp / $num_pos_pp;
	    if ($fdr_pp_decoy < 0) {
		print "DDS: Uh oh!";
	    }
	    $num_corr_pp_decoy = $num_pos_pp - $num_incorr_pp;
	}
	my $PP_uncert = 0;
	if ($num_pos_pp > 0 && $fdr_pp*(1-$fdr_pp)/$num_pos_pp > 0){
	    $PP_uncert = sqrt($fdr_pp*(1-$fdr_pp)/$num_pos_pp);
	}
	my $PP_decoy_uncert = 0;
	if ($num_pos_pp > 0 && $fdr_pp_decoy*(1-$fdr_pp_decoy)/$num_pos_pp > 0) {
	    $PP_decoy_uncert = sqrt($fdr_pp_decoy*(1-$fdr_pp_decoy)/$num_pos_pp);

	}
	my $IP_uncert = 0;
	if ($num_pos_ip > 0 && $has_iprophet == 1 && $fdr_ip*(1-$fdr_ip)/$num_pos_ip > 0){
	    $IP_uncert = sqrt($fdr_ip*(1-$fdr_ip)/$num_pos_ip);
	}
	my $IP_decoy_uncert = 0;
	if ($num_pos_ip > 0 && $has_iprophet == 1 &&  $fdr_ip_decoy*(1-$fdr_ip_decoy)/$num_pos_ip > 0) {
	    $IP_decoy_uncert = sqrt($fdr_ip_decoy*(1-$fdr_ip_decoy)/$num_pos_ip);

	}

# HENRY -- instead of placing unreal FDR values to enforce monotonicity, in my opinion
# it is better to just skip the data point. The fact that FDR values can increase with higher prob cutoff
# is due to randomness in decoy placement; it's a sign that our bin size is too small.
# Skipping data points is equivalent to combining bins, until this abnormality goes away.
# The plot will then "interpolate" at the missing data points, basically "smoothing out" the randomness of decoy placement

	unless ($fdr_pp_decoy > $fdr_pp_decoy_t or $fdr_ip_decoy > $fdr_ip_decoy_t or ($num_corr_pp_decoy_t > 0 and $num_corr_pp_decoy > $num_corr_pp_decoy_t) or ($num_corr_ip_decoy_t > 0 and $num_corr_ip_decoy > $num_corr_ip_decoy_t)) {

	    print PLOTOUT $prob_cutoff, $fdr_ip, $fdr_ip_decoy, $fdr_pp, $fdr_pp_decoy, $num_corr_ip, $num_corr_ip_decoy, $num_corr_pp, $num_corr_pp_decoy,  $PP_uncert, $PP_decoy_uncert, $IP_uncert, $IP_decoy_uncert;
	    
	    $fdr_ip_decoy_t = $fdr_ip_decoy; 
	    $fdr_pp_decoy_t = $fdr_pp_decoy;
	    $num_corr_pp_decoy_t = $num_corr_pp_decoy;
	    $num_corr_ip_decoy_t = $num_corr_ip_decoy;
	}


    }  
    
    close (PLOTOUT);
    
    open(GPOUT, ">$infile_pfx\_FDR.gp") or die "Cannot open $infile_pfx" . "_FDR.gp";
    
    print GPOUT "set terminal png";
    print GPOUT "set output \"$infile_pfx" . "_FDR.png\"";
    print GPOUT "set size 1.0, 1.0";
    print GPOUT "set border 1";
    print GPOUT "set xtics border nomirror out";
    print GPOUT "set ytics border nomirror out";
    print GPOUT "set xzeroaxis linetype -1 linewidth 1.0";
    print GPOUT "set yzeroaxis linetype -1 linewidth 1.0";
    print GPOUT "set key bottom right";
    print GPOUT "set xlabel \"Decoy-estimated Error Rate\"";
    print GPOUT "set ylabel \"Prophet-predicted Error Rate\"";
#print GPOUT "set xlabel \"Probability Cutoff\"";
#print GPOUT "set ylabel \"FDR\"";
    print GPOUT "set mxtics";
    print GPOUT "set mytics";
    print GPOUT "set origin 0.0,0.0";
    
    print GPOUT "set yrange [0.0:]";
# print GPOUT, "set label \"$label\" at $x,$y front";
    
# print GPOUT "plot [0.0 to 1.0] \"$infile_pfx" . "_FDR.tsv\" using 1:2 title \"iProphet_PRED\" with line lc $pos_color , \"$infile_pfx" . "_FDR.tsv\" using 1:3 title \"iProphet_DECOY\" with line lc $neg_color , \"$infile_pfx" . "_FDR.tsv\" using 1:4 title \"PeptideProphet_PRED\" with line lc 4 , \"$infile_pfx" . "_FDR.tsv\" using 1:5 title \"PeptideProphet_DECOY\" with line lc 5";
    
    if ($has_iprophet == 1) {
	
	print GPOUT "plot \"$infile_pfx" . "_FDR.tsv\" using 3:2:13 title \"iProphet\" with xerrorbars  lc $pos_color , \"$infile_pfx" . "_FDR.tsv\" using 5:4:11 title \"PeptideProphet\" with xerrorbars  lc $neg_color , x notitle with line lt 0 lc -1";
	print GPOUT "set output \"$infile_pfx" . "_FDR_10pc.png\"";
	print GPOUT "plot [0.0 to 0.1] \"$infile_pfx" . "_FDR.tsv\" using 3:2:13 title \"iProphet\" with xerrorbars  lc $pos_color , \"$infile_pfx" . "_FDR.tsv\" using 5:4:11 title \"PeptideProphet\" with xerrorbars lc $neg_color , x notitle with line lt 0 lc -1";
	print GPOUT "set output \"$infile_pfx" . "_ROC.png\"";
	print GPOUT "set xlabel \"Decoy Estimated Error Rate\"";
	print GPOUT "set ylabel \"Decoy Estimaed Number of Corrects\"";
	print GPOUT "set key bottom right";	
	
#print GPOUT "plot [0.0:0.05] \"$infile_pfx" . "_FDR.tsv\" using 3:7:13 title \"iProphet(DECOY)\" with xerrorbars lc 4, \"$infile_pfx" . "_FDR.tsv\" using 5:9:11 title \"PeptideProphet(DECOY)\" with xerrorbars lc 5";


#	print GPOUT "plot [0.0 to 0.1]  \"$infile_pfx" . "_FDR.tsv\" using 3:7 title \"iProphet\" with lines lc $pos_color,  \"$infile_pfx" . "_FDR.tsv\" using 5:9 title \"PeptideProphet\" with lines lc $neg_color";

	print GPOUT "plot [0.0 to 0.1]  \"$infile_pfx" . "_FDR.tsv\" using 3:7:13 title \"iProphet(DECOY)\" with xerrorbars lc 4, \"$infile_pfx" . "_FDR.tsv\" using 2:6 title \"iProphet(PREDICTED)\" with lines lc $pos_color, \"$infile_pfx" . "_FDR.tsv\" using 5:9:11 title \"PeptideProphet(DECOY)\" with xerrorbars lc 5, \"$infile_pfx" . "_FDR.tsv\" using 4:8 title \"PeptideProphet(PREDICTED)\" with lines lc $neg_color";


	if ($options{"Make_comp_plots"}) {

	    print GPOUT "set output \"$infile_pfx" . "_CORR.png\"";
	    print GPOUT "set xlabel \"DECOY Estimated Count\"";
	    print GPOUT "set ylabel \"Prophet Estimated Count\"";
	    print GPOUT "set yrange [4000:10000]";
	    print GPOUT "plot [4000 to 10000] \"$infile_pfx" . "_FDR.tsv\" using 7:6 title \"iProphet\" with lines  lc $pos_color , \"$infile_pfx" . "_FDR.tsv\" using 9:8 title \"PeptideProphet\" with lines  lc $neg_color , x notitle with line lt 0 lc -1";
	}

    } else {
	print GPOUT "plot \"$infile_pfx" . "_FDR.tsv\" using 5:4:11 title \"PeptideProphet\" with xerrorbars  lc $neg_color , x notitle with line lt 0 lc -1";
	print GPOUT "set output \"$infile_pfx" . "_FDR_10pc.png\"";
	print GPOUT "plot [0.0 to 0.1] \"$infile_pfx" . "_FDR.tsv\" using 5:4:11 title \"PeptideProphet\" with xerrorbars  lc $neg_color , x notitle with line lt 0 lc -1";
	print GPOUT "set output \"$infile_pfx" . "_ROC.png\"";
	print GPOUT "set xlabel \"FDR\"";
	print GPOUT "set ylabel \"Number of Correct Hits\"";
	print GPOUT "plot [0.0 to 0.1] \"$infile_pfx" . "_FDR.tsv\" using 5:9:11 title \"PeptideProphet(DECOY)\" with xerrorbars lc 5,  \"$infile_pfx" . "_FDR.tsv\" using 4:8 title \"PeptideProphet(PREDICTED)\" with lines lc $neg_color ";
#	print GPOUT "plot [0.0 to 0.05] \"$infile_pfx" . "_FDR.tsv\" using 4:8 title \"PeptideProphet(PREDICTED)\" with line lc $neg_color , \"$infile_pfx" . "_FDR.tsv\" using 5:9 title \"PeptideProphet(DECOY)\" with line lc 5";	

    }
    
    close(GPOUT);
    
    `$gnuplot_cmd $infile_pfx\_FDR\.gp`;
    
    unless ($options{"keep_intermediate_files"}) {
	`rm -f $infile_pfx\_FDR\.gp $infile_pfx\_FDR\.tsv`;
    }
}

if ($options{"Make_comp_plots"}) {

    open(PROUT, ">$infile_pfx\_PPPROB\.tsv") or die "Cannot open file $infile_pfx\_PPPROB\.tsv";
    
    my $wdw = $prob_window;
    my $cnt = 0;
    my $pp_pr = 0;
    my $dec_count = 0;
    
    my @keys = sort {$pp_peps{$b} <=> $pp_peps{$a}} keys %pp_peps;
    
    for ($cnt = 0; $cnt <= $#keys; $cnt++)  {
	my $key = $keys[$cnt];
	
	$pp_pr += $pp_peps{$key};
	$dec_count +=  $pp_dec{$key};
	
	if ($cnt + 1 < $wdw)  {
	    next;
	}
	elsif ($cnt + 1 > $wdw) {
	    $pp_pr -= $pp_peps{$keys[$cnt-$wdw]};
	    $dec_count -=  $pp_dec{$keys[$cnt-$wdw]};
	}
	
	my $dec_prob = ( $wdw - $dec_count / $pp_decoy_ratio ) / $wdw;
	
	print PROUT $pp_pr/$wdw . "\t" . $dec_prob ;
	
	
    }
    close(PROUT);
    
    open(GPOUT, ">$infile_pfx\_PPPROB.gp") or die "Cannot open $infile_pfx" . "_PPPROB.gp";
    
    print GPOUT "set terminal png";
    print GPOUT "set output \"$infile_pfx" . "_PPPROB.png\"";
    print GPOUT "set size 1.0, 1.0";
    print GPOUT "set border 1";
    print GPOUT "set key bottom right";   
    print GPOUT "set xtics border nomirror out";
    print GPOUT "set ytics border nomirror out";
    print GPOUT "set xzeroaxis linetype -1 linewidth 1.0";
    print GPOUT "set yzeroaxis linetype -1 linewidth 1.0";
    print GPOUT "set xlabel \"Decoy-estimated Prob\"";
    print GPOUT "set ylabel \"Prophet-predicted Prob\"";
#print GPOUT "set xlabel \"Probability Cutoff\"";
#print GPOUT "set ylabel \"FDR\"";
    print GPOUT "set mxtics";
    print GPOUT "set mytics";
    print GPOUT "set origin 0.0,0.0";
    
    print GPOUT "set yrange [0.0:]";
    print GPOUT "set xrange [0.0:]";
    
    print GPOUT "plot \"$infile_pfx" . "_PPPROB.tsv\" using 2:1 title \"PeptideProphet\" with line  lt 1 lc $pos_color  , x notitle with line lt 0 lc -1";
    
    
    close(GPOUT);
    `$gnuplot_cmd $infile_pfx\_PPPROB\.gp`;
    
    
    
    if ($has_iprophet) {
	
	open(PROUT, ">$infile_pfx\_IPPROB\.tsv") or die "Cannot open file $infile_pfx\_IPPROB\.tsv";
	
	$wdw = $prob_window;
	$cnt = 0;
	my $ip_pr = 0;
	$dec_count = 0;
	@keys = sort {$ip_peps{$b} <=> $ip_peps{$a}} keys %ip_peps;
	
	
	for ($cnt = 0; $cnt <= $#keys; $cnt++)  {
	    my $key = $keys[$cnt];
	    
	    $ip_pr += $ip_peps{$key};
	    $dec_count +=  $ip_dec{$key};
	    
	    if ($cnt+1 < $wdw)  {
		next;
		
	    }
	    elsif ($cnt+1 > $wdw) {
		$ip_pr -= $ip_peps{$keys[$cnt-$wdw]};
		$dec_count -=  $ip_dec{$keys[$cnt-$wdw]};
	    }
	    
	    my $dec_prob = ( $wdw - $dec_count / $ip_decoy_ratio ) / $wdw;
	    
	    print PROUT $ip_pr/$wdw . "\t" . $dec_prob;
	    
	    
	}
	close(PROUT);
	
	open(GPOUT, ">$infile_pfx\_IPPROB.gp") or die "Cannot open $infile_pfx" . "_IPPROB.gp";
	
	print GPOUT "set terminal png";
	print GPOUT "set output \"$infile_pfx" . "_IPPROB.png\"";
	print GPOUT "set size 1.0, 1.0";
	print GPOUT "set border 1";
	print GPOUT "set key bottom right";   
	print GPOUT "set xtics border nomirror out";
	print GPOUT "set ytics border nomirror out";
	print GPOUT "set xzeroaxis linetype -1 linewidth 1.0";
	print GPOUT "set yzeroaxis linetype -1 linewidth 1.0";
	print GPOUT "set xlabel \"Decoy-estimated Prob\"";
	print GPOUT "set ylabel \"Prophet-predicted Prob\"";
#print GPOUT "set xlabel \"Probability Cutoff\"";
#print GPOUT "set ylabel \"FDR\"";
	print GPOUT "set mxtics";
	print GPOUT "set mytics";
	print GPOUT "set origin 0.0,0.0";
	
	print GPOUT "set yrange [0.0:]";
	print GPOUT "set xrange [0.0:]";
	
	print GPOUT "plot \"$infile_pfx" . "_IPPROB.tsv\" using 2:1 title \"iProphet\" with line lt 1 lc $pos_color  , \"$infile_pfx" . "_PPPROB.tsv\" using 2:1 title \"PeptideProphet\" with line  lt 1 lc $neg_color  , x notitle with line lt 0 lc -1";
	
	close(GPOUT);
	`$gnuplot_cmd $infile_pfx\_IPPROB\.gp`;

	unless ($options{"keep_intermediate_files"}) {
	    `rm -f $infile_pfx\_IPPROB\.gp $infile_pfx\_IPPROB\.tsv`;
	}
	
    }
    
    
    unless ($options{"keep_intermediate_files"}) {
	`rm -f $infile_pfx\_PPPROB\.gp $infile_pfx\_PPPROB\.tsv`;
    }

}

if ($options{Protein_ROC} and ($num_decoy > 0  and ($decoy_ratio > 0 or $pp_prob_array[0] > 0))) {

    print "Processing proteins ...";
    # process proteins 

    my $num_proteins_pp = 0;
    my $num_proteins_decoy_pp = 0;
    my $num_proteins_ip = 0;
    my $num_proteins_decoy_ip = 0;    

    my @pp_proteins_array_uniq = ();
    my @pp_proteins_array_uniq_decoy = ();
    my @ip_proteins_array_uniq = ();
    my @ip_proteins_array_uniq_decoy = ();
    for (my $pp = 0; $pp <= 10000; $pp++) {
	$pp_proteins_array_uniq[$pp] = 0;
	$pp_proteins_array_uniq_decoy[$pp] = 0;
	$ip_proteins_array_uniq[$pp] = 0;
	$ip_proteins_array_uniq_decoy[$pp] = 0;
    }

    foreach my $pr (keys(%pp_prots)) {
	$pp_proteins_array_uniq[int($pp_prots{$pr} * 10000)]++;
	$num_proteins_pp++;
	if ($pr =~ /^$decoy_string/) {
	    $num_proteins_decoy_pp++;
	    $pp_proteins_array_uniq_decoy[int($pp_prots{$pr} * 10000)]++;
	}
    }
    
    print " => Found $num_proteins_pp uniquely-mapped proteins ($num_proteins_decoy_pp decoys).";

    if ($has_iprophet) {

	foreach my $pr (keys(%ip_prots)) {
	    $ip_proteins_array_uniq[int($ip_prots{$pr} * 10000)]++;
	    $num_proteins_ip++;
	    if ($pr =~ /^$decoy_string/) {
		$num_proteins_decoy_ip++;
		$ip_proteins_array_uniq_decoy[int($ip_prots{$pr} * 10000)]++;
	    }
	}
    }
    
    open(PLOTOUT, ">$infile_pfx" . "_PROT.tsv") or die "Cannot open $infile_pfx" . "_PROT.tsv";


    my $fdr_proteins_pp_t = 1.0;
    my $fdr_proteins_ip_t = 1.0;
    my $num_corr_proteins_pp_t = $num_proteins_pp;
    my $num_corr_proteins_ip_t = $num_proteins_ip;

    for (my $prob_cutoff = 0; $prob_cutoff <= 1; $prob_cutoff += 0.0001) {
	my $pp_num_proteins_below_prob = 0;
	my $ip_num_proteins_below_prob = 0;
	my $pp_num_proteins_decoy_below_prob = 0;
	my $ip_num_proteins_decoy_below_prob = 0;

	for (my $pp = 0; $pp < int($prob_cutoff * 10000); $pp++) {
	    $pp_num_proteins_below_prob += $pp_proteins_array_uniq[$pp];
	    $pp_num_proteins_decoy_below_prob += $pp_proteins_array_uniq_decoy[$pp];
	    $ip_num_proteins_below_prob += $ip_proteins_array_uniq[$pp];
	    $ip_num_proteins_decoy_below_prob += $ip_proteins_array_uniq_decoy[$pp];
	}

	my $num_pos_proteins_ip = $num_proteins_ip - $ip_num_proteins_below_prob;
	my $num_pos_proteins_pp = $num_proteins_pp - $pp_num_proteins_below_prob;


	my $fdr_proteins_ip = 0.0;
	my $num_corr_proteins_ip = $num_pos_proteins_ip;
	my $num_incorr_proteins_ip = 0;
	if ($has_iprophet == 1 and $num_pos_proteins_ip > 0) {
	    $num_incorr_proteins_ip = ($num_proteins_decoy_ip - $ip_num_proteins_decoy_below_prob) / $ip_decoy_ratio;
	    $fdr_proteins_ip = $num_incorr_proteins_ip / $num_pos_proteins_ip;
	    $num_corr_proteins_ip = $num_pos_proteins_ip - $num_incorr_proteins_ip;
	}

	my $fdr_proteins_pp = 0.0;
	my $num_corr_proteins_pp = $num_pos_proteins_pp;
	my $num_incorr_proteins_pp = 0;
	if ($num_pos_proteins_pp > 0) {
	    $num_incorr_proteins_pp = ($num_proteins_decoy_pp - $pp_num_proteins_decoy_below_prob) / $pp_decoy_ratio;
	    $fdr_proteins_pp = $num_incorr_proteins_pp / $num_pos_proteins_pp;
	    $num_corr_proteins_pp = $num_pos_proteins_pp - $num_incorr_proteins_pp;
	}

	# print $num_proteins_pp, $num_proteins_decoy_pp, $prob_cutoff, $num_pos_proteins_pp, $num_incorr_proteins_pp, $num_corr_proteins_pp, $fdr_proteins_pp;

	# $fdr_proteins_ip = $fdr_proteins_ip_t if ($fdr_proteins_ip > $fdr_proteins_ip_t);
	# $fdr_proteins_pp = $fdr_proteins_pp_t if ($fdr_proteins_pp > $fdr_proteins_pp_t);
	# $num_corr_proteins_ip = $num_corr_proteins_ip_t if ($num_corr_proteins_ip_t > 0 and $num_corr_proteins_ip > $num_corr_proteins_ip_t);
	# $num_corr_proteins_pp = $num_corr_proteins_pp_t if ($num_corr_proteins_pp_t > 0 and $num_corr_proteins_pp > $num_corr_proteins_pp_t);

	unless ($fdr_proteins_ip > $fdr_proteins_ip_t or $fdr_proteins_pp > $fdr_proteins_pp_t) { #  or $num_corr_proteins_ip > $num_corr_proteins_ip_t or $num_corr_proteins_pp > $num_corr_proteins_pp_t) { 
	    print PLOTOUT $prob_cutoff, $fdr_proteins_ip, $fdr_proteins_pp, $num_pos_proteins_ip, $num_pos_proteins_pp;
		

	    $fdr_proteins_ip_t = $fdr_proteins_ip;
	    $fdr_proteins_pp_t = $fdr_proteins_pp;
	    $num_corr_proteins_pp_t = $num_corr_proteins_pp;
	    $num_corr_proteins_ip_t = $num_corr_proteins_ip;

	}

    }
    close(PLOTOUT);

    open(GPOUT, ">$infile_pfx\_PROT.gp") or die "Cannot open $infile_pfx" . "_PROT.gp";

    if ($has_iprophet == 1) {
	print GPOUT "set terminal png";
	print GPOUT "set output \"$infile_pfx" . "_PROT.png\"";
	print GPOUT "set xlabel \"Decoy-estimated Protein FDR\"";
	print GPOUT "set ylabel \"Number of Uniquely Mapped Proteins\"";
	print GPOUT "set key bottom right";	
	print GPOUT "set yrange [0:]";
	print GPOUT "plot [0:0.1] \"$infile_pfx" . "_PROT.tsv\" using 2:4 title \"iProphet\" with line lc $pos_color , \"$infile_pfx" . "_PROT.tsv\" using 3:5 title \"PeptideProphet\" with line lc $neg_color";

	
    } else {
	print GPOUT "set terminal png";
	print GPOUT "set output \"$infile_pfx" . "_PROT.png\"";
	print GPOUT "set xlabel \"Decoy-estimated Protein FDR\"";
	print GPOUT "set ylabel \"Number of Uniquely Mapped Proteins\"";
	print GPOUT "set key bottom right";
	print GPOUT "set yrange [0:]";
	print GPOUT "plot [0:0.1] \"$infile_pfx" . "_PROT.tsv\" using 3:5 title \"PeptideProphet\" with line lc $neg_color";

    }
    
    close(GPOUT);
    
    `$gnuplot_cmd $infile_pfx\_PROT\.gp`;
    
    unless ($options{"keep_intermediate_files"}) {
	`rm -f $infile_pfx\_PROT\.gp $infile_pfx\_PROT\.tsv`;
    }

}
