#!/usr/bin/perl
#############################################################################
# Program       : tpp_models.pl                                             #
# Author        : Luis Mendoza <lmendoza at systems biology>                #
# Date          : 10.11.12                                                  #
# SVN Info      : $Id: tpp_models.pl 6567 2014-07-18 21:54:41Z real_procopio $
#                                                                           #
# Display TPP analysis info for a set of pepXML and/or protXML files in a   #
# "dashboard' view                                                          #
#                                                                           #
# Copyright (C) 2012-2014 Luis Mendoza                                      #
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
# Luis Mendoza                                                              #
# Institute for Systems Biology                                             #
# 401 Terry Avenue North                                                    #
# Seattle, WA  98109  USA                                                   #
#                                                                           #
#############################################################################
use strict;
use Cwd qw(realpath);
use File::Basename;
use Getopt::Long;
use XML::Parser;
use lib realpath(dirname($0));
use tpplib_perl; # exported TPP lib function points
my $TPPVersionInfo = tpplib_perl::getTPPVersionInfo();

# Where things are
my $js_flot_loc = '/ISB/html/js';

my @files;
my %data;
my %graph_data;
my $debug_info;

my $errors;
$|++; # autoflush

my $USAGE=<<"EOU";
Usage:
   $0 file.xml || directory/

   If an input file is not specified, then all pepXML and protXML files in the input directory are processed.
EOU

my %options;
GetOptions(\%options,'D');

my $DEBUG   = $options{'D'} || '';;

my $inproc  = shift || die "I need an input file or directory!\n\n".$USAGE;

if (-d $inproc) {
    opendir DIR, $inproc || &reportError("BAD_DIR:$!",2);
    @files = map {"$inproc/$_"} grep /\.xml$/i, readdir DIR;
    closedir DIR;

} else {
    push @files, $inproc;
}


for my $infile (@files) {
    print "File: $infile\n";

    %data = {};
    %graph_data = {};
    $errors = 0;

    $debug_info  = "-"x79 ."\n";
    $debug_info .= "TPP DASHBOARD -- started at ". scalar(localtime) ."\n";
    $debug_info .= "-"x79 ."\n";

    my $dfile = $infile;

    unless ($dfile =~ s/.xml$/-MODELS.html/i) {
	&reportError("File name does not end in .xml!  Skipping ($infile)...",0);
	next;
    }

    my $format = &getFileFormat($infile);

    if ($format eq 'protxml') {
	&processProtXML($infile);
	&writeDashboardFile($infile,$dfile);
	print "$debug_info\n\n";

    } elsif ($format eq 'pepxml') {
	&processPepXML($infile);
	&writeDashboardFile($infile,$dfile);
	print "$debug_info\n\n";

    } else {
	&reportError("Unknown file type!",0);
    }

}

## Create dashboard.html file w/frames??

exit;


###################################################
sub writeDashboardFile {
    my $file = shift;
    my $dashboard_file = shift;

    $debug_info .= "--> Trying to write file $dashboard_file\n";
    my $css = &printCSS();

    open(DASHBOARD, ">$dashboard_file");

    print DASHBOARD<<END;
<html>
<head>
<title>TPP Dashboard</title>
$css
<script type="text/javascript" language="javascript" src="$js_flot_loc/jquery.min.js"></script>
<script type="text/javascript" language="javascript" src="$js_flot_loc/jquery.flot.js"></script>
<script type="text/javascript" language="javascript" src="$js_flot_loc/jquery.flot.selection.js"></script>
</head>
<body>

<div class="banner"><b>TPP</b>::Analysis and Models for <u>$file</u>
<br/><br/>
END

    # create empty tabs
    my $pad = '&nbsp;';
    my $table_tag = "<table style='background:#eeeeee' frame='border' rules='all' cellpadding='2'>\n";

    for my $i (0..9) {
	print DASHBOARD "<span>$pad</span>\n<span id='navtab$i'></span>\n";
    }

    print DASHBOARD "</div>\n";
    print DASHBOARD "<div class='opaque'></div>\n";
    print DASHBOARD "<div class='blank'>---</div>\n\n";
    print DASHBOARD "<script type='text/javascript' language='javascript'>addScores();</script>\n";

    $data{"html_table_summary"} .= "</table>\n" if $data{"html_table_summary"};
    $data{"html_table_msruns"} .= "</table>\n" if $data{"html_table_msruns"};
    $data{"html_table_models"} = "$table_tag$data{'html_table_models'}</table>\n" if $data{"html_table_models"};

    $data{"html_table_senserr"} .= "</table>\n" if $data{"html_table_senserr"};
    if ($data{"extra_html_table_error"}) {
	$data{"html_table_senserr"} .= "<span style='width:30px;float:left'>$pad</span>".$data{"extra_html_table_error"};
	$data{"html_table_senserr"} .= "</table>\n";
	$data{"tab_html_table_senserr"} .= "s";
    }

    print DASHBOARD &addTabbedPane(label => 'Models Charts');
    print DASHBOARD<<ENDJS;
    <script type="text/javascript" language="javascript">
	var modelmaxplot = new Array();
        var modelplot = new Array();
	var modelsopts = {
	  legend: {
	    position:'ne'
	    },
	  yaxis: {
	    min: 0
	    },
	  xaxis: {
	    tickSize: 1
	    },
	  selection: {
	    mode: "y"
	    }
	};
     </script>
ENDJS

    if (!keys %graph_data) {
	print DASHBOARD "No models found in input file!\n";
    }

    for my $att (sort keys %graph_data) {
	for my $charge (1..9) {
	    if ($att =~ /js_data_charge${charge}_model_obs/) {

		my $chargetitle = "Model Results for Charge +$charge";
		my $zoombutton  = "<input type='button' id='zoom$charge' value='Zoom' />";
		my $showmodels  = "true";

		if (!$data{"charge_${charge}_has_model"}) {
		    $chargetitle = "No models for Charge +$charge";
		    $zoombutton  = '';
		    $showmodels  = "false";
		    $graph_data{"js_data_charge${charge}_model_pos"} = '';
		    $graph_data{"js_data_charge${charge}_model_neg"} = '';
		}


		print DASHBOARD<<ENDJS;
		<br style='clear:both'>
		<div style='float:left'>
		<div class='tgray'><b>$chargetitle</b>&nbsp;&nbsp;&nbsp;&nbsp;$zoombutton</div>
		<div id='model_$charge' style='width:300px;height:200px;'></div></div>
		<script type="text/javascript" language="javascript">

		scores['plot'] = [];
		if ((scores['charge'] == '$charge') && scores['fval']) {
		    scores['plot'] = {
		      data: [[scores['fval'],0],[scores['fval'],$graph_data{"js_data_charge${charge}_obs_max"}]],
		      label:scores['fval'],
		      color: "#3B5998",
		      lines: {
			show:true
		      },
		      points: {
			show:false
		      }
		    };

		}

		var modeldata$charge = [
		   {
		     data: [ $graph_data{"js_data_charge${charge}_model_obs"} ],
		     label:'Observed',
		     color: "#000000",
		     lines: {
		       lineWidth: 1,
		       show:true
		     },
		     points: {
		       radius:1,
		       shadowSize:0,
		       show:true
		     }
		   },
		   {
		     data: [ $graph_data{"js_data_charge${charge}_model_pos"} ],
		     label:'Pos model',
		     color: "#89A54E",
		     lines: {
		       show:$showmodels
		     },
		     points: {
		       show:false
		     }
		   },
		   {
		     data: [ $graph_data{"js_data_charge${charge}_model_neg"} ],
		     label:'Neg model',
		     color: "#AA4643",
		     lines: {
		       show:$showmodels
		     },
		     points: {
		       show:false
		     }
		   },
		   scores['plot']
					];

	    \$(document).ready(function() {
		modelplot[$charge] = \$.plot(\$("#model_$charge"), modeldata$charge, modelsopts);
		modelmaxplot[$charge] = modelplot[$charge].getAxes().yaxis.max + 1.1*$graph_data{"js_data_charge${charge}_model_pos_max"};
	    });

	    \$("#zoom$charge").click(function () {
		var newy = modelmaxplot[$charge] - modelplot[$charge].getAxes().yaxis.max;
		modelplot[$charge] = \$.plot(\$("#model_$charge"), modeldata$charge,
			\$.extend(true, {}, modelsopts, { legend:{show:false}, yaxis:{max:newy} } ));

	    });
	    </script>
ENDJS

            next unless ($showmodels eq "true");

            for my $att (sort keys %graph_data) {  # another pass, to do sub-models

		if ($att =~ /js_data_($charge.*)_pos$/) {
		    next if ($att =~ /model|fval/);
		    my $gname = "charge$1";

		    my $neg = $att;
		    $neg =~ s/pos$/neg/;

		    my $dispname = uc $gname;
		    $dispname =~ s/charge/+/i;
		    $dispname =~ s/_/ /;

		    my $ticks = '';
		    if ($att =~ /_nmc_/)    { $ticks = 'ticks: [[0,""], [1,"0"], [2, "1-2"], [3, "3+"], [4, ""]]'; }
		    elsif ($att =~ /_ntt_/) { $ticks = 'ticks: [[-1, ""], [0,"0"], [1, "1"], [2, "2"], [3, ""]]';   }

		    print DASHBOARD<<ENDJS;
		    <div style='float:left'>
		    <div class='tgray'><b>$dispname</b> distributions</div>
		    <div id='${gname}_plot' style='width:250px;height:150px;'></div></div>
		    <script type="text/javascript" language="javascript">
		       var ${gname}_pos = $graph_data{$att};
		       var ${gname}_neg = $graph_data{$neg};

		    \$(document).ready(function() {
			\$.plot(\$("#${gname}_plot"), [
			  {
			    data: ${gname}_pos,
			    color: "#89A54E",
			    label:'Pos',
			    bars: {
			      barWidth: .6,
			      align:"center",
			      fill: true,
			      fillColor: "#D0DBB8",
			      horizontal:true,
			      show:true
			      }
			  },
			  {
			    data: ${gname}_neg,
			    color: "#AA4643",
			    label:'Neg',
			    bars: {
			      barWidth: .6,
			      align:"center",
			      fill: true,
			      fillColor: "#DDB5B4",
			      horizontal:true,
			      show:true
			      }
			  }],
			  {
			    legend: {
			      show: false,
			      position:'nw'
			      },
			    yaxis: {
			      tickSize: 1,
			      $ticks
			      },
			    xaxis: {
			      tickSize: 0.2
			      },
			    grid: {
			      markings: [ {
				color:"#000000",
				lineWidth:2,
				xaxis:{from:0, to:0}
			      }]
			    }
			 });
		    });
		    </script>
ENDJS

		}

		if ($att =~ /js_data_($charge)_sens/) {
		    my $sens =  $graph_data{"js_data_${charge}_sens"};
		    my $errr =  $graph_data{"js_data_${charge}_errr"};

		    print DASHBOARD<<ENDJS;
		    <div style='float:left'>
		    <div class='tgray'><b>+$charge Sens/Error Rates</b></div>
		    <div id='senserr_plot$charge' style='width:250px;height:150px;'></div></div>
		    <script type="text/javascript" language="javascript">
		    var sens$charge = [ $sens ];
	            var errr$charge = [ $errr ];

	    \$(document).ready(function() {
		\$.plot(\$("#senserr_plot$charge"), [
		   {
		     data: sens$charge,
		     color: "#89A54E",
		     label:'Sensitivity',
		     lines: {
		       show:true
		     },
		     points: {
		       show:true
		     }
		   },
		   {
		     data: errr$charge,
		     color: "#AA4643",
		     label:'Error',
		     lines: {
		       show:true
		     },
		     points: {
		       show:true
		     }
		   }],
		   {
		     legend: {
		       show:false
		     },
		     yaxis: {
		       min: 0,
		       max: 1.0
		     },
		     xaxis: {
		     },
		   });
	    });
	    </script>
ENDJS
              }

	    }

            }

	}

	if ( ($att =~ /js_data_nsp_pos/) || ($att =~ /js_data_fpkm_pos/) ) {
	    my $gt  = 'Learned NSP distributions';
	    my $pos = "var nsp_pos = $graph_data{js_data_nsp_pos};";
	    my $neg = "var nsp_neg = $graph_data{js_data_nsp_neg};";
	    my $dp  = "		     data: nsp_pos,";
	    my $dn  = "		     data: nsp_neg,";
	    my $pn  = "nsp_plot";

	    if ($att =~ /js_data_fpkm_pos/) {
		$gt  = 'Learned FPKM distributions';
		$pos = "var fpkm_pos = $graph_data{js_data_fpkm_pos};";
		$neg = "var fpkm_neg = $graph_data{js_data_fpkm_neg};";
		$dp  = "		     data: fpkm_pos,";
		$dn  = "		     data: fpkm_neg,";
		$pn  = "fpkm_plot";

	    }

	    print DASHBOARD<<ENDJS;
	    <div style='float:left'>
	    <div class='tgray'><b>$gt</b></div>
	    <div id='$pn' style='width:300px;height:200px;'></div></div>
	    <script type="text/javascript" language="javascript">
	        $pos
		$neg

	    \$(document).ready(function() {
		\$.plot(\$("#$pn"), [
		   {
		     $dp
		     color: "#89A54E",
		     label:'Positive freq',
		     bars: {
		       barWidth: .6,
		       align:"center",
		       fill: true,
		       fillColor: "#D0DBB8",
		       horizontal:true,
		       show:true
		     }
		   },
		   {
		     $dn
		     color: "#AA4643",
		     label:'Negative freq',
		     bars: {
		       barWidth: .6,
		       align:"center",
		       fill: true,
		       fillColor: "#DDB5B4",
		       horizontal:true,
		       show:true
		     }
		   }],
		   {
		     legend: {
		       position:'nw'
		     },
		     yaxis: {
		       min: -0.5,
		       tickSize: 1
		     },
		     xaxis: {
//		       tickSize: 0.1
		     },
                     grid: {
		       markings: [ {
			 color:"#000000",
			 lineWidth:2,
			 xaxis:{from:0, to:0}
		       }]
		     }
		   });
	    });
	    </script>
ENDJS
	}

	if ($att =~ /js_data_pvalue_obs/) {

	    print DASHBOARD<<ENDJS;
	    <br style='clear:both'>
	    <div style='float:left'>
	    <div class='tgray'><b>ASAPRatio pvalue</b></div>
	    <div id='pval_plot' style='width:600px;height:400px;'></div></div>
	    <script type="text/javascript" language="javascript">
	       var pval_obs = [ $graph_data{js_data_pvalue_obs} ];
   	       var pval_mod = [ $graph_data{js_data_pvalue_model} ];

	       \$(document).ready(function() {
		   \$.plot(\$("#pval_plot"), [
		   {
		     data: pval_obs,
		     color: "#000000",
		     label:'Data',
		     lines: {
		       lineWidth: 1,
		       show:true
		     },
		     points: {
		       radius:1,
		       shadowSize:0,
		       show:true
		     }
		   },
		   {
		     data: pval_mod,
		     color: "#89A54E",
		     label:'Fitting',
		     lines: {
		       show:true
		     }
		   }],
		   {
		     legend: {
		       show: true,
		       position:'ne'
		     },
                     grid: {
		       markings: [ {
			 color:"#000000",
			 lineWidth:1,
			 xaxis:{from:0, to:0}
		       }]
		     }
		   });
	      });
	    </script>
ENDJS
	}

	if ($att =~ /js_data_ALL_sens/) {
	    print DASHBOARD<<ENDJS;
	    <div style='float:left'>
	    <div class='tgray'><b>Predicted Sensitivity and Error Rate</b></div>
	    <div id='senserr_plot' style='width:300px;height:200px;'></div></div>
	    <script type="text/javascript" language="javascript">
		var sens = [ $graph_data{js_data_ALL_sens} ];
	        var errr = [ $graph_data{js_data_ALL_errr} ];

	    \$(document).ready(function() {
		\$.plot(\$("#senserr_plot"), [
		   {
		     data: sens,
		     color: "#89A54E",
		     label:'Sensitivity',
		     lines: {
		       show:true
		     },
		     points: {
		       show:true
		     }
		   },
		   {
		     data: errr,
		     color: "#AA4643",
		     label:'Error',
		     lines: {
		       show:true
		     },
		     points: {
		       show:true
		     }
		   }],
		   {
		     legend: {
		       position:'nw'
		     },
		     yaxis: {
		       min: 0,
		       max: 1.0,
		       tickSize: 0.1
		     },
		     xaxis: {
		       tickSize: 0.1
		     },
		   });
	    });
	    </script>
ENDJS

	}

        if ($att =~ /js_data_all_(.*)_posdens$/) {
	    my $gname = "$1";

	    my $neg = $att;
	    $neg =~ s/posdens$/negdens/;

	    my $dispname = uc $gname;
	    $dispname =~ s/_/ /;

	    print DASHBOARD<<ENDJS;
	    <div style='float:left'>
	    <div class='tgray'><b>$dispname</b></div>
	    <div id='${gname}_plot' style='width:250px;height:150px;'></div></div>
	    <script type="text/javascript" language="javascript">
	       var ${gname}_pos = [ $graph_data{$att} ];
	       var ${gname}_neg = [ $graph_data{$neg} ];

	       \$(document).ready(function() {
		   \$.plot(\$("#${gname}_plot"), [
		   {
		     data: ${gname}_pos,
		     color: "#89A54E",
		     label:'Pos Density',
		     lines: {
		       show:true
		     },
		   },
		   {
		     data: ${gname}_neg,
		     color: "#AA4643",
		     label:'Neg Density',
		     lines: {
		       show:true
		     }
		   }],
		   {
		     legend: {
		       show: false,
		       position:'ne'
		     },
		   });
	      });
	    </script>
ENDJS
                 }

    }
    print DASHBOARD "<br style='clear:both'>\n";

    print DASHBOARD
	"<!-- INSERT DECOY VALIDATION HERE -->\n",
	"<!-- END DECOY VALIDATION -->\n";

    print DASHBOARD &closeTabbedPane(selected=>1);

    for my $att (sort keys %data) {

	if ($att =~ /^html_table_/) {
	    print DASHBOARD &addTabbedPane(label => $data{'tab_'.$att});
	    print DASHBOARD "$data{$att}\n\n";
	    print DASHBOARD "<br style='clear:both'>\n";
	    print DASHBOARD &closeTabbedPane();
	}

    }
    print DASHBOARD "<!-- ADD EXTRA TABS HERE -->\n\n";

    $debug_info .= "-"x79 ."\n";
    $debug_info .= "Finished at " . scalar(localtime) . " with $errors errors.\n";
    $debug_info .= "-"x79 ."\n";

    print DASHBOARD &addTabbedPane(label => '(debug)');
    print DASHBOARD "<pre>$debug_info</pre>";
    print DASHBOARD &closeTabbedPane();
    print DASHBOARD "<br><br><hr size='1' noshade><font color='#999999'>TPP Dashboard<br>$TPPVersionInfo</font>\n</body>\n</html>";
    close DASHBOARD;
}

###################################################
sub processPepXML {
    my $file = shift;

    $data{'done_models'} = '';
    $data{'curr_model'} = '';
    $data{'first_param'} = 1;
    $data{'first_mod'} = 1;

    $data{"tab_html_table_summary"} = "Run Options";
    $data{"html_table_summary"} = "<table style='background:#eeeeee' frame='border' rules='all' cellpadding='2'>\n";

    $data{"tab_html_table_models"} = "Learned Models";

    for my $charge (1..9) {
	$graph_data{"js_data_charge${charge}_obs_max"} = 0;
	$graph_data{"js_data_charge${charge}_model_pos_max"} = 0;
    }

    $data{"tab_html_table_msruns"} = "MS Runs";
    $data{"html_table_msruns"} = "<table style='background:#eeeeee' frame='border' rules='all' cellpadding='2'>\n";

    $data{"tab_html_table_senserr"} = "Sens/Error Table";
    $data{"html_table_senserr"} = "<table style='background:#eeeeee;float:left;' frame='border' rules='all' cellpadding='2'>\n<tr class='theader'><td colspan='5'>Predicted Sensitivity and Error Rate</td></tr>\n<tr class='thead'>";
    for my $head (qw(min_prob Sensitivity Error_Rate num_correct num_incorrect)) {
	$data{"html_table_senserr"} .= "<th>$head</th>";
    }
    $data{"html_table_senserr"} .= "</tr>\n";

    #### Set up the XML parser and parse the returned XML
    my $parser = XML::Parser->new(
				  Handlers => {
				      Start => \&start_pepxml_element,
				      End   => \&end_pepxml_element,
#				      Char  => \&pepxml_chars,
				  },
				  ErrorContext => 2 );
    eval { $parser->parsefile( $file ); };
    if($@) {
	if ($@ =~ /^Finished parsing header/) {
	    $debug_info .= "Found end of header\n";
	} else { 
	    die "ERROR_PARSING_PEPXML:$@";
	}
    }

}

sub start_pepxml_element {
    my ($handler, $element, %atts) = @_;

    if ( ($element eq 'search_score') ||
	 ($element eq 'mod_aminoacid_mass') ||
	 ($element eq 'analysis_result') ||
	 ($element eq 'spectrum_query') ||
	 ($element eq 'search_result') ||
	 ($element eq 'search_hit') ||
	 ($element eq 'modification_info') ||
	 ($element eq 'search_score_summary') ||
	 ($element eq 'alternative_protein') ||
	 ($element eq 'intensity') ||
	 ($element eq 'libra_result') ||
	 ($element eq 'isotopic_contributions') ||
	 ($element eq 'interprophet_result') ||
	 ($element eq 'peptideprophet_result') ||
	 ($element eq 'dataset_derivation') ||
	 ($element eq 'msms_pipeline_analysis'))
    {
	    # ignore!
    }

    elsif ($element eq 'parameter') {
	if ($data{'curr_model'} eq 'NONE') {
	    # ignore!
	}

	elsif ($data{'curr_model'} eq 'SUMMARY') {
	    $data{"html_table_summary"} .= "<tr><td>$atts{'name'}</td><td>$atts{'value'}</td></tr>\n";
	}

	elsif ($data{'curr_model'} eq 'search_summary') {
	    if ($data{'first_param'}) {
		$data{"html_table_msruns"} .= "<tr class='thead'><th colspan='2'>Search Parameters</th></tr>\n";
		$data{'first_param'} = 0;
	    }
	    $data{"html_table_msruns"} .= "<tr><td>$atts{'name'}</td><td>$atts{'value'}</td></tr>\n"; # unless ($att eq 'precursor_ion_charge');
	}

	else {
	    $data{"html_table_models"} .= "$atts{'name'} $atts{'value'}, "; # unless ($att eq 'precursor_ion_charge');
	    my $barval = ($data{'curr_posneg'} eq 'neg') ? "-$atts{'value'}" : $atts{'value'};

	    if ($data{'curr_model'} eq 'nmc') {
		my $bin = 0;
		if    ($atts{'name'} eq 'nmc=0')     { $bin = 1; }
		elsif ($atts{'name'} eq '1<=nmc<=2') { $bin = 2; }
		elsif ($atts{'name'} eq 'nmc>=3')    { $bin = 3; }

		$graph_data{"js_data_$data{'curr_charge'}_$data{'curr_model'}_$data{'curr_posneg'}"} .= "[$barval, $bin], ";
	    }
	    else {
		$atts{'name'} =~ s/ntt=//;
		$graph_data{"js_data_$data{'curr_charge'}_$data{'curr_model'}_$data{'curr_posneg'}"} .= "[$barval, $atts{'name'}], ";
	    }
	}
    }

    elsif ($element eq 'msms_run_summary') {
	print " - in ms run: $atts{'base_name'}...\n";

	$data{"html_table_msruns"} .= "<tr class='theader'><th colspan='2'>MS run: $atts{'base_name'}</th></tr>\n";
	for my $att (sort keys %atts) {
	    $data{"html_table_msruns"} .= "<tr><td>$att</td><td>$atts{$att}</td></tr>\n" unless ($att eq 'base_name');
	}
    }

    elsif ($element eq 'search_summary') {
	$data{'curr_model'} = $element;
	$data{"html_table_msruns"} .= "<tr class='thead'><th colspan='2'>$element</th></tr>\n";

	for my $att (sort keys %atts) {
	    $data{"html_table_msruns"} .= "<tr><td>$att</td><td>$atts{$att}</td></tr>\n";
	}
    }

    elsif (($element eq 'search_database') ||
	   ($element eq 'sample_enzyme') ||
	   ($element eq 'enzymatic_search_constraint')
	   ) {
	$data{"html_table_msruns"} .= "<tr class='thead'><th colspan='2'>$element</th></tr>\n";

	for my $att (sort keys %atts) {
	    $data{"html_table_msruns"} .= "<tr><td>$att</td><td>$atts{$att}</td></tr>\n";
	}
    }

    elsif (($element eq 'aminoacid_modification') ||
	   ($element eq 'terminal_modification')) {
	if ($data{'first_mod'}) {
	    $data{"html_table_msruns"} .= "<tr class='thead'><th colspan='2'>Aminoacid and Terminal Modifications</th></tr>\n";
	    $data{'first_mod'} = 0;
	}

	if ($atts{'aminoacid'}) {
	    $data{"html_table_msruns"} .= "<tr><td style='text-align:right'><b>" . delete ($atts{'aminoacid'}) . "&nbsp;&nbsp;&nbsp;</b></td><td>";
	}
	elsif ($atts{'terminus'}) {
	    $data{"html_table_msruns"} .= "<tr><td style='text-align:right'><b>" . delete ($atts{'terminus'}) . "&nbsp;&nbsp;&nbsp;</b></td><td>";
	}
	else {
	    $data{"html_table_msruns"} .= "<tr><td style='text-align:right'>?&nbsp;&nbsp;&nbsp;</td><td>";
	}

	for my $att (sort keys %atts) {
	    $data{"html_table_msruns"} .= "$att: <b>$atts{$att}</b>  ";
	}
	$data{"html_table_msruns"} .= "</td></tr>\n";
    }

    elsif ($element eq 'specificity') {
	$data{"html_table_msruns"} .= "<tr><td>$element</td><td>\n";
	for my $att (sort keys %atts) {
	    $data{"html_table_msruns"} .= "$att: <b>$atts{$att}</b>  ";
	}
	$data{"html_table_msruns"} .= "</td></tr>\n";
    }

    elsif ($element eq 'mixture_model') {
	$data{"html_table_models"} .= "<tr class='theader'><th colspan='2'>+$atts{'precursor_ion_charge'} Models</th></tr>\n";
	$data{"curr_charge"} = $atts{'precursor_ion_charge'} || 'unknown_charge';

	for my $att (sort keys %atts) {
	    $data{"html_table_models"} .= "<tr><td>$att</td><td>$atts{$att}</td></tr>\n" unless ($att eq 'precursor_ion_charge');
	}
    }

    elsif ($element eq 'mixturemodel_distribution') {
	$data{"html_table_models"} .= "<tr class='thead'><th colspan='2'>$atts{'name'}</th></tr>\n";
	if ($atts{'name'} =~ /\[(.*)\]/) {
	    $data{"curr_model"} = $1;
	    $debug_info .= "Found $1 (+$data{'curr_charge'}) model...\n";
	} else {
	    $data{"curr_model"} = 'unknown';
	    $debug_info .= "Could not extract model [name] from: $atts{'name'}\n";
	}
    }

    elsif ($element eq 'mixturemodel') {
	$data{"html_table_models"} .= "<tr class='thead'><th colspan='2'>$atts{'name'}</th></tr>\n";

	$atts{'name'} =~ s/\s//g;
	if ($atts{'name'} =~ /\[(.*)\]/) {
	    $data{"curr_model"} = $1;
	} else {
	    $data{"curr_model"} = $atts{'name'};
	}
	$debug_info .= "Found $atts{'name'} ('+$data{'curr_charge'}') model...\n";

	for my $att (sort keys %atts) {
	    $data{"html_table_models"} .= "<tr><td>$att</td><td>$atts{$att}</td></tr>\n" unless ($att eq 'name');
	}
    }

    elsif ($element eq 'point') {
	if ($data{'done_models'} !~ /$data{"curr_model"}/) {
	    $graph_data{"js_data_all_$data{'curr_model'}_posdens"} .= "[$atts{'value'}, $atts{'pos_dens'}], ";
	    $graph_data{"js_data_all_$data{'curr_model'}_negdens"} .= "[$atts{'value'}, $atts{'neg_dens'}], ";
	} 
    }

    elsif ( ($element eq 'posmodel_distribution') || ($element eq 'negmodel_distribution') ) {
	$data{"curr_posneg"} = substr $element, 0, 3;
	$data{"html_table_models"} .= "<tr><td>$data{'curr_posneg'} model</td><td>";
	for my $att (sort keys %atts) {
	    $data{"html_table_models"} .= "<b>$atts{$att}</b> ";
	}
	$data{"html_table_models"} .= " (";

	$graph_data{"js_data_$data{'curr_charge'}_$data{'curr_model'}_$data{'curr_posneg'}"} = '[';
    }

    elsif ($element eq 'analysis_summary') {
	$data{'curr_model'} = 'SUMMARY';
	$data{"html_table_summary"} .= "<tr class='thead'><th colspan='2'>$atts{'analysis'}</th></tr>\n";

	for my $att (sort keys %atts) {
	    $data{"html_table_summary"} .= "<tr><td>$att</td><td>$atts{$att}</td></tr>\n" unless ($att eq 'analysis');
	}
    }

    elsif ($element eq 'analysis_timestamp') {
	$data{"timestamp_info"} .= "<tr><td>" . delete ($atts{'analysis'}) . "</td><td>" . delete ($atts{'time'});

	for my $att (sort keys %atts) {
	    $data{"timestamp_info"} .= " ($att: <b>$atts{$att}</b>)";
	}

	$data{"timestamp_info"} .= "</td></tr>\n";
    }
    elsif ($element eq 'database_refresh_timestamp') {
	for my $att (sort keys %atts) {
	    $data{"timestamp_info"} .= "<tr><td>$att</td><td>$atts{$att}</td></tr>\n";
	}
    }

    elsif ( ($element eq 'peptideprophet_summary') ||
	    ($element eq 'interact_summary') ||
	    ($element eq 'asapratio_summary') ||
	    ($element eq 'xpressratio_summary') ||
	    ($element eq 'libra_summary') ||
	    ($element eq 'inputfile') ) {

	$atts{$element} = delete $atts{'name'} if ($element eq 'inputfile');

	for my $att (sort keys %atts) {
	    $data{"html_table_summary"} .= "<tr><td>$att</td><td>$atts{$att}</td></tr>\n";
	}
    }

    elsif ($element eq 'roc_error_data') {
	$data{"curr_charge"} = $atts{'charge'} || 'unknown_charge';

	$data{"html_table_senserr"} .= "<tr class='theader'><td colspan='5'>+$atts{'charge'} only</td></tr>\n" unless ($atts{'charge'} eq 'all');
	if ($data{"extra_html_table_error"}) {
	    $data{"extra_html_table_error"} .= "<tr class='theader'><td colspan='4'>+$atts{'charge'} only</td></tr>\n" unless ($atts{'charge'} eq 'all');
	}
    }

    elsif ($element eq 'roc_data_point') {
	$data{"html_table_senserr"} .= "<tr><td>$atts{'min_prob'}</td><td>$atts{'sensitivity'}</td><td>$atts{'error'}</td><td>$atts{'num_corr'}</td><td>$atts{'num_incorr'}</td></tr>\n";

	my $chgstr = $data{"curr_charge"} ? uc $data{'curr_charge'} : 'ALL';

	$graph_data{"js_data_${chgstr}_sens"} .= "[$atts{'min_prob'}, $atts{'sensitivity'}], ";
	$graph_data{"js_data_${chgstr}_errr"} .= "[$atts{'min_prob'}, $atts{'error'}], ";
    }

    elsif ($element eq 'distribution_point') {
	for my $charge (1..9) {
#	    $debug_info .= "---FOUND $charge :: $atts{'obs_'.$charge.'_distr'} \n";
	    if ($atts{"obs_".$charge."_distr"}) {
#		$debug_info .= "---FOUND $atts{'obs_'.$charge.'_distr'} ---\n";
		$graph_data{"js_data_charge${charge}_model_obs"} .= "[$atts{'fvalue'}, ".$atts{"obs_${charge}_distr"}.'], ';
		$graph_data{"js_data_charge${charge}_model_pos"} .= "[$atts{'fvalue'}, ".$atts{"model_${charge}_pos_distr"}.'], ';
		$graph_data{"js_data_charge${charge}_model_neg"} .= "[$atts{'fvalue'}, ".$atts{"model_${charge}_neg_distr"}.'], ';

		$data{"charge_${charge}_has_model"} += $atts{"model_${charge}_pos_distr"} + $atts{"model_${charge}_neg_distr"};

		if ($graph_data{"js_data_charge${charge}_obs_max"} < $atts{"obs_${charge}_distr"}) {
		    $graph_data{"js_data_charge${charge}_obs_max"} = $atts{"obs_${charge}_distr"};
		}

		if ($graph_data{"js_data_charge${charge}_model_pos_max"} < $atts{"model_${charge}_pos_distr"}) {
		    $graph_data{"js_data_charge${charge}_model_pos_max"} = $atts{"model_${charge}_pos_distr"};
		}
	    }
	}
    }

    elsif ($element eq 'error_point') {
	if (!$data{"extra_html_table_error"}) {
	    $data{"extra_html_table_error"} = "<table style='background:#eeeeee;float:left;' frame='border' rules='all' cellpadding='2'>\n<tr class='theader'><td colspan='4'>Error Table</td></tr>\n<tr class='thead'>";
	    for my $head (qw(Error_Rate min_prob num_correct num_incorrect)) {
		$data{"extra_html_table_error"} .= "<th>$head</th>";
	    }
	    $data{"extra_html_table_error"} .= "</tr>\n";
	}

	$data{"extra_html_table_error"} .= "<tr><td>$atts{'error'}</td><td>$atts{'min_prob'}</td><td>$atts{'num_corr'}</td><td>$atts{'num_incorr'}</td></tr>\n";
    }

    elsif ($element eq 'fragment_masses') {
	$data{"html_table_summary"} .= "<tr><td>$element</td><td>";
	for my $att (sort keys %atts) {
	    $data{"html_table_summary"} .= "$att: <b>$atts{$att}</b>  ";
	}
	$data{"html_table_summary"} .= "</td></tr>\n";

    }

    elsif ($element eq 'contributing_channel') {
	$data{"html_table_summary"} .= "<tr><td>isotopic contribution <b>$atts{channel}</b></td><td>";
    }
    elsif ($element eq 'affected_channel') {
	$data{"html_table_summary"} .= "<b>$atts{correction}</b>-->$atts{channel}&nbsp;&nbsp;&nbsp;";
    }

    else {
	if ($DEBUG) {
	    my $name = $atts{'name'} ? "($atts{'name'})" : '';
	    $debug_info .= "Did not handle tag '$element' $name\n";
	}
    }

}

sub end_pepxml_element {
    my ($handler, $element) = @_;

    if ( ($element eq 'posmodel_distribution') || ($element eq 'negmodel_distribution') ) {
	substr $data{"html_table_models"}, -2, 2, ")</td></tr>\n";
#	$data{"html_table_models"} .= ")</td></tr>\n";

	$graph_data{"js_data_$data{'curr_charge'}_$data{'curr_model'}_$data{'curr_posneg'}"} .= ']';
	$data{"curr_posneg"} = '';
    }

    elsif ($element eq 'mixturemodel') {
	$data{'done_models'} .= "$data{'curr_model'} ";
    }

    elsif ($element eq 'msms_run_summary') {
	if ($data{'timestamp_info'}) {
	    $data{"html_table_msruns"} .= "<tr class='thead'><th colspan='2'>Analysis Timestamps</th></tr>\n";
	    $data{'html_table_msruns'} .= delete $data{'timestamp_info'};
	}
    }

    elsif ($element eq 'search_summary') {
	$data{'curr_model'} = "NONE";
	$data{'first_param'} = 1;
	$data{'first_mod'} = 1;
    }

    elsif ( $element eq 'analysis_summary') {
	$data{'curr_model'} = 'NONE';
    }

    elsif ($element eq 'contributing_channel') {
	$data{"html_table_summary"} .= "</td></tr>\n";
    }

}

###################################################
sub processProtXML {
    my $file = shift;

    $data{"tab_html_table_summary"} = "Run Options";
    $data{"html_table_summary"} = "<table style='background:#eeeeee' frame='border' rules='all' cellpadding='2'>\n";

    $data{"tab_html_table_senserr"} = "Sens/Error Table";
    $data{"html_table_senserr"} = "<table style='background:#eeeeee' frame='border' rules='all' cellpadding='2'>\n<tr class='theader'><td colspan='5'>Predicted Sensitivity and Error Rate</td></tr>\n<tr class='thead'>";

    for my $head (qw(min_prob Sensitivity Error_Rate num_correct num_incorrect)) {
	$data{"html_table_senserr"} .= "<th>$head</th>";
    }
    $data{"html_table_senserr"} .= "</tr>\n";

    #### Set up the XML parser and parse the returned XML
    my $parser = XML::Parser->new(
				  Handlers => {
				      Start => \&start_protxml_element,
				      End   => \&end_protxml_element,
#				      Char  => \&protxml_chars,
				  },
				  ErrorContext => 2 );
    eval { $parser->parsefile( $file ); };
    if($@) {
	if ($@ =~ /^Finished parsing header/) {
	    $debug_info .= "Found end of header\n";
	} else { 
	    die "ERROR_PARSING_PROTXML:$@";
	}
    }

}

sub start_protxml_element {
    my ($handler, $element, %atts) = @_;

    if ($element eq 'protein_group') {
	die 'Finished parsing header';  # we're done, so stop reading the file
    }

    elsif ( ($element eq 'isotopic_contributions') ||
	    ($element eq 'dataset_derivation') ||
	    ($element eq 'protein_summary'))
    {
	# ignore
    }

    elsif ( ($element eq 'protein_summary_header') ||
	    ($element eq 'program_details') ||
	    ($element eq 'proteinprophet_details') ) {

	$data{"html_table_summary"} .= "<tr class='thead'><th colspan='2'>$element</th></tr>\n";

	for my $att (sort keys %atts) {
	    $data{"html_table_summary"} .= "<tr><td>$att</td><td>$atts{$att}</td></tr>\n";
	}
    }

    elsif ($element eq 'analysis_summary') {
	$data{"html_table_summary"} .= "<tr class='thead'><th colspan='2'>$atts{'analysis'}</th></tr>\n";

	for my $att (sort keys %atts) {
	    $data{"html_table_summary"} .= "<tr><td>$att</td><td>$atts{$att}</td></tr>\n" unless ($att eq 'analysis');
	}
    }

    elsif ( ($element eq 'ASAP_pvalue_analysis_summary') ||
	    ($element eq 'ASAP_prot_analysis_summary') ||
	    ($element eq 'libra_summary') ||
	    ($element eq 'XPress_analysis_summary') ) {

	for my $att (sort keys %atts) {
	    if ($att eq 'background_ratio_mean') {
		$data{"html_table_summary"} .= "<tr><td>$att</td><td>$atts{$att} (linear ratio: ". sprintf( "%.3f", 10. ** $atts{$att}) .")</td></tr>\n";
	    }
	    else {
		$data{"html_table_summary"} .= "<tr><td>$att</td><td>$atts{$att}</td></tr>\n";
	    }
	}
    }

    elsif ($element eq 'nsp_information') {

	$data{"tab_html_table_NSP"} .= "Learned Models";
	$data{"html_table_NSP"} = "<table style='background:#eeeeee;float:left' frame='border' rules='all' cellpadding='2'>\n<tr class='theader'><td colspan='5'>Learned Number of Sibling Peptide Distributions</td></tr>\n";
	$data{"html_table_NSP"} .= "<tr><td colspan='5'>Neighboring bin smoothing: $atts{'neighboring_bin_smoothing'}</td></tr>\n" if $atts{'neighboring_bin_smoothing'};

	$data{"html_table_NSP"} .= "<tr class='thead'>";
	for my $head (qw(bin_number nsp_range pos_freq neg_freq pos/neg_ratio)) {
	    $data{"html_table_NSP"} .= "<th>$head</th>";
	}
	$data{"html_table_NSP"} .= "</tr>\n";

	$graph_data{"js_data_nsp_pos"} = '[';
	$graph_data{"js_data_nsp_neg"} = '[';
    }

    elsif ($element eq 'nsp_distribution') {
	$data{"html_table_NSP"} .= "<tr id='nsp_bin_$atts{'bin_no'}'><td>$atts{'bin_no'}</td>";

	# range
	$data{"html_table_NSP"} .= "<td>";
	$data{"html_table_NSP"} .= $atts{'nsp_lower_bound_excl'} ? "$atts{'nsp_lower_bound_excl'} &lt; " : "$atts{'nsp_lower_bound_incl'} &lt;= ";
	$data{"html_table_NSP"} .= "nsp ";
	$data{"html_table_NSP"} .= $atts{'nsp_upper_bound_excl'} ? "&lt; $atts{'nsp_upper_bound_excl'}" : "&lt;= $atts{'nsp_upper_bound_incl'}";
	$data{"html_table_NSP"} .= "</td>";

	$data{"html_table_NSP"} .= "<td>$atts{'pos_freq'}</td>";
	$data{"html_table_NSP"} .= "<td>$atts{'neg_freq'}</td>";

	# ratio
	$data{"html_table_NSP"} .= "<td>$atts{'pos_to_neg_ratio'}";
	$data{"html_table_NSP"} .= $atts{'alt_pos_to_neg_ratio'} ? " ($atts{'alt_pos_to_neg_ratio'})" : '';
	$data{"html_table_NSP"} .= "</td>";

	$data{"html_table_NSP"} .= "</tr>\n";

# use these to do vertical bars
#	$data{"js_data_nsp_pos"} .= "[$atts{'bin_no'}, $atts{'pos_freq'}], ";
#	$data{"js_data_nsp_neg"} .= "[$atts{'bin_no'}, -$atts{'neg_freq'}], ";

	# horizontal bars:
	$graph_data{"js_data_nsp_pos"} .= "[$atts{'pos_freq'}, $atts{'bin_no'}], ";
	$graph_data{"js_data_nsp_neg"} .= "[-$atts{'neg_freq'}, $atts{'bin_no'}], ";
    }

    elsif ($element eq 'fpkm_information') {
	$data{"html_table_NSP"} .= "<br style='clear:both'><br/>\n";
	$data{"html_table_NSP"} .= "<table style='background:#eeeeee;float:left' frame='border' rules='all' cellpadding='2'>\n<tr class='theader'><td colspan='5'>Learned FPKM Distributions</td></tr>\n";
	$data{"html_table_NSP"} .= "<tr><td colspan='5'>Neighboring bin smoothing: $atts{'neighboring_bin_smoothing'}</td></tr>\n" if $atts{'neighboring_bin_smoothing'};

	$data{"html_table_NSP"} .= "<tr class='thead'>";
	for my $head (qw(bin_number nsp_range pos_freq neg_freq pos/neg_ratio)) {
	    $data{"html_table_NSP"} .= "<th>$head</th>";
	}
	$data{"html_table_NSP"} .= "</tr>\n";

	$graph_data{"js_data_fpkm_pos"} = '[';
	$graph_data{"js_data_fpkm_neg"} = '[';
    }

    elsif ($element eq 'fpkm_distribution') {
	$data{"html_table_NSP"} .= "<tr><td>$atts{'bin_no'}</td>";

	# range
	$data{"html_table_NSP"} .= "<td>";
	$data{"html_table_NSP"} .= $atts{'fpkm_lower_bound_excl'} ? "$atts{'fpkm_lower_bound_excl'} &lt; " : "$atts{'fpkm_lower_bound_incl'} &lt;= ";
	$data{"html_table_NSP"} .= "fpkm ";
	$data{"html_table_NSP"} .= $atts{'fpkm_upper_bound_excl'} ? "&lt; $atts{'fpkm_upper_bound_excl'}" : "&lt;= $atts{'fpkm_upper_bound_incl'}";
	$data{"html_table_NSP"} .= "</td>";

	$data{"html_table_NSP"} .= "<td>$atts{'pos_freq'}</td>";
	$data{"html_table_NSP"} .= "<td>$atts{'neg_freq'}</td>";

	# ratio
	$data{"html_table_NSP"} .= "<td>$atts{'pos_to_neg_ratio'}";
	$data{"html_table_NSP"} .= $atts{'alt_pos_to_neg_ratio'} ? " ($atts{'alt_pos_to_neg_ratio'})" : '';
	$data{"html_table_NSP"} .= "</td>";

	$data{"html_table_NSP"} .= "</tr>\n";

# use these to do vertical bars
#	$data{"js_data_fpkm_pos"} .= "[$atts{'bin_no'}, $atts{'pos_freq'}], ";
#	$data{"js_data_fpkm_neg"} .= "[$atts{'bin_no'}, -$atts{'neg_freq'}], ";

	# horizontal bars:
	$graph_data{"js_data_fpkm_pos"} .= "[$atts{'pos_freq'}, $atts{'bin_no'}], ";
	$graph_data{"js_data_fpkm_neg"} .= "[-$atts{'neg_freq'}, $atts{'bin_no'}], ";
    }

    elsif ($element eq 'protein_summary_data_filter') {
	$data{"html_table_senserr"} .= "<tr><td>$atts{'min_probability'}</td><td>$atts{'sensitivity'}</td><td>$atts{'false_positive_error_rate'}</td><td>$atts{'predicted_num_correct'}</td><td>$atts{'predicted_num_incorrect'}</td></tr>\n";

	$graph_data{"js_data_ALL_sens"} .= "[$atts{'min_probability'}, $atts{'sensitivity'}], ";
	$graph_data{"js_data_ALL_errr"} .= "[$atts{'min_probability'}, $atts{'false_positive_error_rate'}], ";
    }

    elsif ($element eq 'ASAP_pvalue_analysis_model') {
	$data{"curr_model"} = 'ASAPpvalue';
    }
    elsif ($element eq 'point') {
	if ($data{'done_models'} !~ /$data{"curr_model"}/) {
	    $graph_data{"js_data_pvalue_obs"}   .= "[$atts{'logratio'}, ".$atts{"obs_distr"}.'], ';
	    $graph_data{"js_data_pvalue_model"} .= "[$atts{'logratio'}, ".$atts{"model_distr"}.'], ';
	}
    }

    elsif ($element eq 'fragment_masses') {
	$data{"html_table_summary"} .= "<tr><td>$element</td><td>";
	for my $att (sort keys %atts) {
	    $data{"html_table_summary"} .= "$att: <b>$atts{$att}</b>  ";
	}
	$data{"html_table_summary"} .= "</td></tr>\n";

    }

    elsif ($element eq 'contributing_channel') {
	$data{"html_table_summary"} .= "<tr><td>isotopic contribution <b>$atts{channel}</b></td><td>";
    }
    elsif ($element eq 'affected_channel') {
	$data{"html_table_summary"} .= "<b>$atts{correction}</b>-->$atts{channel}&nbsp;&nbsp;&nbsp;";
    }

    else {
	if ($DEBUG) {
	    $debug_info .= "Did not handle tag '$element'\n";
	}
    }

}

sub end_protxml_element {
    my ($handler, $element) = @_;

    if ( ($element eq 'protein_summary_header') ||
	 ($element eq 'program_details') ||
	 ($element eq 'proteinprophet_details') ) {

#	$data{"html_table_$element"} .= "</table>\n";
    }

    elsif ($element eq 'nsp_information') {
	$data{"html_table_NSP"} .= "</table>\n\n";

	$data{"html_table_NSP"} .=<<"END_JS";
	<SCRIPT LANGUAGE="JavaScript" TYPE="text/javascript">
	    <!--
	    var URLnspbin = getQueryVariable("NSPbin");
	    var URLnspval = getQueryVariable("NSPval");
	    var URLpepseq = getQueryVariable("NSPpep");
	    var URLspchrg = getQueryVariable("NSPchg");
	    var URLprot   = getQueryVariable("NSPprt");
	    var URLnprots = getQueryVariable("NSPind");

	    if (URLpepseq+URLspchrg+URLprot+URLnspval) {
		document.write("<table style='background:#eeeeee;position:relative;left:50;' frame='border' rules='all' cellpadding='2'>");
		document.write("<tr class='thead'><th colspan='2'>Selected Peptide Ion Information</td></tr>");
		if (URLpepseq) document.write("<tr><td>Peptide</td><td><b>" + URLpepseq + "</b></td></tr>");
		if (URLspchrg) document.write("<tr><td>Charge</td><td><b>+" + URLspchrg + "</b></td></tr>");
		if (URLprot) {
		    document.write("<tr><td>Protein</td><td><b>" + URLprot);
		    if (URLnprots != "0") document.write("+" + URLnprots);
		    document.write("</b></td></tr>");
		}
		if (URLnspval) document.write("<tr><td>Calculated NSP</td><td><b>" + URLnspval + "</b></td></tr>");
		document.write("</table>");
	    }

	    if (URLnspbin) {
		document.getElementById("nsp_bin_"+URLnspbin).style.backgroundColor = "#eeaa00";

		var mytab = 'navtab' + (navtabs.length - 1);
		\$(document).ready(function() {
		    document.onload = showPane(mytab);
		});
	    }
	// -->
	</SCRIPT>
END_JS

	$graph_data{"js_data_nsp_pos"} .= '];';
	$graph_data{"js_data_nsp_neg"} .= '];';

    }

    elsif ($element eq 'fpkm_information') {
	$data{"html_table_NSP"} .= "</table>\n\n";
	$graph_data{"js_data_fpkm_pos"} .= '];';
	$graph_data{"js_data_fpkm_neg"} .= '];';

    }

    elsif ($element eq 'ASAP_pvalue_analysis_model') {
	$data{'done_models'} .= "ASAPpvalue ";
    }

    elsif ($element eq 'contributing_channel') {
	$data{"html_table_summary"} .= "</td></tr>\n";
    }

}


###################################################
sub printCSS {
    my $css=<<END;
<style type="text/css">
.hidden {display:none}
.visible{display:block}

table   {border-collapse:collapse; border-color:#000000; font-size: 10pt;}
body    {padding:0px;margin:0px;background-color:#ffffff;font-family: Helvetica, sans-serif; }
.blank  {height:80;background-color:white}
.opaque {position:fixed;width:100%;z-index:9990;height:90;filter:alpha(opacity=85);opacity:.85;background-color:white}
.banner {position:fixed;width:100%;z-index:9995;border-bottom: 5px solid white;padding:0px;margin:0px;font-size:20px;color:#ffffff; background-color:#3B5998}
.thead  {padding:0px;margin:0px;color:#ffffff; background-color:#3B5998}
.theader{font-size:12px;text-align:center;font-weight:bold;background-color:#eeaa00;}
.tgray {border-top: 1px solid black; background: #eeeeee;}

.navtab {
    padding:0px;margin:0px;
    font-size:12px;
    color: #ffffff;
    background:#3B5998;
}
.navtabON {
    padding:6px;margin:0px;
    font-size:12px;
    border-top: 2px solid white;
    border-right: 2px solid white;
    border-left: 2px solid white;
    border-bottom: 0px solid white;
    font-weight: bold;
  background: #ffffff;
  color: #3B5998;
}
.navtab a {
    text-decoration:none;
    color:#ffffff;
  }
.navtab a:hover {
  background:#5B79B8;
  color:#ffffff;
}
.navtabON a {
  text-decoration:none;
  background: #ffffff;
  color: #3B5998;
}



</style>
END

  # add javascript function to create and switch tabs
  $css .=<<'  END_JS';
  <SCRIPT LANGUAGE="JavaScript" TYPE="text/javascript">
  <!--
  var navtabs = new Array();
  var scores = new Object();

  function addPane(tabname) {
      // to close a tabbed pane, simply close a </div> tag
      var index = navtabs.push(tabname) - 1;

      document.getElementById("navtab"+index).innerHTML = "<a href=\"javascript:showPane('navtab" + index + "');\">&nbsp;" + tabname + "&nbsp;</a>";
      document.getElementById("navtab"+index).className = "navtab";

      document.write("<div class=\"hidden\" id=\"navtab"+index+"_content\"><br/>");
  }

  function showPane(divId) {
      var x;
      for (x in navtabs) {
	  var elId = "navtab"+x;
	  document.getElementById(elId).className = "navtab";

	  if (document.getElementById(elId+"_content"))
	      document.getElementById(elId+"_content").className="hidden";
      }

      if (document.getElementById(divId+"_content"))
	  document.getElementById(divId+"_content").className="visible";
      document.getElementById(divId).className="navtabON";
  }

  // taken from http://www.activsoftware.com/
  function getQueryVariable(variable) {
      var query = window.location.search.substring(1);
      var vars = query.split("&");
      for (var i=0;i<vars.length;i++) {
	  var pair = vars[i].split("=");
	  if (pair[0] == variable) {
	      return pair[1];
	  }
      } 
      return false;
  }

  function addScores() {
      var URLspectr = getQueryVariable("Spectrum");
      var URLscores = getQueryVariable("Scores");
      var URLprob   = getQueryVariable("Prob");

      if (URLspectr+URLscores+URLprob) {
	  var chg;

	  addPane('Scores');
	  document.write("<table style='background:#eeeeee' frame='border' rules='all' cellpadding='2'>");
	  document.write("<tr class='thead'><th colspan='2'>Search Hit Info and Scores</th></tr>");

	  if (URLspectr) {
	      document.write("<tr><td><b>Spectrum</b></td><td>" + URLspectr + "</td></tr>");
	      scores['charge'] = URLspectr.charAt(URLspectr.length-1);  // assume single digit for now...
	  }

	  if (URLscores) {
//	      document.write("<tr><td><b>Scores</b></td><td>" + decodeURIComponent(URLscores) + "</td></tr>");

	      var vars = decodeURIComponent(URLscores).split(" ");
	      for (var i=0;i<vars.length;i++) {
		  var pair = vars[i].split(":");
		  document.write("<tr><td><b>" + pair[0] + "</b></td><td>" + pair[1] + "</td></tr>");

		  if (pair[0] == 'fval') {
		      scores['fval'] = pair[1];
		  }
	      }
	  }

	  if (URLprob) document.write("<tr class='theader'><td><b>Prob</b></td><td>" + URLprob + "</td></tr>");

	  document.write("</table>");
	  document.write("</div>");
      }
  }
  // -->
  </SCRIPT>
  END_JS

    return $css;
}


###############################################################################
# addTabbedPane: adds javascript code to add tab entry and opens corresponding div
#
# returns javascript in a scalar
# 
###############################################################################
sub addTabbedPane {
  my %args = @_;

  $args{label} ||= 'tab';
  my $pad = '&nbsp;';

  my $dhtml =<<"  END_JS";
  <SCRIPT LANGUAGE="JavaScript" TYPE="text/javascript">
  <!--
      addPane('$args{label}');
  // -->
  </SCRIPT>
  END_JS

  return $dhtml;

}

###############################################################################
# closeTabbedPane: close div; add hr; select if requested
#
# returns html in a scalar
# 
###############################################################################
sub closeTabbedPane {
  my %args = @_;

  $args{selected} ||= '';

  my $dhtml = "</div>\n";

  if ($args{selected}) {
      $dhtml .=<<"      END_DHTML";
      <SCRIPT LANGUAGE="JavaScript" TYPE="text/javascript">
       <!--
          var mytab = 'navtab' + (navtabs.length - 1);
          document.onload = showPane(mytab);
      // -->
       </SCRIPT>
      END_DHTML
  }

  return $dhtml;

}


sub getFileFormat {
    my $file = shift || die "No file to find format!!";
    my $fmt = 'UNKNOWN';

    $debug_info .= "File $file is ";
    open(FILE, "$file") || return $fmt;

    my $lc = 0;
    while (<FILE>){
        if (/msms_pipeline_analysis/) {
            # only consider versions 1.XXX for now
            $fmt = 'pepxml' if ($_ =~ /pepXML_v1.*\.xsd/);
            last;
        }
        if (/protein_summary/) {
            $fmt = 'protxml';
            last;
	}
        last if ($lc++ > 5);
    }
    close(FILE);

    $debug_info .= "$fmt\n";
    return $fmt;
}

sub printWithDots {
    my $stringLeft  = shift || '';
    my $stringRight = shift || '';
    my $eol = shift() ? "\n" : '';
    my $dots = 38 - length($stringLeft);
    $dots = 3 if ($dots < 3);

    return $stringLeft . "."x$dots . $stringRight . $eol;
}

sub reportError {
    my $errstring = shift;
    my $fatal = shift || 0;

    print "ERROR [" . ++$errors. "]: $errstring\n";
    exit($fatal) if $fatal;
}
