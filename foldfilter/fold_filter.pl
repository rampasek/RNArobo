#! /usr/bin/perl -w

#
#  $Id: fold_filter.pl,v 1.13 2011-12-16 22:47:51 laci Exp $
#
#  Project      : RNA motif searching in genomic sequences
#  Description  : a script enabling to filter RNArobo/RNAbob output by testing
#                 the secondary structure stability by RNAcofold, RNAfold, or
#                 DotKnot
# 
#  Author       : L. Rampasek <rampasek@gmail.com>
#  Institutions : Comenius University in Bratislava,
#                 University of California, Irvine
#

#-------------------------------------------------------------------------------
#TODO:
# - tweak species name parsing
# - add comments to procedures in Filter module
# - test accross platforms
#
#-------------------------------------------------------------------------------

use strict;
use warnings; 
use Filter;
use Getopt::Long;
Getopt::Long::config("no_ignore_case");
use POSIX qw/ceil/;

#set dumper for dumping nested data structures
require Dumpvalue;
my $dump = new Dumpvalue::
			    tick        => q("),
			    compactDump => 1,  # comment these two lines out
			    veryCompact => 1,  # if you want a bigger dump
			;

my $USAGE = "
USAGE: perl $0 [OPTIONS] <descriptor_file> <rnarobo_output_file> <config_file>\n

  A script enabling to filter RNArobo/RNAbob output by testing
the secondary structure stability by RNAcofold, RNAfold, or DotKnot, and further
by culstering (by Uclust) and filtering according to a region of interest.
  You need to have ViennaRNA package installed! (if you use version 2.0.0+, use 
the option --vienna2)

Available options:
--score     show overall score for every sequence

--full      show full score listing for every sequence
  
--anotate   anotate brackets by letters a..z, A..Z according to 
            the level of depth (works only in combination with --full option)
              
--cutoff <num>[%]   set cutoff score for the sequences on the output; if a 
            number <num> is given, it is treated as the minimal absolute score
            threshold; else if a precentage <num>% is given, then top <num>% of
            seqences is taken
  
--seqs <file> store the secondary structure filter output to the external <file>
              e.i. store results of the first part of the pipeline; cutoff
              threshold is applied (if set)
        
--only1part   execute only the first part of the pipeline (folding and scoring
              according to submotifs, and sequence cutoff),  then halt 
              
              
  
--continue <file>    use only culster filter on the sequences from <file> (i.e.
            the second part of the filtering pipeline), <file> must be in fasta
            format with space separated elements; the cutoff threshold cannot be
            applied on these sequences
            
--clusters <file>   store the result of UClust clustering into <file>
            
--minsize <num>     set threshold for minimum number of sequences in a cluster
                    default value = 8

--minspecies <num>  set threshold for minimum number of species in a cluster
                    default value = 3
              
              
              
-p          show progress during computational steps
".
#--vienna2   use this option for compatibility with ViennaRNA version 2.0.0+
"
-h, --help  show this help


Example usage:
  perl $0 sample.des sample.out filter.ftr --full --anotate
  perl $0 sample.des sample.out filter.ftr --out seq.txt --full
        --anotate --cutoff 47%
  perl $0 sample.des filter.ftr -p --continue seq.txt
        --clusters clust.txt --minsize 10 --minspecies 2
";

#parse command line arguments
die $USAGE unless @ARGV>=3;

our $tmpDir = "/tmp/";
$Filter::tmpDir = $tmpDir;

my $showScores = 0;
my $showCompleteScores = 0;
my $showAnotation = 0;
my $showHelp = 0;
my $dumpFileName = undef;
my $fastaFileName = undef;
my $cutoffScore = undef;
my $cutoffIsPercentage = 0;
my $showProgress = 0;
my $clusterResFAFileName = undef;
my $useViennaRNA2 = 0;
my $executeOnlyFirstPart = 0;
#thresholds for cluster filtering, default: 8 for count, 3 for species
my $sizeThreshold = 8;
my $spiecesThreshold = 3;
GetOptions( 'score' => \$showScores,
            'full' => \$showCompleteScores,
            'anotate' => \$showAnotation,
            'p' => \$showProgress,
            'cutoff=s' => \$cutoffScore,
            'seqs=s' => \$dumpFileName,
            'only1part' => \$executeOnlyFirstPart,
            'continue=s' => \$fastaFileName,
            'clusters=s' => \$clusterResFAFileName,
            'minsize=i' => \$sizeThreshold,
            'minspecies=i' => \$spiecesThreshold,
            #'vienna2' => \$useViennaRNA2,
            'help' => \$showHelp ) 
      || die $USAGE;
      
if ($showHelp) { die $USAGE; }

if (defined $cutoffScore and $cutoffScore !~ /^(\d+)(%?)$/ ) { die $USAGE; }
if (defined $2 and $2 eq "%"){
    $cutoffScore = $1;
    $cutoffIsPercentage = 1;
    die "The cutoff percentage must be an integer from 0 to 100" 
      if $cutoffScore >100 or $cutoffScore <0;
}

#find out which version of ViennaRNA package is installed ( >= 2.0.0 or <)
#if the return value is != 0, then ViennaRNA is not installed or the version is < 2.0.0
$Filter::useViennaRNA2 = $useViennaRNA2;
my $viennaVersion = `RNAcofold --version 2> /dev/null`;
if ( $? == 0 ){ 
    $useViennaRNA2 = 1;
    $Filter::useViennaRNA2 = $useViennaRNA2;
}
    
$Filter::showProgress = $showProgress;
if ($showProgress){
    print STDERR "-------- Starting Fold Filter --------\n\n";
    print STDERR "used settings:\n";
    printOptBool("'score'", $showScores);
    printOptBool("'full'", $showCompleteScores);
    printOptBool("'anotate'", $showAnotation);
    printOptBool("'p'", $showProgress);
    printOptVal("'cutoff'", $cutoffScore.("%" x $cutoffIsPercentage));
    printOptVal("'seqs'", $dumpFileName);
    printOptBool("'only1part'", $executeOnlyFirstPart);
    printOptVal("'continue'", $fastaFileName);
    printOptVal("'clusters'", $clusterResFAFileName);
    printOptVal("'minsize'", $sizeThreshold);
    printOptVal("'minspecies'", $spiecesThreshold);
    printOptBool("'vienna2'", $useViennaRNA2);
    print STDERR "\nProgress is shown only for the main computational procedures\n";
}

sub printOptBool {
  my $comment = shift @_;
  my $val = shift @_;
  if($val) {$val = "true";} else {$val = "false";}
  print STDERR "  $comment => $val\n"; 
}

sub printOptVal {
  my $comment = shift @_;
  my $val = shift @_;
  if(!defined $val) {$val = "undef";}
  print STDERR "  $comment => $val\n"; 
}

#print join ("|", @ARGV), "\n";
my ($descFilename, $dataFilename, $fltFilename );
if (scalar @ARGV == 3){
    ($descFilename, $dataFilename, $fltFilename ) = @ARGV;
} elsif (scalar @ARGV == 2 and defined $fastaFileName){
    ($descFilename, $fltFilename ) = @ARGV;
} else {
    die $USAGE;
}

sub trim {
    my($str) = shift =~ m!^\s*(.+?)\s*$!i;
    defined $str ? return $str : return '';
} 

########################################################### parse motif map
#read motif map from descriptor
my @motifMap = @{ Filter::parseMotifMap( $descFilename ) };
#print join ("|", @motifMap), "|\njeeh\n";
 
########################################################### parse filter config
my @result = Filter::parseFilterConfig( $fltFilename, \@motifMap );
my @filter = @{ $result[0] };
my $roi = $result[1];
#$dump->dumpValue(\@filter);

######################## check consistency of the filter config
Filter::checkFilterConfig( \@filter );


my @data = ();

if (defined $fastaFileName){
########################################################### parse given fasta file
  @data = @{ Filter::parseFastaFile($fastaFileName, scalar @motifMap) };
    
} else {

########################################################### parse RNArobo output
@data = @{ Filter::parseRNAroboOutput($dataFilename, scalar @motifMap) };

########################################################### score every sequence
my @seqScores = @{ Filter::computeScores(\@motifMap, \@data, \@filter) };

#sort sequences according to the optained score
@seqScores = sort {$b->[1] <=> $a->[1]} @seqScores;


################################### print out outputs & apply the cutoff threshold
##determine output stream
if (defined $dumpFileName){
    open (DUMP, ">", $dumpFileName) or die "Problem opening dump file $dumpFileName: $!";
    $showCompleteScores = 1;
} else {
    open (DUMP, ">&STDOUT") or die "Can't write to STDOUT: $!";
}

## print & cutoff
my $counter = 0;
my $prevScore = -1;
my @filteredData = ();
foreach my $ref (@seqScores){
    #stop if the sequence score is below cutoff threshold
    if (defined $cutoffScore){
        if ($cutoffIsPercentage == 0){
            last if ($ref->[1] < $cutoffScore);
            
        #take only top $cutoffScore percent of the sequences (if specified)    
        } elsif ($cutoffIsPercentage == 1){ 
            if ($counter > ceil(@seqScores*$cutoffScore/100) ){
                last if $ref->[1] != $prevScore; #print also those on the edge with the same score
            }
        }
    }
    ++$counter;
    $prevScore = $ref->[1];
    push @filteredData, $data[ $ref->[0] -1 ];
    push @filteredData, $data[ $ref->[0] ];
    
    #sequence name
    print DUMP $data[ $ref->[0] -1 ]."\n";
    
    #sequence itself
    my @seq = split " ", $data[ $ref->[0] ];
    my $separator = ($showCompleteScores) ? " " : "";
    print DUMP join($separator, @seq), "\n";
    
    #scores
    if ($showCompleteScores){
        foreach my $submotif (@filter){
            #reference to intermediate score values
            my $submotifScores = $ref->[2]->{$submotif->{'name'}};
            
            my $strmap;
            my @map = split " ", $submotifScores->{'map'};
            my $strstruct;
            my @struct = split " ", $submotifScores->{'structure'};
            my @seq = split " ", $submotifScores->{'sequence'};
            
            #if nested, mark the middle by '&'
            if ($submotifScores->{'type'} eq 'nested'){
                @struct[@struct/2-1] .= " &";
                @seq[@seq/2-1] .= " &";
            }
            
            #indent elements' names according to their seqence/structure length 
            for (my $i=0; $i < @struct; ++$i){
                $strmap .= $map[$i]." ";
                $strstruct .= $struct[$i]." ";
                while ( length($strstruct) > length($strmap) ){ $strmap .= " ";}
            }
            
            printf DUMP "%-10s: %s\n            %s   mfe = %.2f kcal/mol\n",
                $submotif->{'name'},
                join (" ", @seq),  #sequence
                join (" ", @struct),  #structure
                $submotifScores->{'min_en'};  #min. free energy
            
            #show bracket depth anotation if the option is set
            if ($showAnotation){
                my @anotatedSeq = split "", (join (" ", @struct));
                my $level1 = 0;
                my $level2 = 0;
                my @chars = ( 'a'..'z', 'A'..'Z' );
                my @nums = ( '0'..'9' );
                
                for (my $i = 0; $i < scalar @anotatedSeq; ++$i){
                    if ($anotatedSeq[$i] eq "(" ){
                        $anotatedSeq[$i] = $chars[$level1 % scalar @chars];
                        ++$level1;
                    } elsif ($anotatedSeq[$i] eq "[" ){
                        $anotatedSeq[$i] = $nums[$level2 % scalar @nums];
                        ++$level2;
                    } elsif (($anotatedSeq[$i] eq ")" )){
                        --$level1;
                        $anotatedSeq[$i] = $chars[$level1 % scalar @chars];
                    } elsif (($anotatedSeq[$i] eq "]" )){
                        --$level2;
                        $anotatedSeq[$i] = $nums[$level2 % scalar @nums];
                    }
                }
                
                printf DUMP "            %s\n", join ("", @anotatedSeq);
            }
            
            
            printf DUMP "            %s\n",
                $strmap;
                
            my $tmp_type = $submotifScores->{'type'};
            $tmp_type .= ", MFE*10 per nucleotide" if $submotif->{'per_nucleotide'};
            printf DUMP "            type = %s\n", $tmp_type;
                
            printf DUMP "            tail length = %s\n",
                $submotif->{'tail'} if $submotif->{'tail'} ne "";
            
            printf DUMP "            unweighted score = %.2f\n",
                abs($submotifScores->{'unw_score'});
                
            printf DUMP "            weighted score   = %.2f\n",
                abs($submotifScores->{'score'});
                
            print DUMP "--------------------------------------------\n";
        }
    }
    
    printf (DUMP "Overall sequence score: %.2f\n", $ref->[1]) if $showScores or $showCompleteScores;
    printf DUMP "=================================================\n\n" if $showCompleteScores;
}
close (DUMP);
@data = @filteredData;

print STDERR "end of the first part\n" if $showProgress;
}


################# process the region of interest -> RNAfold, fake AA code, cluster, filter
if (defined $roi and !$executeOnlyFirstPart){
    ##extract ROIs from all the sequences
    my @rois = @{ Filter::extractROIs(\@motifMap, \@data, $roi) };
    
    ##cluster the ROIs
    my $deleteClusterFile = 0;
    if (!defined $clusterResFAFileName){
        $clusterResFAFileName = $tmpDir."_tmp_clust_res.fa";
        $deleteClusterFile = 1;
    }
    Filter::clusterROIs(\@rois, \@data, $clusterResFAFileName);
    
    ##filer out insignificant clusters
    my @filteredHeaders = @{ Filter::filterClusters($clusterResFAFileName, $sizeThreshold, $spiecesThreshold) };

#    $dump->dumpValue(\@filteredHeaders);
    unlink($clusterResFAFileName) if $deleteClusterFile;

    #print out FINAL results    
#    print "dokopy: ", scalar @filteredHeaders, "\n";
    foreach my $clusterRef (@filteredHeaders){ 
#      print scalar @{$clusterRef}, "\n";
      foreach my $header (@{$clusterRef}){
        $header =~ s/^(\d+)\s?/>/;
        my $id = $1;
        
        print $header, "\n", $data[ $id+1 ],"\n";
      }
    }
    
    print STDERR "end of the second part\n" if $showProgress;
} else {
    #print "ROI was not defined!\n";
}

