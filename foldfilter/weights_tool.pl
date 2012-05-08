#! /usr/bin/perl -w

#
#  $Id: weights_tool.pl,v 1.3 2011-11-02 20:45:12 laci Exp $
#
#  Project      : RNA motif searching in genomic sequences
#  Description  : a script helping to estimate suitable weights of submotifs 
#                 used for seqence scoring in fold_filter
# 
#  Author       : L. Rampasek <rampasek@gmail.com>
#  Institutions : Comenius University in Bratislava,
#                 University of California, Irvine
#

use strict;
use warnings; 
use Filter;
use Getopt::Long;
Getopt::Long::config("no_ignore_case");
use POSIX qw/ceil/;
use Scalar::Util qw(looks_like_number);


#set dumper for dumping nested data structures
require Dumpvalue;
my $dump = new Dumpvalue::
			    tick        => q("),
			    compactDump => 1,  # comment these two lines out
			    veryCompact => 1,  # if you want a bigger dump
			;
			
my $COMMANDS = "Commands:
  reweight <num>:<num>: ... :<num>
      set new ratio of submotifs' coefficients; the order of submotifs in the
      command is equal to their order in listing

  set coef <sm_name> <new_coef>
      set multiplication coefficient for the submotif <sm_name> to be <new_coef>
      
  export <num>% <file_name>
      export sequences which score is equal or larger than <num>% of maximum
      score; save them to the file <file_name>
 
  draw
      draw a histogram of scores according to current coefficients
      
  set range <left> <right>
      set the histogram range to be form <left>% to <right>% of max score;
      by default it is from 0% to 100%
       
  set width <num>
      set the width of histograms (number of columns); default width is 50

  set height <num>
      set the max height of histogram (number of rows); default height is 20
      
  help
      print this list of commands
   
  quit, exit, or Ctrl+D
      exit the program
";

my $USAGE = "
USAGE: perl $0 [OPTIONS] <fold_filter_seq_file>

  A script helping to estimate suitable weights of submotifs used for seqence
  scoring in fold_filter.

Available options:

  -h, --help  show this help

Example usage:
  perl $0 scored_seq.txt
  
$COMMANDS
";

our $tmpDir = "/tmp/";
$Filter::tmpDir = $tmpDir;

my $showHelp = 0;
GetOptions( 
            'help' => \$showHelp ) 
      || die $USAGE;
      
if ($showHelp) { die $USAGE; }

#parse command line arguments
#print join ("|", @ARGV), "\n";
my $seqFileName;
if (scalar @ARGV == 1){
    $seqFileName = shift @ARGV;
} else {
    die $USAGE;
}

################################################################################
sub parseSeqFile {
  my $seqFileName = shift @_;
  my $numOfElements = shift @_;
  
  my @data = ();
  my @scores = ();
  my @smWeights= ();
  open (IN, "<", $seqFileName) or die "Problem opening sequence file $seqFileName: $!";
  
  while (my $line = <IN>){
      $line = Filter::trim($line);
      #if it is a header of a sequence, then read & store it
      if ($line =~ /^>/){
          #save the header
          push @data, $line;
          
          #read & save the sequence
          $line = <IN>;
          $line = Filter::trim($line);
          push @data, $line;
          
          ### parse all submotifs and their scores
          my @sm = ();
          while ($line = <IN>){
              last if $line =~ /^====*/;
              #parse a submotif
              if ($line =~ /(sm[^:\s]+)/gi){
                  my $smName = $1;
                  $smWeights[@sm] = [$smName] if !defined $smWeights[@sm];
                  die("Inconsistency in submotifs listing") if $smWeights[@sm]->[0] ne $smName;
                  
                  my $smUnweighted = -1;
                  my $smScore = -1;
                  
                  #parse scores
                  while ($line = <IN>){
                      last if $line =~ /^----*/;
                      if ($line =~ /^\s*unweighted\s+score\s+=\s+(\d+\.?\d*)/gi){
                          $smUnweighted = $1;
                      }
                      if ($line =~ /^\s*weighted\s+score\s+=\s+(\d+\.?\d*)/gi){
                          $smScore = $1;
                      }
                  }
                  
                  #save the values
                  if ($smUnweighted == -1 or $smScore == -1){
                      die("All scores for submotif $smName were not found.");
                  }
                  if ($smUnweighted > 0){
                    my $smWeight = $smScore*1. / $smUnweighted;
                    $smWeights[@sm]->[1] += $smWeight;
                  } else { $smWeights[@sm]->[1] += $smWeights[@sm]->[1] / @scores; }
                                    
#                  print "$smName: $smScore, $smUnweighted, $smWeight\n";
                  push @sm, [$smScore, $smUnweighted];
              }
          }
          push @scores, [getSeqScore(\@sm), $#data, \@sm];
      }
  }
  close (IN);
  
  #average the weights
  foreach (@smWeights){
      $_->[1] = $_->[1]*1. / @scores;
  }
  
#  $dump->dumpValue(\@data);
#  $dump->dumpValue(\@smWeights);
#  $dump->dumpValue(\@scores);
  return (\@data, \@smWeights, \@scores);
}

sub getSeqScore {
    my @sm = @{ shift @_ };
    
    my $score = 0;
    foreach (@sm){
        $score += $_->[0];
    }
    
    return $score;
}

sub recalculateSeqScore {
    my @sm = @{ shift @_ };
    my @weights = @{ shift @_ };
    
    foreach (0..$#sm){
        $sm[$_]->[0] = $sm[$_]->[1] * $weights[$_]->[1];
    }
}

sub recalculateAllScores {
    my @scores = @{ shift @_ };
    my $weights = shift @_;
    
    foreach (@scores){
        recalculateSeqScore($_->[2], $weights);
        $_->[0] = getSeqScore($_->[2]);
    }
}

sub computeDistribution {
    my @scores = @{ shift @_ };
    my $numBins = shift @_;
    my $scale = shift @_;
    my $left = shift @_;
    my $right = shift @_;
    
    #if we want N bins, they will have indices 0..(N-1)
    $numBins -= 1;
    
    my @bins = ();
    $bins[$_] = 0 foreach (0..$numBins);
    
    #find max value
    my $maxVal = -1;
    foreach (@scores){
        $maxVal = $_->[0] if $maxVal < $_->[0];
    }
#    print $maxVal,"\n";

    #divide seqences to the bins
    foreach (@scores){
        my $perc = ($_->[0]*100.0)/$maxVal;
        next if $perc < $left or $perc > $right;
        $bins[ int( ($perc-$left)/($right-$left) * $numBins ) ]++;
    }
    
    #scale bins' occupancy relative to total number of sequences
#    foreach (@bins){
#        $_ = int ($_ / @scores *$scale);
#    }
    
    #scale bins' occupancy relative to the most occupied bin
    #find the mostpopulated bin
    my $maxBin = 1;
    foreach (@bins){
        $maxBin = $_ if $maxBin < $_;
    }
    
    foreach (@bins){
        $_ = int ($_ / $maxBin *$scale);
    }
    
#    $dump->dumpValue(\@bins);
    return \@bins;
}

sub drawHistogram {
    my @distrib = @{ shift @_ };
    my $left = shift @_;
    my $right = shift @_;
    
    #find max value
    my $max = 0;
    foreach (@distrib){
        $max = $_ if $max < $_;
    }
    
    for (my $i = $max; $i >= 0; $i--){
        foreach (@distrib){
            print "^" if $i == $_ ;
            print "|" if $i < $_ ;
            print " " if $i > $_ ;
        }
        print "\n";
    }
    
    my $prevPercInt = $left;
    foreach (1..@distrib){
        if($_ % 5 == 0){
            my $perc = $left+ (($_*1.0)/@distrib)*($right-$left);
            my $percInt = int($perc);
            
            if($percInt != $prevPercInt){
              printf "%4.0f%s", $percInt , "%";
            } else {
              my $tmp = sprintf "%.2f", ($perc-$percInt);
              my @chars = split "", $tmp;
              $chars[0] = ' ';
              
              if ($chars[$#chars] eq '0') {
                  print " ", (join "", @chars[0..$#chars-1]), "%";
              } else {
                  print "", (join "", @chars[0..$#chars]), "%";
              }
            }
            
            $prevPercInt = $percInt;
        }
    } print "\n";
}

sub printCoefsAndWeights {
    my @smWeights = @{ shift @_ };
    
    my $cmin = 999999999;
    foreach (@smWeights){
        $cmin = $_->[1] if $_->[1] < $cmin;
    }
    
    foreach (@smWeights){
        printf "%12s -> %2.5f  (%2.2f)\n", $_->[0], $_->[1], $_->[1]*1. / $cmin;
    }
    print "\n";
    
}

################################################################################

my @result = parseSeqFile($seqFileName);
my @data = @{ $result[0] };
my @smWeights = @{ $result[1] };
my @scores = @{ $result[2] };

my $height = 20;
my $width = 50;
my $left = 0;
my $right = 100;
my $distrib = computeDistribution(\@scores, $width, $height, $left, $right);

drawHistogram($distrib, $left, $right);
print "Used coeficients (ratio):\n";
printCoefsAndWeights(\@smWeights);

#$dump->dumpValue(\@smWeights);
#$dump->dumpValue(\@scores);

my $command;
while (my $command = <>){
    my @par = split ' ', $command;
    if( $par[0] eq "set"){  ####################################################
        if( @par < 3){ 
            print "Not enough parameters for the command 'set'\n";
            next;
        }
        if( $par[1] eq "coef"){ ######==========================================
            if (@par != 4){ 
                print "The command 'set coef' expects two more parameters - submotif name and numeric value\n";
                next;
            }
            if (!looks_like_number($par[3])){
                print "The given parameter \"".$par[3]."\" is not a number.\n";
                next;
            }
            
            my $done = 0;
            foreach (@smWeights){
                if($_->[0] eq $par[2]){
                    $_->[1] = $par[3];
                    $done = 1;
                    last;
                }
            }
            recalculateAllScores(\@scores, \@smWeights);
            printCoefsAndWeights(\@smWeights);
            
            print "Unknown submotif \"".$par[2]."\"\n" if (!$done);
        } elsif( $par[1] eq "width"){ ######====================================
            if (looks_like_number($par[2])){
              $width = $par[2];
            } else {
              print "The given parameter \"".$par[2]."\" is not a number."
            }
        } elsif( $par[1] eq "height"){ ######===================================
            if (looks_like_number($par[2])){
              $height = $par[2];
            } else {
              print "The given parameter \"".$par[2]."\" is not a number."
            }
        } elsif( $par[1] eq "range"){ ######====================================
            if (looks_like_number($par[2]) and $par[2]>=0 and $par[2]<=100
                and looks_like_number($par[3]) and $par[3]>=0 and $par[3]<=100
                and $par[2]!=$par[3]){
              $left = $par[2];
              $right = $par[3];
            } else {
              print "The given parameters for 'set range' are not admissible."
            }
        }else { ######=========================================================
            print "Unknown parameter \"".$par[1]."\" for 'set' command; type 'help' for the list of commands\n";
        }
        
    } elsif ($par[0] eq "draw"){ ###############################################
        $distrib = computeDistribution(\@scores, $width, $height, $left, $right);
        drawHistogram($distrib, $left, $right);
        print "Used coeficients (ratio):\n";
        printCoefsAndWeights(\@smWeights);
        
    } elsif ($par[0] eq "reweight"){ ###########################################
        if (!defined $par[1]){
            print "The command 'reweight' expects one more parameter\n";
             next;
        }
        
        my @newWeights = split ':', (join "", @par[1..$#par]);
        if (@newWeights != @smWeights){
            print "In command 'reweight', weights must be specified for every submotif\n";
            print "The command 'reweight' expects weight ratio specification in format <num>:<num>: ... :<num>\n";
            next;
        }
        
        my $wsum = 0;
        foreach (@newWeights){
            if( !looks_like_number($_) ) {
                print "In command 'reweight', weights must be a numeric value\n";
                next;
            }
            $wsum += $_;
        }
        foreach (0..$#newWeights){
            $newWeights[$_] *= 1 / $wsum;
            $newWeights[$_] *= 10;
            $smWeights[$_]->[1] = $newWeights[$_];
        }
        recalculateAllScores(\@scores, \@smWeights);
        printCoefsAndWeights(\@smWeights);
        
    } elsif ($par[0] eq "export"){ #############################################
        if( @par != 3){ 
            print "Not enough parameters for the command 'export'\n";
            next;
        }
        
        my $percentTH = undef;
        if( $par[1] =~ /([\d,\.]+)%/ ){
            $percentTH = $1  if looks_like_number($1);
            $percentTH = undef if $percentTH < 0 or $percentTH > 100;
        }
        if (! defined $percentTH){
            print "In command 'export': \"".$par[1]."\" is not a valid percentage value\n";
            next;
        }
        
        my $exportFileName = $par[2];
        if (! open (OUT, ">", $exportFileName)) {
          print "Problem opening the file for export '$exportFileName': $!";
          next;
        }
        
        @scores = sort {$b->[0] <=> $a->[0]} @scores;
        my $maxScore = $scores[0]->[0];
        my $threshold = $maxScore * ($percentTH/100.);
        for (my $i = 0; $i < @scores and $scores[$i]->[0] >= $threshold; $i++){
            print OUT $data[$scores[$i]->[1]-1], "\n";
            print OUT $data[$scores[$i]->[1]], "\n";
        }
        close (OUT);
        print "  export complete\n";	
        
    } elsif ($par[0] eq "help"){ ###############################################
        print $COMMANDS;
        
    } elsif ($par[0] eq "exit" or $par[0] eq "quit"){ ##########################
        exit;
        
    } else {  ##################################################################
        print "Unknown command \"".$par[0]."\"; type 'help' for the list of commands\n";
    }
}


