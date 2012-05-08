#! /usr/bin/perl -w

#
#  $Id: parse_rnarobo.pl,v 1.8 2011-09-15 19:05:32 laci Exp $
#
#  Project      : RNA motif searching in genomic sequences
#  Description  : a script to convert RNAbob/RNArobo output into FASTA format
# 
#  Authors      : T. Vinar, B. Brejova, R. Jimenez, L. Rampasek
#  Institutions : Comenius University in Bratislava,
#                 University of California Irvine
#

use strict;
use Getopt::Long;

my $USAGE = "
USAGE: perl $0 [OPTIONS] [-i <input_file>] <which_elements> \n

Parse RNAbob/RNArobo output stored in <input_file> or given to the standard
input to fasta format containing the found regions.
If run without any <which_elements> then it outputs all elements.
Otherwise, it selects only those elements specified on the input
(elements are numbered 1,2,... in the order in which they appear
on the first line of the descriptor).

Available options:
  -i <file>   specify path to the input file
  
  -c          number the individual sequences according to their order in input
   
  -f          output full sequence names as thay appear in input
  
  -s          print individual sequence elements separated by space 
  
  -n <number> consider only sequences that contain up to <number> \"N\"s,
              if no <number> is given, than it is set to zero
  
  -a          insert \"&\" into the middle of the elements' list
              WARNING: the program DOESNT check whether the elements 
                       are nested helices
                       
  -h, --help  show this help
  
Example usage:
  perl $0 -f -n 5 -i data.txt 1 3 4 6 > data.fa
  perl $0 -a -c 1 3 < sample.out > sample.fa
  
";

#parse command line arguments
#die $USAGE unless @ARGV>=1;

#parse commandline parameters - input filename, "insert ampersand", ...
my $filename = "";
my $showCounter = 0;
my $showFullName = 0;
my $separateElems = 0;
my $nThreshold = -1;
my $insertAmp = 0;
my $showHelp = 0;
GetOptions( 'i=s' => \$filename,
            'c' => \$showCounter,
            'f' => \$showFullName,
            's' => \$separateElems,
            'n:i' => \$nThreshold,
            'a' => \$insertAmp,
            'help' => \$showHelp ) 
      || die $USAGE;
      
#list of elements
my @which = @ARGV;

if ($showHelp) { die $USAGE; }

#check that all elements of @which are positive integers
foreach my $num (@which) {
    die "Wrong which_element '$num'\n\n$USAGE" unless $num=~/^[0-9]+$/ && $num>0;
}

our $file;
open $file, "<$filename" or die "Cannot open $filename: $!" unless !$filename;

my $active = 0;
my $counter = 1;
while(my $line = get_line() ) {
    chomp $line;
    if ($active) {
        if ($line) {
            #skip this line if it doesnt look like new occurrence
            next unless $line =~ /^\d/;
        
            #read another line with sequence
            my $content = get_line();
            chomp $content;
            
            #if $nThreshold is set, apply the "N" filter
            if ($nThreshold > -1) {
                my $num = $content =~ tr/[n,N]//;
                next if ($num > $nThreshold);
            }
            
            #parse line and print fasta header
            my $countPrefix = ($showCounter) ? ($counter++)." " : "";
            print "> ", $countPrefix, make_header($line), "\n";
            
            my @parts = split '\|',$content;
              
            #first part should be empty
            die "Wrong format" unless $parts[0] eq "";

            #select parts of the output of interest
            my @subparts;
            if (@which == 0) {
                @subparts = @parts;
                shift @subparts;
            } else {
                foreach my $col (@which) {
                die "Part $col not found" unless exists $parts[$col];
                    push @subparts, $parts[$col];
                }
            }
       
            #prints selected parts
            my $separator = ($separateElems) ? " " : ""; 
            if ($insertAmp) {
                #die "Odd number of elements in the list" unless @subparts%2==0;
                print join($separator, @subparts[0 .. @subparts/2 - 1]),
                      $separator, "&", $separator,
                      join($separator, @subparts[@subparts/2 .. @subparts-1]), "\n";
            } else {
                print join($separator, @subparts), "\n";
            }
        } else {
            # empty line
            $active = 0;
        }
          
    } else {
        #start of occurrences, end of file header
        if ($line =~ /^seq-f/) {
            # next line should contain dashed on top of the table
            $line = get_line();
            die "Unexpected end of line" unless defined $line;
            die "Wrong format" unless $line =~ /^[- \t]*\s*$/;
            $active = 1;
        }
    }
}

#erase trailing whitespaces from the string
sub trim($){
	  my $string = shift;
	  chomp $string;
	  $string =~ s/^\s+//;
	  $string =~ s/\s+$//;
	  return $string;
}

#get next unempty line from the input
sub get_line {
    my $str = 0;
    do {
        $str = ($file) ? <$file> : <STDIN>;
        $str = trim($str) if defined $str;
    } while(defined $str and $str eq "");
    
    return $str;
}

sub make_header {
    my ($line) = @_;
    
    #if the option "show full sequence name" is turned on
    if ($showFullName) { return $line; }
    
    #split line into three parts at whitespace
    my @parts = split " ",$line;
    die "Wrong format" unless @parts>=3;

    #if header contains gi numbers, take only the number
    my $header;
    if ($parts[2] =~ /^gi\|[^\|]*\|[^\|]*\|([^\|]*)\|/) {
        $header = $1;
    } else {
        $header = $parts[2];
    }
    
    my $moreInfo = "";
    $moreInfo .= (exists $parts[3]) ? " ".$parts[3] : "";
    $moreInfo .= (exists $parts[4]) ? " ".$parts[4] : "";
    
    return $parts[0]."-".$parts[1]." ".$header.$moreInfo;
}

