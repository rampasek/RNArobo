#
#  $Id: Filter.pm,v 1.4 2011-12-16 22:47:51 laci Exp $
#
#  Project      : RNA motif searching in genomic sequences
#  Description  : a procedural Perl module for fold and cluster filtering
# 
#  Author       : L. Rampasek <rampasek@gmail.com>
#  Institutions : Comenius University in Bratislava,
#                 University of California, Irvine
#

  package Filter;
 
  use strict;
  use warnings;
 
  our $VERSION = '1.2';
 
  use base 'Exporter';
  
  our @EXPORT = qw(runRNAfold runRNAcofold);
 
our $tmpDir = "/tmp/";
our $showProgress = 0;
our $useViennaRNA2 = 0;

################################################################################
sub trim {
    my($str) = shift =~ m!^\s*(.+?)\s*$!i;
    defined $str ? return $str : return '';
} 

################################################################################
sub runRNAfold {
    my $seq = shift @_;
    
#    my ($structure, $min_en) = RNA::fold($roiSeq);
    my $noPS = "-noPS";
    $noPS = "--noPS" if $useViennaRNA2;
    my @viennaResult = split " ", `echo $seq| RNAfold $noPS`;
    if ( $? != 0 ){
        die "Calling RNAfold failed: $!\n".
            "! CHECK IF ViennaRNA IS PROPERLY INSTALLED !\n";
    }
    
    my $structure = $viennaResult[1];
    my $min_en = $viennaResult[2];
    $min_en .= $viennaResult[3] if defined $viennaResult[3];
    $min_en =~ s/[\)\(]//gi;
#    print $structure, '|', $min_en, "\n";

    return ($structure, $min_en);
}

################################################################################
sub runRNAcofold {
    my $seq = quotemeta( shift @_ );
    
#    my ($structure, $min_en) = RNA::fold($roiSeq);
    my $noPS = "-noPS";
    $noPS = "--noPS" if $useViennaRNA2;
    my @viennaResult = split " ", `echo $seq| RNAcofold $noPS`;
    if ( $? != 0 ){
        die "Calling RNAcofold failed: $!\n  input passed: $seq\n".
            "! CHECK IF ViennaRNA IS PROPERLY INSTALLED !\n";
    }
    
    my $structure = $viennaResult[1];
    $structure =~ s/[\&]//g;
    my $min_en = $viennaResult[2];
    $min_en .= $viennaResult[3] if defined $viennaResult[3];
    $min_en =~ s/[\)\(]//g;
#    print $structure, '|', $min_en, "\n";

    return ($structure, $min_en);
}

################################################################################
sub parseMotifMap {
  my $descFilename = shift @_;
  my @motifMap = ();
  
  my $descFile;
  open $descFile, "<$descFilename" or die "Cannot open descriptor $descFilename: $!";
  
  while (scalar @motifMap == 0  and  my $line = <$descFile> ){
      $line = lc trim($line);
#      print $line."|\n";
      
      @motifMap = (split /\s+/, $line) unless $line =~ /^\#/;
  }
#  print join ("|", @motifMap), "|\n";
  die "No motif map found in descriptor file $descFilename!" if (scalar @motifMap == 0);

  #check the motif map
  my $okMap = 1;
  for (my $i = 0; $i < scalar @motifMap; ++$i){
#      print $motifMap[$i], "\n";
      
      #check the format
      $motifMap[$i] =~ /([shr])(\d+)(\'?)/;
      
      if (not (defined $1 and defined $2)){
          $okMap = 0;
          last;
      }
      
      #if it is a helix, check for its second strand
      if ($1 eq "h"){
          my $has2strand = 0;
          
          if ($3 eq "'"){ 
              foreach my $element(@motifMap[0 .. $i-1]){
                  $has2strand = 1 if $element eq $1.$2;
              }
          } else {
              foreach my $element(@motifMap[$i+1 .. @motifMap-1]){
                  $has2strand = 1 if $element eq $1.$2."'";
              }
          }
          
          $okMap = 0 if $has2strand == 0;
#          print "Has second strand: ".$has2strand."\n";
      }
      
#      print "$okMap|$1|$2|$3|\n";
  }
  close $descFile;
  
  die "Incorrect descriptor map \"@motifMap\" in $descFilename!" if not $okMap;
#  print join ("|", @motifMap), "|\n\n";
  return \@motifMap;
}

################################################################################
sub parseFilterConfig {
  my $fltFilename = shift @_;
  my @motifMap = @{ shift @_ };
  
  my @filter = ();
  my $roi = undef;
  
  my $fltFile;
  open $fltFile, "<$fltFilename" or die "Cannot open filter config $fltFilename: $!";
  
  while (my $line = <$fltFile> ){
      $line = trim($line);
      next if $line =~ /^\#/;
#     print $line."|\n";
      
      my @words = (split /\s+/, $line);
      next if @words == 0;
      
      #parse specifications of a submotif
      if ($words[0] =~ /^sm/i ){
#         print "\n-------------------------------------\n";
#         print join ("|", @words), "|\n";
          
          die "Invalid expression \"$words[0]\" in the filter config file" 
            unless $words[0] =~ /(^sm[^\s:%]+)(:\d+)?(%\d+)?(pn)?$/i;
          
          die "Not enough parameters for sumbotif ".$words[0]
            unless scalar @words >= 2;
          
          my $w = "";
          if (defined $2){
              $w = substr $2, 1, length($2);
          }
          
          my $tail = "";
          if (defined $3){
              $tail = substr $3, 1, length($3);
          }
          
          my $perNucleotide = 0;
          if (defined $4){
              $perNucleotide = 1;
              die "In \"$words[0]\": only one of '%' and 'pn' can be set"
                if (defined $3);
          }
          
#         print "---> ".$1.$w."\n";
          my %submotif = (
                "name" => $1,
                "w" => $w,
                "score" => "",
                "per_nucleotide" => $perNucleotide,
                "tail" => $tail,
                "unw_score" => "",
                "elements" => [],
            );
          
          foreach my $element (@words[1..@words-1]){
              die "Invalid expression \"$element\" in the filter config file" 
                  unless $element =~ /^[shr]\d+$/i;
#                 unless $element =~ /^([shr]\d+):?(\d+)?$/i;
              
#             my $w = (defined $2) ? $2 : "";
#             print "=> ".$1.'|'.$w."\n";

              die "Element \"$element\" was not defined in the decriptor"
                unless (grep /^($element)$/i, @motifMap);
          
              my $w = "";
              my %tmp = (
                    "name" => $element,
                    "w" => $w,
                    "score" => "",
                    "unw_score" => "",
                );
              push @{ $submotif{'elements'} }, \%tmp;
          }   
          
          push @filter, \%submotif;
      
      #parse specification of the region of interest    
      } elsif ($words[0] =~ /^roi$/i ){
          die "Region of interest must be exactly one element."
            unless scalar @words == 2;
          
          die "Invalid expression \"".$words[1]."\" in the filter config file" 
            unless $words[1] =~ /^[s]\d+$/i;
          
          $roi = $words[1];
#         print $roi."\n";
          
          die "Element \"$roi\" was not defined in the decriptor"
            unless (grep /^($roi)$/i, @motifMap);

      #else error
      } else {
          die "Invalid expression \"".$words[0]."\" in the filter config file";
      }
      
  }
  close $fltFile;
  
#  $dump->dumpValue(\@filter);
  return (\@filter, $roi);
}

################################################################################
sub checkFilterConfig {
  my @filter = @{ shift @_ };
  
  ## check and recalculate weights of submotifs
  #check
  my $def = 1;
  my $wsum = 0;
  foreach my $submotif (@filter){
      $def = 0 if $submotif->{'w'} eq "";
      $wsum += $submotif->{'w'} if $submotif->{'w'} ne "";
  }
  if ($def == 0 ){
      foreach my $submotif (@filter){
        $def = 1 if $submotif->{'w'} ne "";
      }
      
      die "Submotif weights must be all either defined or all undefined" 
        unless $def == 0;
  }
  #recalculate
  foreach my $submotif (@filter){
      if ($def) {
          $submotif->{'w'} *= 1 / $wsum;
          $submotif->{'w'} *= 10;
      } else {
          $submotif->{'w'} =  10 / scalar @filter;
      }
#      print $submotif->{'name'}," weight => ", $submotif->{'w'}, "\n";
  }

  ## check elements of individual submotifs
#  foreach my $submotif (@filter){
#      my $def = 1;
#      my $wsum = 0;
#      foreach my $element (@{$submotif->{'elements'}}){
#          $def = 0 if $element->{'w'} eq "";
#          $wsum += $element->{'w'} if $element->{'w'} ne "";
#          die "Element \"$element->{'name'}\" is not defined in the provided descriptor" 
#            unless grep $_ eq $element->{'name'}, @motifMap
#      }
#      if ($def == 0 ){
#          foreach my $element (@{$submotif->{'elements'}}){
#            $def = 1 if $element->{'w'} ne "";
#          }
#          
#          die "Element weights in \"$submotif->{'name'}\" must be all either defined or all undefined" 
#            unless $def == 0;
#      }
#      
#      foreach my $element (@{$submotif->{'elements'}}){
#          if ($def) {
#              $element->{'w'} *= 1 / $wsum;
#              $element->{'w'} *= 10;
#          } else {
#              $element->{'w'} =  10 / scalar @{$submotif->{'elements'}};
#          }
#  #        print $element->{'name'}," weight => ", $element->{'w'}, "\n";
#      }
#  }
}

################################################################################
sub parseRNAroboOutput {
  my $dataFilename = shift @_;
  my $numOfElements = shift @_;
  
  my @data = ();
  #convert RNArobo output to FASTA format with space separated elements 
  #call an external script to do the job, and split the output to individual lines
  my $fasta = `perl parse_rnarobo.pl -s -c -f -i $dataFilename`;
  @data = split /\n/, $fasta;
  #get rid of trailing whitespaces
  for (my $i = 0; $i < @data; ++$i) {
      $data[$i] = trim($data[$i]);
    
      if ( $i % 2 == 1){
      my @seq = split " ", $data[$i];
          die "sequence \"".$data[$i-1]."\" has wrong number of elements"
            if scalar @seq != $numOfElements;
      }
  }
  
#  $dump->dumpValue(\@data);
  return \@data;
}

################################################################################
sub parseFastaFile {
  my $fastaFileName = shift @_;
  my $numOfElements = shift @_;
  
  my @data = ();
  open (IN, "<", $fastaFileName) or die "Problem opening dump file $fastaFileName: $!";
  
  while (my $line = <IN>){
      $line = trim($line);
      #if it is a header of a sequence, then read & store it
      if ($line =~ /^>/){
          #push the header
          push @data, $line;
          
          #read, check, push the sequence
          $line = <IN>;
          $line = trim($line);
          my @seq = split " ", $line;
          die "sequence \"".$data[@data-1]."\" has wrong number of elements"
            if scalar @seq != $numOfElements;
          push @data, $line;
      }
  }
  close (IN);
#  $dump->dumpValue(\@data);
  return \@data;
}

################################################################################
sub computeScores {
  my @motifMap = @{ shift @_ };
  my @data = @{ shift @_ };
  my @filter = @{ shift @_ };

  my @seqScores = ();
  
  my $oldProgress = 0;
  
  ## evaluate every sequence
  for (my $id = 1; $id < @data; $id+=2) {
#     print $data[$id-1], "\n";    
      
      if ($showProgress) {
          my $prog = int((($id*1.0) / (@data-1))*10000)/100.0 ;
          print STDERR "running score evaluation: $prog%   \r" if $oldProgress != $prog;
          $oldProgress = $prog;
      }

      
      my @seq = split " ", $data[$id];
      my $seqScore = 0;
      my %subScores = ();
      
      #process all submotifs
      foreach my $submotif (@filter){
#         print "\n----------------\n";
          $subScores{$submotif->{'name'}} = {};
          my $submotifScore = 0;
          my @elements;
          my @map;
          my @subSeq;
          
          #retrieve list of submotif's elements
          @elements = @{$submotif->{'elements'}};
          
          #retrieve sequence and element map of the submotif
          for (my $i = 0; $i < @motifMap; ++$i){
              foreach my $element (@elements){
                  if ($motifMap[$i] =~ /^$element->{name}$/ or 
                      $motifMap[$i] =~ /^$element->{name}\'$/) {
                      
                      push (@subSeq, $seq[$i]);
                      push (@map, $motifMap[$i]);
                  }                       
              }
          }
          
          ###evaluate score of the submotif
          #check whether the submotif consists of nested helices
          my $nested = 0;
          if( @map % 2 == 0 ){
              $nested = 1;
              for (my $i = 0; $i < @map/2; ++$i){
                  $nested = 0 if ($map[$i]."'" ne  $map[@map-$i-1]);
              }
          }
          
          #check whether the submotif contains a pseudoknot
          my $pseudoknotted = 0;
          my @stack;
          foreach my $el (@map){
              if ($el =~ /^[hr]\d+$/i ){
                  push (@stack, $el."'");
              } elsif ($el =~ /^[hr]\d+\'$/i ){
                  $pseudoknotted = 1  if $el ne pop @stack;
              }
          }
          
          #calculate structure and MFE
          my ($structure, $min_en);
          if ($nested){ #if nested, calculate the fold by RNAcofold
              my $tmpseq = (join "", @subSeq[0..@subSeq/2-1] )."&"
                          .(join "", @subSeq[@subSeq/2..@subSeq-1] );
#              print $tmpseq."\n". (join "", @subSeq)."\n";
              ($structure, $min_en) = runRNAcofold($tmpseq);
#              print $structure."\n";
          } elsif($pseudoknotted) { #if contains a pseudoknot use a pseudoknot folder
              my $tmpseq = join "", @subSeq;
              $tmpseq =~ tr/T/U/;
#              print $tmpseq, "\n";
              my $tmpFileName = $tmpDir."_tmp_ff.fa";
              
              chdir "dotknot/";
              open (TMP_FILE, ">", $tmpFileName) || die "Problem opening temp file $tmpFileName for DotKnot run: $!";
              print TMP_FILE ">tmp\n$tmpseq\n"; 
              close(TMP_FILE);
              
              #check that the modified ViennaRNA-1.8.5 is there for DotKnot
              foreach my $RNAprogram ('RNAfold', 'RNAeval'){
                  `echo '' | ViennaRNA-1.8.5/Progs/$RNAprogram`;
                  if ( $? != 0 ){
                    die "Calling $RNAprogram for DotKnot failed: $!\n!".
                    "CHECK IF ViennaRNA-1.8.5 in DotKnot installation IS PROPERLY INSTALLED !\n";
                  }
              }
              
              #call DotKnot
              my $res = `python dotknot.py $tmpFileName -klg`;
              if ( $? != 0 ){
                die "Calling DotKnot failed: $!\n! CHECK IF DotKnot IS PROPERLY INSTALLED !\n";
              }
              unlink($tmpFileName);
              chdir "../";

              #parse out the pseudoknot with the highest MFE          
              if ($res =~ /best.+energy:\s?(\d+)\s(\d+)\s(-?\d+.?\d*)\s?[ACGU]+\s?([\.\(\)\[\]]+)/i){
#                  print $1, ' ', $2, ' ', $3, ' ', $4,"\n";
                  $min_en = $3;
                  $structure = $4;
                  $structure = '.' x ($1 - 1) . $structure;
                  $structure .= '.' x (length($tmpseq)-$2);
#                  print "result: $structure, $min_en \n";
                  
              } elsif ($res =~ /best.+energy:\s?no/i){ #if no pseudoknot was found
#                  print "no pseudoknot detected\n";
                  ($structure, $min_en) = runRNAfold(join "", @subSeq);
              }
          } else { #else put it to RNAfold 
              ($structure, $min_en) = runRNAfold(join "", @subSeq);
          }
          
          #convert structure to array of particular parts
          my @structure = ();
          my $tmplength = 0;
          foreach my $subsq (@subSeq){
              push @structure, substr($structure, $tmplength, length($subsq));
              $tmplength += length($subsq);
          }
          
          #unweighted score
          $submotif->{'unw_score'} = -1 * $min_en;
          #if scoring by MFE per nucelotide
          $submotif->{'unw_score'} /= 0.1 * length($structure) 
            if ($submotif->{'per_nucleotide'});
          
          #if scoring is based on tail base pairs
          if ($submotif->{'tail'} ne ""){
              my @tmpStruct = split ("", $structure);
              my @pairPos = ();
              my @stack1 = ();
              my @stack2 = ();
              for(my $j=0; $j < scalar @tmpStruct; ++$j){
                  if ($tmpStruct[$j] eq "("){
                      push @stack1, $j; 
                  } elsif ($tmpStruct[$j] eq "["){
                      push @stack2, $j; 
                      
                  } elsif ($tmpStruct[$j] eq ")"){
                      my $top = pop @stack1;
                      die "Unexpected error" unless defined $top; 
                      $pairPos[$j] = $top;
                      $pairPos[$top] = $j;
                  } elsif ($tmpStruct[$j] eq "]"){
                      my $top = pop @stack2;
                      die "Unexpected error" unless defined $top; 
                      $pairPos[$j] = $top;
                      $pairPos[$top] = $j;
                  }
              }
#             print $structure, "\n";
#             $dump->dumpValue(\@pairPos);
              
              my $count = 0;
              my $bound = ($submotif->{'tail'} > length($structure)/2) 
                          ? length($structure)/2
                          : $submotif->{'tail'};
              for (my $j=0; $j < $bound; ++$j){
                  if (defined $pairPos[$j] and 
                      $pairPos[$j] >= length($structure) - $bound)
                  {
                      ++$count;
                  }
              }
              
              $submotif->{'unw_score'} = ($count / $bound) *100;
          }
          
          #weighted score
          $submotif->{'score'} = $submotif->{'w'} * $submotif->{'unw_score'};
          
          #store intermediate scores for logging
          $subScores{$submotif->{'name'}}->{'unw_score'} = $submotif->{'unw_score'};
          $subScores{$submotif->{'name'}}->{'score'} = $submotif->{'score'};
          $subScores{$submotif->{'name'}}->{'sequence'} = join " ", @subSeq;
          $subScores{$submotif->{'name'}}->{'map'} = join " ", @map;
          $subScores{$submotif->{'name'}}->{'structure'} = join " ", @structure;
          $subScores{$submotif->{'name'}}->{'min_en'} = $min_en;
          $subScores{$submotif->{'name'}}->{'type'} = "standard";
          $subScores{$submotif->{'name'}}->{'type'} = "nested" if $nested;
          $subScores{$submotif->{'name'}}->{'type'} = "pseudoknotted" if $pseudoknotted;
          
#         print "unweighted: ", $submotif->{'name'}, " score= ", $submotif->{'unw_score'}, "\n";
#         print "weighted: ", $submotif->{'name'}, " score= ", $submotif->{'score'}, "\n";
          
          $seqScore += $submotif->{'score'};
      }
      
      #add the score to @seqScores
      push @seqScores, [$id, $seqScore, \%subScores];
      
#     print $data[$id-1], " score is ", $seqScore, "\n";
  }
  print STDERR "running score evaluation: 100.00%\n" if $showProgress;
  
  return \@seqScores;
}

################################################################################
sub extractROIs {
  my @motifMap = @{ shift @_ };
  my @data = @{ shift @_ };
  my $roi = shift @_;

  my $oldProgress = 0;
      
  my @rois = ();
  for (my $id = 1; $id < scalar @data; $id+=2) {
      
      if ($showProgress) {
          my $prog = int((($id*1.0) / (@data-1))*10000)/100.0 ;
          print STDERR "running ROIs folding: $prog%   \r" if $oldProgress != $prog;
          $oldProgress = $prog;
      }
  
  
      #get the sequence of ROI
      my @seq = split /\s+/i, $data[$id];
      die("Unexpected error") if scalar @seq != scalar @motifMap;
      
      my $roiPos = -1;
      for (my $i=0; $i < scalar @motifMap; $i++){
          if ($motifMap[$i] eq $roi) {
              $roiPos = $i;
          }
      }
      die "Unexpected error, ROI not found in MotifMap" if $roiPos == -1;
      my $roiSeq = $seq[$roiPos];

      #run RNAfold on ROI
      my ($structure, $min_en) = runRNAfold($roiSeq);
      
      #create fake amino acid code
      my @roiAAseqAr = split "", $roiSeq;
      my @structureAr = split "", $structure;
      for(my $i=0; $i < length($structure); ++$i){
          if ($structureAr[$i] eq '('){
              $roiAAseqAr[$i] = "K";
          } elsif ($structureAr[$i] eq ')'){
              $roiAAseqAr[$i] = "R";
          }
      }
      
      my $roiAAseq = (join "", @roiAAseqAr);
      $roiAAseq =~ tr/U/T/;
      $roiAAseq =~ s/[K]+/K/gi;
      $roiAAseq =~ s/[R]+/R/gi;
#        print $roiAAseq, "\n";

      #add to list of extracted ROIs
      $rois[$id/2] = [$id, $roiAAseq];
  }
  print STDERR "running ROIs folding: 100.00%\n" if $showProgress;
  
  return \@rois;
}


################################################################################
sub clusterROIs{
  my @rois = @{ shift @_ };
  my @data = @{ shift @_ };
  my $resFAFileName = shift @_;
  
  my $aaFileName = $tmpDir."_tmp_aa.fa";
  my $resUCFileName = $tmpDir."_tmp_clust_res.uc";
  
  my $os = lc trim `uname -s`;
  my $platform = "32";
  $platform = "64" if `uname -p` =~ /64/g;
  my $uclustBinary = "uclust/uclustq1.2.22_i86".$os.$platform;
    
  #sort sequences according to their length
  @rois = sort {length($b->[1]) <=> length($a->[1])} @rois;

#    $dump->dumpValue(\@data);
#    $dump->dumpValue(\@rois);
  
  open (AA_FILE, ">", $aaFileName) or die "Problem opening temp file $aaFileName for UClust run: $!";
  foreach (@rois){
    #sequence name
      my $id = $_->[0] -1;
      my $header = $data[ $id ];
      $header =~ s/^>/>$id/;
      print AA_FILE $header."\n";
    #fake AA code
      print AA_FILE $_->[1]."\n";
  }
  close(AA_FILE);
  
  system("./$uclustBinary --quiet --input $aaFileName --uc $resUCFileName --id 0.60 --maxaccepts 0 --maxrejects 0");
  if ( $? != 0 ){
      die "Calling UClust failed: $!\n".
          "! CHECK IF UClust IS PROPERLY INSTALLED !\n";
  }
  
  system("./$uclustBinary --quiet --uc2fasta $resUCFileName --input $aaFileName --output $resFAFileName");
  if ( $? != 0 ){
      die "Calling UClust failed: $!\n".
          "! CHECK IF UClust IS PROPERLY INSTALLED !\n";
  }
  
  unlink($aaFileName);
  unlink($resUCFileName);

}

################################################################################
sub filterClusters {
  my $resFAFileName = shift @_;
  my $sizeThreshold = shift @_;
  my $spiecesThreshold = shift @_;
  
  my @filteredHeaders = ();
  
  my @tmpHeaders = ();
  my $clusterSize = -1;
  my $clusterID = -1;
  
  my $oldProgress = 0;
  my $counter = 0;   
  
  open (FA_FILE, "<", $resFAFileName) || die "Problem opening temp file $resFAFileName: $!";
  my @clusters = <FA_FILE>;
  foreach my $loopVar (0 .. $#clusters){
      my $line = $clusters[$loopVar];
      if( $line =~ /^\s*>(\d+)\|[^\|]+\|\s*(.+)$/g){
#          print 'M: ', $1, "_||_", $2, "\n";
          my $cluser = $1;
          my $header = $2;

          if ($showProgress) {
              ++$counter;
              my $prog = int((($counter*1.0) / (@clusters/2))*10000)/100.0;
              print STDERR "running cluster-filtering: $prog%   \r" if $oldProgress != $prog;
              $oldProgress = $prog;

          }

          #if it is the same cluster as previous
          if ($clusterID == $cluser) {
              push @tmpHeaders, $header;
              ++$clusterSize;
          } 
          
          #if it is a new cluster or we are at the end of array
          if ($clusterID != $cluser or $loopVar == $#clusters-1) {
              if ($clusterSize >= $sizeThreshold){
                  ##check for enought variety of species
                  my %hash = ();
                  foreach (@tmpHeaders){
                      my @pieces = split " ", $_;
                      my $hit = quotemeta($pieces[4]);
                      
                      foreach my $i (4..8){
                          last if !defined $pieces[$i];
                          if ($pieces[$i] =~ /^\d|^gi|^ref|^emb|^gb|^NC|^AA|^[A-Z]{1,4}\d/i){
                              $hit = quotemeta($pieces[$i]);
                              $hit = quotemeta($pieces[$i+1]) if defined $pieces[$i+1];
                              last;
                          }
                      }
                      $hit = lc $hit;
                      $hit =~ s/\W//g;
                      
                      $hash{$hit} = 1;	                    
# 	                    print $hit."\n";
                  }
#                    print scalar keys %hash, "\n";
                  if (scalar keys %hash >= $spiecesThreshold){
                      push @filteredHeaders, [@tmpHeaders];
                  }
              }
              
              #sepcial case = if it is the last cluser and is of size 1
              if ($clusterID != $cluser and $loopVar == $#clusters-1){ 
                  push @filteredHeaders, [$header] 
                    if ($spiecesThreshold <= 1) and ($sizeThreshold <= 1);
              }
              
              #start a new cluster
              $clusterID = $cluser;
              $clusterSize = 1;
              @tmpHeaders = ($header);
          }
      }
  }
  close(FA_FILE);
  print STDERR "running cluster-filtering: 100.00%\n" if $showProgress;
  
  #sort resulting clusters according to their size (from largest to smallest)
  @filteredHeaders = sort {scalar @{$b} <=> scalar @{$a}} @filteredHeaders;
  
  return \@filteredHeaders;
}

################################################################################


1;
