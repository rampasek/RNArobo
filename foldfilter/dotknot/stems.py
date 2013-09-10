"""
    DotKnot: software for RNA pseudoknot prediction
    Copyright (C) 2010-2011 Jana Sperschneider	

    This program is free software; you can redistribute it and/or modify  
    it under the terms of the GNU General Public License as published by  
    the Free Software Foundation; either version 2 of the License, or     
    (at your option) any later version.    
"""

import os
import sys
import re
from subprocess import Popen, PIPE

  
def RNAfold_dotplot(seq, RNAFOLD_PATH, FILE_ID):
    """ Function: RNAfold_dotplot()

        Purpose:  Given an input RNA sequence, obtain stack probabilities from RNAfold -p.
                  Probability threshold E-11 set by hand in RNAfold.c.
                  For each base pair with an assigned stack probability, store it in a
                  dictionary: Base pair -> Stack probability
                  
        Input:    RNA sequence
        
        Return:   A dictionary of stack probabilities. 
    """  
    try:
        infile = open(FILE_ID + "_input.fasta", "w")
        print >> infile, '>', FILE_ID
        print >> infile, seq
        infile.close()
    except IOError:
        print "Error: could not write to disk!"
        sys.exit(1)
	
    command = RNAFOLD_PATH + "RNAfold -p2 -d2 -noLP < " + FILE_ID + "_input.fasta"
    result = Popen(command, shell=True, stdout = PIPE)
    mfe, err = result.communicate()
    
    if not mfe:
        print "Check if RNAfold works."
        sys.exit(0)   
        
    try:
        f = open(FILE_ID + "_dp2.ps",'r')  
    except IOError:
        print "Error: the dot plot file 'dot2.ps' created by RNAfold could not be opened!"
        sys.exit(1)
        
    basepairs = {}
    
    for line in f:  
        if "lbox" in line and not "/lbox" in line:
            info = line.split()     
            if int(info[0]) < len(seq):
                base_pair = int(info[0]), int(info[1])
                basepairs[base_pair] = float(info[2])*float(info[2])        
      
    return basepairs


def find_stems(basepairs):
    """ Function: find_stems()

        Purpose:  Given the dictionary of stack probabilities, assemble regular stems.
                  For a given stem start position (i,j), find stacked base pairs
                  (i+1,j-1), (i+2,j-2), ... where absolute percentage increase/decrease
                  between stacked base pairs is < 66.7%. A stem needs to have at least
                  three base pairs. Delete duplicate stems from dictionary to keep only
                  stems with maximum length. This guarantees that each stem can be stored
                  using an unique key (i,j), where i is the start position and j is the
                  end position of the stem in the sequence. For each stem, calculate
                  confidence indicator and store it in dictionary: Stem -> (Length, Confidence)
                  
        Input:    A dictionary of stack probabilities. 

        Return:   A dictionary of stems.
    """
    stem_list = []
    
    for bp, stack_prob_ij in sorted(basepairs.items()):
        stem = []
        stem.append((bp, stack_prob_ij))   # Stacked base pair (i,j) and (i+1,j-1)
        count = 1
        
        stacked_pair = True
        
        while stacked_pair:                # Look for stacked base pairs (i+1,j-1), (i+2,j-2), ...           
            bp_next = (bp[0] + count, bp[1] - count)            
            if bp_next in basepairs:
                stack_prob_next = basepairs[bp_next]
                minimum = min(stack_prob_ij, stack_prob_next)
                maximum = max(stack_prob_ij, stack_prob_next)

                if (abs((maximum - minimum)/maximum))*100 < 66.7: # Calculate absolute percentage increase/decrease
                    stem.append((bp_next, stack_prob_next))                
                    count = count + 1
                    stack_prob_ij = stack_prob_next                    
                else:                    
                    stem.append((bp_next, 0.0))
                    stacked_pair = False
                    
            else:   # The last stacked base pair is discovered, append with stack probability = 0.0                                
                triple = bp_next, 0.0
                stem.append(triple)
                stacked_pair = False
                                
        # Stem has format [((75, 96), 0.0011170568217599998), ((76, 95), 0.00089892632041000005), ((77, 94), 0.0)]    
        if len(stem) >= 3:            
            prob = sum([item[1] for item in stem])
            confidence = prob/(len(stem)-1)  # Calculate confidence indicator                                       
            key = stem[0][0][0], stem[0][0][1]    
            values = len(stem), confidence
            stem_list.append((key, values))  
                      
    stem_dic = dict(stem_list)               # Store all stems in dictionary
    all_stems = stem_dic.items()

    # Delete duplicate stems from dictionary to keep only stems with maximum length, delete (57,65): 3 if (56,66): 4 exists
    for (start, end), (length, confidence) in all_stems:
        duplicate = start + 1, end - 1         
        if duplicate in stem_dic and stem_dic[duplicate][0] == length - 1:  
            del stem_dic[duplicate]

    return stem_dic              


def evaluation_stems(stem_dic, seq, RNAFOLD_PATH, FILE_ID):
    """ Function: evaluation_stems()

        Purpose:  Function for evaluation of local stem energies using RNAeval.
                  Create dot-bracket notation and input file to RNAeval. Call
                  RNAeval with and without hairpin loop entropies. Scan output
                  files and extract free stacking energies and free energies.
                  Update stem dictionary: Stem -> (Length, Confidence, Stacking Energy, Free Energy)
                  
        Input:    A dictionary of stems and the sequence. 

        Return:   An updated dictionary of stems.
    """    
    try:
        stem_structure = file(FILE_ID + "_stem_structure.txt",'w')
    except IOError:
        print "Error: file could not be written to disk."
        sys.exit(1)
  
    for stem, (length_stem, confidence) in sorted(stem_dic.items()):   # Iterate like this because of shortened stem keys
        start, end = stem[0], stem[1]
        
        # Create dot-bracket notation for each stem
        string_list = ["(" for i in xrange(length_stem)] 
        string_list = string_list + [":" for i in xrange(end - start + 1 - 2*length_stem)]
        string_list = string_list + [")" for i in xrange(length_stem)]      
        
        # Watch out for dangling ends
        if start == 1:
            if end != len(seq):     # Dangling end on the right only
                local_sequence = seq[start-1:end+1]                                
                string_list.append(":")                
            else:                   # No dangling ends
                local_sequence = seq[start-1:end]
        else:
            if end != len(seq):     # Dangling end on both ends
                local_sequence = seq[start-2:end+1]
                string_list.insert(0,":")
                string_list.append(":")         
            else:                   # Dangling end on the left only
                local_sequence = seq[start-2:end]
                string_list.insert(0,":")
                
        structure = "".join(string_list)
        header = "> " + str(start) + " " + str(end) + " " + str(length_stem) + " " + str(confidence) + "\n"
        stem_structure.write(header)
        stem_structure.write(local_sequence)
        stem_structure.write("\n")
        stem_structure.write(structure)
        stem_structure.write("\n")
        
    stem_structure.close()

    RNAeval_files(RNAFOLD_PATH, FILE_ID)
    stem_dic = update_stem_dictionary(stem_dic, FILE_ID)
    
    return stem_dic

  
def RNAeval_files(RNAFOLD_PATH, FILE_ID):
    """ Function: RNAeval_files()

        Purpose:  Call RNAeval with and without hairpin loop entropies and
                  write to files.               
    """  
    try: 
        # Evaluate with stacking only
        command = RNAFOLD_PATH + "RNAeval -d2 -P stacking_only.par < " + FILE_ID + "_stem_structure.txt > " + FILE_ID + "_stacking_energy.txt"

        result = Popen(command, shell=True, stdout = PIPE)
        eval_stacking, err = result.communicate()         
        
        # Evaluate with loop entropies
        command = RNAFOLD_PATH + "RNAeval -d2 < " + FILE_ID + "_stem_structure.txt > " + FILE_ID + "_loops_energy.txt"

        result = Popen(command, shell=True, stdout = PIPE)
        eval_loops, err = result.communicate()         

    except:
        print "Error: RNAeval could not be called!"
        
    return


def update_stem_dictionary(stem_dic, FILE_ID):
    """ Function: update_stem_dictionary()

        Purpose:  Scan files produced by RNAeval. Update stem dictionary:
                  Stem -> (Length, Confidence, Stacking Energy, Free Energy)
                  
        Input:    A dictionary of stems. 

        Return:   An updated dictionary of stems.
    """    
    output = file(FILE_ID + "_stacking_energy.txt",'r')
    stem_list = [i for i in output]
    output.close()    

    output = file(FILE_ID + "_loops_energy.txt",'r')
    stem_loops_list = [i for i in output]
    output.close()

    # The energy evaluation result is stored in two lists now    
    for i in xrange(0, len(stem_list), 3):
        
        info_stack = stem_list[i+2].split()              
        float_in_lists = [re.findall(r"-?\d+\.?\d*", s) for s in info_stack]
        for list in float_in_lists:
            if list:
                energy_stack = list[0]
            
        info_loops = stem_loops_list[i+2].split()
        float_in_lists = [re.findall(r"-?\d+\.?\d*", s) for s in info_loops]
        for list in float_in_lists:
            if list:
                energy_loops = list[0]    
                
        info = stem_list[i].split()            
        key = int(info[1]), int(info[2])
        key_short = int(info[1]), int(info[2]), int(info[3])
            
        if key in stem_dic:         # Case for evaluation of regular stems
            values = stem_dic[key][0], stem_dic[key][1], float(energy_stack), float(energy_loops)
            stem_dic[key] = values

        if key_short in stem_dic:   # Case for evaluation of shortened stems
            values = stem_dic[key_short][0], stem_dic[key_short][1], float(energy_stack), float(energy_loops)
            stem_dic[key_short] = values
            
    return stem_dic


def update_confidence(stem_dic, stems_shortened_dic):
    """ Function: update_confidence()

        Purpose:  Replace entry 0.0 with confidence indicator in shortened stem dictionary.              
                  
        Input:    A dictionary of shortened stems. 

        Return:   An updated dictionary of shortened stems.
    """    
    stem_list = stems_shortened_dic.items()

    for stem, values in stem_list:
        if (stem[0], stem[1]) in stem_dic:
            confidence = stem_dic[(stem[0], stem[1])][1]
            stems_shortened_dic[stem] = values[0], confidence, values[2], values[3]

    return stems_shortened_dic


def filter_stems(stem_dic, CUTOFF_STACK, CUTOFF_LOOP):    
    """ Function: filter_stems()

        Purpose:  Filtering of stems with high free (stacking) energy.              
                  
        Input:    A dictionary of stems. 

        Return:   An filtered dictionary of stems.
    """
    dictionary = stem_dic.copy()
    stem_list = dictionary.items()

    for stem, values in stem_list:
        if values[2] >= CUTOFF_STACK or values[3] >= CUTOFF_LOOP:
            del dictionary[stem]          

    return dictionary


def filter_stems_prob(stem_dic, CUTOFF_PROB):
    """ Function: filter_stems_prob()

        Purpose:  Create new dictionary which stores stems with high probability.              
                  
        Input:    A dictionary of stems. 

        Return:   An new dictionary of stems.
    """    
    stems_ib = stem_dic.copy()
    stem_list = stems_ib.items()

    for stem, values in stem_list:
        if values[1] < CUTOFF_PROB:   
            del stems_ib[stem]

    return stems_ib

