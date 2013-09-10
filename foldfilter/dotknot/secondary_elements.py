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
import functions
import re
from subprocess import Popen, PIPE 


def internal_mwis(stems_ib):
    """ Function: internal_mwis()

        Purpose:  Scan the list of stems. For each outer stem, a maximum weight
                  independent set (MWIS) calculation with confidence weights on
                  the set of all possible inner stems is performed. This
                  information is stored in a dictionary: Outer stem -> [Inner stems]
                  
        Input:    A dictionary of stems. 
        
        Return:   A dictionary of secondary structure information.
    """    
    structures_dic, candidate_list = {}, []    

    for stem, values in stems_ib.items():        
        entry = (stem[0], stem[1], values[0], values[1], values[1])
        candidate_list.append(entry)
        
    candidate_list.sort()    
    sorted_endpoint_list = functions.create_sorted_endpointlist(candidate_list)

    for endpoint in sorted_endpoint_list:   # Scan sorted endpoints list                
        if endpoint[1] == 'r':              # If a right endpoint is scanned
            sorted_endpoint_list_recursive, nested = [], []
            index = endpoint[3]
            interval = candidate_list[index - 1]
            
            nested = find_nested(interval, candidate_list)

            if nested:                      # MWIS on the set of inner stems
                sorted_endpoint_list_recursive = functions.create_sorted_endpointlist(nested)                
                result = functions.MWIS(nested, sorted_endpoint_list_recursive)
                confidence = interval[4] + sum([element[4] for element in result])  # Confidence sum

                # Store updated confidence for outer stem
                candidate_list[index-1] = (interval[0], interval[1], interval[2], interval[3], confidence)
                stem = interval[0], interval[1], interval[2]
                
                # Store inner structure elements in dictionary
                structures_dic[stem] = result

            else:                
                stem = interval[0], interval[1], interval[2]                
                structures_dic[stem] = []
    
    return structures_dic


def find_nested(interval, candidate_list):
    """ Function: find_nested()

        Purpose:  Find all inner stems, given a right endpoint of an outer stem
                  Flexible inner stem finding, allow certain base pair overlap.
                  
        Input:    An outer interval and list of candidate inner intervals. 
        
        Return:   Inner intervals.
    """
    result = []
    
    interval_left = interval[0] + interval[2]
    interval_right = interval[1] - interval[2]
    stemlength_outer = interval[2]
    
    for compare_interval in candidate_list:
        stemlength_inner = compare_interval[2]
        if compare_interval[0] >= interval_left:                              
            if compare_interval[1] <= interval_right:
                result.append(compare_interval)
                
        if stemlength_inner > 3 and stemlength_outer > 3:   # Overlap of 1 bp
            if compare_interval[0] >= interval_left - 1:                                                          
                if compare_interval[1] <= interval_right + 1:
                    result.append(compare_interval)
                    
        if stemlength_inner > 4 and stemlength_outer > 4:   # Overlap of 2 bp
            if compare_interval[0] >= interval_left - 2:                                                           
                if compare_interval[1] <= interval_right + 2:
                    result.append(compare_interval)
                    
    result = list(set(result))          # Remove duplicates from list
    result.sort()
    
    return result


def evaluation_secondary_structures(structures_dic, seq, RNAFOLD_PATH, FILE_ID):    
    """ Function: evaluation_secondary_structures()

        Purpose:  Find all inner stems, given a right endpoint of an outer stem
                  Flexible inner stem finding, allow certain base pair overlap.
                  
        Input:    An outer interval and list of candidate inner intervals. 
        
        Return:   Inner intervals.
    """
    try:
        bulges_internal = file(FILE_ID + "_bulge_internal_structures.txt",'w')
        multiloops = file(FILE_ID + "_multiloop_structures.txt",'w')
    except IOError:
        print "Error: could not write files to disk!"
        sys.exit(1)

    for stem in sorted(structures_dic):
        
        ml = False
        stem_list_recursive = []

        start, end, length_stem = stem[0], stem[1], stem[2]
        string = "> " + str(start) + " " + str(end)+ " " + str(length_stem) + "\n"
        local_sequence = seq[start-1:end]        
        
        dot_bracket = ['(' for i in xrange(length_stem)]
        dot_bracket = dot_bracket + ['.' for i in xrange(len(local_sequence) - 2*length_stem)]
        dot_bracket = dot_bracket + [')' for i in xrange(length_stem)]
                  
        inner_stems = structures_dic[stem]        

        if inner_stems:            
            if len(inner_stems) > 1:
                ml = True            

            for stem_internal in inner_stems:    # Look up all internal structure elements
                stem_list_recursive.append(stem_internal)
                nested_recursive = structures_dic[stem_internal[0], stem_internal[1], stem_internal[2]]

                if len(nested_recursive) > 1:
                    ml = True
                    
                stem_list_recursive, ml = find_recursive(structures_dic, stem_list_recursive, stem_internal,ml)
                
            stem_list_recursive.sort()
            
            local_sequence, structure = dot_bracket_notation(start, end, seq, dot_bracket, stem_list_recursive)
            
            if ml == True:                 # If it is a multiloop                
                multiloops.write(string)                                        
                multiloops.write(str(local_sequence))
                multiloops.write('\n')                     
                multiloops.write(str(structure))
                multiloops.write('\n')
                
            if ml == False:                 # If it is not a multiloop
                bulges_internal.write(string)                               
                bulges_internal.write(str(local_sequence))
                bulges_internal.write('\n')                     
                bulges_internal.write(str(structure))
                bulges_internal.write('\n')
                
    bulges_internal.close()
    multiloops.close()
    
    RNAeval_files(RNAFOLD_PATH, FILE_ID)
    bulges_internal, multiloops = create_stem_dictionaries(FILE_ID)
                                
    return bulges_internal, multiloops


def find_recursive(structures_dic, stem_list_recursive, stem_internal, ml):
    """ Function: find_recursive()

        Purpose:  Find all inner stems in a recursive fashion. Remember if outer
                  stem will form a multiloop structure by counting the number of
                  inner stems. 
                  
        Input:    Information on inner stems and multiloop marker. 
        
        Return:   Inner intervals and multiloop marker.
    """    
    stem_internal_list = structures_dic[stem_internal[0], stem_internal[1], stem_internal[2]]

    if stem_internal_list:
        if len(stem_internal_list) > 1:
            ml = True
            
        for stem in stem_internal_list:            
            stem_list_recursive.append(stem)
            stem_list_recursive, ml = find_recursive(structures_dic, stem_list_recursive, stem,ml)
            
    return stem_list_recursive, ml


def dot_bracket_notation(start, end, seq, dot_bracket, stem_list_recursive):
    """ Function: dot_bracket_notation()

        Purpose:  Insert recursive stems into dot-bracket notation.
                  
        Input:    List of recursive stems and outer stem dot bracket structure. 
        
        Return:   Local sequence and dot bracket structure for secondary structure element.
    """    
    for stem in stem_list_recursive:
        stem_start = stem[0]
        stem_end = stem[1]
        length_stem = stem[2]
        
        for i in xrange(length_stem):
            if dot_bracket[stem_start - start + i] == '.':
                if dot_bracket[stem_end - start - i] == '.':
                    dot_bracket[stem_start - start + i] = '('
                    dot_bracket[stem_end - start - i] = ')'
                    
    if start == 1:
        if end != len(seq):     # Dangling end on the right only
            local_sequence = seq[start-1:end+1]
            dot_bracket.append(":")     
        else:                   # No dangling ends
            local_sequence = seq[start-1:end]
    else:
        if end != len(seq):     # Dangling end on both ends
            local_sequence = seq[start-2:end+1]
            dot_bracket.insert(0,":")
            dot_bracket.append(":")         
        else:                   # Dangling end on the left only
            local_sequence = seq[start-2:end]
            dot_bracket.insert(0,":")                    

    structure = "".join(dot_bracket)
    
    return local_sequence, structure

   
def RNAeval_files(RNAFOLD_PATH, FILE_ID):
    """ Function: RNAeval_files()

        Purpose:  Call RNAeval with hairpin loop entropies and write to file.
                  Only stacking energies does not work with RNAeval. Handle
                  this later. 
    """     
    try: 
        # Evaluate with stacking only
        command = RNAFOLD_PATH + "RNAeval < " + FILE_ID + "_multiloop_structures.txt > " + FILE_ID + "_multiloops_energy.txt"
        
        result = Popen(command, shell=True, stdout = PIPE)
        eval_stacking, err = result.communicate()               
        
        # Evaluate with loop entropies
        command = RNAFOLD_PATH + "RNAeval < " + FILE_ID + "_bulge_internal_structures.txt > " + FILE_ID + "_bulge_internal_energy.txt"
        
        result = Popen(command, shell=True, stdout = PIPE)
        eval_loops, err = result.communicate()           

    except:
        print "Error: RNAeval could not be called!"
        sys.exit(1)
        
    return


# Dictionary which stores hairpin loop entropies for bulge/internal loop energy evaluation.
# RNAfold does not return the correct stacking free energies without hairpin loop entropies.
hairpin_dic = {
 0: 0.0 ,   1: 0.0 , 2: 0.0,  3: 4.1,  4: 4.9,  5: 4.4,  6: 4.7,  7: 5.0,  8: 5.1,  9: 5.2, 10: 5.3,
11: 5.4 , 12: 5.5 , 13: 5.6, 14: 5.7, 15: 5.8, 16: 5.8, 17: 5.9, 18: 5.9, 19: 6.0, 20: 6.1, 21: 6.1,
22: 6.2 , 23: 6.2 , 24: 6.3, 25: 6.3, 26: 6.3, 27: 6.4, 28: 6.4, 29: 6.5, 30: 6.5
}

def create_stem_dictionaries(FILE_ID):
    """ Function: create_stem_dictionaries()

        Purpose:  Scan files produced by RNAeval. Create stem dictionaries:
                  Stem -> (Length, Dot-Bracket Notation, Free Energy)              

        Return:   Two dictionary of bulge/internal loop structures and multiloops.
    """        
    output = file(FILE_ID + "_bulge_internal_energy.txt",'r')
    stem_list = [i for i in output]
    output.close()
    
    bulge_internal_dic = {}     # Bulge and internal loop dictionary
    
    for i in xrange(0, len(stem_list), 3):
        
        info_stack = stem_list[i+2].split()
        dot_bracket = info_stack[0]
        
        float_in_lists = [re.findall(r"-?\d+\.?\d*", s) for s in info_stack]
        for list in float_in_lists:
            if list:
                free_energy = list[0]        
          
        info = stem_list[i].split()
        key = int(info[1]), int(info[2])
        length = int(info[3])

        # Calculate stacking free energy by substracting inner hairpin loop
        left = dot_bracket.rfind('(') - 1     # Position of last bracket '(' 
        right = dot_bracket.find(')') - 1     # Position of first bracket ')'     
        if key[0] == 1:                       # Case of dangling end 
            left = left + 1
            right = right + 1
        hairpin_loop = len(dot_bracket[left:right])-1
        if hairpin_loop > 30:
            hairpin_loop = 30
        entropy = hairpin_dic[hairpin_loop] 
        stack_energy = float(free_energy) - entropy

        bulge_internal_dic[key] = length, dot_bracket, float(free_energy), float(stack_energy)
        
    output = file(FILE_ID + "_multiloops_energy.txt",'r')
    stem_list = [i for i in output]
    output.close()
    
    multiloops = {}         # Multiloop dictionary
    
    for i in range(0, len(stem_list), 3):
        
        info_stack = stem_list[i+2].split()
        dot_bracket = info_stack[0]        

        float_in_lists = [re.findall(r"-?\d+\.?\d*", s) for s in info_stack]
        for list in float_in_lists:
            if list:
                free_energy = list[0]
                
        info = stem_list[i].split()
        key = int(info[1]), int(info[2])
        length = int(info[3])
        multiloops[key] = length, dot_bracket, float(free_energy)

    return bulge_internal_dic, multiloops


def filter_stems(secondary_dic, cutoff_loops):
    """ Function: filter_stems()

        Purpose:  Filtering of secondary structure elements with
                  high free (stacking) energy.              
                  
        Input:    A dictionary of secondary structure elements. 

        Return:   An filtered dictionary of secondary structure elements.
    """
    stem_list = secondary_dic.items()
    
    for stem, values in stem_list:
        if values[2] >= cutoff_loops:
            del secondary_dic[stem]

    return secondary_dic    
