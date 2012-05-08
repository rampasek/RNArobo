"""
    DotKnot: software for RNA pseudoknot prediction
    Copyright (C) 2010-2011 Jana Sperschneider	

    This program is free software; you can redistribute it and/or modify  
    it under the terms of the GNU General Public License as published by  
    the Free Software Foundation; either version 2 of the License, or     
    (at your option) any later version.
"""

import functions
import pk_tools

L1_LOWER = 0
L1_UPPER = 400
L2_LOWER = 0 
L2_UPPER = 400
L3_LOWER = 0
L3_UPPER = 400

MAX_LENGTH = 400
    
    
def loops_fulfilled(l1, l2, l3):
    """ Function: loops_fulfilled()

        Purpose:  Check if loop length requirements are fulfilled
                  
        Input:    Looplengths
        
        Return:   True or False
    """    
    loops_fulfilled = False
    
    if l1 >= L1_LOWER and l1 < L1_UPPER:
        if l2 >= L2_LOWER and l2 < L2_UPPER:                               
            if l3 >= L3_LOWER and l3 < L3_UPPER:
                loops_fulfilled = True
    
    return loops_fulfilled


def add_pk(key, l1, l2, l3, all_pseudoknots):
    """ Function: add_pk()

        Purpose:  Add simple H-type pseudoknot constructed with two stems to
                  first and second pseudoknot dictionary. Check that loop 
                  length requirements are fulfilled.
                  
        Input:    A pseudoknot and a dictionary of H-type pseudoknots. 
        
        Return:   An updated dictionary of H-type pseudoknots.
    """
    if loops_fulfilled(l1, l2, l3):        
        if (l1 > 3 or l3 > 3) and (l1 + l2 + l3) < MAX_LENGTH:        
            all_pseudoknots[key] = 0.0      # Default energy 0.0
            
    return all_pseudoknots 

               
def resolve_overlap(stem_start, stem_end, stemlength, loop_overlap, loop_test):
    """ Function: resolve_overlap()

        Purpose:  Handle base pair overlap for core H-type pseudoknot.
                  
        Input:    Pseudoknot stems and loop lengths. 
        
        Return:   Updated pseudoknot stems and loop lengths. 
    """    
    # Case that base pair overlap of 1 bp occurs at loop
    if loop_overlap == -1 and stemlength >= 4 and loop_test >= 0:          
        stemlength = stemlength - 1       # Cut stem if possible
        loop_overlap = loop_overlap + 1
        loop_test = loop_test + 1
        stem_shortended = stem_start, stem_end, stemlength

    # Case that base pair overlap of 2 bp occurs at loop
    elif loop_overlap == -2 and stemlength >= 5 and loop_test >= -1:          
        stemlength = stemlength - 2       # Cut stem if possible
        loop_overlap = loop_overlap + 2
        loop_test = loop_test + 2
        stem_shortended = stem_start, stem_end, stemlength 

    else:
        stem_shortended = None
    
    return stem_shortended, loop_overlap, loop_test


def build_pseudoknots(stem_dic):
    """ Function: build_pseudoknots()

        Purpose:  Construct core H-type pseudoknots with regular stems.
                  Allow certain base pair overlaps.
                  
        Input:    A dictionary of regular stems. 
        
        Return:   Core pseudoknot dictionary and dictionary for shortened stems.
    """
    all_pseudoknots = {}    
    
    pk_dic, stems_shortened = {}, {}
    
    stem_dic_list = stem_dic.items()
    stem_dic_list.sort()
    
    for index, ((i, j), S1_values) in enumerate(stem_dic_list):
        for (k, l), S2_values in stem_dic_list[index:]:                
            
            if (l - i + 1) >= 16 and (l - i + 1) < 400:

                stemlength1 = S1_values[0]            
                stemlength2 = S2_values[0]
                l1 = k - (i + stemlength1)            
                l2 = (j - stemlength1 + 1) - (k + stemlength2)
                l3 = (l - stemlength2) - j                  
                s1_shortended, s2_shortended = None, None
                
                # Case that no base pair overlap occurs
                key = i, j, stemlength1, k, l, stemlength2
                all_pseudoknots = add_pk(key, l1, l2, l3, all_pseudoknots)
                
                # Loop L1 is shorter than 3 nt. Loop L3 has to be longer than 3 nt
                # to accommodate the kissing interaction. An overlap of 1 bp (2 bps)
                # can occur at loops L1 or L2 and stem S1 is shortened by 1 bp (2 bps).
                if l1 <= 3 and l3 > 3:
                    if (l1 == -1 or l1 == -2) and l2 >= -1:                       
                        s1_shortended, l1, l2 = resolve_overlap(i, j, stemlength1, l1, l2)                                                   
                    if (l2 == -1 or l2 == -2) and l1 >= -1:
                        s1_shortended, l2, l1 = resolve_overlap(i, j, stemlength1, l2, l1)

                    if s1_shortended:
                        key = s1_shortended[0], s1_shortended[1], s1_shortended[2], k, l, stemlength2 
                        stems_shortened[s1_shortended] = s1_shortended[2], 0.0
                        all_pseudoknots = add_pk(key, l1, l2, l3, all_pseudoknots)                                
                            
                # Loop L3 is shorter than 3 nt. Loop L1 has to be longer than 3 nt 
                # to accommodate the kissing interaction. An overlap of 1 bp (2 bps)
                # can occur at loops L2 or L3 and stem S2 is shortened by 1 bp (2 bps).
                if l1 > 3 and l3 <= 3:
                    if (l2 == -1 or l2 == -2) and l3 >= -1:
                        s2_shortended, l2, l3 = resolve_overlap(k, l, stemlength2, l2, l3)                           
                    if (l3 == -1 or l3 == -2) and l2 >= -1:                        
                        s2_shortended, l3, l2 = resolve_overlap(k, l, stemlength2, l3, l2)

                    if s2_shortended:
                        key = i, j, stemlength1, s2_shortended[0], s2_shortended[1], s2_shortended[2]
                        stems_shortened[s2_shortended] = s2_shortended[2], 0.0
                        all_pseudoknots = add_pk(key, l1, l2, l3, all_pseudoknots)
                            
                # Loops L1 and L3 are both >= 3 nt and therefore it is not clear where the
                # kissing interaction might form. An overlap of 1 bp (2 bps) can occur at loop
                # L2 and in such a case, the longer stem S1 or S2 is shortened by 1 bp (2 bps).
                if l1 >= 3 and l3 >= 3:
                    # Case that base pair overlap of 1 bp occurs at L2
                    overlap = l2
                    if overlap == -1 or overlap == -2:
                        if stemlength1 >= stemlength2 and stemlength1 >= 3 - overlap:      
                            stemlength1 = stemlength1 + overlap       # Cut longer stem S1
                            l1 = l1 - overlap
                            l2 = l2 - overlap
                            s1_shortended = i, j, stemlength1   # Add to stem dictionary
                            stems_shortened[s1_shortended] = s1_shortended[2], 0.0                            
                            key = i, j, stemlength1, k, l, stemlength2                                                                     
                            all_pseudoknots = add_pk(key, l1, l2, l3, all_pseudoknots)
                            
                        if stemlength2 > stemlength1 and stemlength2 >= 3 - overlap:      
                            stemlength2 = stemlength2 + overlap       # Cut longer stem S2
                            l2 = l2 - overlap
                            l3 = l3 - overlap                            
                            s2_shortended = k, l, stemlength2   # Add to stem dictionary
                            stems_shortened[s2_shortended] = s2_shortended[2], 0.0                        
                            key = i, j, stemlength1, k, l, stemlength2                        
                            all_pseudoknots = add_pk(key, l1, l2, l3, all_pseudoknots)
    
    return stems_shortened, all_pseudoknots


def filter_khp_dictionary(all_pseudoknots, stems_shortened_dic, matrix_stems):
    """ Function: filter_khp_dictionary()

        Purpose:  For a pseudoknot dictionary, delete those pseudoknots
                  where the shortened stem has been filtered because of
                  high free energy. Construct two new dictionaries which
                  store pseudoknots by their first and second stem. 
                  
        Input:    Pseudoknot dictionary.
        
        Return:   Filtered pseudoknot dictionary and two dictionaries which
                  store pseudoknots by first and second stem.
    """
    pseudoknots = all_pseudoknots.items()
    
    for pk, default in pseudoknots:
        stem1_key = pk[0], pk[1]
        stemlength1 = pk[2]
        stem1_short = pk[0], pk[1], pk[2]

        stem2_key = pk[3], pk[4]
        stemlength2 = pk[5]
        stem2_short = pk[3], pk[4], pk[5]
                    
        if matrix_stems[stem1_key][0] != stemlength1:
            if stem1_short not in stems_shortened_dic:  # Shortened stem was filtered because of high free energy
                if pk in all_pseudoknots:
                    del all_pseudoknots[pk]                    
                
        if matrix_stems[stem2_key][0] != stemlength2:
            if stem2_short not in stems_shortened_dic:  # Shortened stem was filtered because of high free energy                                    
                if pk in all_pseudoknots:
                    del all_pseudoknots[pk]                                    

    # Construct new dictionaries for first and second pseudoknot stems
    pseudoknot_first, pseudoknot_second = {}, {}
    
    for pk in all_pseudoknots:        
        stem1 = pk[0:2]         # Store without stem length because of overlap
        stem2 = pk[3:5]         # Store without stem length because of overlap

        if stem1 in pseudoknot_first:
            entries = pseudoknot_first[stem1]
            entries.append(pk)
            pseudoknot_first[stem1] = entries
        else:
            entries = []
            entries.append(pk)            
            pseudoknot_first[stem1] = entries

        if stem2 in pseudoknot_second:
            entries = pseudoknot_second[stem2]
            entries.append(pk)
            pseudoknot_second[stem2] = entries
        else:
            entries = []
            entries.append(pk)
            pseudoknot_second[stem2] = entries
    
    return all_pseudoknots, pseudoknot_first, pseudoknot_second 


def kissing_hairpins(seq, pseudoknot_second, pseudoknot_first, all_stems, array_traceback, INIT, UNPAIRED_NT, UNPAIRED_NT_L3):
    """ Function: kissing_hairpins()

        Purpose:  Assemble kissing hairpin candidates from core pseudoknots.
                  Search for recursive secondary structure elements and filter
                  kissing hairpins with low probabilities and high free energies. 
                  
        Input:    Pseudoknot dictionary and secondary structure element dictionaries.
        
        Return:   Filtered kissing hairpin dictionary.
    """        
    best_khps = {}

    # pseudoknot_second[key2] stores all pks which have key2 as second stem
    # pseudoknot_first[key1] stores all pks which have key1 as first stem
    for key2 in sorted(pseudoknot_second):
        if key2 in pseudoknot_first:
            pks_2 = pseudoknot_second[key2]
            pks_1 = pseudoknot_first[key2]
            for pk2 in pks_2:
                for pk1 in pks_1:                    
                    if pk2[1] < pk1[3]:     # j < m, otherwise it is a triple helix
                        
                        stem1 = pk2[0], pk2[1], pk2[2]                                                
                        stem2 = pk2[3], pk2[4], min(pk2[5],pk1[2])                        
                        stem3 = pk1[3], pk1[4], pk1[5]
                        length = stem3[1] - stem1[0] + 1                                              

                        if length < MAX_LENGTH:    # Length restriction
                                                    
                            prob1 = all_stems[stem1][1]
                            prob2 = all_stems[stem2][1]
                            prob3 = all_stems[stem3][1]
                            
                            # Set a probability threshold 
                            if prob1 + prob2 + prob3 > 0.001: 

                                stemlength1 = stem1[2]
                                stemlength2 = stem2[2] 
                                stemlength3 = stem3[2]
                            
                                l1 = stem2[0] - (stem1[0] + stemlength1)
                                l2 = (stem1[1] - stemlength1 + 1) - (stem2[0] + stemlength2)
                                l3 = stem3[0] - stem1[1] - 1
                                l4 = (stem2[1] - stemlength2 + 1) - (stem3[0] + stemlength3)
                                l5 = (stem3[1] - stemlength3) - stem2[1]                    
                                                            
                                # Loop length requirement             
                                if l1 >= 0 and l2 >= 0 and l3 >= 0 and l4 >= 0 and l5 >= 0:

                                    stack_energy1 = all_stems[stem1][2]
                                    stack_energy2 = all_stems[stem2][2]
                                    stack_energy3 = all_stems[stem3][2]                            
                                                                         
                                    kissing_hairpin = stem1, stemlength1, stem2, stemlength2, stem3, stemlength3                                                                        

                                    L1_key = stem1[0] + stemlength1, stem2[0] - 1        
                                    L2_key = stem2[0] + stemlength2, stem1[1] - stemlength1
                                    L3_key = stem1[1] + 1, stem3[0] - 1
                                    L4_key = stem3[0] + stemlength3, stem2[1] - stemlength2
                                    L5_key = stem2[1] + 1, stem3[1] - stemlength3
                                    
                                    energy_l1, result1, effective_l1 = array_traceback[L1_key[0]][L1_key[1]]
                                    energy_l2, result2, effective_l2 = array_traceback[L2_key[0]][L2_key[1]]
                                    energy_l3, result3, effective_l3 = array_traceback[L3_key[0]][L3_key[1]]
                                    energy_l4, result4, effective_l4 = array_traceback[L4_key[0]][L4_key[1]]
                                    energy_l5, result5, effective_l5 = array_traceback[L5_key[0]][L5_key[1]]                                    

                                    # Estimate kissing hairpin energy                                                  
                                    free_energy = stack_energy1 + stack_energy2 + stack_energy3 + energy_l1 + energy_l2 + energy_l3 + energy_l4 + energy_l5 
                                    entropy = UNPAIRED_NT * (effective_l1 + effective_l2 + effective_l4 + effective_l5) + UNPAIRED_NT_L3 * effective_l3
                                    estimated_energy = free_energy + entropy + INIT                                                                                                                                            
                                    
                                    # Indication of stability                                    
                                    if (estimated_energy/length) <= -0.25:              # Filtering step   
                                        
                                        # Leave at least one base unpaired in loops L1 or L2 and loops L4 and L5
                                        if (effective_l1 > 0 or effective_l2 > 0) and (effective_l4 > 0 or effective_l5 > 0):    
                                                                                                                                                                                                                                                                                                                             
                                            interval = kissing_hairpin[0][0], kissing_hairpin[4][1]
                                            # Keep best kissing hairpin for interval 
                                            if interval not in best_khps:
                                                best_khps[interval] = kissing_hairpin, estimated_energy, result1, result2, result3, result4, result5
                                            else:
                                                if best_khps[interval][1] > estimated_energy:
                                                    best_khps[interval] = kissing_hairpin, estimated_energy, result1, result2, result3, result4, result5
    
    return best_khps



