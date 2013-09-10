"""
    DotKnot: software for RNA pseudoknot prediction
    Copyright (C) 2010-2011 Jana Sperschneider	

    This program is free software; you can redistribute it and/or modify  
    it under the terms of the GNU General Public License as published by  
    the Free Software Foundation; either version 2 of the License, or     
    (at your option) any later version.
"""

import functions

def output_best_ratio(pk_not_filtered, best_khps, NUMBER_OF_BEST_PKS):    
    """ Function: output_best_ratio()

        Purpose:  Scan the dictionary of pseudoknots and kissing
                  hairpins and output the x best pseudoknots
                  in terms of energy to length ratio.
                  
        Input:    Number x and pseudoknot/kissing hairpin dictionary. 
        
        Return:   Best x pseudoknots. 
    """ 
    candidates = pk_not_filtered.copy()
    best_pks = []
    
    if best_khps:
        for khp in best_khps:
            stems = best_khps[khp][0]
            result = best_khps[khp][1:]
            candidates[khp] = result    
    
    for i in xrange(NUMBER_OF_BEST_PKS):
        minimum = 0.0
        if candidates:    
            for pk, values in candidates.items():
                length_pk = pk[1] - pk[0] + 1
                energy_pk = values[0]
                ratio = energy_pk/length_pk                
                if ratio < minimum:
                    minimum = ratio
                    best_pk = pk
            best_pks.append(best_pk)
            del candidates[best_pk]                   
        
    return best_pks


def output_best_energy(pk_not_filtered, best_khps, NUMBER_OF_BEST_PKS):    
    """ Function: output_best_energy()

        Purpose:  Scan the dictionary of pseudoknots and kissing
                  hairpins and output the x best pseudoknots in
                  terms of free energy.
                  
        Input:    Number x and pseudoknot/kissing hairpin dictionary. 
        
        Return:   Best x pseudoknots. 
    """ 
    candidates = pk_not_filtered.copy()
    best_pks = []
    
    if best_khps:
        for khp in best_khps:
            stems = best_khps[khp][0]
            result = best_khps[khp][1:]
            candidates[khp] = result
    
    for i in xrange(NUMBER_OF_BEST_PKS):
        minimum = 0.0
        if candidates:    
            for pk, values in candidates.items():
                energy_pk = values[0]
                if energy_pk < minimum:
                    minimum = energy_pk
                    best_pk = pk
            best_pks.append(best_pk)
            del candidates[best_pk]                   
                
    return best_pks


def dot_bracket(seq, best_pks, stem_dic, stems_shortened_dic, bulges_internal, multiloops, pk_not_filtered, best_khps):
    """ Function: dot_bracket()

        Purpose:  Obtain dot-bracket notation for set of pseudoknots,
                  with kissing hairpins. 
                  
        Input:    Dictionaries of structure elements. 
        
        Return:   Pseudoknots in dot-bracket notation.
    """
    pseudoknot_list = []    

    for pk in best_pks:
        # Kissing hairpin 
        if len(pk) == 2:            
            kissing_hairpin = best_khps[pk][0]
            i, j, stemlength1 = kissing_hairpin[0][0], kissing_hairpin[0][1], kissing_hairpin[1]
            k, l, stemlength2 = kissing_hairpin[2][0], kissing_hairpin[2][1], kissing_hairpin[3]
            m, n, stemlength3 = kissing_hairpin[4][0], kissing_hairpin[4][1], kissing_hairpin[5]            
            energy = round(best_khps[pk][1],2)
            recursive_loop1 = best_khps[pk][2]
            recursive_loop2 = best_khps[pk][3]
            recursive_loop3 = best_khps[pk][4]
            recursive_loop4 = best_khps[pk][5]
            recursive_loop5 = best_khps[pk][6]            
                                                   
            L1_start = kissing_hairpin[0][0] + kissing_hairpin[1]
            L1_end = kissing_hairpin[2][0] - 1
            l1 = L1_end - L1_start + 1        
            L2_start = kissing_hairpin[2][0] + kissing_hairpin[3]
            L2_end = kissing_hairpin[0][1] - kissing_hairpin[1]
            l2 = L2_end - L2_start + 1      
            L3_start = kissing_hairpin[0][1] + 1
            L3_end = kissing_hairpin[4][0] - 1
            l3 = L3_end - L3_start + 1             
            L4_start = kissing_hairpin[4][0] + kissing_hairpin[5]
            L4_end = kissing_hairpin[2][1] - kissing_hairpin[3]
            l4 = L4_end - L4_start + 1       
            L5_start = kissing_hairpin[2][1] + 1
            L5_end = kissing_hairpin[4][1] - kissing_hairpin[5]
            l5 = L5_end - L5_start + 1

            pk_seq = seq[int(i-1):int(n)]

            pk_structure = ['(' for x in xrange(stemlength1)]
            pk_structure = pk_structure + ['.' for x in xrange(l1)]
            pk_structure = pk_structure + ['[' for x in xrange(stemlength2)]
            pk_structure = pk_structure + ['.' for x in xrange(l2)]
            pk_structure = pk_structure + [')' for x in xrange(stemlength1)]
            pk_structure = pk_structure + ['.' for x in xrange(l3)]
            pk_structure = pk_structure + ['(' for x in xrange(stemlength3)]
            pk_structure = pk_structure + ['.' for x in xrange(l4)]
            pk_structure = pk_structure + [']' for x in xrange(stemlength2)]            
            pk_structure = pk_structure + ['.' for x in xrange(l5)]
            pk_structure = pk_structure + [')' for x in xrange(stemlength3)]

            # Now add recursive structure elements
            if recursive_loop1:
                pk_structure = add_recursive_elements(i, recursive_loop1, pk_structure, stem_dic, bulges_internal, multiloops)
            if recursive_loop2:
                pk_structure = add_recursive_elements(i, recursive_loop2, pk_structure, stem_dic, bulges_internal, multiloops)
            if recursive_loop3:
                pk_structure = add_recursive_elements(i, recursive_loop3, pk_structure, stem_dic, bulges_internal, multiloops)
            if recursive_loop4:
                pk_structure = add_recursive_elements(i, recursive_loop4, pk_structure, stem_dic, bulges_internal, multiloops)
            if recursive_loop5:
                pk_structure = add_recursive_elements(i, recursive_loop5, pk_structure, stem_dic, bulges_internal, multiloops)

            pk_structure = ''.join(pk_structure)                                        
            pseudoknot = [i, n, energy, 'khp', pk_seq, pk_structure]             
            pseudoknot_list.append(pseudoknot)

        else:
            # H-type pseudoknots
            i, j, stemlength1 = pk[2], pk[3], pk[4]        
            k, l, stemlength2 = pk[5], pk[6], pk[7]
            marker = pk[8]
            key = i, l, i, j, stemlength1, k, l, stemlength2, marker
            energy = round(pk_not_filtered[key][0],2)
            recursive_loop1 = pk_not_filtered[key][1]
            recursive_loop2 = pk_not_filtered[key][2]
            recursive_loop3 = pk_not_filtered[key][3]
            looplength1 = k - (i + stemlength1)            
            looplength2 = (j - stemlength1 + 1) - (k + stemlength2)
            looplength3 = (l - stemlength2) - j 
            pk_seq = seq[int(i-1):int(l)]
                
            if marker == 'r':            
                # Now assemble the core pseudoknot structure with regular stems
                pk_structure = ['(' for x in xrange(stemlength1)]
                pk_structure = pk_structure + ['.' for x in xrange(looplength1)]
                pk_structure = pk_structure + ['[' for x in xrange(stemlength2)]
                pk_structure = pk_structure + ['.' for x in xrange(looplength2)]
                pk_structure = pk_structure + [')' for x in xrange(stemlength1)]
                pk_structure = pk_structure + ['.' for x in xrange(looplength3)]
                pk_structure = pk_structure + [']' for x in xrange(stemlength2)]

            else:                            
                # Case 1: stem S1 is interrupted
                if marker == 'iS1':
                    stem1 = i, j
                    structure_stem1 = bulges_internal[stem1][1]            
                    # Delete dangling ends ':'
                    structure_stem1 = structure_stem1.replace(':','')

                    pk_structure = list(structure_stem1)
                    start = k - i               
                    for x in xrange(stemlength2):
                        pk_structure = pk_structure[0:start] + list('[') + pk_structure[start+1:]
                        start = start + 1
                    pk_structure = pk_structure + ['.' for x in xrange(looplength3)]
                    pk_structure = pk_structure + [']' for x in xrange(stemlength2)]
                    
                # Case 2: stem S2 is interrupted, change brackets '(' to '[' and ')' to ']'
                if marker == 'iS2':
                    stem2 = k, l
                    structure_stem2 = bulges_internal[stem2][1]                    
                    # Delete dangling ends ':'
                    structure_stem2 = structure_stem2.replace(':','')
                    structure_stem2 = structure_stem2.replace('(','[')
                    structure_stem2 = structure_stem2.replace(')',']')
                    
                    pk_structure = list(structure_stem2)
                    pk_structure = ['.' for x in xrange(looplength1)] + pk_structure 
                    pk_structure = ['(' for x in xrange(stemlength1)] + pk_structure 
                    end = j - i
                    for x in xrange(stemlength1):
                        pk_structure = pk_structure[0:end] + list(')') + pk_structure[end+1:]
                        end = end - 1                
                    
            # Now add recursive structure elements
            if recursive_loop1:
                pk_structure = add_recursive_elements(i, recursive_loop1, pk_structure, stem_dic, bulges_internal, multiloops)
                
            if recursive_loop2:
                pk_structure = add_recursive_elements(i, recursive_loop2, pk_structure, stem_dic, bulges_internal, multiloops)

            if recursive_loop3:
                pk_structure = add_recursive_elements(i, recursive_loop3, pk_structure, stem_dic, bulges_internal, multiloops)

            pk_structure = ''.join(pk_structure)
            pseudoknot = [i, l, energy, pk_seq, pk_structure]             
            pseudoknot_list.append(pseudoknot)                                          

    return pseudoknot_list


def add_recursive_elements(i, recursive_loop, pk_structure, stem_dic, bulges_internal, multiloops):
    """ Function: add_recursive_elements()

        Purpose:  For a core pseudoknot, add recursive secondary structure
                  elements for a given loop L1, L2 or L3.
                  
        Input:    Pseudoknot and elements for a loop. 
        
        Return:   Pseudoknots in dot-bracket notation.
    """

    # Calculate MWIS to avoid overlapping base pairs
    if recursive_loop:      
        candidate_list = []
        for item in recursive_loop:
            # Weights need to be positive
            element = (item[0], item[1], -1.0*item[2], -1.0*item[2], -1.0*item[2], item[4])
            candidate_list.append(element)
        
        candidate_list.sort()        
        sorted_endpoint_list = functions.create_sorted_endpointlist(candidate_list)
        result = functions.MWIS(candidate_list, sorted_endpoint_list)

        mwis_set = []
        overlapping_set = []
        
        for item in result:
            item_format = (item[0], item[1], -1.0*item[2], -1.0*item[2], item[5])
            mwis_set.append(item_format)
        for item in recursive_loop:
            if item not in mwis_set:
                overlapping_set.append(item)
    
    for element in mwis_set:
        start = element[0] - i
        end = element[1] - i
        
        if element[4] == 'hp':
            element_length = stem_dic[element[0],element[1]][0] + 1
            for counter in xrange(1,element_length):
                pk_structure = pk_structure[0:start] + list('(') + pk_structure[start+1:]
                pk_structure = pk_structure[0:end] + list(')') + pk_structure[end+1:]
                start = start + 1
                end = end - 1
                
        elif element[4] == 'ib':
            structure_ib = bulges_internal[element[0],element[1]][1]                    
            structure_ib = structure_ib.replace(':','') # Cut off dangling ends
            pk_structure = pk_structure[0:start] + list(structure_ib) + pk_structure[end+1:]
            
        elif element[4] == 'ml':
            structure_ml = multiloops[element[0],element[1]][1]                    
            structure_ml = structure_ml.replace(':','') # Cut off dangling ends
            pk_structure = pk_structure[0:start] + list(structure_ml) + pk_structure[end+1:]

    for element in overlapping_set:
        start = element[0] - i + 1
        end = element[1] - i - 1
        
        if element[4] == 'hp':
            element_length = stem_dic[element[0],element[1]][0] 
            for counter in xrange(1,element_length):
                pk_structure = pk_structure[0:start] + list('(') + pk_structure[start+1:]
                pk_structure = pk_structure[0:end] + list(')') + pk_structure[end+1:]
                start = start + 1
                end = end - 1
                
        elif element[4] == 'ib':
            structure_ib = bulges_internal[element[0],element[1]][1]                    
            structure_ib = structure_ib.replace(':','') # Cut off dangling ends
            structure_ib = structure_ib[1:-1]
            pk_structure = pk_structure[0:start] + list(structure_ib) + pk_structure[end+1:]
            
        elif element[4] == 'ml':
            structure_ml = multiloops[element[0],element[1]][1]                    
            structure_ml = structure_ml.replace(':','') # Cut off dangling ends
            structure_ml = structure_ml[1:-1]
            pk_structure = pk_structure[0:start] + list(structure_ml) + pk_structure[end+1:]

    return pk_structure
