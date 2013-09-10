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


def method(stem_dic_mwis, pk_recursive_dic, bulges_internal, multiloops, best_khps):
    """ Function: method()

        Purpose:  Maximum weight independent set (MWIS) calculation using the
                  set of secondary structure elements, pseudoknots and kissing
                  hairpins. Hairpin loops may contain inner structure elements. 
                  
        Input:    Dictionaries of structure elements. 
        
        Return:   Structure elements in the MWIS.
    """ 
    crossing_structures, secondary_structures, mwis_dic, candidate_list = {}, {}, {}, []
    
    for stem, values in stem_dic_mwis.items():
        if values[3] < 0.0:
            element = (stem[0], stem[1], values[0], values[1], -1*round(values[3], 2), "hp")
            candidate_list.append(element)
            
    for pk_stem, pk_energy in pk_recursive_dic.items():
        element = (pk_stem[0], pk_stem[1], pk_stem[4], pk_stem[7], -1*round(pk_energy[0], 2), "pk", pk_stem[2], pk_stem[3], pk_stem[4], pk_stem[5], pk_stem[6], pk_stem[7], pk_stem[8])
        candidate_list.append(element)
        
    for stem, values in bulges_internal.items():
        element = (stem[0], stem[1], values[0], values[1], -1*round(values[2], 2), "ib")
        candidate_list.append(element)
        
    for stem, values in multiloops.items():
        element = (stem[0], stem[1], values[0], values[1], -1*round(values[2], 2), "ml")
        candidate_list.append(element)

    for stem, values in best_khps.items():
        element = (stem[0], stem[1], values[1], 0.0, -1*round(values[1], 2), "khp")
        candidate_list.append(element)

    if candidate_list:
        candidate_list.sort()
        sorted_endpoint_list = functions.create_sorted_endpointlist(candidate_list)    
        
        for endpoint in sorted_endpoint_list:     # Scan sorted endpoints list         
            if endpoint[1] == 'r':                # If a right endpoint is scanned              
                outer_interval = candidate_list[endpoint[3] - 1]            
                if outer_interval[5] == 'hp':
                    nested, only_hp_ib_ml = find_nested(outer_interval, candidate_list)
                    # MWIS on the set of nested structure elements           
                    if nested and only_hp_ib_ml == False:         
                        endpoint_list_recursive = functions.create_sorted_endpointlist(nested)
                        result = functions.MWIS(nested, endpoint_list_recursive)
                        # Free energy sum
                        energy = outer_interval[4]
                        for element in result:              
                            energy = energy + element[4]                      
                        # Store updated free energy for outer stem
                        candidate_list[endpoint[3] - 1] = (outer_interval[0], outer_interval[1], outer_interval[2], outer_interval[3], energy, outer_interval[5])                    
                        stem = outer_interval[0], outer_interval[1], outer_interval[2]
                        # Store inner structure elements in dictionary
                        mwis_dic[stem] = result
                        
        # Main MWIS calculation
        sorted_endpoint_list_recursive = functions.create_sorted_endpointlist(candidate_list)
        result = functions.MWIS(candidate_list, sorted_endpoint_list_recursive)    

        # Free energy sum        
        energy = sum([item[4] for item in result])

        # Search for detected pseudoknots and kissing hairpins
        for element in result:
            if element[5] == 'khp' or element[5] == 'pk':
                crossing_structures[element] = element[4]
            if element[5] == 'hp' or element[5] == 'ib' or element[5] == 'ml':
                secondary_structures[element] = element[4]
            if element[5] == 'hp':  # Hairpin loop can have nested elements
                crossing_structures, secondary_structures = print_recursive(element, mwis_dic, crossing_structures, secondary_structures)
                
    return mwis_dic, crossing_structures, secondary_structures


def find_nested(interval, candidate_list):
    """ Function: find_nested()

        Purpose:  Function for finding nested intervals, given an outer interval.
                  Outer stem is a hairpin loop. Watch out for special case of
                  inner nested pseudoknot. 
                  
        Input:    Outer interval and set of candidate intervals. 
        
        Return:   Nested intervals and marker which indicates whether pseudoknots
                  are in the interval set. 
    """  
    only_hp_ib_ml = True
    result = []
    
    interval_left = interval[0] + interval[2]  
    interval_right = interval[1] - interval[2]
    
    for compare_interval in candidate_list:
        if compare_interval[0] >= interval_left:                              
            if compare_interval[1] <= interval_right:
                # Special case for nested pseudoknot
                # Add one base to each side as a safeguard
                # {{{.(((..[[[.)))....]]].}}}
                if compare_interval[5] == 'pk':
                    if compare_interval[0] > interval_left:
                        if compare_interval[1] < interval_right:
                            result.append(compare_interval)
                            only_hp_ib_ml = False
                else:
                    result.append(compare_interval)
                    if compare_interval[5] == 'khp':
                        only_hp_ib_ml = False
                        
    return result, only_hp_ib_ml


def print_recursive(element, mwis_dic, crossing_structures, secondary_structures):
    """ Function: print_recursive()

        Purpose:  Given the MWIS result, look for nested pseudoknots, kissing 
                  hairpins and secondary structures recursively. Store free 
                  energy information. 
                  
        Input:    Outer structure and set of candidate structures. 
        
        Return:   Nested structures with free energies.  
    """ 
    structure = element[0], element[1], element[2]

    if structure in mwis_dic:
        stem_internal_list = mwis_dic[element[0], element[1], element[2]]

        for item in stem_internal_list:
            if item[5] == 'khp' or item[5] == 'pk': 
                crossing_structures[item] = item[4]
            if item[5] == 'hp' or item[5] == 'ib' or item[5] == 'ml': 
                secondary_structures[item] = item[4]
            if item[5] == 'hp':
                substructure = item[0], item[1], item[2]
                if substructure in mwis_dic:
                    crossing_structures, secondary_structures = print_recursive(substructure, mwis_dic, crossing_structures, secondary_structures)
                    
    return crossing_structures, secondary_structures 


def pk_khp_structures(seq, crossing_structures, best_khps, stem_dic, bulges_internal, multiloops, pk_recursive_dic, pk_dic_ib):
    """ Function: pk_khp_structures()

        Purpose:  Obtain dot-bracket notation for set of pseudoknots,
                  with kissing hairpins. 
                  
        Input:    Dictionaries of structure elements. 
        
        Return:   Pseudoknots and kissing hairpins in dot-bracket notation.
    """ 
    pseudoknot_list = []    

    for pk in crossing_structures:

        if pk[5] == 'khp':            
            i, n = pk[0], pk[1]
            energy = pk[2]
            for khp, khp_value in best_khps.items():
                if khp_value[0][0][0] == i and khp_value[0][4][1] == n and khp_value[1] == energy:
                    kissing_hairpin = khp_value[0]
                    value = khp_value[1:]
            i, j, stemlength1 = kissing_hairpin[0][0], kissing_hairpin[0][1], kissing_hairpin[1]
            k, l, stemlength2 = kissing_hairpin[2][0], kissing_hairpin[2][1], kissing_hairpin[3]
            m, n, stemlength3 = kissing_hairpin[4][0], kissing_hairpin[4][1], kissing_hairpin[5]
            # Find the recursive structure elements            
            recursive_loop1 = value[1]
            recursive_loop2 = value[2]
            recursive_loop3 = value[3]
            recursive_loop4 = value[4]
            recursive_loop5 = value[5]

            l1 = kissing_hairpin[2][0] - (kissing_hairpin[0][0] + kissing_hairpin[1])      
            l2 = (kissing_hairpin[0][1] - kissing_hairpin[1]) - (kissing_hairpin[2][0] + kissing_hairpin[3]) + 1
            l3 = kissing_hairpin[4][0] - (kissing_hairpin[0][1] + 1)
            l4 = (kissing_hairpin[2][1] - kissing_hairpin[3]) - (kissing_hairpin[4][0] + kissing_hairpin[5]) + 1
            l5 = (kissing_hairpin[4][1] - kissing_hairpin[5]) - kissing_hairpin[2][1]
            
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
                pk_structure = pk_tools.add_recursive_elements(i, recursive_loop1, pk_structure, stem_dic, bulges_internal, multiloops)
                
            if recursive_loop2:
                pk_structure = pk_tools.add_recursive_elements(i, recursive_loop2, pk_structure, stem_dic, bulges_internal, multiloops)

            if recursive_loop3:
                pk_structure = pk_tools.add_recursive_elements(i, recursive_loop3, pk_structure, stem_dic, bulges_internal, multiloops)

            if recursive_loop4:
                pk_structure = pk_tools.add_recursive_elements(i, recursive_loop4, pk_structure, stem_dic, bulges_internal, multiloops)

            if recursive_loop5:
                pk_structure = pk_tools.add_recursive_elements(i, recursive_loop5, pk_structure, stem_dic, bulges_internal, multiloops)

            pk_structure = ''.join(pk_structure)
            pseudoknot = [pk[0], pk[1], -1*float(pk[4]), pk[5], pk_seq, pk_structure]             
            pseudoknot_list.append(pseudoknot)

        if pk[5] == 'pk':
            i, j, stemlength1 = pk[6], pk[7], pk[8]        
            k, l, stemlength2 = pk[9], pk[10], pk[11]
            marker = pk[12]                               
            key = i, l, i, j, stemlength1, k, l, stemlength2, marker
            energy = round(pk_recursive_dic[key][0],2)
            recursive_loop1 = pk_recursive_dic[key][1]
            recursive_loop2 = pk_recursive_dic[key][2]
            recursive_loop3 = pk_recursive_dic[key][3]
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
                pk_structure = pk_tools.add_recursive_elements(i, recursive_loop1, pk_structure, stem_dic, bulges_internal, multiloops)
                
            if recursive_loop2:
                pk_structure = pk_tools.add_recursive_elements(i, recursive_loop2, pk_structure, stem_dic, bulges_internal, multiloops)

            if recursive_loop3:
                pk_structure = pk_tools.add_recursive_elements(i, recursive_loop3, pk_structure, stem_dic, bulges_internal, multiloops)

            pk_structure = ''.join(pk_structure)
            pseudoknot = [pk[0], pk[1], -1*float(pk[4]), pk[5], pk_seq, pk_structure]             
            pseudoknot_list.append(pseudoknot)
        
    pseudoknot_list.sort()                                     

    return pseudoknot_list    

 
def assemble_global_structure(seq, secondary_structures, pseudoknot_list):
    """ Function: assemble_global_structure()

        Purpose:  Assemble the global structure including pseudoknots.
                  
        Input:    List of predicted pseudoknots and secondary structures. 
        
        Return:   Global structure in dot-bracket notation.
    """
    predicted_global_structure = ['.' for i in xrange(len(seq))]

    # Insert secondary structure elements
    if secondary_structures:        
        for element, energy in sorted(secondary_structures.items()):

            # Insert detected hairpin loops
            if element[5] == 'hp':  
                structure = ['.' for i in xrange(element[1] - element[0] + 1)]
                for i in xrange(element[2]):    
                    structure[i] = '('
                    structure[len(structure) - i - 1] = ')'

                predicted_global_structure[element[0]-1:element[1]] = structure                

            # Insert detected bulge loop, internal loop and multiloop structures
            else:               
                structure = element[3].replace(':','')                        
                predicted_global_structure[element[0]-1:element[1]] = list(structure)
              
    # Insert detected pseudoknots and kissing hairpins
    if pseudoknot_list:
        for pk in pseudoknot_list:                
            predicted_global_structure[pk[0]-1:pk[1]] = list(pk[5])

    predicted_global_structure = ''.join(predicted_global_structure)        
      
    return predicted_global_structure

