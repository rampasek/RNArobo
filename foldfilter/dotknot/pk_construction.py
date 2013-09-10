"""
    DotKnot: software for RNA pseudoknot prediction
    Copyright (C) 2010-2011 Jana Sperschneider	

    This program is free software; you can redistribute it and/or modify  
    it under the terms of the GNU General Public License as published by  
    the Free Software Foundation; either version 2 of the License, or     
    (at your option) any later version.
""" 

L1_LOWER = 1
L1_UPPER = 100

L2_LOWER = 0 
L2_UPPER = 50

L3_LOWER = 2 
L3_UPPER = 100


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
    
    
def build_pseudoknots(stem_dic):
    """ Function: build_pseudoknots()

        Purpose:  Assemble H-type pseudoknot from two stems. Allow
                  certain base pair overlap and store shortened stems
                  in a new stem dictionary. Use a different key for
                  shortened stems, i.e. (i, j, Length) -> Length, 0.0
                  
        Input:    Dictionary of stems. 
        
        Return:   A dictionary of H-type pseudoknots and a dictionary of 
                  shortened stems. 
    """    
    pk_dic, stems_shortened = {}, {}
    stem_dic_list = stem_dic.items()
    stem_dic_list.sort() 
                         
    for index, ((i,j), S1_values) in enumerate(stem_dic_list):
        for (k,l), S2_values in stem_dic_list[index:]:      

            #i, j = S1
            #k, l = S2 
      
            if (l - i + 1) > 15 and (l - i + 1) < 400:
            
                stemlength1 = S1_values[0]         
                stemlength2 = S2_values[0]
                l1 = k - (i + stemlength1)            
                l2 = (j - stemlength1 + 1) - (k + stemlength2)
                l3 = (l - stemlength2) - j             
                
                if loops_fulfilled(l1, l2, l3):
                    pk_dic = add_pseudoknot(i, j, k, l, stemlength1, stemlength2, l1, l2, l3, pk_dic)
                      
                # Case that overlap of 1 bp or 2 bps occurs at L2
                if (l2 == -1 or l2 == -2) and (stemlength1 > 3 or stemlength2 > 3):
        
                    marker, l1, l2, l3, stemlength1, stemlength2 = overlap(l1, l2, l3, stemlength1, stemlength2)

                    if marker:           
                        if loops_fulfilled(l1, l2, l3):
                            pk_dic = add_pseudoknot(i, j, k, l, stemlength1, stemlength2, l1, l2, l3, pk_dic)            

                            # Add shortened stems to dictionary                            
                            s1_shortened = i, j, stemlength1     
                            key = i, j
                            original_length = stem_dic[key][0]                            
                            if stemlength1 != original_length:
                                stems_shortened[s1_shortened] = stemlength1, 0.0                                         
                                
                            s2_shortened = k, l, stemlength2
                            key = k, l
                            original_length = stem_dic[key][0]      
                            if stemlength2 != original_length:
                                stems_shortened[s2_shortened] = stemlength2, 0.0    

    return stems_shortened, pk_dic


def add_pseudoknot(i, j, k, l, stemlength1, stemlength2, l1, l2, l3, pk_dic):
    """ Function: add_pseudoknot()

        Purpose:  Add simple H-type pseudoknot constructed with two stems to
                  H-type pseudoknot dictionary. Use default free energy of zero.
                  
        Input:    A pseudoknot and a dictionary of H-type pseudoknots. 
        
        Return:   An updated dictionary of H-type pseudoknots.
    """       
    key = i, l, i, j, stemlength1, k, l, stemlength2, 'r'                                 
    pk_dic[key] = 0.0
  
    return pk_dic


def overlap(l1, l2, l3, stemlength1, stemlength2):
    """ Function: overlap()

        Purpose:  Handle case that overlap of 1 bp or 2 bps occurs at L2.
                  A marker is used to indicate whether overlap can be resolved. 
                  
        Input:    Loop and stem lengths of a pseudoknot 
        
        Return:   Marker and updated loop and stem lengths. 
    """
    marker = True
    
    while l2 != 0:
        if l1 == 0:                
            if stemlength1 > 3:         # Cut S1
                stemlength1, l1, l2 = stemlength1 - 1, l1 + 1, l2 + 1                                    
            else:
                l2, marker = 0, False   # It is not possible to resolve the overlap
          
        elif l3 == 1: 
            if stemlength2 > 3:         # Cut S2
                stemlength2, l2, l3 = stemlength2 - 1, l2 + 1, l3 + 1                                  
            else:
                l2, marker = 0, False   # It is not possible to resolve the overlap                          

        elif l1 >= 1 and l3 >= 2:
            if stemlength1 > 3:                 # Cut S1
                stemlength1, l1, l2 = stemlength1 - 1, l1 + 1, l2 + 1                   
            elif stemlength2 > 3:               # Cut S2
                stemlength2, l2, l3 = stemlength2 - 1, l2 + 1, l3 + 1    
            else:
                l2, marker = 0, False   # It is not possible to resolve the overlap
          
        else:
            l2, marker = 0, False       # It is not possible to resolve the overlap          

    return marker, l1, l2, l3, stemlength1, stemlength2
    
    
def build_pseudoknots_ib(bulges_internal, stem_dic):
    """ Function: build_pseudoknots_ib()

        Purpose:  Assemble H-type pseudoknot from two stems, where one of 
                  the stems can be interrupted. Allow certain base pair 
                  overlap and store shortened stems in a new stem dictionary. 
                  Use a different key, i.e. (i, j, Length) -> Length, 0.0
                  
        Input:    Dictionaries of stems and interrupted stems. 
        
        Return:   A dictionary of interrupted H-type pseudoknots and a dictionary of 
                  shortened stems. 
    """
    pk_dic_ib = {}    
    stems_shortened_ib = {}

    for stem_ib, values_ib in bulges_internal.items():
        structure = values_ib[1]            # :((((.(((((.....))))).))))                    
        left = structure.rfind('(') - 1     # Position of last bracket '(' 
        right = structure.find(')') - 1     # Position of first bracket ')'     

        if stem_ib[0] == 1:                 # Case of dangling end 
            left = left + 1
            right = right + 1

        stack_energy_ib = values_ib[3]      # bulges_internal[key] = length, dot_bracket, float(free_energy), float(stack_energy)
      
        left_index = stem_ib[0] + left        
        right_index = stem_ib[0] + right - 1
        stemlength_left  = left_index - stem_ib[0] + 1
        stemlength_right = stem_ib[1] - right_index 

        for stem, values in stem_dic.items():                

            if stack_energy_ib + values[2] < -12.0:
                length = values[0]
                energy_sum = values_ib[2] + values[2]
                            
                """ First Case, combine s_ib with normal stem s
                    (((...((((.xxx...))))...)))........xxx
                            L1    L2             L3
                """           
                if stem[0] >= left_index:                    
                    l1 = stem[0] - (stem_ib[0] + left) - 1
                    l2 = (stem_ib[0] + right) - (stem[0] + length)
                    l3 = (stem[1] - length) - stem_ib[1] 
                    stem_length_eff = structure.rfind(')') - right + 1                                   
                    
                    if loops_fulfilled(l1, l2, l3):                            
                        key = stem_ib[0], stem[1], stem_ib[0], stem_ib[1], values_ib[0], stem[0], stem[1], length, 'iS1'                                                                                                   
                        pk_dic_ib[key] = energy_sum, values_ib[2], l1, l2, l3, stack_energy_ib, left, right, values[3], stemlength_left, stemlength_right 

                    """ Case that overlap of 1 bp occurs at loop L2
                        (((...((((......))))...)))
                        ............xxxxx............xxxxx
                    """                                     
                    if l2 == -1 and length > 3:   # Allow overlap of 1 nt                          
                        length = length - 1         # Cut regular stem S2
                        l2 = l2 + 1
                        l3 = l3 + 1
                        if loops_fulfilled(l1, l2, l3):                                                                        
                            key = stem_ib[0], stem[1], stem_ib[0], stem_ib[1], values_ib[0], stem[0], stem[1], length, 'iS1'
                            pk_dic_ib[key] = energy_sum, values_ib[2], l1, l2, l3, stack_energy_ib, left, right, values[3], stemlength_left, stemlength_right                            
                            s2_shortended = stem[0], stem[1], length
                            stems_shortened_ib[s2_shortended] = length, 0.0  
                                            
                """ Second Case, combine normal stem s with s_ib
                    xxx........(((...((((.xxx...))))...)))
                          L1             L2   L3
                """
                if stem[0] < stem_ib[0]:                    
                    l1 = stem_ib[0] - (stem[0] + length)
                    l2 = (stem[1] - length) - (stem_ib[0] + left)
                    l3 = (stem_ib[0] + right) - stem[1] - 1 
                    stem_length_eff = left

                    if loops_fulfilled(l1, l2, l3):                                                                                                                                        
                        key = stem[0], stem_ib[1], stem[0], stem[1], length, stem_ib[0], stem_ib[1], values_ib[0], 'iS2'                                
                        pk_dic_ib[key] = energy_sum, values_ib[2], l1, l2, l3, stack_energy_ib, left, right, values[3], stemlength_left, stemlength_right
            
                    """ Case that overlap of 1 bp occurs at loop L2
                        ((((............))))
                        .......xxx..xxxxx........xxxx..xxxxx
                    """
                    if l2 == -1 and length > 3:   # Allow overlap of 1 nt                                
                        length = length - 1         # Cut regular stem S1                                    
                        l1 = l1 + 1
                        l2 = l2 + 1                                    
                        if loops_fulfilled(l1, l2, l3):                                                                    
                            key = stem[0], stem_ib[1], stem[0], stem[1], length, stem_ib[0], stem_ib[1], values_ib[0], 'iS2'
                            pk_dic_ib[key] = energy_sum, values_ib[2], l1, l2, l3, stack_energy_ib, left, right, values[3], stemlength_left, stemlength_right                            
                            s1_shortended = stem[0], stem[1], length
                            stems_shortened_ib[s1_shortended] = length, 0.0

    return stems_shortened_ib, pk_dic_ib

                
def pk_dic_scan(pk_dic, stems_dic, stems_shortened_dic):
    """ Function: pk_dic_scan()

        Purpose:  Task: Scan pseudoknot dictionary and create three new dictionaries,
                  according to interhelix loop size L2. Watch out for pseudoknot candidates
                  where one of the shortened stems S1 or S2 may have been deleted because of 
                  high free energy in the stem filtering step.
                  
        Input:    Dictionaries of pseudoknots, stems and shortened stems. 
        
        Return:   Three dictionaries separating the pseudoknot types. 
    """  
    pk_dic_cc06 = {}
    pk_dic_cc09 = {}
    pk_dic_longpk = {}

    for pk_stem in pk_dic:
  
        i, j, k, l = pk_stem[2], pk_stem[3], pk_stem[5], pk_stem[6]        
    
        stem1 = i, j
        length_in_pk = pk_stem[4]
        exists_stem1 = check_if_stem_exists(stem1, length_in_pk, stems_dic, stems_shortened_dic)
        
        stem2 = k, l  
        length_in_pk = pk_stem[7]  
        exists_stem2 = check_if_stem_exists(stem2, length_in_pk, stems_dic, stems_shortened_dic)    
                                   
        l1 = k - (i + pk_stem[4])
        l2 = (j - pk_stem[4] + 1) - (k + pk_stem[7])            
        l3 = (l - pk_stem[7]) - j 
            
        if exists_stem1 and exists_stem2:
            if l2 <= 1:
                pk_dic_cc06[pk_stem] = 0.0
            else:
                if l2 < 7: 
                    pk_dic_cc09[pk_stem] = 0.0
                else:
                    pk_dic_longpk[pk_stem] = 0.0
          
    return pk_dic_cc06, pk_dic_cc09, pk_dic_longpk


def check_if_stem_exists(stem, length_in_pk, stems_dic, stems_shortened_dic):
    """ Function: check_if_stem_exists()

        Purpose:  Check if a shortened stem was deleted in the filtering step.
                  
        Input:    Stem and its length, dictionaries of stems and shortened stems. 
        
        Return:   True or False. 
    """  
    exists = False

    # Get original stem length and compare to length in pseudoknot
    stemlength = stems_dic[stem][0]      

    if stemlength != length_in_pk:  # It is a shortened stem            
        stemlength = length_in_pk
        key_short = stem[0], stem[1], length_in_pk
      
        if key_short in stems_shortened_dic:
            exists = True
        # If not found, the shortened stem was deleted before because of high free energy
        else:
            exists = False
      
    else:                           # It is not a shortened stem and thus exists by default
        exists = True

    return exists
  
