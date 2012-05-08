""" 
    DotKnot: software for RNA pseudoknot prediction
    Copyright (C) 2010-2011 Jana Sperschneider	

    This program is free software; you can redistribute it and/or modify  
    it under the terms of the GNU General Public License as published by  
    the Free Software Foundation; either version 2 of the License, or     
    (at your option) any later version.
"""

import functions
import mwis_2d


# Dictionary for coaxial stacking
stack_dic ={
('AU','AU') : -0.9 , ('AU','CG') : -2.2 , ('AU','GC') : -2.1 , ('AU','GU') : -0.6 , ('AU','UA') : -1.1 , ('AU','UG') : -1.4 ,
('CG','AU') : -2.1 , ('CG','CG') : -3.3 , ('CG','GC') : -2.4 , ('CG','GU') : -1.4 , ('CG','UA') : -2.1 , ('CG','UG') : -2.1 ,
('GC','AU') : -2.4 , ('GC','CG') : -3.4 , ('GC','GC') : -3.3 , ('GC','GU') : -1.5 , ('GC','UA') : -2.2 , ('GC','UG') : -2.5 ,
('GU','AU') : -1.3 , ('GU','CG') : -2.5 , ('GU','GC') : -2.1 , ('GU','GU') : -0.5 , ('GU','UA') : -1.4 , ('GU','UG') :  1.3 ,
('UA','AU') : -1.3 , ('UA','CG') : -2.4 , ('UA','GC') : -2.1 , ('UA','GU') : -1.0 , ('UA','UA') : -0.9 , ('UA','UG') : -1.3 ,
('UG','AU') : -1.0 , ('UG','CG') : -1.5 , ('UG','GC') : -1.4 , ('UG','GU') :  0.3 , ('UG','UA') : -0.6 , ('UG','UG') : -0.5 }


def evaluate_pk_with_IB(pk_dic_ib, stem_dic, stems_shortened_dic, INIT, PENALTY, seq):
    """ Function: evaluate_pk_with_IB()

        Purpose:  Calculate pseudoknot free energies where one of the stems is
                  interrupted under energy model LongPK.              
                  
        Input:    Dictionary with pseudoknots where one of the stems is
                  interrupted.
        
        Return:   Dictionary with pseudoknots and associated free energy. 

    """  
    pk_dic_ib_result = {}
    
    for pk, values in pk_dic_ib.items():                    

        stack_energy_ib = values[5]                  
        energy_stems = calculate_stem_energy_sum(pk, values, stack_energy_ib, stem_dic, stems_shortened_dic)
        
        l1, l2, l3 = values[2], values[3], values[4]        
        
        left, right = values[6], values[7]        

        if l2 == 0:          # Coaxial stacking
            cs = cs_calculation(seq, pk, left, right)             
            cs = 0.75*cs      # Weighting parameter
        else:
            cs = 0.0
            
        pk_energy = energy_stems + PENALTY*(l1 + l2 + l3) + INIT + cs                 
        
        # If S1 is interrupted and > 10 bp: L3 >= 6 nt
        if pk[8] == 'iS1' and (values[9] >= 10 or values[10] >= 10):
            if l3 < 5:
                pk_energy = 100.0
        # If S2 is interrupted and > 10 bp: L1 >= 2 nt
        if pk[8] == 'iS2' and (values[9] >= 10 or values[10] >= 10):
            if l1 < 2:
                pk_energy = 100.0
        
        if pk_energy < 0.0:
            if pk_energy < values[1] and pk_energy < values[8]:            
                pk_dic_ib_result[pk] = pk_energy, l1, l2, l3, energy_stems, left, right, pk[8]
          
    return pk_dic_ib_result


def cs_calculation(seq, pk, left, right):
    """ Function: cs_calculation()

        Purpose:  Calculate coaxial stacking free energies. 
                  
        Input:    RNA sequence and two coaxially stacked stems. 
        
        Return:   Coaxial stacking free energy.         
    """
    
    # Example:
    # L2 == 0:
    # AACcUUCACCAAUUagGUUCAAAuAAGUGGU
    # ((((:::[[[[.[[[))))::::]]].]]]]
    # AACCUUCcCCAAUUagGUUCAAAuAAGuGGU
    # [[[[.[[[...((((]]].]]]]....))))
    
    if pk[8] == 'iS2':    # S1 is regular stem                
        pair1 = seq[((pk[2] - 1) + pk[4] - 1)] + seq[pk[3] - pk[4]] 
        pair2 = seq[pk[5] + left - 1] + seq[pk[5] + right - 1]
    
    if pk[8] == 'iS1':    # S2 is regular stem                
        pair2 = seq[((pk[5] - 1) + pk[7] - 1)] + seq[pk[6] - pk[7]]
        pair1 = seq[pk[2] + left - 1] + seq[pk[2] + right - 1]
    
    stack = pair1, pair2          

    if stack in stack_dic:
        coaxial_stacking = stack_dic[stack]
    else:
        coaxial_stacking = 0.0

    return coaxial_stacking


def calculate_stem_energy_sum(pk, values, stack_energy_ib, stem_dic, stems_shortened_dic):
    """ Function: evaluate_pk_with_IB()

        Purpose:  Calculate pseudoknot free energies where one of the stems is
                  interrupted under energy model LongPK.              
                  
        Input:    Dictionary with pseudoknots where one of the stems is
                  interrupted.
        
        Return:   Dictionary with pseudoknots and associated free energy. 

    """  
    if pk[8] == 'iS1':    # S1 is an interrupted stem
        s2 = pk[5], pk[6], pk[7]    
        if s2 in stems_shortened_dic:   # S2 is a shortened stem      
            stack_energy2 = stems_shortened_dic[s2][2]
            energy_stems = stack_energy2 + stack_energy_ib
        else:
            s2 = pk[5], pk[6]            
            if pk[7] == stem_dic[s2][0]:
                energy_stems = stem_dic[s2][2] + stack_energy_ib
            else:                       # S2 is a shortened stem and was filtered 
                energy_stems = 100.0                    
    
    if pk[8] == 'iS2':    # S2 is an interrupted stem     
        s1 = pk[2], pk[3], pk[4]      
        if s1 in stems_shortened_dic:   # S1 is a shortened stem        
            stack_energy1 = stems_shortened_dic[s1][2]
            energy_stems = stack_energy1 + stack_energy_ib        
        else:
            s1 = pk[2], pk[3]
            if pk[4] == stem_dic[s1][0]:
                energy_stems = stem_dic[s1][2] + stack_energy_ib
            else:                       # S1 is a shortened stem and was filtered
                energy_stems = 100.0  
    
    return energy_stems

  
def recursive_pk(pk_dic_ib, array_traceback, array_traceback_positive):
    """ Function: recursive_pk()

        Purpose:  Given core pseudoknot dictionary where one of the stems is interrupted
                  and secondary structure dictionaries, calculate set of recursive
                  elements in the three loops with MWIS calculation. 
                  
        Input:    Core pseudoknot dictionary and secondary structure dictionaries.
        
        Return:   Recursive pseudoknot dictionary.
    """ 
    for pk, values_pk in pk_dic_ib.items():
        
        result1, result2, result3 = [], [], []
        
        l1, l2, l3 = values_pk[1], values_pk[2], values_pk[3]
        left, right = values_pk[5], values_pk[6]        
        marker = pk[8]        
         
        if marker == 'iS1':
            # First Case, combine s_ib with normal stem s
            # (((...((((.xxx...))))...)))........xxx        
            loop1_start = pk[0] + left + 1
            loop1_end = pk[5] - 1            
            loop2_start = pk[5] + pk[7]
            loop2_end = pk[0] + right - 1           
            loop3_start = pk[3] + 1
            loop3_end = pk[6] - pk[7]
            
        if marker == 'iS2':
            # Second Case, combine normal stem s with s_ib
            # xxx........(((...((((.xxx...))))...)))    
            loop1_start = pk[2] + pk[4] 
            loop1_end = pk[5] - 1            
            loop2_start = pk[5] + left + 1
            loop2_end = pk[3] - pk[4]           
            loop3_start = pk[3] + 1
            loop3_end = pk[5] + right -1

        L1_key = loop1_start, loop1_end   
        L2_key = loop2_start, loop2_end   
        L3_key = loop3_start, loop3_end        
        
        energy_l1, result1, effective_l1 = array_traceback[L1_key[0]][L1_key[1]]        
        energy_l2, result2, effective_l2 = array_traceback[L2_key[0]][L2_key[1]]
        energy_l3, result3, effective_l3 = array_traceback[L3_key[0]][L3_key[1]]
                    
        # Check for sterically infeasible configurations in loops L1 and L3
        if l1 >= 9 and result1:
            energy_l1, result1, effective_l1 = functions.steric_constraints(result1, L1_key, array_traceback)
        if l3 >= 9 and result3:
            energy_l3, result3, effective_l3 = functions.steric_constraints(result3, L3_key, array_traceback)       
        
        if not result1 and l1 >= 9:
            energy_l1, result1, effective_l1 = array_traceback_positive[L1_key[0]][L1_key[1]]      
            if result1:
                energy_l1, result1, effective_l1 = functions.steric_constraints(result1, L1_key, array_traceback_positive)
        if not result3 and l3 >= 9:
            energy_l3, result3, effective_l3 = array_traceback_positive[L3_key[0]][L3_key[1]]
            if result3:
                energy_l3, result3, effective_l3 = functions.steric_constraints(result3, L3_key, array_traceback_positive)     

        pseudoknot_energy = pk_dic_ib[pk]
        pk_dic_ib[pk] = pseudoknot_energy, result1, result2, result3
        
    return pk_dic_ib

 
def re_evaluate_pk_with_IB(pk_dic_ib, INIT, PENALTY):
    """ Function: re_evaluate_pk_with_IB()

        Purpose:  Calculate pseudoknot free energies under energy model LongPK
                  where one of the stems is interrupted. Take into account
                  recursive structure elements which can occur in loops L1 and L3.               
                  
        Input:    Dictionary with pseudoknots where one of the stems is
                  interrupted.
        
        Return:   Dictionary with pseudoknots and associated free energy. 

    """ 
    for pk_stem, values in pk_dic_ib.items():

        if values[1] or values[2] or values[3]:
            list1 = values[1]
            list2 = values[2]
            list3 = values[3]                   
            l1 = values[0][1]
            l2 = values[0][2]     
            l3 = values[0][3]
            entropy_l1, energy_l1, entropy_l2, energy_l2, entropy_l3, energy_l3 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0            
            energy = values[0][4]         # Store stem energies for S1 and S2
            
            # Calculate effective loop lengths
            if list1:
                energy_l1 = sum([item[3] for item in list1])          # Add free energies
                effective_l1 = functions.effective_length(list1, l1)  # Calculate effective loop length    
            else:
                effective_l1 = l1
                
            if list2:
                energy_l2 = sum([item[3] for item in list2])          # Add free energies
                effective_l2 = functions.effective_length(list2, l2)  # Calculate effective loop length                      
            else:
                effective_l2 = l2
                
            if list3:
                energy_l3 = sum([item[3] for item in list3])          # Add free energies
                effective_l3 = functions.effective_length(list3, l3)  # Calculate effective loop length                    
            else:
                effective_l3 = l3
                
            looplength = effective_l1 + effective_l2 + effective_l3            
            entropy = PENALTY*(looplength)
            
            pk_energy = energy + entropy + energy_l1 + energy_l2 + energy_l3 + INIT

            # Only include internal loop energies if this leads to more stable pseudoknots
            if pk_energy < values[0][0]:
                pk_dic_ib[pk_stem] = pk_energy, values[1], values[2], values[3]
            else:
                energy = values[0][0]
                pk_dic_ib[pk_stem] = energy, [], [], []                                           
        else:
            energy = values[0][0]
            pk_dic_ib[pk_stem] = energy, values[1], values[2], values[3]

    return pk_dic_ib
