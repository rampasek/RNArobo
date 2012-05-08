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


def recursive_pk(pk_core_dic, seq, array_traceback, array_traceback_positive):   
    """ Function: recursive_pk()

        Purpose:  Given core pseudoknot dictionary and secondary structure dictionaries,
                  calculate set of recursive elements in the three loops with MWIS
                  calculation. 
                  
        Input:    Core pseudoknot dictionary and secondary structure dictionaries.
        
        Return:   Recursive pseudoknot dictionary.
    """     
    for pk in pk_core_dic:      

        result1, result2, result3 = [], [], []
        
        loop1_start = pk[2] + pk[4] 
        loop1_end = pk[5] - 1       
        l1 = loop1_end - loop1_start + 1        
        loop2_start = pk[5] + pk[7] 
        loop2_end = pk[3] - pk[4]
        l2 = loop2_end - loop2_start + 1           
        loop3_start = pk[3] + 1     
        loop3_end = pk[6] - pk[7]
        l3 = loop3_end - loop3_start + 1        
        
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
                
        if not result2 and l2 >= 9:
            energy_l2, result2, effective_l2 = array_traceback_positive[L2_key[0]][L2_key[1]]                  
                
        if not result3 and l3 >= 9:
            energy_l3, result3, effective_l3 = array_traceback_positive[L3_key[0]][L3_key[1]]
            if result3:
                energy_l3, result3, effective_l3 = functions.steric_constraints(result3, L3_key, array_traceback_positive)                                                           

        pseudoknot_energy = pk_core_dic[pk]
        pk_core_dic[pk] = pseudoknot_energy, result1, result2, result3

    return pk_core_dic


def pk_dic_scan_recursive(pk_core_dic):
    """ Function: pk_dic_scan_recursive()

        Purpose:  Scan pseudoknot dictionary and create three new dictionaries,
                  according to interhelix loop size L2. 
                  
        Input:    Dictionary of pseudoknots. 
        
        Return:   Three dictionaries separating the pseudoknot types. 
    """    
    pk_dic_cc06 = {}
    pk_dic_cc09 = {}
    pk_dic_longpk = {}
  
    for pk_stem in pk_core_dic:

        i, j, k, l = pk_stem[2], pk_stem[3], pk_stem[5], pk_stem[6]        
        stemlength1 = pk_stem[4]
        stemlength2 = pk_stem[7]
                                   
        l1 = k - (i + stemlength1)            
        l2 = (j - stemlength1 + 1) - (k + stemlength2)
        l3 = (l - stemlength2) - j 

        if l2 <= 1:            
            pk_dic_cc06[pk_stem] = pk_core_dic[pk_stem]
        else:
            if l2 < 7: 
                pk_dic_cc09[pk_stem] = pk_core_dic[pk_stem]
            else:
                pk_dic_longpk[pk_stem] = pk_core_dic[pk_stem]
                
    return pk_dic_cc06, pk_dic_cc09, pk_dic_longpk


def pk_filter(pk_recursive_dic):
    """ Function: pk_filter()

        Purpose:  Filter pseudoknots with high free energy. 
                  
        Input:    Dictionary of pseudoknots. 
        
        Return:   Dictionary of filtered pseudoknots. 
    """    
    for pk, values in pk_recursive_dic.items():
        
        if values[0] >= -5.25:
            del pk_recursive_dic[pk]
        else:
            length = pk[1] - pk[0] + 1
            # Length-normalized free energy
            if (values[0]/length) > -0.25:
                del pk_recursive_dic[pk]                

    return pk_recursive_dic

