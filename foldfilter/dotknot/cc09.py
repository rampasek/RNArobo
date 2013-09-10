"""
    DotKnot: software for RNA pseudoknot prediction
    Copyright (C) 2010-2011 Jana Sperschneider	

    This program is free software; you can redistribute it and/or modify  
    it under the terms of the GNU General Public License as published by  
    the Free Software Foundation; either version 2 of the License, or     
    (at your option) any later version.
"""

import os
import math
import functions

  
def dic_caochen09(pk_dic_cc09, stem_dic):
    """ Function: dic_caochen09()

        Purpose:  Energy evaluation for pseudoknots with energy parameters CC09.
                  Note that Cao & Chen (2009) save pseudoknot in slightly 
                  different format (s1,s2,l1,l3,l2) due to loop naming. 
                  Note that no shortened stems occur here. 
                  
        Input:    Dictionary with pseudoknots where 2 <= L2 <= 6.
        
        Return:   Dictionary with pseudoknots and associated free energy. 
    """
    pk_dic_cc09_result = {}
    entropies_dic, entropies_dic_L1, entropies_dic_L3 = {}, {}, {}    

    # Store entropy parameters in dictionary in format (s1, s2, l1, l3, l2)    
    entropies = file("CaoChenParameters.txt",'r')
    for line in entropies:
        i = line.split()
        quintet = int(i[0][3:]), int(i[1][3:]), int(i[2][3:]), int(i[3][3:]), int(i[4][3:])
        entropies_dic[quintet] = float(i[5][3:])

    # For long loops L1 > 7 nt AND L3 <= 7 nt 
    # Store entropy parameters in dictionary in format (s1, s2, 'long', l3, l2)
    entropies = file("CaoChenParameters_L1.txt",'r')
    for line in entropies:
        i = line.split()
        a_b = float(i[4][7:]), float(i[5][7:])
        quintet = int(i[0][3:]), int(i[1][3:]), 'long', int(i[2][3:]), int(i[3][3:])
        entropies_dic_L1[quintet] = a_b

    # For long loops L3 > 7 nt 
    # Store entropy parameters in dictionary in format (s1, s2, l1, 'long', l2)
    entropies = file("CaoChenParameters_L3.txt",'r')
    for line in entropies:
        i = line.split()
        a_b = float(i[4][7:]), float(i[5][7:])
        quintet = int(i[0][3:]), int(i[1][3:]), int(i[2][3:]), 'long', int(i[3][3:])
        entropies_dic_L3[quintet] = a_b
      
    for pk_stem in pk_dic_cc09:

        i, j, k, l = pk_stem[2], pk_stem[3], pk_stem[5], pk_stem[6]

        stem1 = i, j        
        stemlength1 = pk_stem[4]
        stack_s1 = stem_dic[stem1][2]
        energy_s1 = stem_dic[stem1][3]        
            
        stem2 = k, l                
        stemlength2 = pk_stem[7]  
        stack_s2 = stem_dic[stem2][2]
        energy_s2 = stem_dic[stem2][3]    
                                
        l1 = k - (i + stemlength1)
        l2 = (j - stemlength1 + 1) - (k + stemlength2)            
        l3 = (l - stemlength2) - j

        if stemlength1 > 10:
            stemlength1 = 10
        if stemlength2 > 10:
            stemlength2 = 10

        # In case configuration is not defined in virtual bond model
        entropy = 0.0
        
        if l1 <= 7 and l3 >= 1 and l3 <= 7:     
            quintet = stemlength1, stemlength2, l1, l3, l2
            if quintet in entropies_dic:
                entropy = entropies_dic[quintet]
            
        elif l1 <= 7 and l3 > 7:                 
            quintet = stemlength1, stemlength2, l1, 'long', l2
            if quintet in entropies_dic_L3:
                a_b = entropies_dic_L3[quintet]
                entropy = a_b[0] * math.log(l3) + a_b[1]
                
        elif l1 > 7 and l3 >= 1 and l3 <= 7:     
            quintet = stemlength1, stemlength2, 'long', l3, l2    
            if quintet in entropies_dic_L1:
                a_b = entropies_dic_L1[quintet]
                entropy = a_b[0] * math.log(l1) + a_b[1]
                
        elif l1 > 7 and l3 > 7:                 
            quintet = stemlength1, stemlength2, l1, 'long', l2
            if quintet in entropies_dic_L3:
                a_b = entropies_dic_L3[quintet]
                entropy = a_b[0] * math.log(l3) + a_b[1]

            
        if entropy:
            pk_energy = stack_s1 + stack_s2 - (0.62 * entropy)
            if pk_energy < 0.0: 
                pk_dic_cc09_result[pk_stem] = pk_energy, stack_s1, stack_s2, l1, l3, l2, 0.62 * entropy            
            
    return pk_dic_cc09_result, entropies_dic, entropies_dic_L1, entropies_dic_L3

    
def pk_energy_reevaluation_09(pk_dic_cc09, entropies_dic, entropies_dic_L1, entropies_dic_L3):
    """ Function: pk_energy_reevaluation_09()

        Purpose:  Energy evaluation for pseudoknots with energy parameters CC09.
                  Take into account recursive structure elements which can 
                  occur in loops L1 and L3.  
                  
        Input:    Dictionary with pseudoknots where 2 <= L2 <= 6.
        
        Return:   Dictionary with recursive pseudoknots and associated free energy. 
    """      
    pk_dic_cc09_result = {}
  
    for pk_stem, values in pk_dic_cc09.items():
        entropy = 0.0
        list1 = values[1]
        list3 = values[3]
    
        if list1 or list3:
            i, j, k, l = pk_stem[2], pk_stem[3], pk_stem[5], pk_stem[6]
            stemlength1 = pk_stem[4]                        
            stemlength2 = pk_stem[7]  
                                                                  
            l1 = k - (i + stemlength1)
            l2 = (j - stemlength1 + 1) - (k + stemlength2)             
            l3 = (l - stemlength2) - j            
      
            energy_l1, energy_l3 = 0.0, 0.0
      
            if stemlength1 > 10:
                stemlength1 = 10
            if stemlength2 > 10:
                stemlength2 = 10

            effective_l1 = functions.effective_length(list1, l1) # Calculate effective loop length                
            effective_l3 = functions.effective_length(list3, l3) # Calculate effective loop length       

            if list1:   
                energy_l1 = sum([item[3] for item in list1])    # Add free energies        
                                              
            if list3:   
                energy_l3 = sum([item[3] for item in list3])    # Add free energies
                             
            if effective_l1 <= 7 and effective_l3 <= 7:   
                quintet = stemlength1, stemlength2, effective_l1, effective_l3, l2
                if quintet in entropies_dic:
                    entropy = entropies_dic[quintet]
                  
            elif effective_l1 <= 7 and effective_l3 > 7:          
                quintet = stemlength1, stemlength2, effective_l1, 'long', l2
                if quintet in entropies_dic_L3:
                    a_b = entropies_dic_L3[quintet]
                    entropy = a_b[0] * math.log(effective_l3) + a_b[1]
                
            elif effective_l1 > 7 and effective_l3 <= 7:             
                quintet = stemlength1, stemlength2, 'long', effective_l3, l2                   
                if quintet in entropies_dic_L1:
                    a_b = entropies_dic_L1[quintet]
                    entropy = a_b[0] * math.log(effective_l1) + a_b[1]     
                                       
            elif effective_l1 > 7 and effective_l3 > 7:   
                quintet = stemlength1, stemlength2, effective_l1, 'long', l2
                if quintet in entropies_dic_L3:
                    a_b = entropies_dic_L3[quintet]
                    entropy = a_b[0] * math.log(l3) + a_b[1]
            else:
                entropy = 0.0
        
            stack_s1 = values[0][1]
            stack_s2 = values[0][2]   
            
            pk_energy = stack_s1 + stack_s2 - (0.62 * entropy) + energy_l1 + energy_l3 
            
            if entropy:
                if pk_energy < 0.0:
                    if pk_energy < values[0][0]:                 
                        pk_dic_cc09_result[pk_stem] = pk_energy, values[1], values[2], values[3]
                    else:
                        pk_dic_cc09_result[pk_stem] = values[0][0], [], [], []                            
        else:
            pk_energy = pk_dic_cc09[pk_stem][0][0]
            pk_dic_cc09_result[pk_stem] = pk_energy, values[1], values[2], values[3]                
        
    return pk_dic_cc09_result
