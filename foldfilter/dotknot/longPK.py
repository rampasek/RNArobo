""" DotKnot: software for RNA pseudoknot prediction
    Copyright (C) 2010-2011 Jana Sperschneider	

    This program is free software; you can redistribute it and/or modify  
    it under the terms of the GNU General Public License as published by  
    the Free Software Foundation; either version 2 of the License, or     
    (at your option) any later version.
"""

import functions


def dic_longpks(pk_dic_longpk, stem_dic, INIT, PENALTY):
    """ Function: dic_longpks()

        Purpose:  Calculate pseudoknot free energies under energy model LongPK.
                  Note that no shortened stems will occur here. 
                  
        Input:    Dictionary with pseudoknots where L2 >= 7.
        
        Return:   Dictionary with pseudoknots and associated free energy. 

    """    
    pk_dic_longpk_result = {}    

    for pk_stem in pk_dic_longpk:

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

        looplength = l1 + l2 + l3        
      
        entropy = PENALTY*(looplength)
      
        pk_energy = stack_s1 + stack_s2 + entropy + INIT
      
        if pk_energy < 0.0:            
            pk_dic_longpk_result[pk_stem] = pk_energy, stack_s1, stack_s2, l1, l2, l3, entropy, looplength
                
    return pk_dic_longpk_result


def pk_energy_reevaluation_long(pk_dic_longpk, INIT, PENALTY):
    """ Function: pk_energy_reevaluation_long()

        Purpose:  Calculate pseudoknot free energies under energy model LongPK. 
                  Take into account recursive structure elements which can 
                  occur in loops L1, L2 and L3.               
                  
        Input:    Dictionary with recursive pseudoknots where L2 >= 7.
        
        Return:   Dictionary with recursive pseudoknots and associated free energy. 

    """  
    for pk_stem, values in pk_dic_longpk.items():

        list1 = values[1]
        list2 = values[2]
        list3 = values[3]
        
        if list1 or list2 or list3:
            i, j, k, l = pk_stem[2], pk_stem[3], pk_stem[5], pk_stem[6]
            
            stem1 = i, j        
            stemlength1 = pk_stem[4]
            stem2 = k, l                
            stemlength2 = pk_stem[7]  
            stack_s1 = values[0][1]
            stack_s2 = values[0][2]
            
            l1 = values[0][3]
            l2 = values[0][4]     
            l3 = values[0][5]                                                
            
            entropy_l1, energy_l1 = 0.0, 0.0
            entropy_l2, energy_l2 = 0.0, 0.0
            entropy_l3, energy_l3 = 0.0, 0.0            
            
            energy = stack_s1 + stack_s2    # Store stem energies for S1 plus S2

            # Calculate effective loop lengths
            effective_looplength1 = functions.effective_length(list1, l1)
            effective_looplength2 = functions.effective_length(list2, l2)
            effective_looplength3 = functions.effective_length(list3, l3)
                                    
            if list1:      
                energy_l1 = sum([item[3] for item in list1])      # Add free energies          
                      
            if list2:            
                energy_l2 = sum([item[3] for item in list2])      # Add free energies
                                       
            if list3:       
                energy_l3 = sum([item[3] for item in list3])      # Add free energies    
            
            # Calculate free energy of pseudoknot
            entropy = PENALTY * (effective_looplength1 + effective_looplength2 + effective_looplength3)

            pk_energy = stack_s1 + stack_s2 + entropy + energy_l1 + energy_l2 + energy_l3 + INIT
            
            if pk_energy < values[0][0]:
                pk_dic_longpk[pk_stem] = pk_energy, values[1], values[2], values[3]
            else:
                energy = values[0][0]
                pk_dic_longpk[pk_stem] = energy, [], [], []                                           

        else:
            energy = values[0][0]
            pk_dic_longpk[pk_stem] = energy, values[1], values[2], values[3]
              
    return pk_dic_longpk
