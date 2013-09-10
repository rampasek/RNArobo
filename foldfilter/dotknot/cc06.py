"""
    DotKnot: software for RNA pseudoknot prediction
    Copyright (C) 2010-2011 Jana Sperschneider	

    This program is free software; you can redistribute it and/or modify  
    it under the terms of the GNU General Public License as published by  
    the Free Software Foundation; either version 2 of the License, or     
    (at your option) any later version.
"""

import math
import functions
    
# Loop entropy for loop L1 and stem S2 (deep groove)
# First entry is loop size, second entry is stem size    
loop1_dic_cc ={
(1,3): 100.0, (1,4) : 4.4 , (1,5) : 2.3 , (1,6) : 2.3 , (1,7) : 2.3 , (1,8) : 4.4 , (1,9) : 5.5 , (1,10) : 6.9 , (1,11) : 8.7 , (1,12) : 9.8 ,
(2,3) : 6.4 , (2,4) : 4.4 , (2,5) : 4.4 , (2,6) : 4.4 , (2,7) : 4.4 , (2,8) : 4.4 , (2,9) : 5.5 , (2,10) : 6.9 , (2,11) : 8.7 , (2,12) : 9.8 , 
(3,3) : 6.4 , (3,4) : 4.5 , (3,5) : 4.6 , (3,6) : 4.8 , (3,7) : 5.0 , (3,8) : 5.2 , (3,9) : 5.5 , (3,10) : 6.9 , (3,11) : 8.7 , (3,12) : 9.8 ,
(4,3) : 6.4 , (4,4) : 5.4 , (4,5) : 5.7 , (4,6) : 5.8 , (4,7) : 5.9 , (4,8) : 5.7 , (4,9) : 6.4 , (4,10) : 6.9 , (4,11) : 8.7 , (4,12) : 9.8 ,
(5,3) : 6.6 , (5,4) : 5.6 , (5,5) : 6.0 , (5,6) : 6.0 , (5,7) : 6.2 , (5,8) : 6.4 , (5,9) : 6.7 , (5,10) : 7.5 , (5,11) : 8.7 , (5,12) : 9.8 ,
(6,3) : 6.6 , (6,4) : 6.0 , (6,5) : 6.5 , (6,6) : 6.5 , (6,7) : 6.8 , (6,8) : 6.7 , (6,9) : 7.2 , (6,10) : 7.7 , (6,11) : 8.8 , (6,12) : 9.2 ,
(7,3) : 6.8 , (7,4) : 6.3 , (7,5) : 6.9 , (7,6) : 6.8 , (7,7) : 7.0 , (7,8) : 7.1 , (7,9) : 7.5 , (7,10) : 8.1 , (7,11) : 8.9 , (7,12) : 9.5 ,
(8,3) : 6.9 , (8,4) : 6.6 , (8,5) : 7.2 , (8,6) : 7.1 , (8,7) : 7.3 , (8,8) : 7.3 , (8,9) : 7.9 , (8,10) : 8.3 , (8,11) : 9.1 , (8,12) : 9.6 ,
(9,3) : 7.1 , (9,4) : 6.9 , (9,5) : 7.5 , (9,6) : 7.4 , (9,7) : 7.6 , (9,8) : 7.5 , (9,9) : 8.1 , (9,10) : 8.6 , (9,11) : 9.2 , (9,12) : 9.7 ,
(10,3): 7.3 , (10,4): 7.1 , (10,5): 7.8 , (10,6): 7.6 , (10,7): 7.8 , (10,8): 7.7 , (10,9): 8.3 , (10,10): 8.8 , (10,11): 9.3 , (10,12): 9.8 ,
(11,3): 7.5 , (11,4): 7.3 , (11,5): 8.0 , (11,6): 7.8 , (11,7): 8.0 , (11,8): 7.9 , (11,9): 8.5 , (11,10): 8.9 , (11,11): 9.3 , (11,12): 9.8 }

# Loop entropy for loop L2 (L3) and stem S1 (shallow groove)
# First entry is loop size, second entry is stem size
loop3_dic_cc ={
(1,3) : 100.0 , (1,4) : 100.0 , (1,5) : 100.0 , (1,6) : 100.0 , (1,7) : 100.0 , (1,8) : 100.0 , (1,9) : 100.0 , (1,10) : 100.0 , (1,11) : 100.0 , (1,12) : 100.0 ,
(2,3) : 6.5 , (2,4) : 9.2 , (2,5) : 9.8 , (2,6) : 11.9 , (2,7) : 12.4 , (2,8) : 12.1 , (2,9) : 13.7 , (2,10) : 13.7 , (2,11) : 15.9 , (2,12) : 18.7 ,
(3,3) : 6.5 , (3,4) : 9.2 , (3,5) : 9.8 , (3,6) : 11.9 , (3,7) : 12.4 , (3,8) : 12.1 , (3,9) : 13.7 , (3,10) : 13.7 , (3,11) : 15.9 , (3,12) : 18.7 ,
(4,3) : 6.5 , (4,4) : 9.2 , (4,5) : 9.8 , (4,6) : 11.9 , (4,7) : 12.4 , (4,8) : 12.1 , (4,9) : 13.7 , (4,10) : 13.7 , (4,11) : 15.9 , (4,12) : 18.7 ,
(5,3) : 6.6 , (5,4) : 9.2 , (5,5) : 9.8 , (5,6) : 11.9 , (5,7) : 12.4 , (5,8) : 12.1 , (5,9) : 13.7 , (5,10) : 13.7 , (5,11) : 15.9 , (5,12) : 18.7 ,
(6,3) : 6.7 , (6,4) : 8.9 , (6,5) : 9.8 , (6,6) : 11.9 , (6,7) : 12.4 , (6,8) : 12.1 , (6,9) : 13.7 , (6,10) : 13.7 , (6,11) : 15.9 , (6,12) : 18.7 ,
(7,3) : 6.9 , (7,4) : 8.9 , (7,5) : 9.1 , (7,6) : 11.9 , (7,7) : 12.4 , (7,8) : 12.1 , (7,9) : 13.7 , (7,10) : 13.7 , (7,11) : 15.9 , (7,12) : 18.7 ,
(8,3) : 7.1 , (8,4) : 8.9 , (8,5) : 8.9 , (8,6) : 11.0 , (8,7) : 12.4 , (8,8) : 12.1 , (8,9) : 13.7 , (8,10) : 13.7 , (8,11) : 15.9 , (8,12) : 18.7 ,
(9,3) : 7.2 , (9,4) : 9.0 , (9,5) : 8.8 , (9,6) : 10.4 , (9,7) : 11.4 , (9,8) : 11.6 , (9,9) : 13.7 , (9,10) : 12.7 , (9,11) : 15.9 , (9,12) : 18.7 ,
(10,3): 7.4 , (10,4): 9.0 , (10,5): 8.8 , (10,6): 10.1 , (10,7): 11.0 , (10,8): 11.4 , (10,9): 12.6 , (10,10): 12.2 , (10,11): 14.1 , (10,12): 15.8 ,
(11,3): 7.6 , (11,4): 9.1 , (11,5): 8.8 , (11,6) : 9.9 , (11,7): 10.7 , (11,8): 11.2 , (11,9): 12.0 , (11,10): 11.8 , (11,11): 13.0 , (11,12): 14.2 ,
(12,3): 7.7 , (12,4): 9.2 , (12,5): 8.8 , (12,6) : 9.8 , (12,7): 10.5 , (12,8): 11.1 , (12,9): 11.5 , (12,10): 11.5 , (12,11): 12.4 , (12,12): 13.2 }

# l_min dictionary for given stem length S1
# stemlength1 : (l_min, a, b, c)
lmin_S1 ={
2:  (4,0.95,1.84,-0.67),  3: (2,0.32,1.92,-3.9),   4: (3,1.77,1.82,-5.76), 5: (4,3.99,1.55,-5.86),
6:  (4,7.73,1.29,-12.67), 7: (5,8.38,1.16,-11.45), 8: (5,4.52,1.61,-7.58), 9: (6,9.05,1.15,-11.45),
10: (6,4.77,1.68,-6.78), 11: (9,2.74,2.05,1.38),   12: (9,4.69,1.8,-1.11) }

# l_min dictionary for given stem length S2
# stemlength2 : (l_min, a, b, c)
lmin_S2 ={
2:  (4,0.12,1.96,0.52),   3:  (2,0.39,1.92,-3.89),  4: (1,-2.14,2.15,-2.09), 5: (1,-2.22,2.11,-2.25),
6:  (1,-2.4,2.18,-2.33),  7:  (1,-2.61,2.21,-2.32), 8: (2,-1.17,2.03,-1.96), 9: (2,-1.66,2.09,-1.98),
10: (2,-1.43,2.09,-2.93), 11: (5,-0.14,2.06,0.15), 12: (5,0.77,1.84,-0.65) }

# Dictionary for coaxial stacking
stack_dic ={
('AU','AU') : -0.9 , ('AU','CG') : -2.2 , ('AU','GC') : -2.1 , ('AU','GU') : -0.6 , ('AU','UA') : -1.1 , ('AU','UG') : -1.4 ,
('CG','AU') : -2.1 , ('CG','CG') : -3.3 , ('CG','GC') : -2.4 , ('CG','GU') : -1.4 , ('CG','UA') : -2.1 , ('CG','UG') : -2.1 ,
('GC','AU') : -2.4 , ('GC','CG') : -3.4 , ('GC','GC') : -3.3 , ('GC','GU') : -1.5 , ('GC','UA') : -2.2 , ('GC','UG') : -2.5 ,
('GU','AU') : -1.3 , ('GU','CG') : -2.5 , ('GU','GC') : -2.1 , ('GU','GU') : -0.5 , ('GU','UA') : -1.4 , ('GU','UG') :  1.3 ,
('UA','AU') : -1.3 , ('UA','CG') : -2.4 , ('UA','GC') : -2.1 , ('UA','GU') : -1.0 , ('UA','UA') : -0.9 , ('UA','UG') : -1.3 ,
('UG','AU') : -1.0 , ('UG','CG') : -1.5 , ('UG','GC') : -1.4 , ('UG','GU') :  0.3 , ('UG','UA') : -0.6 , ('UG','UG') : -0.5 }
  

def dic_caochen06(pk_dic_cc06, stem_dic, stems_shortened_dic, seq):    
    """ Function: dic_caochen06()

        Purpose:  Calculate pseudoknot free energies under energy model CC06.              
                  
        Input:    Dictionary with pseudoknots where L2 <= 1.
        
        Return:   Dictionary with pseudoknots and associated free energy. 

    """
    pk_dic_cc06_result = {}
    entropy_l1, entropy_l3 = 0.0, 0.0

    for pk_stem in pk_dic_cc06:
        
        i, j, k, l = pk_stem[2], pk_stem[3], pk_stem[5], pk_stem[6]
      
        stem1 = i, j
        stemlength1 = pk_stem[4]
        stem1_short = i, j, stemlength1

        threshold_s1, stack_s1, energy_s1 = get_values(stem1, stem1_short, stems_shortened_dic, stem_dic)
                 
        stem2 = k, l
        stemlength2 = pk_stem[7]
        stem2_short = k, l, stemlength2
      
        threshold_s2, stack_s2, energy_s2 = get_values(stem2, stem2_short, stems_shortened_dic, stem_dic)    

        l1 = k - (i + stemlength1)
        l2 = (j - stemlength1 + 1) - (k + stemlength2)            
        l3 = (l - stemlength2) - j
      
        if stemlength1 > 12:
            stemlength1 = 12
        if stemlength2 > 12:
            stemlength2 = 12
            
        # Coaxial stacking
        if l2 == 0 or l2 == 1:
            coaxial_stacking = cs_calculation(seq, stemlength1, stemlength2, i, j, k, l) 
            # Weighting parameter              
            coaxial_stacking = 0.75 * coaxial_stacking    

        # Calculate loop entropy for L1 dependent on stem S2     
        loop1_stem2 = l1, stemlength2    
                 
        if loop1_stem2 in loop1_dic_cc:
            entropy_l1 = 0.62 * loop1_dic_cc[loop1_stem2] 
        else:
            ln_w_coil = 2.14 * l1 + 0.10
            fitting = lmin_S2[stemlength2]
            l_min = fitting[0]
            ln_w = fitting[1] * math.log(l1 - l_min + 1) + fitting[2] * (l1 - l_min + 1) + fitting[3]
            entropy_l1 = 0.62*(ln_w_coil - ln_w)

        # Calculate loop entropy for L3 dependent on stem S1
        loop3_stem1 = l3, stemlength1

        if loop3_stem1 in loop3_dic_cc:
            entropy_l3 = 0.62 * loop3_dic_cc[loop3_stem1]
        else:
            ln_w_coil = 2.14 * l3 + 0.10
            fitting = lmin_S1[stemlength1]
            l_min = fitting[0]
            ln_w = fitting[1] * math.log(l3 - l_min + 1) + fitting[2] * (l3 - l_min + 1) + fitting[3]
            entropy_l3 = 0.62*(ln_w_coil - ln_w)

        # Calculate free energy for pseudoknot
        pk_energy = stack_s1 + stack_s2 + entropy_l1 + entropy_l3 + 1.3 + coaxial_stacking

        if pk_energy < 0.0:
            pk_dic_cc06_result[pk_stem] = pk_energy, stack_s1, stack_s2, entropy_l1, 0.0, entropy_l3, coaxial_stacking              

    return pk_dic_cc06_result  


def cs_calculation(seq, stemlength1, stemlength2, i, j, k, l):
    """ Function: cs_calculation()

        Purpose:  Calculate coaxial stacking free energies. 
                  
        Input:    RNA sequence and two coaxially stacked stems. 
        
        Return:   Coaxial stacking free energy.         
    """
    # Example:
    
    # L2 == 0:
    # ACGGaUUGUguCCGUAAUcACA
    # (((((.[[[[)))))...]]]]
    # Stack is AU,GC = -2.10
    
    # L2 == 1:
    # ACGGaUUGUgAuCCGUAAUcACA
    # (((((.[[[[:)))))...]]]]
    # Stack is AU,GC = -2.10    
     
    pair1 = seq[((i - 1) + stemlength1 - 1)] + seq[j - stemlength1]
    pair2 = seq[((k - 1) + stemlength2 - 1)] + seq[l - stemlength2]

    stack = pair1, pair2

    if stack in stack_dic:
        coaxial_stacking = stack_dic[stack]
    else:
        coaxial_stacking = 0.0

    return coaxial_stacking


def get_values(stem, stem_short, stems_shortened_dic, stem_dic):
    """ Function: get_values()

        Purpose:  For a given stem, look up values in dictionary.     
                  Watch out for shortened stems.         
                  
        Input:    Stem and stem dictionaries.
        
        Return:   Values for stem in dictionary.

    """    
    if stem_short in stems_shortened_dic:    # Look up whether S1 is a shortened stem
        stack = stems_shortened_dic[stem_short][2]
        energy = stems_shortened_dic[stem_short][3]
    else:     
        stack = stem_dic[stem][2]
        energy = stem_dic[stem][3] 
  
    threshold = stem_dic[stem][3] 
  
    return threshold, stack, energy

 
def pk_energy_reevaluation_06(pk_dic_cc06):
    """ Function: pk_energy_reevaluation_06()

        Purpose:  Calculate pseudoknot free energies under energy model CC06. 
                  Take into account recursive structure elements which can 
                  occur in loops L1 and L3.               
                  
        Input:    Dictionary with recursive pseudoknots where L2 <= 1.
        
        Return:   Dictionary with recursive pseudoknots and associated free energy. 

    """ 
    for pk_stem, values in pk_dic_cc06.items():
    
        i, j, k, l = pk_stem[2], pk_stem[3], pk_stem[5], pk_stem[6]  
        stemlength1 = pk_stem[4]
        stemlength2 = pk_stem[7]
      
        l1 = k - (i + stemlength1)
        l2 = (j - stemlength1 + 1) - (k + stemlength2) 
        l3 = (l - stemlength2) - j
    
        list1 = values[1]
        list3 = values[3]
      
        if list1:   
            energy_l1 = sum([item[3] for item in list1])          # Add free energies
            effective_l1 = functions.effective_length(list1, l1)  # Calculate effective loop length            
        
            if stemlength2 > 12:
                stemlength2 = 12                            
        
            # Calculate loop entropy for L1 dependent on stem S2    
            loop1_stem2 = effective_l1, stemlength2
        
            if loop1_stem2 in loop1_dic_cc:
                entropy_l1 = 0.62 * loop1_dic_cc[loop1_stem2]
            else:
                if effective_l1 > 0:
                    ln_w_coil = 2.14 * effective_l1 + 0.10
                    fitting = lmin_S2[stemlength2]
                    l_min = fitting[0]                    
                    ln_w = fitting[1] * math.log(effective_l1 - l_min + 1) + fitting[2] * (effective_l1 - l_min + 1) + fitting[3]
                    entropy_l1 = 0.62*(ln_w_coil - ln_w)
        else:
            energy_l1 = 0.0
            entropy_l1 = values[0][3]
          
        if list3:   
            energy_l3 = sum([item[3] for item in list3])          # Add free energies
            effective_l3 = functions.effective_length(list3, l3)  # Calculate effective loop length    
                
            if stemlength1 > 12:
                stemlength1 = 12   
        
            # Calculate loop entropy for L3 dependent on stem S1
            loop3_stem1 = effective_l3, stemlength1
        
            if loop3_stem1 in loop3_dic_cc:
                entropy_l3 = 0.62 * loop3_dic_cc[loop3_stem1]
            else:
                if effective_l3 > 1:
                    ln_w_coil = 2.14 * effective_l3 + 0.10
                    fitting = lmin_S1[stemlength1]
                    l_min = fitting[0]
                    ln_w = fitting[1] * math.log(effective_l3 - l_min + 1) + fitting[2] * (effective_l3 - l_min + 1) + fitting[3]
                    entropy_l3 = 0.62*(ln_w_coil - ln_w)    
        else:
            energy_l3 = 0.0
            entropy_l3 = values[0][5]          

        stack_s1 = values[0][1]
        stack_s2 = values[0][2]
      
        energy = stack_s1 + stack_s2 + entropy_l1 + energy_l1 + entropy_l3 + energy_l3 + 1.3 + values[0][6]
      
        if energy < values[0][0]:     # If recursive elements stabilize pseudoknot
            pk_dic_cc06[pk_stem] = energy, values[1], values[2], values[3]
        else:
            energy = values[0][0]
            pk_dic_cc06[pk_stem] = energy, [],[],[]
            
    return pk_dic_cc06
