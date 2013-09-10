"""
    DotKnot: software for RNA pseudoknot prediction
    Copyright (C) 2010-2011 Jana Sperschneider	

    This program is free software; you can redistribute it and/or modify  
    it under the terms of the GNU General Public License as published by  
    the Free Software Foundation; either version 2 of the License, or     
    (at your option) any later version.
"""

             
def effective_length(list_of_elements, looplength):
    """ Function: effective_length()

        Purpose:  For a given pseudoknot loop with recursive secondary structure
                  elements, calculate effective looplength. That is the number
                  of unpaired nucleotides outside the recursive helices plus the
                  number of internal helices.
                  
        Input:    Recursive elements and looplength. 
        
        Return:   Effective looplength.
    """    
    effective_looplength = looplength 
  
    if list_of_elements:                         
    # Subtract number of nucleotides in recursive elements
        for item in list_of_elements:
            effective_looplength = effective_looplength - (item[1] - item[0] + 1)                    
            effective_looplength = effective_looplength + 1     # Plus number of helices                           
                
    return effective_looplength

     
def MWIS(interval_set, sorted_endpointlist):
    """ Function: MWIS()

        Purpose:  Calculation of maximum weight independent set. 
                  Step 1, initialization
                  Step 2, scan sorted endpoints list
                  Step 3, traceback step
                  
        Input:    Set of intervals sorted by right endpoints. 
        
        Return:   Maximum weight independent set.
    """ 
    value, temp_max, Smax1, last_interval = [], 0.0, [], 0
    
    for j in xrange(len(sorted_endpointlist)):    # Step 1
        value.insert(j, 0.0)

    for endpoint in sorted_endpointlist:          # Step 2
        if endpoint[1] == 'l':                    # If left endpoint is scanned
            c = endpoint[3] - 1
            value[c] = temp_max + endpoint[2]
            
        if endpoint[1] == 'r':                    # If right endpoint is scanned            
            c = endpoint[3] - 1                            
            if value[c] > temp_max:
                temp_max = value[c]
                last_interval = c
                
    Smax1.insert(0, interval_set[last_interval])    
    temp_max = temp_max - interval_set[last_interval][4]
    
    for j in xrange(last_interval-1,-1,-1):       # Step 3
        if round(value[j], 2) == round(temp_max, 2):            
            if interval_set[j][1] < interval_set[last_interval][0]:                
                Smax1.append(interval_set[j])
                temp_max = temp_max - interval_set[j][4]
                last_interval = j

    return Smax1


def create_sorted_endpointlist(intervals):
    """ Function: create_sorted_endpointlist()

        Purpose:  Create sorted endpoints list for a set of intervals.
                  If two intervals share the same endpoint where one is
                  a right and one is a left endpoint, scan left endpoint
                  before right endpoint.
                  
        Input:    Set of intervals. 
        
        Return:   Set of intervals sorted by right endpoints.
    """  
    sorted_list= []
  
    for index, interval in enumerate(intervals):
        firstpoint = (interval[0], 'l', float(interval[4]), index + 1)
        sorted_list.append(firstpoint)
        secondpoint = (interval[1], 'r', float(interval[4]), index + 1)
        sorted_list.append(secondpoint)  
  
    sorted_list = sorted(sorted_list)
        
    return sorted_list

 
def steric_constraints(result, key, array_traceback):
    """ Function: steric_constraints()

        Purpose:  While looking for recursive secondary structure
                  elements in the loops, check for sterically
                  infeasible configurations. Add one base at
                  either side of the loop region.
                  
        Input:    Set of secondary structure elements, loop start
                  and end and dynamic programming matrix.
        
        Return:   Set of sterically feasible secondary structure elements.
    """ 
    energy, result, effective = array_traceback[key[0]][key[1]]      

    for item in result:
    
        if (item[0], item[1]) == key:
        
            left_result = array_traceback[key[0]+1][key[1]]      
            right_result = array_traceback[key[0]][key[1]-1]                   
        
            if left_result[0] < right_result[0]:
                energy, result, effective = left_result                
            else:
                energy, result, effective = right_result
            
    return energy, result, effective


