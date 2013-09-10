"""
    DotKnot: software for RNA pseudoknot prediction
    Copyright (C) 2010-2011 Jana Sperschneider	

    This program is free software; you can redistribute it and/or modify  
    it under the terms of the GNU General Public License as published by  
    the Free Software Foundation; either version 2 of the License, or     
    (at your option) any later version.
"""

import functions


def DP_MWIS(seq, elements, MAX_LOOP_UPPER):
    """ Function: DP_MWIS()

        Purpose:  To pre-compute MWIS for all loop intervals via dynamic programming.
                  Given set of structure elements and a limit on the upper loop length, 
                  calculate f(i,j) for each loop interval [i,j] where i < j using the
                  following recursion:
                  f(i,i) = 0
                               {  f(i,j-1)
                  f(i,j) = max {
                               {  f(i,k-1) + w([k:j]) for all intervals which end at j

                  For each loop interval [i,j], f(i,j) holds the MWIS of non-overlapping
                  intervals which start after or at i and end before or at j.
                  
                  Implementation of dynamic programming matrix as array to avoid
                  expensive dictionary lookup operations and reduce memory requirements.              
                  
        Input:    Set of structure elements and upper loop length limit. 
        
        Return:   Dynamic programming matrix.
    """     
    # Create a dictionary which stores for each right endpoint a
    # list of intervals which share this right endpoint.     
    endpoints = {}
    for item in elements:
        right_endpoint = item[1]
        if right_endpoint in endpoints:
            values = endpoints[right_endpoint]
            values.append(item)            
        else:
            values = []
            values.append(item)
        endpoints[right_endpoint] = values
            
    # Init dynamic programming matrix, i.e. two-dimensional array.
    array = [[(0.0, None) for i in xrange(len(seq) + 1)] for j in xrange(len(seq) + 1)]                 
    
    for i in xrange(1, len(seq)):
        for j in xrange(i + 1, len(seq)):

            if j - i + 1 <= MAX_LOOP_UPPER and j - i + 1 >= 9:

                previous_value = array[i][j-1][0]
                temp_min = array[i][j-1][0]
                added_interval = None
                
                # At least one interval ends at endpoint j            
                if j in endpoints:
                    # All intervals which end at j
                    list_of_endpoints = endpoints[j]                                    
                    for interval in list_of_endpoints:
                        # All intervals [k:j] where i <= k
                        if interval[0] >= i:                        
                            local_weight = array[i][interval[0]-1][0] + interval[3]                                                      
                            if local_weight <= temp_min:    # This allows for free energy weight of 0.0 kcal/mol.
                                temp_min = local_weight
                                added_interval = interval

                # Lower score when adding an interval which ends at j            
                if temp_min <= previous_value and added_interval:                                                     
                    array[i][j] = temp_min, added_interval
                    
                # Lower score when no interval which ends at j is added
                else:
                    array[i][j] = previous_value, None                

    """mwis_2d_file = file("mwis_2d_file.txt",'w')
    for index, item in enumerate(array):        
        for index2, item2 in enumerate(array[index]):
            #if item2[0] != 0.0 and item2[1] != [] and index < index2:
            mwis_2d_file.write(str(index))
            mwis_2d_file.write(" ")            
            mwis_2d_file.write(str(index2))
            mwis_2d_file.write(" result ")
            mwis_2d_file.write(str(item2))
            mwis_2d_file.write("\n")"""
    
    return array


def traceback_mwis_2d(seq, array, MAX_LOOP_UPPER):    
    """ Function: traceback_mwis_2d()

        Purpose:  Traceback.              
                  
        Input:    Dynamic programming matrix.
        
        Return:   Result for each loop interval plus effective looplength. 
    """
    # Traceback for given loop interval  
    array_traceback = [[(0.0, None) for i in xrange(len(seq) + 1)] for j in xrange(len(seq) + 1)]                 
    
    for i in xrange(1, len(seq)):
        for j in xrange(i-1, len(seq)):       

            if j - i + 1 <= MAX_LOOP_UPPER:
                
                energy = array[i][j][0]
                loop_start = i
                loop_end = j
                elements = []
                
                #if energy != 0.0:
                energy = array[loop_start][loop_end][0]
                element = array[loop_start][loop_end][1]
                new_loop_end = loop_end

                if element:
                    elements.append(element)    
                    energy = array[loop_start][loop_end][0] - element[3]
                    new_loop_end = element[0]    
                    
                while round(energy, 2) != 0.0:
                    element = array[loop_start][new_loop_end][1]
                    if element:
                        elements.append(element)    
                        energy = array[loop_start][new_loop_end][0] - element[3]
                        new_loop_end = element[0]
                    elif new_loop_end > loop_start:
                        new_loop_end = new_loop_end - 1
                    else:
                        break

                # Calculate effective loop length
                looplength = loop_end - loop_start + 1
                effective_looplength = functions.effective_length(elements, looplength)                 

                #else:
                #    effective_looplength = j - i + 1
                
                array_traceback[i][j] = array[i][j][0], elements, effective_looplength


    """mwis_2d_file = file("mwis_2d_file.txt",'w')
    for index, item in enumerate(array_traceback):        
        for index2, item2 in enumerate(array_traceback[index]):
            #if item2[0] != 0.0 and item2[1] != [] and index < index2:
            mwis_2d_file.write(str(index))
            mwis_2d_file.write(" ")            
            mwis_2d_file.write(str(index2))
            mwis_2d_file.write(" result ")
            mwis_2d_file.write(str(item2))
            mwis_2d_file.write("\n")"""

    return array_traceback


def DP_MWIS_positive_weights(seq, elements, MAX_LOOP_UPPER):
    """ Function: DP_MWIS_positive_weights()

        Purpose:  To pre-compute MWIS for all loop intervals via dynamic programming.
                  In this version, use positive weights.             
                  
        Input:    Set of structure elements and upper loop length limit. 
        
        Return:   Dynamic programming matrix.
    """ 
    endpoints = {}
    for item in elements:
        right_endpoint = item[1]
        if right_endpoint in endpoints:
            values = endpoints[right_endpoint]
            values.append(item)            
        else:
            values = []
            values.append(item)
        endpoints[right_endpoint] = values
            
    # Base case: f(i,i) = {}, 0
    # Initialize array instead of dictionary to save space      
    array = [[(0.0, None)  for i in xrange(len(seq) + 1)] for j in xrange(len(seq) + 1)]                 
    
    for i in xrange(1, len(seq)):
        for j in xrange(i + 1, len(seq)):

            if j - i + 1 <= MAX_LOOP_UPPER and j - i + 1 >= 9:

                previous_value = array[i][j-1][0]
                temp_min = array[i][j-1][0]
                            
                # At least one interval ends at endpoint j            
                if j in endpoints:

                    # All intervals which end at j
                    list_of_endpoints = endpoints[j]                
                    
                    for interval in list_of_endpoints:                     
                        if interval[0] >= i:                        
                            local_weight = array[i][interval[0]-1][0] + interval[2]                  
                                    
                            if local_weight < temp_min:
                                temp_min = local_weight
                                added_interval = interval

                # Lower score when adding an interval which ends at j            
                if temp_min < previous_value:                                                     
                    array[i][j] = temp_min, added_interval
                    
                # Lower score when no interval which ends at j is added
                else:
                    array[i][j] = previous_value, None
    
    return array


def traceback_mwis_2d_positive(seq, array, MAX_LOOP_UPPER):    
    """ Function: traceback_mwis_2d_positive()

        Purpose:  Traceback for positive weights.              
                  
        Input:    Dynamic programming matrix.
        
        Return:   Result for each loop interval plus effective looplength. 
    """
    # Traceback for given loop interval
    rows = len(seq)        
    columns = len(seq)    
    array_traceback = [[(0.0, None, None)  for i in xrange(columns + 1)] for j in xrange(rows + 1)]                 
    
    for i in xrange(1, len(seq)):
        for j in xrange(i-1, len(seq)):       

            if j - i + 1 <= MAX_LOOP_UPPER:
                
                energy = array[i][j][0]
                loop_start = i
                loop_end = j
                elements = []
                
                if energy != 0.0:
                    energy = array[loop_start][loop_end][0]
                    element = array[loop_start][loop_end][1]
                    new_loop_end = loop_end

                    if element:
                        elements.append(element)    
                        energy = array[loop_start][loop_end][0] - element[2]
                        new_loop_end = element[0]    
                        
                    while round(energy, 2) != 0.0:
                        element = array[loop_start][new_loop_end][1]
                        if element:
                            elements.append(element)    
                            energy = array[loop_start][new_loop_end][0] - element[2]
                            new_loop_end = element[0]
                        elif new_loop_end > loop_start:
                            new_loop_end = new_loop_end - 1
                        else:
                            break

                    # Calculate effective loop length
                    looplength = loop_end - loop_start + 1
                    effective_looplength = functions.effective_length(elements, looplength)                 

                else:                    
                    effective_looplength = j - i + 1
                
                array_traceback[i][j] = array[i][j][0], elements, effective_looplength
                
                    
    # Now write real, positive free energy weights to array
    for i in xrange(1, len(seq)):
        for j in xrange(i-1, len(seq)):
            elements = array_traceback[i][j][1]            
            if elements:                
                energy = sum([item[3] for item in elements])                
                array_traceback[i][j] = energy, elements, array_traceback[i][j][2]
                    
    return array_traceback

    

