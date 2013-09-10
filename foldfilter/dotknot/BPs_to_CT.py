"""
    DotKnot: software for RNA pseudoknot prediction
    Copyright (C) 2010-2011 Jana Sperschneider	

    This program is free software; you can redistribute it and/or modify  
    it under the terms of the GNU General Public License as published by  
    the Free Software Foundation; either version 2 of the License, or     
    (at your option) any later version.
"""

import os
import sys
import math  


def db_to_ct_brackets(structure):
    """ Function: db_to_ct_brackets()

        Purpose:  Convert dot-bracket notation to CT file notation.
                  
        Input:    Structure given as dot-bracket notation, this may include
                  pseudoknots with two types of brackets, round brackets or
                  square brackets.
        
        Return:   List of base pairs, which will be written to CT file later.
    """    
    stack_round = []
    stack_square = []

    base_pairs = []
    
    """ Scan the structure and put different types of opening brackets on different stacks,
        i.e. in three different lists. """
    for index, element in enumerate(structure):
        if element == '(':
            stack_round.append(index + 1)
        if element == ')':
            bp = (stack_round.pop(), (index + 1))
            bp_rev = (bp[1], bp[0])
            base_pairs.append(bp)
            base_pairs.append(bp_rev)
            
        if element == '[':
            stack_square.append(index + 1)
        if element == ']':
            bp = (stack_square.pop(), (index + 1))
            bp_rev = (bp[1], bp[0])
            base_pairs.append(bp)
            base_pairs.append(bp_rev)        
            
    return base_pairs


def write_CT_file(sequence, base_pairs, output_ct):
    """ Function: write_CT_file()

        Purpose:  Use list of base pairs derived from dot-bracket notation
                  to create CT file.
                  
        Input:    Desired file name, sequence and structure given as list of base pairs. 
        
        Return:   Create CT file.
    """  
    output = file(output_ct,'w')

    output.write(str(len(sequence)))
    output.write('\n')    
        
    for index, base in enumerate(sequence):

        found = False
        
        for (b1, b2) in base_pairs:
            if b1 == index + 1:
                line = (index + 1, base, index, index + 2, b2, index + 1)
                found = True
                
        if not found:
            line = (index + 1, base, index, index + 2, 0, index + 1)

        for element in line:
            output.write(str(element))
            output.write(' ')
        output.write('\n')
        

    output.close()

    return
        
