""" 
    DotKnot: software for RNA pseudoknot prediction
    Copyright (C) 2010-2011 Jana Sperschneider	

    This program is free software; you can redistribute it and/or modify  
    it under the terms of the GNU General Public License as published by  
    the Free Software Foundation; either version 2 of the License, or     
    (at your option) any later version.
"""

import sys
import os
import getopt


def usage():
    """ Function: usage()

        Purpose:  Print helpful information for the user.               
    """
    print 
    print "Usage for DotKnot: ", 
    print "python dotknot.py <input_file> [-k][-l][-g]"
    print "[-k] include kissing hairpins"
    print "[-l] show best local pseudoknots"
    print "[-g] show predicted global structure"
    sys.exit(0)    

    return


def process_inputfile(commandline):
    """ Function: process_inputfile()

        Purpose:  Read fasta file(s) from user input.
        
        Input:    Command line.
    
        Return:   Sequence(s) and sequence ID(s).                        
    """  
    list_of_ids_seqs = []              # Store ids and sequences in list
    input_list = []
      
    if len(commandline) > 1:
        if commandline[1] == '-h':
            usage()
        else:
            try:
                fasta_file = commandline[1]
                f = open(fasta_file,'U')

                for line in f:                                  
                    if line:
                        list_of_ids_seqs.append(line)                                                
                
                for index, line in enumerate(list_of_ids_seqs):
                    if line[0] == '>':
                        identifier = line
                        complete_seq = ''

                        # Now check if the following lines are all in sequence format
                        for following_line in list_of_ids_seqs[index + 1:]:

                            if '>' in following_line:
                                break
                            elif following_line == '\n':
                                break
                            else:                                                                
                                following_line = following_line.strip()
                                following_line = following_line.upper()
                                marker = True                                
                                for base in following_line:
                                    if base in 'ACGUTacgut':
                                        pass
                                    else:
                                        marker = False
                                if marker:
                                    complete_seq = complete_seq + following_line
                                else:
                                    print "Invalid characters in sequence. Make sure to use only A, C, G, U, T."
                                    sys.exit(0)                                    

                        if len(complete_seq) == 0:
                            print "Sequence is empty."
                            sys.exit(0)
                            
                        if complete_seq:
                            complete_seq = complete_seq.strip()
                            complete_seq = complete_seq.upper()
                            complete_seq = complete_seq.replace('T','U')                                                    

                        result = identifier.strip(), complete_seq.strip()
                        input_list.append(result)                                                          
                                       
            except IOError:
                print
                print "Specify input file!"
                usage()            
                
    else:
        usage()        
    
    return input_list


def process_arguments(commandline):
    """ Function: process_arguments()

        Purpose:  Read arguments given by user.
        
        Input:    Command line.
    
        Return:   Optional arguments.                       
    """ 
    khp, local, global_structure = False, False, False

    try:
        optlist, arguments = getopt.getopt(commandline[2:], 'hklg')
        for opt in optlist:
            if opt[0] == '-h':
                usage()
            elif opt[0] == '-k':            
                khp = True
                print "Include kissing hairpins."
            elif opt[0] == '-l':
                local = True
                print "Show best local pseudoknots."              
            elif opt[0] == '-g':
                global_structure = True
                print "Show predicted global structure."
            else:
                print "Error while parsing options."
                sys.exit(0)
        print
    except getopt.GetoptError:
        print "Error while parsing options."
        usage()

    return khp, local, global_structure

