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
import user
import stems
import secondary_elements
import pk_construction
import cc06
import cc09
import longPK
import pk_interrupted
import pk_recursive
import mwis
import mwis_2d
import kissing_hairpins
import pk_tools
import BPs_to_CT

# Path to RNAfold program which uses low probability threshold
RNAFOLD_PATH = "./ViennaRNA-1.8.5/Progs/"

def main(identifier, seq, khp, local, global_structure):
    #---------------------------------------------------------------------------------------------------------                
    # Assign ID to the sequence, so that file names will be unique, hopefully, if DotKnot is run in parallel
    FILE_ID = identifier.strip()
    FILE_ID = FILE_ID.replace('>','')
    FILE_ID = FILE_ID.split()
    if len(FILE_ID) > 1:
        FILE_ID = FILE_ID[0]
    else:
        FILE_ID = str(FILE_ID[0])
    FILE_ID = FILE_ID.replace('.','')
    FILE_ID = FILE_ID.replace('_','')
    FILE_ID = FILE_ID[:12]
    #---------------------------------------------------------------------------------------------------------                
    print identifier
    print seq
    print "Sequence length: ", len(seq)
    print "Using " + FILE_ID + " as file identifier"
    print
    print "DotKnot is running..."
    print "Predicting pseudoknots..."
    #---------------------------------------------------------------------------------------------------------                
    # Construct stems from stack probabilities and store in stem dictionary
    basepairs = stems.RNAfold_dotplot(seq, RNAFOLD_PATH, FILE_ID)    
    stem_dic = stems.find_stems(basepairs)
    stem_dic = stems.evaluation_stems(stem_dic, seq, RNAFOLD_PATH, FILE_ID)    
    CUTOFF_STACK = 0.0 
    CUTOFF_LOOP = 4.0
    stem_dic = stems.filter_stems(stem_dic, CUTOFF_STACK, CUTOFF_LOOP)
    #---------------------------------------------------------------------------------------------------------
    # Construct stems with bulge loops, internal loops and multiloops
    CUTOFF_PROB =  0.001
    stems_ib = stems.filter_stems_prob(stem_dic, CUTOFF_PROB)
    structures_dic = secondary_elements.internal_mwis(stems_ib)
    bulges_internal, multiloops = secondary_elements.evaluation_secondary_structures(structures_dic, seq, RNAFOLD_PATH, FILE_ID)
    CUTOFF_IB_ML = 0.0
    bulges_internal = secondary_elements.filter_stems(bulges_internal, CUTOFF_IB_ML)
    multiloops = secondary_elements.filter_stems(multiloops, CUTOFF_IB_ML) 
    #---------------------------------------------------------------------------------------------------------
    # Construct H-type core pseudoknots
    stems_shortened, pk_dic = pk_construction.build_pseudoknots(stem_dic)
    # Construct pseudoknots with one interrupted stem
    stems_shortened_ib, pk_dic_ib = pk_construction.build_pseudoknots_ib(bulges_internal, stem_dic)      
    #---------------------------------------------------------------------------------------------------------
    # Re-evaluate stem energies to account for shortended stems in pseudoknots
    stems_shortened.update(stems_shortened_ib)
    if stems_shortened:
        stems_shortened_dic = stems.evaluation_stems(stems_shortened, seq, RNAFOLD_PATH, FILE_ID)
        stems_shortened_dic = stems.update_confidence(stem_dic, stems_shortened_dic)
        stems_shortened_dic = stems.filter_stems(stems_shortened_dic, CUTOFF_STACK, CUTOFF_LOOP)
    else:
        stems_shortened_dic = {}    
    #---------------------------------------------------------------------------------------------------------
    # Construct three different pseudoknot dictionaries
    pk_dic_cc06, pk_dic_cc09, pk_dic_longpk = pk_construction.pk_dic_scan(pk_dic, stem_dic, stems_shortened_dic)
    # Heuristic energy parameters
    INIT = 7.0 
    PENALTY = 0.1
    # Evaluate energies for pseudoknots with loop L2 <= 1nt
    pk_dic_cc06_result = cc06.dic_caochen06(pk_dic_cc06, stem_dic, stems_shortened_dic, seq)
    # Evaluate energies for pseudoknots with loop L2 >= 2nt and L2 <= 6nt
    pk_dic_cc09_result, entropies_dic, entropies_dic_L1, entropies_dic_L3 = cc09.dic_caochen09(pk_dic_cc09, stem_dic)
    # Evaluate energies for pseudoknots with loop L2 >= 7nt
    pk_dic_longpk_result = longPK.dic_longpks(pk_dic_longpk, stem_dic, INIT, PENALTY)
    # Now we form a pseudoknot dictionary out of the three results
    pk_core_dic = {}
    pk_core_dic.update(pk_dic_cc06_result)
    pk_core_dic.update(pk_dic_cc09_result)
    pk_core_dic.update(pk_dic_longpk_result)
    #---------------------------------------------------------------------------------------------------------
    # Create list of all secondary structure elements for dynamic programming MWIS calculation
    # Format is (start, end, weight, free energy, marker)
    elements = [(stem[0], stem[1], values[3], values[3], "hp") for stem, values in stem_dic.items() if values[3] <= 0.0]
    elements = elements + [(stem[0], stem[1], values[2], values[2], "ib") for stem, values in bulges_internal.items()]
    elements = elements + [(stem[0], stem[1], values[2], values[2], "ml") for stem, values in multiloops.items()]    
    elements.sort()
    MAX_LOOP_UPPER = 400
    array_DP = mwis_2d.DP_MWIS(seq, elements, MAX_LOOP_UPPER)
    array_traceback = mwis_2d.traceback_mwis_2d(seq, array_DP, MAX_LOOP_UPPER)
    MAX_LOOP_UPPER = 100      
    elements_positive = [(stem[0], stem[1], -1.0/round(values[3],2), values[3], "hp") for stem, values in stem_dic.items() if values[3] > 0.0]
    elements_positive.sort()
    array_DP_positive = mwis_2d.DP_MWIS_positive_weights(seq, elements_positive, MAX_LOOP_UPPER)
    array_traceback_positive = mwis_2d.traceback_mwis_2d_positive(seq, array_DP_positive, MAX_LOOP_UPPER)   
    #---------------------------------------------------------------------------------------------------------
    # Now, look for recursive elements in the pseudoknot loops (MWIS) and add them to dictionary
    pk_core_dic = pk_recursive.recursive_pk(pk_core_dic, seq, array_traceback, array_traceback_positive)
    # Construct three different pseudoknot dictionaries for energy re-evaluation because of the loop entropies
    pk_dic_cc06, pk_dic_cc09, pk_dic_longpk = pk_recursive.pk_dic_scan_recursive(pk_core_dic)  
    pk_dic_cc06_result = cc06.pk_energy_reevaluation_06(pk_dic_cc06)
    pk_dic_cc09_result = cc09.pk_energy_reevaluation_09(pk_dic_cc09, entropies_dic, entropies_dic_L1, entropies_dic_L3)
    pk_dic_longpk_result = longPK.pk_energy_reevaluation_long(pk_dic_longpk, INIT, PENALTY) 
    pk_recursive_dic = {}    
    pk_recursive_dic.update(pk_dic_cc06_result)
    pk_recursive_dic.update(pk_dic_cc09_result)
    pk_recursive_dic.update(pk_dic_longpk_result)
    #---------------------------------------------------------------------------------------------------------
    # Evaluate energies for pseudoknots with stems s_ib
    pk_dic_ib = pk_interrupted.evaluate_pk_with_IB(pk_dic_ib, stem_dic, stems_shortened_dic, INIT, PENALTY, seq)
    pk_dic_ib = pk_interrupted.recursive_pk(pk_dic_ib, array_traceback, array_traceback_positive)
    pk_dic_ib = pk_interrupted.re_evaluate_pk_with_IB(pk_dic_ib, INIT, PENALTY)
    pk_dic_ib = pk_recursive.pk_filter(pk_dic_ib)
    # Pseudoknots constructed with stems s_ib are stored in the pseudoknot dictionary
    pk_recursive_dic.update(pk_dic_ib)
    #---------------------------------------------------------------------------------------------------------
    # Store these pseudoknots for best local pseudoknot search later
    pk_not_filtered = pk_recursive_dic.copy()
    #---------------------------------------------------------------------------------------------------------
    # Filter all pseudoknots with energy < 5.25 and normalized energy < -0.25    
    pk_recursive_dic = pk_recursive.pk_filter(pk_recursive_dic) 
    #---------------------------------------------------------------------------------------------------------    
    # If user wants to predict kissing hairpins
    if khp:   
        print "Predicting kissing hairpins..."      
        kissing_hairpin_dic = {}          
        #---------------------------------------------------------------------------------------------------------
        # Kissing hairpin parameters
        CUTTOFF_STACK_KHP = -5.0 
        CUTTOFF_LOOP_KHP = 2.0       
        INIT = 9.0 
        UNPAIRED_NT = 0.5
        UNPAIRED_NT_L3 = 0.0    
        #---------------------------------------------------------------------------------------------------------    
        stem_dic_kissing = stems.filter_stems(stem_dic, CUTTOFF_STACK_KHP, CUTTOFF_LOOP_KHP)
        #---------------------------------------------------------------------------------------------------------
        stems_shortened, all_pseudoknots = kissing_hairpins.build_pseudoknots(stem_dic_kissing)    
        #---------------------------------------------------------------------------------------------------------
        # Re-evaluate stem energies to account for shortended stems
        if stems_shortened:
          stems_shortened_dic = stems.evaluation_stems(stems_shortened, seq, RNAFOLD_PATH, FILE_ID)
          stems_shortened_dic = stems.update_confidence(stem_dic_kissing, stems_shortened_dic)
          stems_shortened_dic = stems.filter_stems(stems_shortened_dic, CUTTOFF_STACK_KHP, CUTTOFF_LOOP_KHP)  
        else:
          stems_shortened_dic = {}    
        #---------------------------------------------------------------------------------------------------------
        # Filter pseudoknot candidates where one of the stems has been shortened and filtered
        all_pseudoknots, pseudoknot_first, pseudoknot_second = kissing_hairpins.filter_khp_dictionary(all_pseudoknots, stems_shortened_dic, stem_dic_kissing)    
        #---------------------------------------------------------------------------------------------------------     
        # Merge shortenend and regular stems into one dictionary to disallow unsuccessful lookups
        all_stems = stems_shortened_dic.copy()
        stems_temp = dict([((stem[0], stem[1], values[0]), values) for stem, values in stem_dic_kissing.items()])
        all_stems.update(stems_temp) 
        #---------------------------------------------------------------------------------------------------------                
        best_khps = kissing_hairpins.kissing_hairpins(seq, pseudoknot_second, pseudoknot_first, all_stems, array_traceback, INIT, UNPAIRED_NT, UNPAIRED_NT_L3)                        
    #---------------------------------------------------------------------------------------------------------
    # The user wants to predict only pseudoknots, no kissing hairpins
    if not khp:
        best_khps = {}   
    #---------------------------------------------------------------------------------------------------------
    # Global MWIS calculation
    stem_dic_mwis = stems.filter_stems_prob(stem_dic, CUTOFF_PROB)    
    mwis_dic, crossing_structures, secondary_structures = mwis.method(stem_dic_mwis, pk_recursive_dic, bulges_internal, multiloops, best_khps) 
    pseudoknot_list = mwis.pk_khp_structures(seq, crossing_structures, best_khps, stem_dic, bulges_internal, multiloops, pk_recursive_dic, pk_dic_ib)     
    #---------------------------------------------------------------------------------------------------------
    # Display predicted pseudoknots and kissing hairpins    
    print
    if pseudoknot_list:
        pseudoknot_list.sort()
        if khp:
            print "Detected pseudoknots and kissing hairpins:"
        else:
            print "Detected pseudoknots:"
        for pk in pseudoknot_list:
            print pk[0], pk[1], pk[2]
            print pk[4]
            print pk[5]
            print
    else:
        if khp:
            print "No pseudoknots or kissing hairpins were detected."
        else:
            print "No pseudoknots were detected."
    #---------------------------------------------------------------------------------------------------------              
    # The user wants to see near-optimal pseudoknots 
    if local:
        # Here, return the best five near-optimal pseudoknots
        NUMBER_OF_BEST_PKS = 5    
        best_pks_ratio = pk_tools.output_best_ratio(pk_not_filtered, best_khps, NUMBER_OF_BEST_PKS)
        best_pks_ratio = pk_tools.dot_bracket(seq, best_pks_ratio, stem_dic, stems_shortened_dic, bulges_internal, multiloops, pk_not_filtered, best_khps)
        best_pks_energy = pk_tools.output_best_energy(pk_not_filtered, best_khps, NUMBER_OF_BEST_PKS)
        best_pks_energy = pk_tools.dot_bracket(seq, best_pks_energy, stem_dic, stems_shortened_dic, bulges_internal, multiloops, pk_not_filtered, best_khps)                
        print
        if khp:
            print "Best", NUMBER_OF_BEST_PKS, "pseudoknots and kissing hairpins in terms of energy to length ratio:"
        else:
            print "Best", NUMBER_OF_BEST_PKS, "pseudoknots in terms of energy to length ratio:"        
        if not best_pks_ratio:
            print "None"
        for item in best_pks_ratio:
            if item[3] == 'khp':
                print item[0], item[1], item[2]
                print item[4]
                print item[5]
            else:                
                print item[0], item[1], item[2]
                print item[3]
                print item[4]
        print
        if khp:
            print "Best", NUMBER_OF_BEST_PKS, "pseudoknots and kissing hairpins in terms of free energy:"
        else:
            print "Best", NUMBER_OF_BEST_PKS, "pseudoknots in terms of free energy:"
        if not best_pks_energy:
            print "None"
        for item in best_pks_energy:
            if item[3] == 'khp':
                print item[0], item[1], item[2]
                print item[4]
                print item[5]
            else:                  
                print item[0], item[1], item[2]
                print item[3]
                print item[4]
        print
    #---------------------------------------------------------------------------------------------------------
    # Assemble global structure if desired by the user
    if global_structure:
        predicted_global_structure = mwis.assemble_global_structure(seq, secondary_structures, pseudoknot_list)
        print "Predicted global structure"
        print seq
        print predicted_global_structure

        # Create CT file
        base_pairs = BPs_to_CT.db_to_ct_brackets(predicted_global_structure)
        output_ct = FILE_ID + '.ct'
        BPs_to_CT.write_CT_file(seq, base_pairs, output_ct)
        
    #---------------------------------------------------------------------------------------------------------
    try:
        os.remove(FILE_ID + '_stem_structure.txt')
        os.remove(FILE_ID + '_stacking_energy.txt')
        os.remove(FILE_ID + '_loops_energy.txt')
        os.remove(FILE_ID + '_input.fasta')
        os.remove(FILE_ID + '_bulge_internal_structures.txt')
        os.remove(FILE_ID + '_bulge_internal_energy.txt')
        os.remove(FILE_ID + '_multiloops_energy.txt')
        os.remove(FILE_ID + '_multiloop_structures.txt')
        os.remove(FILE_ID + '_ss.ps')
        os.remove(FILE_ID + '_dp.ps')
        os.remove(FILE_ID + '_dp2.ps')
    except:
        print
        print "Error deleting the files DotKnot created."
        
    return


#--------------------------------------------------------------------------------------------------------- 
path_exists = os.access(RNAFOLD_PATH, os.F_OK)
if not path_exists:
    print "Path to Vienna RNA package does not exist!"
    print "Check path: ", RNAFOLD_PATH
    sys.exit(0)
#--------------------------------------------------------------------------------------------------------- 
commandline = sys.argv
input = user.process_inputfile(commandline)
khp, local, global_structure = user.process_arguments(commandline)

for identifier, seq in input:    
    main(identifier, seq, khp, local, global_structure)
   
