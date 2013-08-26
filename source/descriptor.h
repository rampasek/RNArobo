/*
 * Project      : RNA motif searching in genomic sequences
 * Description  : the header file for descriptor.cpp
 *
 * Author       : Ladislav Rampasek <rampasek@gmail.com>
 * Institution  : Comenius University in Bratislava
 *
 */

#ifndef DESCRIPTOR_H
#define DESCRIPTOR_H

#include "types.h"

// the class of an descriptor
class Descriptor{
    public:
        vector<SSE> sses;           // the array where all SSEs are stored, indexed from 1 (0-th position is grabage)
        vector<int> motif;          // the list of all SSEs, i-th number is index to the array "sses",
                                    //   if the number is negative it signs for the second strand
                                    //   of the helix which has index "-number" in "sses"
        map<string, int> transl;    // a transl. dictionary from original names (h1,s5,...) to indices in "sses"
        vector<int> predef_srch_order;// the ordered list of SSEs to be searched (indices to the vector "sses") predefined by a user
        vector<int> pknot_levels;   // UNUSED

        Descriptor():initialized(false){};
        Descriptor(ifstream &fin);
        bool is_initialized(){ return initialized;}
        string error_str(){ return err_str;}
        string search_order_to_str(vector<int> order);
        int get_max_motif_length();
        string get_dotbracket_notation();

    private:
        bool initialized;
        string err_str;

        bool parse_desc_map(string line);
        bool parse_desc_properties(string line);
        bool parse_desc_order(string line);
        bool expand_wildcards(string &s);
        int  check_consistency();
        bool has_no_duplicates(vector<int> vec);
        void compute_inf_contents();
        void compute_auxiliary_stats();
        void compile_pattern(SSE &se, bool isFwdPattern);
        void compute_pknot_levels();
};

#endif
