/*
 * Project      : RNA motif searching in genomic sequences
 * Description  : the header file for search.cpp
 *
 * Author       : Ladislav Rampasek <rampasek@gmail.com>
 * Institution  : Comenius University in Bratislava
 *
 */

#ifndef SEARCH_H
#define SEARCH_H

#include <sstream>
#include <algorithm>
#include <vector>
#include <list>

#include "generalfuncs.h"
#include "descriptor.h"
#include "orderer.h"

using namespace std;

class Simple_Search{
    public:
        bool have_solution;
        Descriptor *desc;
        Orderer *orderer;           // pointer to an instance of Orderer class - manages the element search ordering
        vector<intervals> solutions;

        Simple_Search(Descriptor &dsc, Orderer &ord);
        void search(string &seq);
        void find_motif(int ind, string &seq, intervals &grid);
        list<interval_pair> get_domain(intervals &grid, string &seq, SSE &se);
        interval_pair get_motif_element_domain(intervals &grid, string &seq, int index_in_motif);
        list<interval> get_next_match(SSE &se);
        void get_bndm_ss_matches(SSE &se, string &seq, interval &begin_reg, interval &end_reg);
        void run_fwd_ss_filter(SSE &se, string &seq, interval &begin_reg, interval &end_reg);
        void get_naive_ss_matches(SSE &se, string &seq, interval &begin_reg, interval &end_reg);
        void get_simple_ss_matches(SSE &se, string &seq, interval &begin_reg, interval &end_reg);
        void get_ss_matches(SSE &se, string &seq, interval &begin_reg, interval &end_reg);
        void run_fwddp_ss(SSE &se, string &seq, interval &begin_reg);                       //DEPRECATED
        void trace_bckdp_ss(SSE &se, string &seq, interval &begin_reg, interval &end_reg);  //DEPRECATED
        void get_h_matches(SSE &se, string &seq, list<interval_pair> &domain);
        void run_fwddp_h(SSE &se, string &seq, int strand1_begin, int strand2_end);         //DEPRECATED
        void trace_bckdp_h(SSE &se, string &seq, int lower_bound1, int upper_bound1, int end,
                           int begin, int lower_bound2, int upper_bound2);                  //DEPRECATED
        void set_grid(intervals &grid, SSE &se, list<interval> &match);
        void reset_grid(intervals &grid, SSE &se);

        string solution_to_str(unsigned int ind, string &seq, string &separator);
        string solution_to_dotbracket(unsigned int ind, string &separator);
};

#endif
