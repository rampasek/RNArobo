/*
 * $Id: search.h,v 1.7 2012-04-16 19:27:26 laci Exp $
 *
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
#include "regression_based_orderer.h"

using namespace std;

class Simple_Search{
    public:
        bool have_solution;
        Descriptor *desc;
        RegressionBasedOrderer *orderer;           // pointer to an instance of Orderer class - manages the element search ordering
        vector<intervals> solutions;

        Simple_Search(Descriptor &dsc, RegressionBasedOrderer &ord);
        void search(string &seq);
        void find_motif(int ind, string &seq, intervals &grid);
        list<interval_pair> get_domain(intervals &grid, string &seq, SSE &se);
        interval_pair get_motif_element_domain(intervals &grid, string &seq, int index_in_motif);
        list<interval> get_next_match(SSE &se);
        void flood_ss(SSE &se, string &seq, int index_in_seq);
        void get_ss_matches(SSE &se, string &seq, int lower_bound, int upper_bound, int end);
        void flood_h(SSE &se, string &seq, int strand1_begin, int strand2_end);
        void get_h_matches(SSE &se, string &seq, int lower_bound1, int upper_bound1, int end,
                              int begin, int lower_bound2, int upper_bound2);
        void set_grid(intervals &grid, SSE &se, list<interval> &match);
        void reset_grid(intervals &grid, SSE &se);

        string solution_to_str(unsigned int ind, string &seq, string &separator);
        string solution_to_dotbracket(unsigned int ind, string &separator);
};

#endif