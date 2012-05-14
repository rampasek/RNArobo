/*
 * $Id: types.h,v 1.7 2012-04-14 19:10:25 laci Exp $
 *
 * Project      : RNA motif searching in genomic sequences
 * Description  : definition of data types used throughout this project
 *
 * Author       : Ladislav Rampasek <rampasek@gmail.com>
 * Institution  : Comenius University in Bratislava
 *
 */

#ifndef TYPES_H
#define TYPES_H

#include <vector>
#include <queue>
#include <map>
#include <set>
#include "matrix.h"

using namespace std;

typedef long long ll;
typedef unsigned long long ull;
typedef pair<int, int> interval;
typedef vector< interval > intervals;
typedef pair<interval, interval> interval_pair;

// the structure to represent an Secondary Structure Element
struct SSE{
    int id;
    bool is_helix;              // TRUE if the SSE is a helix, FALSE if a single strand
    interval size_range;        // the interval of eligible length of the SSE
    string pattern;             // primary structure restrictions for the SSE
    string complement;          // restrictions for a complement of the "pattern" (for helix)
    int num_insertions;         // the number of allowed insertions in the primary structure
    char allowed_insertion;     // allowed insertions coded in IUPAC notation
    int num_mismatches;         // the number of tolerated mismatches in the primary structure
    int num_mispairings;        // the number of tolerated mispaired positions in the helix
    interval strand_dist;       // the number of nucleotides between strands of the helix
    string transf_matrix;       // transformation matrix for strand complementarity
    
    double infContent;          //estimated information content (relative entropy) of the SSE

    Matrix table;                   // n-dimensional sparse matrix for dynamic programming
    set< pair<int, int> > occurrences;  // whether on given position(s) ends a match
    queue<interval> match_buffer;   // buffer for unprocessed matches (in a domain)
    //map< pair<int,int>, set< pair<int,int> > > h_beginnings_cache; // cache of traceback if helical element
    //map< int, set<int> > ss_beginnings_cache;   // cache of traceback if single strand element
};

// a structure to keep empirical data gathered for a K-elements reorder
struct TupleStats{
    vector<int> tuple; // vector of elements that form this tuple, values are indices to Descriptor::sses
    double heuristicScore;
    vector<double> ICscores;
    vector<int> DFscores;
    
    vector<double> sampledOpsPerWindow;
    vector<unsigned long long> memOps;
    unsigned long long basesScanned;
};

#endif