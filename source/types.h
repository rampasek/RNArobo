/*
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
#include <stdint.h>
#include <emmintrin.h>
#include "matrix.h"

using namespace std;

typedef long long ll;
typedef unsigned long long ull;
typedef pair<int, int> interval;
typedef vector< interval > intervals;
typedef pair<interval, interval> interval_pair;

static void aligned_free(void* p){
    if (!p) return;
    free((void*)((uintptr_t)p-((uint16_t*)p)[-1]));
}

// the structure to represent an Secondary Structure Element
struct SSE{
    int id;
    string name;
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
    string stripped_pattern;    // as 'pattern' but stripped of leading/trailing '*'
    string stripped_complement; // as 'complement' but stripped of leading/trailing '*'
    interval num_wc_padding;    // number of leading/trailing '*' in original 'pattern'
    
    double infContent;          //estimated information content (relative entropy) of the SSE
    bool recordOps;             //if true, measure the elapsed time spend by search for this element
    ull ops_counter;            //time/ops counter for DDEO
    
    Matrix table;                   // n-dimensional sparse matrix for dynamic programming
    set<interval> occurrences;  // whether on given position(s) ends a match
    queue<interval> match_buffer;   // buffer for unprocessed matches (in a domain)
    //map< pair<int,int>, set< pair<int,int> > > h_beginnings_cache; // cache of traceback if helical element
    //map< int, set<int> > ss_beginnings_cache;   // cache of traceback if single strand element

    __m128i *maskv;             //precomputed "pattern" table for BNDM
    uint8_t used[256];
    
    SSE() : ops_counter(0ULL), maskv(NULL) {}
    ~SSE() { aligned_free(maskv); }
    
    ull getOpsCount() { return ops_counter; }
    void incOpsCounter() { ++ops_counter; }
    void incOpsCounter(ull c) { ops_counter+=c; }
    void resetOpsCounter() { ops_counter=0; }
};

// a structure to keep empirical data gathered for a K-elements reorder
struct TupleStats{
    vector<int> tuple; // vector of elements that form this tuple, values are indices to Descriptor::sses
    double heuristicScore;
    
    vector<double> sampledOpsPerWindow;
    vector<unsigned long long> memOps;
};

#endif
