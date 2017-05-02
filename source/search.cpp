/*
 * Project      : RNA motif searching in genomic sequences
 * Description  : the search core implementation
 *
 * Author       : Ladislav Rampasek <rampasek@gmail.com>
 * Institution  : Comenius University in Bratislava
 *
 */

#include <iostream>
#include <sstream>
#include <cassert>
#include <algorithm>
#include <vector>
#include <array>
#include <list>
#include <stack>
#include <queue>
#include <set>
#include <stdint.h>
#include <emmintrin.h>

#if defined(_WIN32)
    #include <Windows.h>
#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
    #include <unistd.h>
    #include <sys/resource.h>
    #include <sys/times.h>
    #include <time.h>
#else
    #error "Unable to define getCPUTime() for an unknown OS."
#endif

#include "search.h"
#include "generalfuncs.h"
#include "descriptor.h"

using namespace std;
using namespace GF;

//#define DEBUG
//#define DO_CACHE
#define SKIP

#define MAX_INT 2147483647
#define BEGIN first
#define END second

Simple_Search::Simple_Search(Descriptor &dsc, Orderer &ord){
    desc = &dsc;
    orderer = &ord;
}

void Simple_Search::search(string &seq){
    solutions.clear();
    intervals grid(desc->motif.size(), make_pair(-1,-1));

    //use another order
    orderer->setNewSearchOrder(seq.size());
    
    #ifdef DEBUG
    cout<<"Seq:\n"<<seq<<endl;
    cout<<"Motif:"<<endl;
    for(int i=0;i<desc->motif.size();i++){
        cout<<"  "<<desc->motif[i]<<": ";
        cout<<desc->sses[abs(desc->motif[i])].id<<" "<<desc->sses[abs(desc->motif[i])].is_helix<<' ';
        cout<<desc->sses[abs(desc->motif[i])].pattern<<" "<<desc->sses[abs(desc->motif[i])].complement<<endl;
    }

    cout<<"Order:"<<endl;
    for(int i=0;i<orderer->searchOrder.size();i++){
        cout<<"  "<<orderer->searchOrder[i]<<": ";
        sse tmp = desc->sses[abs(orderer->searchOrder[i])];
        
        cout<<tmp.id<<" "<<tmp.is_helix<<' ';
        cout<<tmp.pattern<<" "<<tmp.complement<<endl;
        cout<<tmp.strand_dist.first<<" "<<tmp.strand_dist.second<<endl;

    }

    cout<<orderer->searchOrder.back()<<" <- last element"<<endl;
    #endif

    // clear & resize all tables
    for(int i=1;i<desc->sses.size();i++){
        desc->sses[i].table.clear();
        desc->sses[i].occurrences.clear();
        assert(desc->sses[i].match_buffer.empty());

        if(desc->sses[i].is_helix) {
            desc->sses[i].table.set_dimensions(7);
            //desc->sses[i].occurrences.set_dimensions(2);
            #ifdef DO_CACHE
                desc->sses[i].h_beginnings_cache.clear();
            #endif
        } else {
            desc->sses[i].table.set_dimensions(5);
            //desc->sses[i].occurrences.set_dimensions(1);
            #ifdef DO_CACHE
                desc->sses[i].ss_beginnings_cache.clear();
            #endif
        }
    }
    
    // start search
    have_solution=false;
    find_motif(0, seq, grid);
}

/*
 * Returns the amount of CPU time used by the current process
 * in the most precise units possible, or 0 if an error occurred.
 *
 * Original author:  David Robert Nadeau
 * Site:    http://NadeauSoftware.com/
 * License: Creative Commons Attribution 3.0 Unported License
 */
unsigned long long getCPUTime(void){
    #if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    FILETIME createTime;
    FILETIME exitTime;
    FILETIME kernelTime;
    FILETIME userTime;
    if ( GetProcessTimes( GetCurrentProcess(),
        &createTime, &exitTime, &kernelTime, &userTime ) != -1 )
    {
        return userTime.dwLowDateTime;
    }
    
    #elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
    /* AIX, BSD, Cygwin, HP-UX, Linux, OSX, and Solaris --------- */
    
    #if defined(_POSIX_TIMERS) && (_POSIX_TIMERS > 0)
    /* Prefer high-res POSIX timers, when available. */
    {
        clockid_t id;
        struct timespec ts;
        #if _POSIX_CPUTIME > 0
        /* Clock ids vary by OS.  Query the id, if possible. */
        if ( clock_getcpuclockid( 0, &id ) == -1 )
            #endif
        #if defined(CLOCK_PROCESS_CPUTIME_ID)
        /* Use known clock id for AIX, Linux, or Solaris. */
        id = CLOCK_PROCESS_CPUTIME_ID;
        #elif defined(CLOCK_VIRTUAL)
        /* Use known clock id for BSD or HP-UX. */
        id = CLOCK_VIRTUAL;
        #else
        id = (clockid_t)-1;
        #endif
        if ( id != (clockid_t)-1 && clock_gettime( id, &ts ) != -1 )
            return ts.tv_sec*1000000000ULL + ts.tv_nsec;
    }
    #endif
    
    #if defined(RUSAGE_SELF)
    {
        struct rusage rusage;
        if ( getrusage( RUSAGE_SELF, &rusage ) != -1 )
            return rusage.ru_utime.tv_sec*1000000ULL + rusage.ru_utime.tv_usec;
    }
    #endif
    
    #if defined(CLOCKS_PER_SEC)
    {
        clock_t cl = clock();
        if ( cl != (clock_t)-1 )
            return cl / (CLOCKS_PER_SEC/1000);
    }
    #endif
    
    #if defined(_SC_CLK_TCK)
    {
        struct tms tms;
        if ( times( &tms ) != (clock_t)-1 )
            return tms.tms_utime;
    }
    #endif
    
    #endif
    
    return 0;      /* Failed. */
}

void Simple_Search::find_motif(int ind, string &seq, intervals &grid){
    SSE& se = desc->sses[ orderer->searchOrder[ind] ]; //sse that is going to be searched for
    list<interval_pair> domain = get_domain(grid, seq, se);

    #ifdef DEBUG
    cout<<ind<<' '<<se.id<<" domain: \n   ";
    cout<<"   1.BEGIN:"<<domain.front().BEGIN.first<<" to "<<domain.front().BEGIN.second<<"\n";
    cout<<"   1.END:"<<domain.front().END.first<<" to "<<domain.front().END.second<<"\n";
    if(domain.size()==2){
        cout<<"   2.BEGIN:"<<domain.back().BEGIN.first<<" to "<<domain.back().BEGIN.second<<"\n";
        cout<<"   2.END:"<<domain.back().END.first<<" to "<<domain.back().END.second<<"\n";
    }
    #endif


    /// find the element occurrences in the domain
    // measure CPU time the search takes (instead of memOps)
    unsigned long long se_start_t = 0;
    if(se.recordOps) se_start_t = getCPUTime(); 
    /// DP for helices
    if(se.is_helix){
        get_h_matches(se, seq, domain);
        
        /* //for run with fwd/bckwd DP approach
        int min_dist = se.strand_dist.first + 2*se.size_range.first -1;
        int max_dist = se.strand_dist.second + 2*se.size_range.second;
        
        for(int s1_b=domain.front().BEGIN.first; s1_b<domain.front().BEGIN.second; s1_b++){
            for(int s2_e=max(domain.back().END.first, s1_b+min_dist); s2_e<min(domain.back().END.second, s1_b+max_dist); s2_e++){
                run_fwddp_h(se, seq, s1_b, s2_e);
                //cout<<"flooding "<<s1_b<<' '<<s2_e<<endl;
            }
        }

        for(int s1_e=domain.front().END.first; s1_e<domain.front().END.second; s1_e++){
            for( int s2_b=max(s1_e+se.strand_dist.first+1, domain.back().BEGIN.first);
                s2_b<min(domain.back().BEGIN.second, s1_e+se.strand_dist.second+2); s2_b++ )
            {
                //printf("tracing <%d, %d) %d; %d <%d, %d)  -> %d\n", domain.front().BEGIN.first, domain.front().BEGIN.second, s1_e, s2_b, domain.back().END.first, domain.back().END.second,se.occurrences.get(2,s1_e+1,s2_b+1));
                trace_bckdp_h(se, seq, domain.front().BEGIN.first, domain.front().BEGIN.second, s1_e,
                              s2_b, domain.back().END.first, domain.back().END.second);
            }
        }*/
        
    /// battery of algorithms for single strand elements
    } else {
        //true if the element pattern has fix-sized core, i.e. wild cards are only as prefix/suffix
        bool fixed_core = (se.size_range.first==(se.size_range.second-se.num_wc_padding.first-se.num_wc_padding.second));

        // BNDM for single strand element with fix-sized core and with NO other wild cards, nor insertions
        if(fixed_core && se.num_insertions==0 && se.stripped_pattern.size()<=128){
            get_bndm_ss_matches(se, seq, domain.front().BEGIN, domain.front().END);
            
        // bit-parallel forward scan prefiltering - then calls DP algorithms to find matches
        } else if(se.stripped_pattern.size()<=128){
            run_fwd_ss_filter(se, seq, domain.front().BEGIN, domain.front().END);
            
        // naive alg. for single strand element with fix-sized core and with NO other wild cards nor insertions
        } else if(fixed_core && se.num_insertions==0){
            get_naive_ss_matches(se, seq, domain.front().BEGIN, domain.front().END);
            
        // DP for single strand element with NO insertions
        } else if(se.num_insertions==0){
            get_simple_ss_matches(se, seq, domain.front().BEGIN, domain.front().END);

        // general single strand element DP
        } else {
            get_ss_matches(se, seq, domain.front().BEGIN, domain.front().END);
            
            /* //for run with fwd/bckwd DP approach
            run_fwddp_ss(se, seq, domain.front().BEGIN);
            //printf("tracing <%d, %d) - <%d, %d)\n", domain.front().BEGIN.first, domain.front().BEGIN.second, domain.front().END.first, domain.front().END.second);
            trace_bckdp_ss(se, seq, domain.front().BEGIN, domain.front().END);
            */
        }
    }
    //cout<<getCPUTime() - se_start_t<<" ";
    if(se.recordOps) se.incOpsCounter(getCPUTime() - se_start_t);
    
    list<interval> match = get_next_match(se);

    while(match.front().first != -1 && !have_solution){
        set_grid(grid, se, match);
        if( se.id == orderer->searchOrder.back()){
            #ifdef DEBUG
                cout<<"SOLUTION!!! ";
                for(int i=0;i<desc->motif.size();i++){
                    cout<<desc->motif[i]<<"("<<grid[i].first<<"-"<<grid[i].second<<") ";
                }cout<<endl;
            #endif
            
            solutions.push_back(grid);
            
            #ifdef SKIP
                have_solution=true;
            #endif
        } else {
            #ifdef DEBUG
                cout<<"DOWN"<<endl;
                SSE& t_se = desc->sses[ orderer->searchOrder[ind+1] ];
                cout<<"looking for "<<t_se.id<<": ";
                for(int i=0;i<desc->motif.size();i++){
                    cout<<desc->motif[i]<<"("<<grid[i].first<<"-"<<grid[i].second<<") ";
                }cout<<endl;
            #endif
            
            find_motif(ind+1, seq, grid);
        }

        #ifdef SKIP
            if( se.id == orderer->searchOrder.front() ) have_solution=false;
        #endif
        if( !have_solution ) match = get_next_match(se);
    }
    #ifdef DEBUG
        cout<<"UP "<<match.front().first<<endl;
    #endif
    #ifdef SKIP
        if(!se.match_buffer.empty()) se.match_buffer=queue<interval>();
    #endif
    reset_grid(grid, se);
}

interval_pair Simple_Search::get_motif_element_domain(intervals &grid, string &seq, int index_in_motif){
    //find left boundaries
    int left_search_bound = -1;
    int left_ncover_bound = -1;
    int min_left_shift = 0;
    int max_left_shift = 0;
    bool fixed_left = false;
    for(int i=index_in_motif-1;i>=0;i--){
        if( grid[i].second != -1){
            left_search_bound = grid[i].second + min_left_shift;
            left_ncover_bound = grid[i].second + max_left_shift;
            fixed_left = true;
            break;
        } else {
            min_left_shift += desc->sses[ abs(desc->motif[i]) ].size_range.first;  //minimum of size of the sse
            max_left_shift += desc->sses[ abs(desc->motif[i]) ].size_range.second; //maximum of size of the sse
        }
    }
    if(left_search_bound == -1) left_search_bound = min_left_shift;

    //find right boundaries
    int right_search_bound = -1;
    int right_ncover_bound = -1;
    int min_right_shift = 0;
    int max_right_shift = 0;
    bool fixed_right = false;
    for(int i=index_in_motif+1;i<desc->motif.size();i++){
        if( grid[i].first != -1){
            right_search_bound = grid[i].first - min_right_shift;
            right_ncover_bound = grid[i].first - max_right_shift;
            fixed_right = true;
            break;
        } else {
            min_right_shift += desc->sses[ abs(desc->motif[i]) ].size_range.first; //minimum of size of the sse
            max_right_shift += desc->sses[ abs(desc->motif[i]) ].size_range.second; //maximum of size of the sse
        }
    }
    if(right_search_bound == -1) right_search_bound = seq.size() - min_right_shift;
    
    //if(right_ncover_bound!=-1 && left_ncover_bound > right_ncover_bound) left_ncover_bound = right_ncover_bound = -1;
    
    //search interval is <left_search_boundary,right_search_boundary)
    //interval to be necessarily covered is <left_ncover_boundary,right_ncover_boundary)
    //if -1 is a boundary, it means there is no boundary
    // but if both boundaries are -1, it means there IS NOT an interval to be necessarily covered
    
    
    ///combine search interval, necessary-to-cover interval and minimal size of the element to obtain the *search domain*
    int lower_begin_bound = left_search_bound;
    //if an occurrence is already fixed at the right side, then this element must begin so that the gap between these two
    // elements could be (at least theoretically) filled in by the elmenets inbetween them
    if(fixed_right) lower_begin_bound = max(left_search_bound,
        right_ncover_bound - desc->sses[ abs(desc->motif[index_in_motif]) ].size_range.second);
    
    //int upper_begin_bound = max(left_search_bound+1, right_search_bound - desc->sses[ abs(desc->motif[index_in_motif]) ].size_range.first);
    int upper_begin_bound = max(left_search_bound+1, right_search_bound+1);
    //if necessary cover interval beginning is defined then a match must start at/before it
    if(fixed_left) upper_begin_bound = min(left_ncover_bound+1, upper_begin_bound);
    
    
    int lower_end_bound = min(right_search_bound-1, left_search_bound-1);
    //if necessary cover interval end is defined then a match must end at/after it
    if(fixed_right) lower_end_bound = max(right_ncover_bound-1, lower_end_bound);
    
    int upper_end_bound = right_search_bound;
    //if an occurrence is already fixed at the left side, then this element must end so that the gap between these two
    // elements could be (at least theoretically) filled in by the elmenets inbetween them
    if(fixed_left) upper_end_bound = min(right_search_bound,
        left_ncover_bound + desc->sses[ abs(desc->motif[index_in_motif]) ].size_range.second);
    
    return make_pair( make_pair(lower_begin_bound,upper_begin_bound), make_pair(lower_end_bound,upper_end_bound) );
}

list<interval_pair> Simple_Search::get_domain(intervals &grid, string &seq, SSE &se){
    list<interval_pair> domain;
    int index_in_motif=-1, index2_in_motif=-1;

    for(int i=0;i<desc->motif.size();i++) {
        if( se.id == desc->motif[i]) index_in_motif=i;
        if( se.id == -desc->motif[i]) index2_in_motif=i; //the second strand of the helix (if sse is a helix)
    }

    if(se.is_helix){
        domain.push_back( get_motif_element_domain(grid, seq, index_in_motif) );
        domain.push_back( get_motif_element_domain(grid, seq, index2_in_motif) );
    } else { //is a single strand sse
        domain.push_back( get_motif_element_domain(grid, seq, index_in_motif) );
    }

    return domain;
}

//get value of p-th bit
inline bool getbit(void *v, int p) {
    return ( ((uint32_t*)v)[p >> 5] & (1 << (p & 31)) ) != 0;
}

//bitwise shift left on __m128i
inline void bsl_m128(__m128i *v){
    *v = _mm_or_si128(_mm_slli_epi64(*v, 1), _mm_srli_epi64(_mm_slli_si128(*v, 8), 63));
}

/* Run BNDM pattern search for @se.pattern in @seq. All occurrences must begin
 * at index within <@begin_reg.first, @begin_reg.second) and end within @end_reg.
 * !!!Works for single strand elements with *NO inner wild cards and NO insertions*!!!
 */
void Simple_Search::get_bndm_ss_matches(SSE &se, string &seq, interval &begin_reg, interval &end_reg){
    int patt_length=se.stripped_pattern.size();
    int seq_length=seq.size();
    interval match;
    set<interval> found_matches;
    
    assert(patt_length <= 128);
    
    ///SEARCH
    if(patt_length==0){ //special case for an all-wild-card pattern
        for(int i=begin_reg.first; i<begin_reg.second+se.num_wc_padding.first; ++i){
            //add all possible matches of leading/trailing wild cards
            for(int prefix_l=0; prefix_l<=se.num_wc_padding.first; ++prefix_l){
                for(int suffix_l=0; suffix_l<=se.num_wc_padding.second; ++suffix_l){
                    match.first = i - prefix_l;
                    match.second = i + suffix_l;
                    //check if the match is inside the search domain
                    if(end_reg.first < match.second && match.second <= end_reg.second
                        && match.first >= begin_reg.first && match.first < begin_reg.second)
                    {
                        if(found_matches.count(match)==0){
                            found_matches.insert( match );
                            se.match_buffer.push( match );
                        }
                        #ifdef DEBUG
                            cout<<se.id<<" has match(bndm) "<<match.first+1<<" to "<<match.second<<" | ";
                            cout<<se.stripped_pattern<<" "<<seq.substr(match.first, match.second-match.first)<<endl;
                        #endif
                        
                        //se.table.incOpsCounter();
                    }
                }
            }
        }
    } else {
        int j, last;
        __m128i newR, oldR;
        __m128i R[se.num_mismatches + 1];
        __m128i tmp, zero = {}, ones = se.maskv[4];
        
        for(int i = begin_reg.first; i < begin_reg.second+se.num_wc_padding.first; i += last){
            if(patt_length - 1 + i >= seq_length) break;

            j = patt_length - 1;
            last = max(1, patt_length - 1);
            
            R[0] = se.maskv[(int)seq[i + j]];
            for(int x = 1; x <= se.num_mismatches; ++x) R[x] = ones;
            newR = ones;
            while(0xFFFF != _mm_movemask_epi8(_mm_cmpeq_epi8(zero, newR))){
                oldR = R[0];
                newR = R[0];
                if (j!=0) {
                    bsl_m128(&newR);
                    newR = _mm_and_si128(newR, se.maskv[(int)seq[i + j - 1]]);
                    R[0] = newR;
                }
                //run "dynamic programming" for number of mismatches
                for(int x = 1; x <= se.num_mismatches; ++x){
                    //align text to pattern at this "level"
                    tmp = R[x];
                    if (j!=0) {
                        bsl_m128(&tmp);
                        tmp = _mm_and_si128(tmp, se.maskv[(int)seq[i + j - 1]]);
                     }
                    //shift left previous "level" - prepare for mismatch
                    bsl_m128(&oldR);
                    
                    //by taking OR allow for mismatch and alignment at the same time
                    newR = _mm_or_si128(tmp, oldR);
                    
                    oldR = R[x];
                    R[x] = newR;
                    
                    //se.table.incOpsCounter();
                }
                
                --j;
                
                //if we have a suffix in the text that is a prefix of the pattern
                //if( ((uint32_t*)&newR)[mask_offset] & patlen_mask ){ //getbit(&newR, patt_length-1)
                if( 0 > (int16_t)_mm_movemask_epi8(newR) ){
                    if(j>0){
                        last = j;
                    } else { //we have a complete match
                        int start=i;
                        int end=i+patt_length;
                        //add also all possible matches of leading/trailing wild cards
                        for(int prefix_l=0; prefix_l<=se.num_wc_padding.first; ++prefix_l){
                            for(int suffix_l=0; suffix_l<=se.num_wc_padding.second; ++suffix_l){
                                match.first = start - prefix_l;
                                match.second = end + suffix_l;
                                //check if the match is inside the search domain
                                if(end_reg.first < match.second && match.second <= end_reg.second
                                    && match.first >= begin_reg.first && match.first < begin_reg.second)
                                {
                                    if(found_matches.count(match)==0){
                                        found_matches.insert( match );
                                        se.match_buffer.push( match );
                                    }
                                    #ifdef DEBUG
                                        cout<<se.id<<" has match(bndm) "<<match.first+1<<" to "<<match.second<<" | ";
                                        cout<<se.stripped_pattern<<" "<<seq.substr(match.first, match.second-match.first)<<endl;
                                    #endif

                                    //se.table.incOpsCounter();
                                }
                            }
                        }
                        break;
                    }
                }

                //se.table.incOpsCounter();
            }
        }
    }
}

inline void subtract_m128(__m128i *a, __m128i *b, __m128i *res){
    *res = _mm_sub_epi64(_mm_sub_epi64(*a, *b), _mm_sub_epi64(_mm_loadl_epi64(a), _mm_loadl_epi64(b)));
}

/* Run bit-parallel forward scan filtering for @se.stripped_pattern in @seq. Uses DP algorithms to
 * find the actual pattern occurrences in filtered regions. Handles all types of allowed errors
 * (i.e. mismatches, wild-cards, insertions).
 */
void Simple_Search::run_fwd_ss_filter(SSE &se, string &seq, interval &begin_reg, interval &end_reg){
    int patt_length=se.stripped_pattern.size();
    interval match;
    set<interval> found_matches;
    
    assert(patt_length <= 128);
    
    if(patt_length==0){
        ///special treatment for an all-wild-card pattern
        for(int i=begin_reg.first; i<begin_reg.second+se.num_wc_padding.first; ++i){
            //add all possible matches of leading/trailing wild cards
            for(int prefix_l=0; prefix_l<=se.num_wc_padding.first; ++prefix_l){
                for(int suffix_l=0; suffix_l<=se.num_wc_padding.second; ++suffix_l){
                    match.first = i - prefix_l;
                    match.second = i + suffix_l;
                    //check if the match is inside the search domain
                    if(end_reg.first < match.second && match.second <= end_reg.second
                        && match.first >= begin_reg.first && match.first < begin_reg.second)
                    {
                        if(found_matches.count(match)==0){
                            found_matches.insert( match );
                            se.match_buffer.push( match );
                        }
                        #ifdef DEBUG
                            cout<<se.id<<" has match(fwd scan) "<<match.first+1<<" to "<<match.second<<" | ";
                            cout<<se.stripped_pattern<<" "<<seq.substr(match.first, match.second-match.first)<<endl;
                        #endif
                        
                        //se.table.incOpsCounter();
                    }
                }
            }
        }
    } else {
        ///run filtering
        queue<int> match_ends;
        int k = se.num_mismatches + se.num_insertions;
        __m128i newR, oldR;
        __m128i R[k + 1];
        __m128i tmp, tmp2, zero = {};
        __m128i I = se.maskv[0], F = se.maskv[1], nF = se.maskv[2], one = se.maskv[3];
        
        R[0] = zero;
        for(int x = 1; x <= k; ++x){
            R[x] = zero;
            /*
            //init to allow insertions and wild cards in front of the pattern - *we don't need that*
            tmp = R[x-1];
            bsl_m128(&tmp);
            R[x] = _mm_or_si128(_mm_or_si128(R[x-1], tmp), one);
            R[x] = _mm_or_si128(R[x-1], one);
            
            tmp = _mm_and_si128(R[x], I);
            subtract_m128(&F, &tmp, &R[x]);
            R[x] = _mm_and_si128(R[x], nF);
            */
        }
        
        //start bit-parallel forward searching
        for(int i = begin_reg.first; i < end_reg.second; ++i){
            oldR = R[0];
            tmp = R[0];
            bsl_m128(&tmp);
            R[0] = _mm_and_si128(_mm_or_si128(tmp, one), se.maskv[(int)seq[i]]);
            
            tmp = _mm_and_si128(R[0], I);
            subtract_m128(&F, &tmp, &tmp2);
            R[0] = _mm_or_si128(R[0], _mm_and_si128(tmp2, nF));
            
            //run "dynamic programming" for number of mismatches+insertions
            for(int x = 1; x <= k; ++x){
                //align text to pattern at this "level"
                tmp = R[x];
                bsl_m128(&tmp);
                tmp = _mm_and_si128(tmp, se.maskv[(int)seq[i]]);
                
                //ORing in oldR = allowing insert; ORing in One = new beginning
                tmp2 = _mm_or_si128(oldR, one);
                //shift left previous "level" = allow mismatch
                bsl_m128(&oldR);

                //allow for alignment, mismatch, insertion at the same time + new beginning
                newR = _mm_or_si128(_mm_or_si128(tmp, oldR), tmp2);
                
                //add epsilon jumps (wild cards)
                tmp = _mm_and_si128(newR, I);
                subtract_m128(&F, &tmp, &tmp2);
                newR = _mm_or_si128(newR, _mm_and_si128(tmp2, nF));

                oldR = R[x];
                R[x] = newR;

                //se.table.incOpsCounter();
            }

            //if we have a match ending at seq[i]
            //if( ((uint32_t*) &R[k])[mask_offset] & patlen_mask){
            if( 0 > (int16_t)_mm_movemask_epi8(R[k]) ){
                if(end_reg.first <= i+1 + se.num_wc_padding.second + se.num_insertions){
                    match_ends.push(i+1);
                }
                #ifdef DEBUG
                    //uint32_t *Z = ((uint32_t*) &R[k]);
                    //cout<<i<<" "<<(Z[0]&8)<<(Z[0]&4)<<(Z[0]&2)<<(Z[0]&1)<<endl;
                    cout<<se.id<<" found ending(fwd scan) at "<<i+1<<endl;
                #endif
            }
        }

        ///check for complete matches in prefiltered positions
        int end, start;
        interval begin_window = make_pair(-1, -1);
        bool keep_running = true;
        while(keep_running){
            if(!match_ends.empty()){
                end = match_ends.front() + 1;
                match_ends.pop();
                start = end - se.size_range.second - 1;
            } else {
                start = end = MAX_INT;
                keep_running = false;
            }
            
            //expand current window if there is an overlap

            if(begin_window.first <= start && start <= begin_window.second){
                begin_window.second = end - se.size_range.first + 1;
            
            //otherwise run search in the current window and create a new one
            } else {
                //run search in the current window
                if(begin_window.first != -1){
                    begin_window.first = max(begin_window.first, begin_reg.first);
                    begin_window.second = min(begin_window.second, begin_reg.second);
                    
                    // DP for single strand element with NO insertions
                    if(se.num_insertions==0){
                        get_simple_ss_matches(se, seq, begin_window, end_reg);
                     
                    // general single strand element DP
                    } else {
                        get_ss_matches(se, seq, begin_window, end_reg);
                    }
                }
                //create a new window
                begin_window.first = start;
                begin_window.second = end - se.size_range.first + 1;
            }
        }
        
    }
}

/* Run naive pattern search for @se.pattern in @seq. All occurrences must begin
 * at index within <@begin_reg.first, @begin_reg.second) and end within @end_reg.
 * !!!Works for single strand elements with *NO inner wild cards or insertions*!!!
 */
void Simple_Search::get_naive_ss_matches(SSE &se, string &seq, interval &begin_reg, interval &end_reg){
    int patt_length=se.stripped_pattern.size();
    int seq_length=seq.size();
    interval match;
    set<interval> found_matches;
    
    for(int i=begin_reg.first; i<begin_reg.second+se.num_wc_padding.first; i++){
        //align seq and stripped_pattern
        int j=0;
        int mm=0;
        if(se.num_mismatches==0){
            while(j<patt_length && i+j<seq_length && fits(seq[i+j], se.stripped_pattern[j])){
                j++;
                //se.table.incOpsCounter();
            }
        } else { //allow mismatches
            while(j<patt_length && i+j<seq_length && mm<=se.num_mismatches){
                if(fits(seq[i+j], se.stripped_pattern[j])){
                    j++;
                } else {
                    j++;
                    mm++;
                }
                //se.table.incOpsCounter();
            }
        }
        //if it is a complete match, put it to the list of occurrences
        if(j==patt_length && mm<=se.num_mismatches){
            int start=i;
            int end=i+j;
            //add also all possible matches of leading/trailing wild cards
            for(int prefix_l=0; prefix_l<=se.num_wc_padding.first; ++prefix_l){
                for(int suffix_l=0; suffix_l<=se.num_wc_padding.second; ++suffix_l){
                    match.first = start - prefix_l;
                    match.second = end + suffix_l;
                    /*if(se.id==4) {
                        cout<<"SDOMAIN "<<begin_reg.first<<"-"<<begin_reg.second<<"   "<<end_reg.first<<"-"<<end_reg.second<<endl;
                        cout<<"MATCH   "<<match.first<<"-"<<match.second<<"  "<<(match.first >= begin_reg.first)<<(end_reg.first < match.second)<<(match.second <= end_reg.second)<<endl;
                    }*/
                    //check if the match is inside the search domain
                    if(end_reg.first < match.second && match.second <= end_reg.second
                        && match.first >= begin_reg.first && match.first < begin_reg.second)
                    {
                        //se.occurrences.insert( make_pair(end, 0) ); //at position 'i+j-1' in seq ends a match
                        if(found_matches.count(match)==0){
                            found_matches.insert( match );
                            se.match_buffer.push( match );
                        }
                        #ifdef DEBUG
                            cout<<se.id<<" has match(n) "<<match.first+1<<" to "<<match.second<<" | ";
                            cout<<se.stripped_pattern<<" "<<seq.substr(match.first, match.second-match.first)<<endl;
                        #endif
                            
                        //se.table.incOpsCounter();
                    }
                }
            }
            
        }
    }
}

/* Run simple DP pattern search for @se.pattern in @seq. All occurrences must begin
 * at index within <@begin_reg.first, @begin_reg.second) and end within @end_reg.
 * !!!Works for single strand elements with *NO insertions*!!!
 */
void Simple_Search::get_simple_ss_matches(SSE &se, string &seq, interval &begin_reg, interval &end_reg){
    int patt_length=se.stripped_pattern.size();
    int seq_length=seq.size();
    interval match;
    set<interval> found_matches;
    
    set < pair<interval,int> > visited;
    queue < pair<interval,int> > frontier;
    for(int pos=begin_reg.first; pos<begin_reg.second+se.num_wc_padding.first; pos++){
        //align seq and pattern
        int i, j, mm;
        visited.clear();    
        
        frontier.push(make_pair(make_pair(pos, 0), 0));
        
        while(!frontier.empty()){
            if(visited.count(frontier.front())!=0){
                frontier.pop();
                continue;
            }
            i=frontier.front().first.first;
            j=frontier.front().first.second;
            mm=frontier.front().second;
            visited.insert(frontier.front());
            frontier.pop();
            
            //se.table.incOpsCounter();
            
            if(j<patt_length && i<seq_length){
                if(fits(seq[i], se.stripped_pattern[j])){ //if matches
                    frontier.push(make_pair(make_pair(i+1, j+1), mm));
                } else if(mm+1<=se.num_mismatches){ //if a mismatch
                    frontier.push(make_pair(make_pair(i+1, j+1), mm+1));
                }
                //se.table.incOpsCounter();
            }
            if(j<patt_length && se.stripped_pattern[j]=='*'){ //skip the wild card
                frontier.push(make_pair(make_pair(i, j+1), mm));
                //se.table.incOpsCounter();
            }
            
            //if it is a complete match, put it to the list of occurrences
            if(j==patt_length && mm<=se.num_mismatches){
                int start=pos;
                int end=i;
                //add also all possible matches of leading/trailing wild cards
                for(int prefix_l=0; prefix_l<=se.num_wc_padding.first; ++prefix_l){
                    for(int suffix_l=0; suffix_l<=se.num_wc_padding.second; ++suffix_l){
                        match.first = start - prefix_l;
                        match.second = end + suffix_l;
                        //check if the match is inside the search domain
                        if(end_reg.first < match.second && match.second <= end_reg.second
                            && match.first >= begin_reg.first && match.first < begin_reg.second)
                        {
                            //se.occurrences.insert( make_pair(end, 0) ); //at position 'pos-1' in seq ends a match
                            if(found_matches.count(match)==0){
                                found_matches.insert( match );
                                se.match_buffer.push( match );
                            }
                            #ifdef DEBUG
                                cout<<se.id<<" has match(s) "<<match.first+1<<" to "<<match.second<<" | ";
                                cout<<se.stripped_pattern<<" "<<seq.substr(match.first, match.second-match.first)<<endl;
                            #endif
                            
                            //se.table.incOpsCounter();
                        }
                    }
                }
                
            }
        }
    }
}

/* Run general DP pattern search for @se.pattern in @seq. All occurrences must begin
 * at index within <@begin_reg.first, @begin_reg.second) and end within @end_reg. Push all 
 * found matches to se.match_buffer.
 */
void Simple_Search::get_ss_matches(SSE &se, string &seq, interval &begin_reg, interval &end_reg){
    int patt_length=se.pattern.size();
    int seq_length=seq.size();
    interval match;
    
    stack< array<unsigned int, 7> > vertex_queue;
    for(int pos=begin_reg.first; pos<begin_reg.second; pos++){
        array<unsigned int, 7> tmp_vertex = {{0, 0, 0, 0, 0, 0, 0}};
        assert(vertex_queue.empty());
        se.table.clear();
        
        //seq is indexed from 0 to seq.size()-1, but we need it now from 1..seq.size()
        tmp_vertex[0] = pos;
        tmp_vertex[4] = 1;
        vertex_queue.push(tmp_vertex);
        
        //start flood
        while(!vertex_queue.empty()){
            //ignore already visited vertices
            while( !vertex_queue.empty() && se.table.get(vertex_queue.top())==true ) vertex_queue.pop();
            if( vertex_queue.empty() ) break;
            
            //get next unvisited vertex to be processed
            tmp_vertex=vertex_queue.top();
            vertex_queue.pop();
            
            /// process the vertex
            //mark the vertex as reachable
            se.table.set(tmp_vertex);
            
            //read coordinates of the vertex
            int i,j,m,n,b, x;
            i = tmp_vertex[0];
            j = tmp_vertex[1];
            m = tmp_vertex[2];
            n = tmp_vertex[3];
            b = tmp_vertex[4];
            
            //if it is a complete match, put it to the list of occurrences
            if(j==patt_length && (b==0 || n==0)){
                match.first = pos;
                match.second = i;
                if(end_reg.first < match.second && match.second <= end_reg.second)                {
                    //se.occurrences.insert(make_pair(i,0)); //at position (i-1) in seq ends a match
                    //cout<<i<<endl;
                    se.match_buffer.push( match );
                    
                    #ifdef DEBUG
                        cout<<se.id<<" has match(fwd) "<<match.first+1<<" to "<<match.second<<" | ";
                        cout<<se.pattern<<" "<<seq.substr(match.first, match.second-match.first)<<endl;
                    #endif
                }
            } else {
                //align next symbol from pattern to text if possible
                if(i<seq_length){
                    //not i+1 and j+1, because seq and se.pattern are indexed from 0
                    x = 1-(int)(fits(seq[i],se.pattern[j]));
                    //cout<<"porovnavam "<<seq[i]<<" k "<<se.pattern[j]<<endl;
                    
                    if(m+x<=se.num_mismatches){
                        tmp_vertex[0] = i+1;
                        tmp_vertex[1] = j+1;
                        tmp_vertex[2] = m+x;
                        tmp_vertex[3] = n;
                        tmp_vertex[4] = 0;
                        vertex_queue.push(tmp_vertex);
                    }
                }
                
                //skip a wild card if possible
                if(se.pattern[j]=='*'){
                    tmp_vertex[0] = i;
                    tmp_vertex[1] = j+1;
                    tmp_vertex[2] = m;
                    tmp_vertex[3] = n;
                    tmp_vertex[4] = b;
                    vertex_queue.push(tmp_vertex);
                }
                
                //do insertion if possible
                if(n<se.num_insertions && b==0 && i<seq_length && fits(seq[i],se.allowed_insertion)){
                    tmp_vertex[0] = i+1;
                    tmp_vertex[1] = j;
                    tmp_vertex[2] = m;
                    tmp_vertex[3] = n+1;
                    tmp_vertex[4] = 1;
                    vertex_queue.push(tmp_vertex);
                }
            }
            
        }
    }
}


/* DEPRECATED => no longer used in the battery of pattern matching algorithms
 * Do flood (breadth first search) from @seq[@index_in_seq]. In other words do forward
 * dynamic programming to compute ends of all possible matches starting at @index_in_seq.
 */
void Simple_Search::run_fwddp_ss(SSE &se, string &seq, interval &begin_reg){
    int patt_length=se.pattern.size();
    int seq_length=seq.size();

    stack< array<unsigned int, 7> > vertex_queue;
    for(int pos=begin_reg.first; pos<begin_reg.second; pos++){
        array<unsigned int, 7> tmp_vertex = {{0, 0, 0, 0, 0, 0, 0}};
        assert(vertex_queue.empty());
        
        //seq is indexed from 0 to seq.size()-1, but we need it now from 1..seq.size()
        tmp_vertex[0] = pos;
        tmp_vertex[4] = 1;
        vertex_queue.push(tmp_vertex);

        //start flood
        while(!vertex_queue.empty()){
            //ignore already visited vertices
            while( !vertex_queue.empty() && se.table.get(vertex_queue.top())==true ) vertex_queue.pop();
            if( vertex_queue.empty() ) break;

            //get next unvisited vertex to be processed
            tmp_vertex=vertex_queue.top();
            vertex_queue.pop();

        /// process the vertex
            //mark the vertex as reachable
            se.table.set(tmp_vertex);

            //read coordinates of the vertex
            int i,j,m,n,b, x;
            i = tmp_vertex[0];
            j = tmp_vertex[1];
            m = tmp_vertex[2];
            n = tmp_vertex[3];
            b = tmp_vertex[4];

            //if it is a complete match, put it to the list of occurrences
            if(j==patt_length && (b==0 || n==0)){
                se.occurrences.insert(make_pair(i,0)); //at position (i-1) in seq ends a match
                //cout<<i<<endl;
                #ifdef DEBUG
                    cout<<se.id<<" has match(fwd) "<<pos+1<<" to "<<i<<" | ";
                    cout<<se.pattern<<" "<<seq.substr(pos, i-pos)<<endl;
                #endif
            }

            //align next symbol from pattern to text if possible
            if(i+1<=seq_length && j+1<=patt_length){
                //not i+1 and j+1, because seq and se.pattern are indexed from 0
                x = 1-(int)(fits(seq[i],se.pattern[j]));
                //cout<<"porovnavam "<<seq[i]<<" k "<<se.pattern[j]<<endl;

                if(m+x<=se.num_mismatches){
                    tmp_vertex[0] = i+1;
                    tmp_vertex[1] = j+1;
                    tmp_vertex[2] = m+x;
                    tmp_vertex[3] = n;
                    tmp_vertex[4] = 0;
                    vertex_queue.push(tmp_vertex);
                }
            }

            //skip a wild card if possible
            if(j+1<=patt_length && se.pattern[j]=='*'){
                tmp_vertex[0] = i;
                tmp_vertex[1] = j+1;
                tmp_vertex[2] = m;
                tmp_vertex[3] = n;
                tmp_vertex[4] = b;
                vertex_queue.push(tmp_vertex);
            }

            //do insertion if possible
            if(b==0 && i+1<=seq_length && fits(seq[i],se.allowed_insertion) && n+1<=se.num_insertions){
                tmp_vertex[0] = i+1;
                tmp_vertex[1] = j;
                tmp_vertex[2] = m;
                tmp_vertex[3] = n+1;
                tmp_vertex[4] = 1;
                vertex_queue.push(tmp_vertex);
            }
        }
    }
}


/* DEPRECATED => no longer used in the battery of pattern matching algorithms 
 * Push into @se.match_buffer all matches of the single strand @se in the sequence @seq such that
 * seq[S..@end] is a correct match and S is in <@lower_bound, @upper_bound).
 * In other words do traceback of dynamic programming.
 */
void Simple_Search::trace_bckdp_ss(SSE &se, string &seq, interval &begin_reg, interval &end_reg){
    set<int> beginnings;
    stack< array<unsigned int, 7> > vertex_queue;
    Matrix visited(5);
    //set< array<unsigned int, 7> > visited;

    for(int end=end_reg.first; end<end_reg.second; end++){
        if(se.occurrences.find(make_pair(end+1, 0)) == se.occurrences.end()) continue;

        beginnings.clear();
        visited.clear();
        assert(vertex_queue.empty());
        
        #ifdef DO_CACHE
            map< int, set<int> >::iterator mapit = se.ss_beginnings_cache.find(end);
            if( mapit != se.ss_beginnings_cache.end() ){ //if in cache
                beginnings = mapit->second;
            } else {
        #endif

        array<unsigned int, 7> tmp_vertex = {{0, 0, 0, 0, 0, 0, 0}};
        
        tmp_vertex[0] = end+1;
        tmp_vertex[1] = se.pattern.size();

        for(int i=0; i<=se.num_insertions; i++){
            tmp_vertex[3] = i;         //set position corresponding to #insertions
            for(int j=0; j<=se.num_mismatches; j++){
                tmp_vertex[2] = j;      //set position corresponding to #mismatches
                vertex_queue.push(tmp_vertex);
                
                //special case for "all wild card patterns"
                //=> if all '*' were skiped when matching, flag 'b' ("prev insert") never got set back to zero
                tmp_vertex[3] = 0;    //set #insertions to zero
                tmp_vertex[4] = 1;    //"prev insert" true
                vertex_queue.push(tmp_vertex);
                tmp_vertex[4] = 0;
            }
        }

        //start traceback
        while(!vertex_queue.empty()){
            //ignore "bad" vertices
            while( !vertex_queue.empty() &&
                    //(visited.find(vertex_queue.top())!=visited.end() || se.table.get(vertex_queue.top())==false)
                    (se.table.get(vertex_queue.top())==false || visited.get(vertex_queue.top())==true)
                ) vertex_queue.pop();
            if( vertex_queue.empty() ) break;

            //get next vertex to be processed
            tmp_vertex=vertex_queue.top();
            vertex_queue.pop();

        /// process the vertex
            visited.set(tmp_vertex);

            //read coordinates of the vertex
            int i,j,m,n,b, x;
            i = tmp_vertex[0];
            j = tmp_vertex[1];
            m = tmp_vertex[2];
            n = tmp_vertex[3];
            b = tmp_vertex[4];

            //if we are at the beginning of the pattern then we have found a beginning of a match
            if(j==0 && b==1){
                beginnings.insert(i); //at position i in seq starts a match
                //cout<<i<<endl;
            }

            //we can have got here by aligning current symbol from pattern to text
            if(i-1>=0 && j-1>=0 && (b==0 || n==0)){
                //not i and j, because seq and se.pattern are indexed from 0
                x = 1-(int)(fits(seq[i-1],se.pattern[j-1]));

                if(m-x >= 0){
                    tmp_vertex[0] = i-1;
                    tmp_vertex[1] = j-1;
                    tmp_vertex[2] = m-x;
                    tmp_vertex[3] = n;
                    tmp_vertex[4] = 0;
                    vertex_queue.push(tmp_vertex);
                    tmp_vertex[4] = 1;
                    vertex_queue.push(tmp_vertex);
                }
            }

            //we can have got here by skipping a wild card
            if(j-1>=0 && se.pattern[j-1]=='*'){
                tmp_vertex[0] = i;
                tmp_vertex[1] = j-1;
                tmp_vertex[2] = m;
                tmp_vertex[3] = n;
                tmp_vertex[4] = b;
                vertex_queue.push(tmp_vertex);
            }

            //we can have got here by doing insertion
            if(b==1 && i-1>=0 && fits(seq[i-1],se.allowed_insertion) && n-1>=0){
                tmp_vertex[0] = i-1;
                tmp_vertex[1] = j;
                tmp_vertex[2] = m;
                tmp_vertex[3] = n-1;
                tmp_vertex[4] = 0;
                vertex_queue.push(tmp_vertex);
            }
        }

        #ifdef DO_CACHE
            //insert into the cache
            se.ss_beginnings_cache[end] = beginnings;
            }
        #endif

        assert(!beginnings.empty());
        for(set<int>::iterator itt=beginnings.begin(); itt!=beginnings.end(); ++itt){
            if(begin_reg.first <= *itt && *itt < begin_reg.second){
                se.match_buffer.push( make_pair(*itt, end+1) );
                #ifdef DEBUG
                    cout<<begin_reg.first<<" "<<begin_reg.second<<endl;
                    cout<<se.id<<" has match "<<*itt+1<<" to "<<end+1<<endl;
                #endif
            }
        }
    }
}

/* Run general DP pattern search for helical element @se in @seq. All occurrences must be within the
 * corrseponding search domain @domain. Push all matches to se.match_buffer.
 */
void Simple_Search::get_h_matches(SSE &se, string &seq, list<interval_pair> &domain){
    int patt_length=se.pattern.size();
    int seq_length=seq.size();
    int min_dist = se.strand_dist.first + 2*se.size_range.first -1;
    int max_dist = se.strand_dist.second + 2*se.size_range.second;
    
    stack< array<unsigned int, 7> > vertex_queue;
    for(int s1_b=domain.front().BEGIN.first; s1_b<domain.front().BEGIN.second; s1_b++)
    {
    for(int s2_e=max(domain.back().END.first, s1_b+min_dist); s2_e<min(domain.back().END.second, s1_b+max_dist); s2_e++)
    {
        array<unsigned int, 7> tmp_vertex = {{0, 0, 0, 0, 0, 0, 0}};
        assert(vertex_queue.empty());
        se.table.clear();
    
        //seq is indexed from 0 to seq.size()-1, but we need it now from 1..seq.size()
        tmp_vertex[0] = s1_b;
        tmp_vertex[1] = s2_e+2;
        tmp_vertex[6] = 1;
        vertex_queue.push(tmp_vertex);
        
        //start flood
        while(!vertex_queue.empty()){
            //ignore already visited vertices
            while( !vertex_queue.empty() && se.table.get(vertex_queue.top())==true ) vertex_queue.pop();
            if( vertex_queue.empty() ) break;
            
            //get next unvisited vertex to be processed
            tmp_vertex=vertex_queue.top();
            vertex_queue.pop();
            
            /// process the vertex
            //mark the vertex as reachable
            se.table.set(tmp_vertex);
            
            //read coordinates of the vertex
            int i,j,k,m,p,n,b, x,y;
            i = tmp_vertex[0];
            j = tmp_vertex[1];
            k = tmp_vertex[2];
            m = tmp_vertex[3];
            p = tmp_vertex[4];
            n = tmp_vertex[5];
            b = tmp_vertex[6];
            
            //if it is a complete match, put it to the list of occurrences
            if(k==patt_length && (b==0 || n==0)){
                //at position i in seq ends a match of 1.strand
                //at position j-1 in seq begins a match of 2.strand
                //cout<<i<<' '<<j-1<<endl;
                int s1_e = i;
                int s1_e_left = domain.front().END.first;
                int s1_e_right = domain.front().END.second;
                int s2_b = j-1;
                int s2_b_left = max(s1_e+se.strand_dist.first, domain.back().BEGIN.first);
                int s2_b_right = min(domain.back().BEGIN.second, s1_e+se.strand_dist.second+1);
                
                //cout<<s1_e_left<<" "<<s1_e_right<<" | "<<s2_b_left<<" "<<s2_b_right<<endl;
                
                if( s1_e_left < s1_e && s1_e <= s1_e_right &&
                    s2_b_left <= s2_b && s2_b < s2_b_right)
                {
                    se.match_buffer.push( make_pair(s1_b, s1_e) );
                    se.match_buffer.push( make_pair(s2_b, s2_e+1) );
                    #ifdef DEBUG
                        cout<<se.id<<" has match "<<s1_b+1<<" to "<<s1_e<<"; "<<s2_b+1<<" to "<<s2_e+1<<" | ";
                        cout<<" "<<seq.substr(s1_b, s1_e-s1_b);
                        cout<<" - "<<seq.substr(s2_b, s2_e+1-s2_b)<<"     "<<seq<<endl;
                    #endif
                }
                
            } else {
                //align next symbol from pattern to text if possible
                if(i<seq_length && j-1>0 && fits(seq[j-2],se.complement[k])){
                    //not i+1 nor j-1 nor k+1, because seq and se.pattern are indexed from 0
                    x = 1-(int)(fits(seq[i],se.pattern[k]));
                    y = 1-(int)(is_complemntary(seq[j-2],seq[i],se.transf_matrix));
                    
                    if(m+x<=se.num_mismatches && p+y<=se.num_mispairings){
                        tmp_vertex[0] = i+1;
                        tmp_vertex[1] = j-1;
                        tmp_vertex[2] = k+1;
                        tmp_vertex[3] = m+x;
                        tmp_vertex[4] = p+y;
                        tmp_vertex[5] = n;
                        tmp_vertex[6] = 0;
                        vertex_queue.push(tmp_vertex);
                    }
                }
                
                //skip a wild card if possible
                if(se.pattern[k]=='*'){
                    tmp_vertex[0] = i;
                    tmp_vertex[1] = j;
                    tmp_vertex[2] = k+1;
                    tmp_vertex[3] = m;
                    tmp_vertex[4] = p;
                    tmp_vertex[5] = n;
                    tmp_vertex[6] = b;
                    vertex_queue.push(tmp_vertex);
                }
                
                //do insertion to 1.strand if possible
                if(n<se.num_insertions && b==0 && i<seq_length && fits(seq[i],se.allowed_insertion)){
                    tmp_vertex[0] = i+1;
                    tmp_vertex[1] = j;
                    tmp_vertex[2] = k;
                    tmp_vertex[3] = m;
                    tmp_vertex[4] = p;
                    tmp_vertex[5] = n+1;
                    tmp_vertex[6] = 1;
                    vertex_queue.push(tmp_vertex);
                }
                
                //do insertion to 2.strand if possible
                if(n<se.num_insertions && b==0 && j-1>0 && fits(seq[j-2],se.allowed_insertion)){
                    tmp_vertex[0] = i;
                    tmp_vertex[1] = j-1;
                    tmp_vertex[2] = k;
                    tmp_vertex[3] = m;
                    tmp_vertex[4] = p;
                    tmp_vertex[5] = n+1;
                    tmp_vertex[6] = 1;
                    vertex_queue.push(tmp_vertex);
                }
            }
            
        }
    }
    }
}


/* DEPRECATED => no longer used in the battery of pattern matching algorithms 
 * Do flood (breadth first search) from @seq[@strand1_begin], @seq[@strand2_end].
 * In other words do forward dynamic programming to compute ends of all possible 1.strand matches
 * starting at @strand1_begin and beginnings of all possible 2.strand matches ending
 * at @strand2_end such that they form a correct helix according to @se.
 */
void Simple_Search::run_fwddp_h(SSE &se, string &seq, int strand1_begin, int strand2_end){
    //vector<int> tmp_vertex(7);
    array<unsigned int, 7> tmp_vertex = {{0, 0, 0, 0, 0, 0, 0}};
    stack< array<unsigned int, 7> > vertex_queue;
    int patt_length=se.pattern.size();
    int seq_length=seq.size();
    
    //seq is indexed from 0 to seq.size()-1, but we need it now from 1..seq.size()
    tmp_vertex[0] = strand1_begin;
    tmp_vertex[1] = strand2_end+2;
    tmp_vertex[6] = 1;
    vertex_queue.push(tmp_vertex);

    //start flood
    while(!vertex_queue.empty()){
        //ignore already visited vertices
        while( !vertex_queue.empty() && se.table.get(vertex_queue.top())==true ) vertex_queue.pop();
        if( vertex_queue.empty() ) break;

        //get next unvisited vertex to be processed
        tmp_vertex=vertex_queue.top();
        vertex_queue.pop();

      /// process the vertex
        //mark the vertex as reachable
        se.table.set(tmp_vertex);

        //read coordinates of the vertex
        int i,j,k,m,p,n,b, x,y;
        i = tmp_vertex[0];
        j = tmp_vertex[1];
        k = tmp_vertex[2];
        m = tmp_vertex[3];
        p = tmp_vertex[4];
        n = tmp_vertex[5];
        b = tmp_vertex[6];

        //if it is a complete match, put it to the list of occurrences
        if(k==patt_length && (b==0 || n==0)){
            //at position (i-1) in seq ends a match of 1.strand
            //at position (j) in seq begins a match of 2.strand
            //se.occurrences.set2(2,i,j);
            se.occurrences.insert(make_pair(i,j));
            //cout<<i-1<<' '<<j-1<<endl;
        }

        //align next symbol from pattern to text if possible
        if(i+1<=seq_length && j-1>0 && k+1<=patt_length && fits(seq[j-2],se.complement[k])){
            //not i+1 nor j-1 nor k+1, because seq and se.pattern are indexed from 0
            x = 1-(int)(fits(seq[i],se.pattern[k]));
            y = 1-(int)(is_complemntary(seq[j-2],seq[i],se.transf_matrix));

            if(m+x<=se.num_mismatches && p+y<=se.num_mispairings){
                tmp_vertex[0] = i+1;
                tmp_vertex[1] = j-1;
                tmp_vertex[2] = k+1;
                tmp_vertex[3] = m+x;
                tmp_vertex[4] = p+y;
                tmp_vertex[5] = n;
                tmp_vertex[6] = 0;
                vertex_queue.push(tmp_vertex);
            }
        }

        //skip a wild card if possible
        if(k+1<=patt_length && se.pattern[k]=='*'){
            tmp_vertex[0] = i;
            tmp_vertex[1] = j;
            tmp_vertex[2] = k+1;
            tmp_vertex[3] = m;
            tmp_vertex[4] = p;
            tmp_vertex[5] = n;
            tmp_vertex[6] = b;
            vertex_queue.push(tmp_vertex);
        }

        //do insertion to 1.strand if possible
        if(n+1<=se.num_insertions && b==0 && i+1<=seq_length && fits(seq[i],se.allowed_insertion)){
            tmp_vertex[0] = i+1;
            tmp_vertex[1] = j;
            tmp_vertex[2] = k;
            tmp_vertex[3] = m;
            tmp_vertex[4] = p;
            tmp_vertex[5] = n+1;
            tmp_vertex[6] = 1;
            vertex_queue.push(tmp_vertex);
        }

        //do insertion to 2.strand if possible
        if(n+1<=se.num_insertions && b==0 && j-1>0 && fits(seq[j-2],se.allowed_insertion)){
            tmp_vertex[0] = i;
            tmp_vertex[1] = j-1;
            tmp_vertex[2] = k;
            tmp_vertex[3] = m;
            tmp_vertex[4] = p;
            tmp_vertex[5] = n+1;
            tmp_vertex[6] = 1;
            vertex_queue.push(tmp_vertex);
        }
    }
}

/* DEPRECATED => no longer used in the battery of pattern matching algorithms
 * Push into @se.match_buffer all matches of the helix @se in the sequence @seq such that
 * seq[S..@end], seq[@begin..E] is a correct match and S is in <@lower_bound1, @upper_bound1)
 * and E is in <@lower_bound2, @upper_bound2). In other words do traceback of dynamic programming.
 */
void Simple_Search::trace_bckdp_h(SSE &se, string &seq, int lower_bound1, int upper_bound1, int end,
                                  int begin, int lower_bound2, int upper_bound2){

    if(se.occurrences.find(make_pair(end+1, begin+1)) == se.occurrences.end()) return;

    int seq_length=seq.size();
    set< pair<int,int> > beginnings;
    #ifdef DO_CACHE
        map< pair<int,int>, set< pair<int,int> > >::iterator mapit = se.h_beginnings_cache.find(make_pair(end,begin));
        if( mapit != se.h_beginnings_cache.end() ){ //if in cache
            beginnings = mapit->second;
        } else {
    #endif

    //vector<int> tmp_vertex(7);
    array<unsigned int, 7> tmp_vertex = {{0, 0, 0, 0, 0, 0, 0}};
    stack< array<unsigned int, 7> > vertex_queue;
    //Matrix visited(7);
    set< array<unsigned int, 7> > visited;

    tmp_vertex[0] = end+1;
    tmp_vertex[1] = begin+1;
    tmp_vertex[2] = se.pattern.size();

    for(int i=0; i<=se.num_insertions; i++){
        tmp_vertex[5] = i;           //set position corresponding to #insertions
        for(int j=0; j<=se.num_mispairings; j++){
            tmp_vertex[4] = j;       //set position corresponding to #mispairings
            for(int k=0; k<=se.num_mismatches; k++){
                tmp_vertex[3] = k;    //set position corresponding to #mismatches
                vertex_queue.push(tmp_vertex);
            }
            
            //special case for "all wild card patterns"
            //=> if all '*' were skiped when matching, flag 'b' ("prev insert") never got set back to zero
            tmp_vertex[5] = 0;    //set #insertions to zero
            tmp_vertex[6] = 1;    //"prev insert" true
            vertex_queue.push(tmp_vertex);
            tmp_vertex[6] = 0;
        }
    }

    //start traceback
    while(!vertex_queue.empty()){
        //ignore "bad" vertices
        while( !vertex_queue.empty() &&
                (se.table.get(vertex_queue.top())==false || visited.find(vertex_queue.top())!=visited.end())
             ) vertex_queue.pop();
        if( vertex_queue.empty() ) break;

        //get next vertex to be processed
        tmp_vertex=vertex_queue.top();
        vertex_queue.pop();

      /// process the vertex
        visited.insert(tmp_vertex);

        //read coordinates of the vertex
        int i,j,k,m,p,n,b, x,y;
        i = tmp_vertex[0];
        j = tmp_vertex[1];
        k = tmp_vertex[2];
        m = tmp_vertex[3];
        p = tmp_vertex[4];
        n = tmp_vertex[5];
        b = tmp_vertex[6];
        //printf("trace cez %d %d %d %d %d %d %d \n",i,j,k,m,p,n,b);

        //if we are at the beginning of the pattern then we have found a begining of helical match
        if(k==0 && b==1){
            //at position i in seq (indexed from 0) starts a 1.strand match
            //at position (j-2) in seq (indexed from 0) ends a 2.strand match
            beginnings.insert( make_pair(i,j-2) );
           // cout<<"+beginning: "<<i<<"  "<<j-2<<endl;
        }

        //we can have got here by aligning current symbols from pattern to text
        if(i-1>=0 && j<=seq_length && k-1>=0 && fits(seq[j-1],se.complement[k-1])){
            //recall that seq, se.pattern and se.complement are indexed from 0
            x = 1-(int)(fits(seq[i-1],se.pattern[k-1]));
            y = 1-(int)(is_complemntary(seq[j-1],seq[i-1],se.transf_matrix));
            //cout<<"is complementary? "<<seq[i-1]<<' '<<seq[j-1]<<"    "<<se.pattern[k-1]<<' '<<se.complement[k-1]<<endl;

            if(m-x >= 0 && p-y >= 0){
                tmp_vertex[0] = i-1;
                tmp_vertex[1] = j+1;
                tmp_vertex[2] = k-1;
                tmp_vertex[3] = m-x;
                tmp_vertex[4] = p-y;
                tmp_vertex[5] = n;
                tmp_vertex[6] = 0;
                vertex_queue.push(tmp_vertex);
                tmp_vertex[6] = 1;
                vertex_queue.push(tmp_vertex);
            }
        }

        //we can have got here by skipping a wild card
        if(k-1>=0 && se.pattern[k-1]=='*'){
            tmp_vertex[0] = i;
            tmp_vertex[1] = j;
            tmp_vertex[2] = k-1;
            tmp_vertex[3] = m;
            tmp_vertex[4] = p;
            tmp_vertex[5] = n;
            tmp_vertex[6] = b;
            vertex_queue.push(tmp_vertex);
        }

        //we can have got here by doing insertion to 1.strand
        if(b==1 && i-1>=0 && fits(seq[i-1],se.allowed_insertion) && n-1>=0){
            tmp_vertex[0] = i-1;
            tmp_vertex[1] = j;
            tmp_vertex[2] = k;
            tmp_vertex[3] = m;
            tmp_vertex[4] = p;
            tmp_vertex[5] = n-1;
            tmp_vertex[6] = 0;
            vertex_queue.push(tmp_vertex);
        }

        //we can have got here by doing insertion to 2.strand
        if(b==1 && j<=seq_length && fits(seq[j-1],se.allowed_insertion) && n-1>=0){
            tmp_vertex[0] = i;
            tmp_vertex[1] = j+1;
            tmp_vertex[2] = k;
            tmp_vertex[3] = m;
            tmp_vertex[4] = p;
            tmp_vertex[5] = n-1;
            tmp_vertex[6] = 0;
            vertex_queue.push(tmp_vertex);
        }
    }

    #ifdef DO_CACHE
        //insert into the cache
        se.h_beginnings_cache[make_pair(end, begin)] = beginnings;
        }
    #endif

    assert(!beginnings.empty());
    for(set< pair<int,int> >::iterator itt=beginnings.begin(); itt!=beginnings.end(); ++itt){
        if( lower_bound1 <= itt->first && itt->first < upper_bound1 &&
            lower_bound2 <= itt->second && itt->second < upper_bound2){
                se.match_buffer.push( make_pair(itt->first, end+1) );
                se.match_buffer.push( make_pair(begin, itt->second+1) );
                #ifdef DEBUG
                cout<<se.id<<" has match "<<itt->first<<" to "<<end+1<<"; "<<begin<<" to "<<itt->second+1<<endl;
                #endif
        }
    }
}

//returns next match from the @se.match_buffer
list<interval> Simple_Search::get_next_match(SSE &se){
    list<interval> new_match;

    /// get next match from match_buffer
    if(se.is_helix){ //if a helix
        if(!se.match_buffer.empty()){
            new_match.push_back( se.match_buffer.front() );
            se.match_buffer.pop();
            new_match.push_back( se.match_buffer.front() );
            se.match_buffer.pop();
        } else {
            new_match.push_back( make_pair(-1,-1) );
            new_match.push_back( make_pair(-1,-1) );
        }

    } else { //if a single strand
        if(!se.match_buffer.empty()){
            new_match.push_back( se.match_buffer.front() );
            se.match_buffer.pop();
        } else {
            new_match.push_back( make_pair(-1,-1) );
        }
    }

    #ifdef DEBUG
        cout<<"match: "<<se.id<<": ";
        cout<<new_match.front().first<<" to "<<new_match.front().second<<": ";
        //if(new_match.front().first != -1)
        //    cout<<seq.substr(new_match.front().first,new_match.front().second-new_match.front().first)<<"   ";
        if(new_match.size()==2){
            cout<<new_match.back().first<<" to "<<new_match.back().second<<": ";
            //if(new_match.front().first != -1)
            //    cout<<seq.substr(new_match.back().first,new_match.back().second-new_match.back().first)<<"   ";
        }
    #endif
    return new_match;
}

void Simple_Search::set_grid(intervals &grid, SSE &se, list<interval> &match){
    for(int i=0;i<desc->motif.size();i++){
        if( desc->motif[i] == se.id) grid[i] = match.front();
        else if( -desc->motif[i] == se.id) grid[i] = match.back();
    }
}

void Simple_Search::reset_grid(intervals &grid, SSE &se){
    for(int i=0;i<desc->motif.size();i++){
        if( desc->motif[i] == se.id) grid[i] = make_pair(-1,-1);
        else if( -desc->motif[i] == se.id) grid[i] = make_pair(-1,-1);
    }
}

string Simple_Search::solution_to_str(unsigned int ind, string &seq, string &separator){
    if(ind >= solutions.size()) return "";

    string str = separator;
    for(int i=0;i<solutions[ind].size();i++){
        str+= seq.substr(solutions[ind][i].first, solutions[ind][i].second-solutions[ind][i].first);
        str+= separator;

        if(i>0) assert(solutions[ind][i].first==solutions[ind][i-1].second);
    }
    return str;
}

string Simple_Search::solution_to_dotbracket(unsigned int ind, string &separator){
    if(ind >= solutions.size()) return "";
    
    string res = separator;
    for(int i=0;i<solutions[ind].size();i++){
        char c;
        int elementid = abs(desc->motif[i]);
       
        if(desc->sses[elementid].is_helix == true){
            if(desc->motif[i]>0){ //if it is the first strand of the helix
                switch(desc->pknot_levels[elementid]){
                  case 0: c = '(';
                          break;
                  case 1: c = '[';
                          break;
                  case 2: c = '{';
                          break;
                  default: c = '-'; die("Unexpected situation in Descriptor::get_dotbracket_notation()");
                }
            } else {        //if it is the second strand of the helix
                switch(desc->pknot_levels[elementid]){
                  case 0: c = ')';
                          break;
                  case 1: c = ']';
                          break;
                  case 2: c = '}';
                          break;
                  default: c = '+'; die("Unexpected situation in Descriptor::get_dotbracket_notation()");
                }
            }
            
        } else {
            c = '.';
        }
            
        for(int j=0;j<solutions[ind][i].second-solutions[ind][i].first;j++){
            res += c;
        }
        res += separator;
    }
    
    return res;
}
