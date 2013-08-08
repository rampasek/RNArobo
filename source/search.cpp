/*
 * Project      : RNA motif searching in genomic sequences
 * Description  : the search core implementation
 *
 * Author       : Ladislav Rampasek <rampasek@gmail.com>
 * Institution  : Comenius University in Bratislava
 *
 */

//#define NDEBUG
#include <assert.h>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <tr1/array>
#include <list>
#include <stack>
#include <queue>
#include <set>

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
    // DP for helices
    if(se.is_helix){
        int min_dist = se.strand_dist.first + 2*se.size_range.first -1;
        int max_dist = se.strand_dist.second + 2*se.size_range.second;

        for(int i=domain.front().BEGIN.first; i<domain.front().BEGIN.second; i++){
            for(int j=max(domain.back().END.first, i+min_dist); j<min(domain.back().END.second, i+max_dist); j++){
                run_fwddp_h(se, seq, i, j);
                //cout<<"flooding "<<i<<' '<<j<<endl;
            }
        }

        for(int i=domain.front().END.first; i<domain.front().END.second; i++){
            for( int j=max(i+se.strand_dist.first+1, domain.back().BEGIN.first);
                 j<min(domain.back().BEGIN.second, i+se.strand_dist.second+2); j++ )
            {
                //printf("tracing <%d, %d) %d; %d <%d, %d)  -> %d\n", domain.front().BEGIN.first, domain.front().BEGIN.second, i, j, domain.back().END.first, domain.back().END.second,se.occurrences.get(2,i+1,j+1));
                trace_bckdp_h(se, seq, domain.front().BEGIN.first, domain.front().BEGIN.second, i,
                              j, domain.back().END.first, domain.back().END.second);
            }
        }
        
    // single strand element with NO wild cards or insertions
    } else if(se.size_range.first==se.size_range.second && se.num_insertions==0){
        get_naive_ss_matches(se, seq, domain.front().BEGIN, domain.front().END);

    // single strand element with NO mismatches or insertions
    } else if(se.num_insertions==0){
        get_simple_ss_matches(se, seq, domain.front().BEGIN, domain.front().END);
      
    // general single strand element DP
    } else {
        run_fwddp_ss(se, seq, domain.front().BEGIN);
        //printf("tracing <%d, %d) - <%d, %d)\n", domain.front().BEGIN.first, domain.front().BEGIN.second, domain.front().END.first, domain.front().END.second);
        trace_bckdp_ss(se, seq, domain.front().BEGIN, domain.front().END);
    }

    list<interval> match = get_next_match(se);

    while(match.front().first != -1 && !have_solution){
        set_grid(grid, se, match);
        if( se.id == orderer->searchOrder.back()){
            #ifdef DEBUG
                cout<<"SOLUTION!!!"<<endl;
            #endif
            solutions.push_back(grid);
            #ifdef SKIP
                have_solution=true;
            #endif
        } else {
            #ifdef DEBUG
                cout<<"DOWN"<<endl;
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
    for(int i=index_in_motif-1;i>=0;i--){
        if( grid[i].second != -1){
            left_search_bound = grid[i].second + min_left_shift;
            left_ncover_bound = grid[i].second + max_left_shift;
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
    for(int i=index_in_motif+1;i<desc->motif.size();i++){
        if( grid[i].first != -1){
            right_search_bound = grid[i].first - min_right_shift;
            right_ncover_bound = grid[i].first - max_right_shift;
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
    //  if -1 is a boundary, it means there is no boundary
    //  but if both boundaries are -1, it means there IS NOT an interval to be necessarily covered


    //combine search interval, necessary cover interval and minimal size of the element to optain domain
    int lower_begin_bound = left_search_bound;
    //int upper_begin_bound = max(left_search_bound+1, right_search_bound - desc->sses[ abs(desc->motif[index_in_motif]) ].size_range.first);
    int upper_begin_bound = max(left_search_bound+1, right_search_bound);
        //if necessary cover interval beginning is defined then a match must start at/before it
        if(left_ncover_bound != -1) upper_begin_bound = min(left_ncover_bound+1, upper_begin_bound);

    //int lower_end_bound = min(right_search_bound-1, left_search_bound + desc->sses[ abs(desc->motif[index_in_motif]) ].size_range.first);
    int lower_end_bound = min(right_search_bound-1, left_search_bound );
        //if necessary cover interval end is defined then a match must end at/after it
        if(right_ncover_bound != -1) lower_end_bound = max(right_ncover_bound-1, lower_end_bound);
    int upper_end_bound = right_search_bound;

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

/* Run naive pattern search for @se.pattern in @seq. All occurrences must begin
 * at index within @begin_reg.first...@begin_reg.second and end within @end_reg.
 * !!!Works for single strand elements with *NO wild cards or insertions*!!!
 */
void Simple_Search::get_naive_ss_matches(SSE &se, string &seq, interval &begin_reg, interval &end_reg){
    int patt_length=se.pattern.size();
    int seq_length=seq.size();
    for(int i=begin_reg.first; i<begin_reg.second; i++){
        //align seq and pattern
        int j=0;
        int mm=0;
        if(se.num_mismatches==0){
            while(j<patt_length && i+j<seq_length && fits(seq[i+j], se.pattern[j])){
                j++;
                se.table.incOpsCounter();
            }
        } else { //allow mismatches
            while(j<patt_length && i+j<seq_length && mm<=se.num_mismatches){
                if(fits(seq[i+j], se.pattern[j])){
                    j++;
                } else {
                    j++;
                    mm++;
                }
                se.table.incOpsCounter();
            }
        }
        //if it is a complete match, put it to the list of occurrences
        if(j==patt_length && mm<=se.num_mismatches){
            int end=i+j;
            if(end_reg.first < end && end <= end_reg.second){
                //se.occurrences.insert( make_pair(end, 0) ); //at position 'i+j-1' in seq ends a match
                se.match_buffer.push( make_pair(i, end) );
                #ifdef DEBUG
                cout<<se.id<<" has match "<<i+1<<" to "<<end<<endl;
                #endif
            }
        }
    }
}

/* Run simple DP pattern search for @se.pattern in @seq. All occurrences must begin
 * at index within @begin_reg.first...@begin_reg.second and end within @end_reg.
 * !!!Works for single strand elements with *NO mismatches or insertions*!!!
 */
void Simple_Search::get_simple_ss_matches(SSE &se, string &seq, interval &begin_reg, interval &end_reg){
    int patt_length=se.pattern.size();
    int seq_length=seq.size();

    set < pair<interval,int> > visited;
    queue < pair<interval,int> > frontier;
    for(int pos=begin_reg.first; pos<begin_reg.second; pos++){
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
            
            se.table.incOpsCounter();
            
            if(j<patt_length && i<seq_length){
                if(fits(seq[i], se.pattern[j])){ //if matches
                    frontier.push(make_pair(make_pair(i+1, j+1), mm));
                } else if(mm+1<=se.num_mismatches){ //if a mismatch
                    frontier.push(make_pair(make_pair(i+1, j+1), mm+1));
                }
                se.table.incOpsCounter();
            }
            if(j<patt_length && se.pattern[j]=='*'){ //skip the wild card
                frontier.push(make_pair(make_pair(i, j+1), mm));
                se.table.incOpsCounter();
            }
            
            //if it is a complete match, put it to the list of occurrences
            if(j==patt_length && mm<=se.num_mismatches){
                int end=i;
                if(end_reg.first < end && end <= end_reg.second){
                    //se.occurrences.insert( make_pair(end, 0) ); //at position 'i-1' in seq ends a match
                    se.match_buffer.push( make_pair(pos, end) );
                    #ifdef DEBUG
                    cout<<se.id<<" has match "<<pos+1<<" to "<<end<<endl;
                    #endif
                }
            }
        }
    }
}

/* Do flood (breadth first search) from @seq[@index_in_seq]. In other words do forward
 * dynamic programming to compute ends of all possible matches starting at @index_in_seq.
 */
void Simple_Search::run_fwddp_ss(SSE &se, string &seq, interval &begin_reg){
    int patt_length=se.pattern.size();
    int seq_length=seq.size();

    stack< tr1::array<unsigned int, 7> > vertex_queue;
    for(int pos=begin_reg.first; pos<begin_reg.second; pos++){
        tr1::array<unsigned int, 7> tmp_vertex = {{0, 0, 0, 0, 0, 0, 0}};
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
            if(j==patt_length && b==0){
                se.occurrences.insert(make_pair(i,0)); //at position (i-1) in seq ends a match
                //cout<<i<<endl;
                #ifdef DEBUG
                    cout<<se.id<<" has match "<<pos+1<<" to "<<i<<" | ";
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


/*  Push into @se.match_buffer all matches of the single strand @se in the sequence @seq such that
 *  seq[S..@end] is a correct match and S is in <@lower_bound, @upper_bound).
 *  In other words do traceback of dynamic programming.
 */
void Simple_Search::trace_bckdp_ss(SSE &se, string &seq, interval &begin_reg, interval &end_reg){
    set<int> beginnings;
    stack< tr1::array<unsigned int, 7> > vertex_queue;
    Matrix visited(5);
    //set< tr1::array<unsigned int, 7> > visited;

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

        tr1::array<unsigned int, 7> tmp_vertex = {{0, 0, 0, 0, 0, 0, 0}};
        
        tmp_vertex[0] = end+1;
        tmp_vertex[1] = se.pattern.size();

        for(int i=0; i<=se.num_insertions; i++){
            tmp_vertex[3] = i;         //set position corresponding to #insertions
            for(int j=0; j<=se.num_mismatches; j++){
                tmp_vertex[2] = j;      //set position corresponding to #mismatches
                vertex_queue.push(tmp_vertex);
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
            if(i-1>=0 && j-1>=0 && b==0){
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

/* Do flood (breadth first search) from @seq[@strand1_begin], @seq[@strand2_end].
    In other words do forward dynamic programming to compute ends of all possible 1.strand matches
    starting at @strand1_begin and beginnings of all possible 2.strand matches ending
    at @strand2_end such that they form a correct helix according to @se.
*/
void Simple_Search::run_fwddp_h(SSE &se, string &seq, int strand1_begin, int strand2_end){
    //vector<int> tmp_vertex(7);
    tr1::array<unsigned int, 7> tmp_vertex = {{0, 0, 0, 0, 0, 0, 0}};
    stack< tr1::array<unsigned int, 7> > vertex_queue;
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
        if(k==patt_length && b==0){
            //at position (i-1) in seq ends a match of 1.strand
            //at position (j-1) in seq ends a match of 2.strand
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


/* Push into @se.match_buffer all matches of the helix @se in the sequence @seq such that
    seq[S..@end], seq[@begin..E] is a correct match and S is in <@lower_bound1, @upper_bound1)
    and E is in <@lower_bound2, @upper_bound2). In other words do traceback of dynamic programming.
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
    tr1::array<unsigned int, 7> tmp_vertex = {{0, 0, 0, 0, 0, 0, 0}};
    stack< tr1::array<unsigned int, 7> > vertex_queue;
    //Matrix visited(7);
    set< tr1::array<unsigned int, 7> > visited;

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
        if(new_match.front().first != -1);
            //cout<<seq.substr(new_match.front().first,new_match.front().second-new_match.front().first)<<"   ";
        if(new_match.size()==2){
            cout<<new_match.back().first<<" to "<<new_match.back().second<<": ";
            if(new_match.front().first != -1);
               // cout<<seq.substr(new_match.back().first,new_match.back().second-new_match.back().first)<<"   ";
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
