/*
 * Project      : RNA motif searching in genomic sequences
 * Description  : the implementation of a class representing a motif descriptor
 *
 * Author       : Ladislav Rampasek <rampasek@gmail.com>
 * Institution  : Comenius University in Bratislava
 *
 */

#include <sstream>
#include <algorithm>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <stack>
#include <iostream>
#include <stdint.h>
#include <emmintrin.h>
#include <cstring>
#include <cassert>

#include "descriptor.h"
#include "types.h"
#include "generalfuncs.h"

#define CHECK_OK -1

using namespace std;
using namespace GF;

Descriptor::Descriptor(ifstream &fin){
    initialized=false;
    err_str="";

    bool desc_ok = true;
    string line;

    // read the map of the descriptor
    desc_ok &= get_valuable_line(fin,line);
    desc_ok &= parse_desc_map(line);
    if(! desc_ok) { err_str="Invalid descriptor topology map"; return; }

    // read properties of sses declared in the map
    while( get_valuable_line(fin,line) ){
        normalize_seq(line.begin(), line.end());
        desc_ok &= parse_desc_properties(line);
        if(! desc_ok) desc_ok = parse_desc_order(line);
        if(! desc_ok) {
            err_str="Invalid descriptor line: "+line+
                "\nInvalid syntax of an secondary structure element or of a reorder command";
            return;
        }
    }

    // check whether all sses are defined and descriptor contains no inconsistencies
    int err_code = check_consistency();
    if(err_code != CHECK_OK) {
        err_str="The descriptor is incomplete or contains a logical error";
        
        if(err_code == -2){
            err_str += " in the reorder command";
        } else {
            vector<int> err_element_id(1);
            err_element_id[0] = err_code;
            string err_element_name = search_order_to_str(err_element_id);
            
            err_str += " in element " + err_element_name;
        }
        
        return;
    }

    // compute auxiliary stats, e.g. sizes of SSEs, distances between strands of helices
    compute_auxiliary_stats();

    //compile patterns that can be searched for by BNDM
    for(int i=1;i<sses.size();i++){
        if(sses[i].is_helix) continue;

        // single strand element with fix-sized core and with NO other wild cards, nor insertions will be search by BNDM;
        // for the other ss elements will be used forward filtering
        bool fixed_core = (sses[i].size_range.first==(sses[i].size_range.second-sses[i].num_wc_padding.first-sses[i].num_wc_padding.second));
        bool is_bckwd = fixed_core && sses[i].num_insertions==0;
        if(sses[i].stripped_pattern.size() <= 128){
            compile_pattern(sses[i], !is_bckwd);
        }
    }

    
    // compute information content of all SSEs
    compute_inf_contents();
    
    //compute_pknot_levels();

    initialized=true;
}

bool Descriptor::parse_desc_map(string line){
    transform(line.begin(), line.end(), line.begin(), ::toupper); //uppercase the line
    istringstream sin(line);

    sses.clear();
    sses.resize(1); //to have this array indexed from 1 !!!
    motif.clear();
    transl.clear();

    string item;
    int open_helices=0;
    while( sin>>item ){
        if( item[item.length()-1] == '\'' ){ //if it is the second strand of a helix
            map<string, int>::iterator it=transl.find( item.substr(0,item.length()-1) );
            if( it==transl.end() ) return false;
            motif.push_back( -(it->second) );
            open_helices--;
        } else {
            SSE new_sse;
            switch (item[0]){
                case 'S': new_sse.is_helix=false;
                          break;
                case 'R':
                case 'H': new_sse.is_helix=true;
                          open_helices++;
                          break;
                default: return false;
            }
            new_sse.id=sses.size();
            sses.push_back(new_sse);
            motif.push_back(sses.size()-1);
            transl[item]=motif.back();
        }
    }

    //for(int i=0;i<motif.size();i++) cout<<motif[i]<<" "; cout<<endl;
    if( !has_no_duplicates(motif) ) return false;
    return open_helices==0;
}

bool Descriptor::parse_desc_properties(string line){
    istringstream sin(line);
    string item;
    char first_char_in_name;

    //get "name" of the sse
    if(! (sin>>item)) return false; //cout<<item<<endl;
    //transform(item.begin(), item.end(), item.begin(), ::toupper);
    map<string, int>::iterator it=transl.find(item);
    if( it==transl.end() ) return false; //if it is unknown
    first_char_in_name=item[0];
    sses[it->second].name=item;
    sses[it->second].recordOps=false;
    
    if( sses[it->second].is_helix){ //if sse is helix
        //get number of tolerated mismatches, mispairings and insertions
        if(! (sin>>item)) return false; //cout<<item<<endl;
        int tmp_num=0;
        for(int i=0;i<item.size();i++) if(item[i]==':') { item[i]=' '; tmp_num++; }
        if( tmp_num != 1 && tmp_num != 2) return false;

        istringstream sin2(item);
        if( !(sin2>>sses[it->second].num_mismatches) || !(sin2>>sses[it->second].num_mispairings) ) return false;
        //by defalut we do not allow insertions
        sses[it->second].num_insertions=0;
        if( tmp_num == 2 && !(sin2>>sses[it->second].num_insertions)) return false;

        //get primary structure restrictions
        if(! (sin>>item)) return false; //cout<<item<<endl;
        tmp_num=0;
        for(int i=0;i<item.size();i++) if(item[i]==':') { item[i]=' '; tmp_num++; }
        if( tmp_num != 1 && tmp_num != 2) return false;

        sin2.clear();
        sin2.str(item); //cout<<item<<endl;
        if( !(sin2>>sses[it->second].pattern) || !(sin2>>sses[it->second].complement) ) {
            //cout<<sses[it->second].pattern<<endl;
            return false;
        }
        //by defalut we allow insertion of all nucleotides
        sses[it->second].allowed_insertion='N';
        if( tmp_num == 2 && !(sin2>>sses[it->second].allowed_insertion)) return false;

        //get transformation matrix (if it is R)
        switch(first_char_in_name){
            case 'R': if(! (sin>>item)) return false; //cout<<item<<endl;
                        if(item.size()!=4) return false;
                        //transform(item.begin(), item.end(), item.begin(), ::toupper);
                        sses[it->second].transf_matrix = item;
                        break;
            case 'H': sses[it->second].transf_matrix = "TGYR";
                        break;
            default: return false; //unexpected value of first_char_in_name
        }
        //cout<<sses[it->second].transf_matrix<<endl;

        if(sin>>item) return false; //no more input is expected in this line

    } else { //if sse is single strand
        //get number of tolerated mismatches and insertions
        if(! (sin>>item)) return false; //cout<<item<<endl;
        int tmp_num=0;
        for(int i=0;i<item.size();i++) if(item[i]==':') { item[i]=' '; tmp_num++; }
        if( tmp_num != 0 && tmp_num != 1) return false;

        istringstream sin2(item);
        if( !(sin2>>sses[it->second].num_mismatches) ) return false;
        //by defalut we do not allow insertions
        sses[it->second].num_insertions=0;
        if( tmp_num == 1 && !(sin2>>sses[it->second].num_insertions)) return false;

        //get primary structure restrictions
        if(! (sin>>item)) return false; //cout<<item<<endl;
        tmp_num=0;
        for(int i=0;i<item.size();i++) if(item[i]==':') { item[i]=' '; tmp_num++; }
        if( tmp_num != 0 && tmp_num != 1) return false;

        sin2.clear();
        sin2.str(item); //cout<<item<<endl;
        if( !(sin2>>sses[it->second].pattern) ) return false;
        //by defalut we allow insertion of all nucleotides
        sses[it->second].allowed_insertion='N';
        if( tmp_num == 1 && !(sin2>>sses[it->second].allowed_insertion)) return false;

        if(sin>>item) return false; //no more input is expected in this line
    }
    return true;
}

bool Descriptor::parse_desc_order(string line){
    //transform(line.begin(), line.end(), line.begin(), ::toupper); //uppercase the line
    istringstream sin(line);

    string item;
    if(! (sin>>item)) return false;
    if(item!="R") return false;

    predef_srch_order.clear();
    map<string, int>::iterator it;
    while( sin>>item ){
        it=transl.find(item);
        if( it==transl.end() ) return false; //if it is unknown
        predef_srch_order.push_back(it->second);
    }

    if( !has_no_duplicates(predef_srch_order) ) return false;
    return true;
}

// returns true if vector contains no duplicates
bool Descriptor::has_no_duplicates(vector<int> vec){
    if(vec.size() < 2) return true;

    sort(vec.begin(), vec.end());
    for(int i=1;i<vec.size();i++) if(vec[i-1]==vec[i]) return false;

    return true;
}

//expand [x] to (*)^x in the given string
bool Descriptor::expand_wildcards(string &s){
    string news;
    news.reserve(s.size());
    int i=0;

    while(i<s.size()){
        if(s[i]!='[') {
            news+=s[i];
            ++i;
            continue;
        }

        ++i;
        int num=0;
        while(i<s.size() && s[i]>='0' && s[i]<='9'){
            num*=10;
            num+=s[i] - (int)'0';
            ++i;
        }
        if(i>=s.size() || s[i]!=']') return false;
        for(int j=0;j<num;j++) news+='*';
        ++i;
    }
    s=news;
    return true;
}

// consistency check + wildcards exansion
int Descriptor::check_consistency(){
    /*
    for(int i=1;i<sses.size();i++){
        cout<<i<<": \n";
        cout<<sses[i].pattern<<" "<<sses[i].complement<<endl;
        cout<<sses[i].num_mismatches<<" "<<sses[i].num_mispairings<<endl<<endl;
    }
    for(int i=0;i<srch_order.size();i++) cout<<srch_order[i]<<" "; cout<<endl;
    */

    for(int i=1;i<sses.size();i++){
        if (sses[i].is_helix){
            //a helix
            if( !expand_wildcards(sses[i].pattern) ) return sses[i].id;
            if( !expand_wildcards(sses[i].complement) ) return sses[i].id;
            if( sses[i].pattern.size() != sses[i].complement.size() ) return sses[i].id;
            if( sses[i].pattern.size() == 0) return sses[i].id;
            if( sses[i].num_mismatches > sses[i].pattern.size()) return sses[i].id;
            if( sses[i].num_mispairings > sses[i].pattern.size()) return sses[i].id;

            reverse(sses[i].complement.begin(), sses[i].complement.end()); //reverse the complement
            for(int j=0;j<sses[i].pattern.size();j++){ //complementarity check
                //if there is a wild card in one of the strands, there must be one also in the other strand
                if((sses[i].pattern[j]=='*') != (sses[i].complement[j]=='*')) return sses[i].id;
                //check complementarity
                if(!is_complemntary(sses[i].pattern[j], sses[i].complement[j], sses[i].transf_matrix)) return sses[i].id;
            }
        } else {
            //a single strand
            if( !expand_wildcards(sses[i].pattern) ) return sses[i].id;
            if( sses[i].pattern.size() == 0) return sses[i].id;
            if( sses[i].num_mismatches > sses[i].pattern.size()) return sses[i].id;
        }
    }

    if( !has_no_duplicates(predef_srch_order) ) return -2;
    return CHECK_OK;
}

void Descriptor::compute_auxiliary_stats(){
    //compute size range of sses
    for(int i=1;i<sses.size();i++){
        int count=0;
        for(int j=0;j<sses[i].pattern.size();j++) if(sses[i].pattern[j]=='*') ++count;

        sses[i].size_range.first = sses[i].pattern.size() - count;
        sses[i].size_range.second = sses[i].pattern.size() + sses[i].num_insertions;
    }
    //compute strand distance of helices
    for(int i=1;i<sses.size();i++){
        if (sses[i].is_helix){
            bool start=false;
            int min_dis=0, max_dis=0;
            for(int j=0;j<motif.size();j++){
                if(motif[j] == -i) break;
                if(start){
                    min_dis+=sses[abs(motif[j])].size_range.first;
                    max_dis+=sses[abs(motif[j])].size_range.second;
                }
                if(motif[j] == i) start=true;
            }
            sses[i].strand_dist.first = min_dis;
            sses[i].strand_dist.second = max_dis;
        }
    }
    //get counts of leading/trailing wild cards and get stripped patterns
    for(int i=1;i<sses.size();i++){
        //get the wild card prefix/suffix lengths
        int num_leading=0, num_trailing=0;
        int j=0, count=0;
        bool b=true;
        while(j<sses[i].pattern.size()){
            if(sses[i].pattern[j]=='*'){
                ++count;
            } else {
                if(b){
                    num_leading = count; //length of wild card prefix
                    b = false;
                }
                count = 0;
            }
            j++;
        }
        num_trailing = count; //length of wild card suffix
        
        sses[i].num_wc_padding = make_pair(num_leading, num_trailing);
            
        //strip the patterns of wild card prefix and suffix
        if (sses[i].is_helix){
            sses[i].stripped_pattern = sses[i].pattern.substr(num_leading, sses[i].pattern.size()-num_trailing-num_leading);
            sses[i].stripped_complement = sses[i].complement.substr(num_leading, sses[i].complement.size()-num_trailing-num_leading);
        } else {
            sses[i].stripped_pattern = sses[i].pattern.substr(num_leading, sses[i].pattern.size()-num_trailing-num_leading);
        }
    }
    
    /*
    for(int i=1;i<sses.size();i++){
        cout<<i<<": \n";
        cout<<"size_range "<<sses[i].size_range.first<<" "<<sses[i].size_range.second<<endl;
        cout<<"strand_dist "<<sses[i].strand_dist.first<<" "<<sses[i].strand_dist.second<<endl;
        cout<<"seq(comp) "<<sses[i].pattern<<" "<<sses[i].complement<<endl;
        cout<<"mismach(pair) "<<sses[i].num_mismatches<<" "<<sses[i].num_mispairings<<endl;
        cout<<"insert. "<<sses[i].num_insertions<<" "<<sses[i].allowed_insertion<<endl;
        cout<<"wc_padding "<<sses[i].num_wc_padding.first<<" "<<sses[i].num_wc_padding.second<<endl;
        cout<<"stripped pattern    "<<sses[i].stripped_pattern<<" vs "<<sses[i].pattern<<endl;
        cout<<"stripped complement "<<sses[i].stripped_complement<<" vs "<<sses[i].complement<<endl;
        cout<<endl;
    }
    */
}

// Alignment must be power of 2
void* aligned_malloc(size_t size, size_t alignment){
    //| unused_data0 | start_address | aligned_data | unused_data1 |
    uintptr_t r = (uintptr_t)malloc(size + --alignment + 2);
    if (!r) return NULL;
    //aligned offset
    uintptr_t o = (r + 2 + alignment) & ~(uintptr_t)alignment;
    //store address of the actual start of the allocated memory
    ((uint16_t*)o)[-1] = (uint16_t)(o-r);
    return (void*)o;
}

inline void setbit(void *v, int p) {
    ((uint32_t*)v)[p >> 5] |= 1 << (p & 31);
}

//compile pattern for Forward filtering resp. Backward search in BNDM
void Descriptor::compile_pattern(SSE &se, bool isFwdPattern){
    int patt_length=se.stripped_pattern.size();
    assert(patt_length <= 128);
    if(patt_length == 0) return;
    int j;
    __m128i zero = {};
    
    se.maskv = (__m128i*)aligned_malloc(256 * sizeof(__m128i), 16);
    __m128i I = {}, F = {}, nF = {}, one = {}, ones = {};
    
    for(int i=0; i<256; ++i) se.used[i] = 0;

    //IUPAC codes
    vector<string> iupac(4);
    iupac[0] = "NWMRDHV";  //these contain A
    iupac[1] = "NSMYBHV";  //these contain C
    iupac[2] = "NSKRBDV";  //these contain G
    iupac[3] = "NWKYBDH";  //these contain T
    string ambig_codes = "NWSMKRYBDHV";
    
    ///PREPROCESSING
    ///precompute se.maskv for forward filtering / backward search in BNDM
    int i = 0;
    while(i < patt_length) {
        //use the upper bits in the mask (align to left), i.e 128 ... (128 - patt_length)
        int pos_in_mask = (isFwdPattern) ? 128 - patt_length + i : 128 - 1 - i;

        //wild cards are admissible ONLY in forward filtering
        //NOTE: wildcards can't be a prefix or suffix of the stripped_pattern
        if(se.stripped_pattern[i] == '*'){
            assert(isFwdPattern);
            int block_end = i;
            while(block_end<patt_length && se.stripped_pattern[block_end]=='*') ++block_end;
            
            setbit(&I, 128 - patt_length + i - 1); //"gap-initial" state
            setbit(&F, 128 - patt_length + block_end); //"gap-final" state
            
            //set that every nucleotide matches a wild card
            for(int k=i; k<block_end; ++k){
                for(int ind=0; ind<4; ++ind){
                    unsigned int nuc = "ACGT"[ind];
                    if (!se.used[nuc]){
                        se.used[nuc] = 1;
                        se.maskv[nuc] = zero;
                    }
                    setbit(&se.maskv[nuc], 128 - patt_length + k);
                }
            }
            i = block_end;
            
        //if it's a IUPAC code for multiple nucleotides - "classes in pattern"
        } else if(ambig_codes.find(se.stripped_pattern[i]) != string::npos){
            //allow this position to match all the nucleotides it codes for
            //(by adding this position to their mask of matching positions)
            for(int ind=0; ind<4; ++ind){
                unsigned int nuc = "ACGT"[ind];
                //if iupac code se.stripped_pattern[i] includes nucelotide nuc
                if(iupac[ind].find(se.stripped_pattern[i]) != string::npos){
                    if (!se.used[nuc]){
                        se.used[nuc] = 1;
                        se.maskv[nuc] = zero;
                    }
                    setbit(&se.maskv[nuc], pos_in_mask);
                }
            }
            ++i;

        //simple character
        } else {
            j = se.stripped_pattern[i];
            if (!se.used[j]){
                se.used[j] = 1;
                se.maskv[j] = zero;
            }
            setbit(&se.maskv[j], pos_in_mask);
            ++i;
        }
    }
    
    ///precompute bitwise negation of F
    for(int i=0; i<4; ++i) ((uint32_t*) &nF)[i] = ~((uint32_t*) &F)[i];
    // store I, F, ~F in se.maskv[0..2]
    assert(!se.used[0] && !se.used[1] && !se.used[2]);
    se.used[0] = 1, se.maskv[0] = I;
    se.used[1] = 1, se.maskv[1] = F;
    se.used[2] = 1, se.maskv[2] = nF; 
    //_mm_storeu_si128(&se.maskv[i], maskv[i]);
    
    ///precompute initial masks
    //set 1 at position that corresponds to the beginning of the pattern in the bit mask
    setbit(&one, (128 - patt_length));
    
    //set to 1^(patt_length)0^(128 - patt_length)
    ones = _mm_set_epi32(0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF); //initialize to all ones
    for(int i = 0; i < (128 - patt_length)/8; ++i) //shift by bytes
        ones = _mm_slli_si128(ones, 1); 
    for(int i = 0; i < (128 - patt_length)%8; ++i) //shift by bits to the proper position
        ones = _mm_or_si128(_mm_slli_epi64(ones, 1), _mm_srli_epi64(_mm_slli_si128(ones, 8), 63));
    
    //store to se.maskv[3, 4]
    assert(!se.used[3] && !se.used[4]);
    se.used[3] = 1, se.maskv[3] = one;
    se.used[4] = 1, se.maskv[4] = ones;
    
    ///"classes in text" - when searched sequence contains ambiguous IUPAC codes
    for(int i=0; i<ambig_codes.size(); ++i){ //initialize to zeros first
        se.used[(unsigned int)ambig_codes[i]] = 1;
        se.maskv[(unsigned int)ambig_codes[i]] = zero;
    }
    //for each ambiguous letter: set its mask to union of positions that match the coded nucelotides
    for(int i=0; i<4; ++i){
        unsigned int nuc = "ACGT"[i];
        for (j=0; j<7; ++j){
            if(se.used[nuc]) se.maskv[(unsigned int)iupac[i][j]] |= se.maskv[nuc];
        }
    }
}

inline double log2(double x){
    return log(x) / log(2);
}

const static double ninf = -1. * INFINITY;
inline double logsum(double x, double y){
    if(x == ninf) return y;
    if(y == ninf) return x;
    double a = max(x, y);
    double b = min(x, y);
    return a + log2(1. + pow(2., b-a) );
}

//compute information content of all elements
void Descriptor::compute_inf_contents(){
    string bases = "ACGT";
    string one = "ACGTU";
    string two = "MRWSYK";
    string three = "VHBD";

    //precompute binomial coefficients in log-space
    int maxlength = 0;
    for(int i = 1; i < sses.size(); ++i) maxlength = max(maxlength, sses[i].size_range.second);
    ++maxlength;

    double logbin[maxlength][maxlength];
    logbin[0][0] = 0;
    for(int n = 1; n < maxlength; n++) {
        logbin[n][0] = 0;
        logbin[n][n] = 0;
        for(int k = 1; k < n; k++) {
            logbin[n][k] = logsum(logbin[n - 1][k], logbin[n - 1][k - 1]);
        }
    }

    //compute IC for each element
    for(int i = 1; i < sses.size(); ++i){
        double logX = 0, entropy_before = 0;
        int N = sses[i].size_range.second;
        
        double insEntropy = 0;
        if(sses[i].allowed_insertion == 'N'){
            insEntropy += log2(4);
        } else if(one.find(sses[i].allowed_insertion) != string::npos){
            insEntropy += log2(1);
        } else if(two.find(sses[i].allowed_insertion) != string::npos){
            insEntropy += log2(2);
        } else if(three.find(sses[i].allowed_insertion) != string::npos){
            insEntropy += log2(3);
        }
        
        ///Helix ---------------------------------------------------------------->
        if (sses[i].is_helix){
            entropy_before = 2 * N * log2(4);
            double meanP = 0, countP = 0;
            double P = 0;
            
            //calculate P
            for(int p1=0; p1<4; ++p1){
                for(int p2=0; p2<4; ++p2){
                    if(is_complemntary(bases[p1], bases[p2], sses[i].transf_matrix)){
                        ++P;
                    }
                }
            }
            
            //1.) count (in log-space) number of sequences matching sses[i] without any distortion or wild card
            for(int j=0; j<sses[i].pattern.size(); j++){
                int Ppos = 0;
                if(sses[i].pattern[j]=='*' || sses[i].complement[j]=='*') continue;
                
                for(int p1=0; p1<4; ++p1){
                    for(int p2=0; p2<4; ++p2){
                        if( fits(bases[p1], sses[i].pattern[j]) &&
                            fits(bases[p2], sses[i].complement[j]) &&
                            is_complemntary(bases[p1], bases[p2], sses[i].transf_matrix)
                        ){
                            ++Ppos;
                        }
                    }
                }
                logX += log2(Ppos);
                meanP += Ppos;
                ++countP;
            }
            meanP /= countP; //average number of nuleotide pairs allowed per non-wildcard pair
            //cout<<"mean P:"<<meanP<<"   P:"<<P<<endl;
            
            //2.) add blocks of wild cards
            int block_size = 0;
            int pos = 0;
            while(pos < sses[i].pattern.size()){
                while(pos < sses[i].pattern.size() && sses[i].pattern[pos] != '*') ++pos;
                while(pos < sses[i].pattern.size() && sses[i].pattern[pos] == '*'){
                    ++pos;
                    ++block_size;
                }
                
                double tmp = ninf;
                for(int j=0; j<=block_size; ++j){
                    tmp = logsum(tmp, j * log2(P) + (block_size-j) * log2(16-P));
                }
                logX += tmp;
                block_size = 0;
            }
            
            //3.) include mismatches
            double tmp = ninf;
            for(int j=1; j<sses[i].num_mismatches; j++){
                tmp = logsum(tmp, logbin[(int)countP][j] + j*log2((P-meanP) / meanP) );
            }
            logX = logsum(logX, logX + tmp);
            
            //4.) include mispairs
            tmp = ninf;
            for(int j=1; j<sses[i].num_mispairings; j++){
                tmp = logsum(tmp, logbin[N][j] + j*log2((16-meanP) / meanP) );
            }
            logX = logsum(logX, logX + tmp);
            
            //5.) insertions
            tmp = ninf;
            for(int j=1; j<sses[i].num_insertions; j++){
                tmp = logsum(tmp, logbin[max(1, 2*N - 4)][j] + j*insEntropy + (2*sses[i].num_insertions - j)*log2(4));
            }
            logX = logsum(logX + 2 * sses[i].num_insertions * log2(4), logX + tmp);
            
        ///Single Strand -------------------------------------------------------->
        } else {
            entropy_before = N * log2(4);
            double meanC = 0, countC = 0;
            
            //1.) count (in log-space) number of sequences matching sses[i] without any distortion or wild card
            for(int j=0; j<sses[i].pattern.size(); j++){
                char c = sses[i].pattern[j];
                if(c == 'N'){
                    logX += log2(4);
                } else if(c == '*'){
                    ; //ignore wild cards for now
                } else if(one.find(c) != string::npos){
                    logX += log2(1);
                    meanC += 1;
                    ++countC;
                } else if(two.find(c) != string::npos){
                    logX += log2(2);
                    meanC += 2;
                    ++countC;
                } else if(three.find(c) != string::npos){
                    logX += log2(3);
                    meanC += 3;
                    ++countC;
                }
            }
            meanC /= countC; //average number of nuleotides allowed per non-"N"/"*" position 
            
            //2.) add blocks of wild cards
            int block_size = 0;
            int pos = 0;
            while(pos < sses[i].pattern.size()){
                while(pos < sses[i].pattern.size() && sses[i].pattern[pos] != '*') ++pos;
                while(pos < sses[i].pattern.size() && sses[i].pattern[pos] == '*'){
                    ++pos;
                    ++block_size;
                }
                logX += block_size * log2(4) + log2(block_size + 1);
                block_size = 0;
            }
            
            //3.) include mismatches
            double tmp = ninf;
            for(int j=1; j<sses[i].num_mismatches; j++){
                tmp = logsum(tmp, logbin[(int)countC][j] + j*log2((4.-meanC) / meanC) );
            }
            logX = logsum(logX, logX + tmp);
            
            //4.) insertions
            tmp = ninf;
            for(int j=1; j<sses[i].num_insertions; j++){
                tmp = logsum(tmp, logbin[max(1, N-2)][j] + j*insEntropy + (sses[i].num_insertions - j)*log2(4));
            }
            logX = logsum(logX + sses[i].num_insertions * log2(4), logX + tmp);
        }
        
        sses[i].infContent = entropy_before - logX;
        
        //debug output
        /*
        for(map<string, int>::iterator it=transl.begin();it!=transl.end();++it){
            if( it->second == sses[i].id ){
                cout<<it->first<<"("<<sses[i].id<<") before = "<<entropy_before<<"  after = "<<logX<<" IC= "<<sses[i].infContent<<endl;
                break;
            }
        }*/
    }
}

//returns the order in which elements of the motif are going to be search in human-readable form
string Descriptor::search_order_to_str(vector<int> order){
    string result;
    for(int i=0;i<order.size();i++){
        for(map<string, int>::iterator it=transl.begin();it!=transl.end();++it){
            if( it->second == order[i] ){
                result+=it->first+' ';
                break;
            }
        }
    }
    return result;
}

//returns the maximal length of the motif specified by the instance of Descriptor
int Descriptor::get_max_motif_length(){
    int result=0;
    for(int i=0;i<motif.size();i++){
        result+=sses[abs(motif[i])].size_range.second;
    }
    return result;
}

//compute which helical elements cause a pseudoknot and of which level, pseudoknots
// of up to level 2 are supported; the result is stored in vector pknot_levels
void Descriptor::compute_pknot_levels(){
    pknot_levels.clear();
    pknot_levels.resize(motif.size());
    
    stack<int> hstack;
    stack<int> temp;
    
    
    //decide which helical elements cause a pseudoknot and of which level
    for(int i=0;i<motif.size();i++){
        int elementid = abs(motif[i]);
        
        if(sses[elementid].is_helix == true){
            if(motif[i]>0){ //if it is the first strand of the helix
                hstack.push(elementid);
            } else {        //if it is the second strand of the helix
                int top = hstack.top();
                hstack.pop();
                
                while(!temp.empty()) temp.pop();
                
                //take out all helices crossing the particular one
                while(top != elementid){ 
                    temp.push(top);
                    top = hstack.top();
                    hstack.pop();
                }
                
                //adjust "pseudoknot level" of the crossing helcies
                while(!temp.empty()){
                    top = temp.top();
                    temp.pop();
                    
                    hstack.push(top);
                                        
                    if(pknot_levels[top] <= pknot_levels[elementid]){
                        pknot_levels[top] = pknot_levels[elementid] + 1;
                    }
                }
            }
            
        }
    }

}

//returns secondary structure representation in dot-bracked notation, pseudoknots
// of up to level 2 are supported
string Descriptor::get_dotbracket_notation(){
    string res;

    //create the dot-bracked anotation according to the determined levels
    for(int i=0;i<motif.size();i++){
        char c;
        int elementid = abs(motif[i]);
       
        if(sses[elementid].is_helix == true){
            if(motif[i]>0){ //if it is the first strand of the helix
                switch(pknot_levels[elementid]){
                  case 0: c = '(';
                          break;
                  case 1: c = '[';
                          break;
                  case 2: c = '{';
                          break;
                  default: c = '-'; die("Unexpected situation + in Descriptor::get_dotbracket_notation()");
                }
            } else {        //if it is the second strand of the helix
                switch(pknot_levels[elementid]){
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
            
        for(int j=0;j<sses[abs(motif[i])].size_range.second;j++){
            res += c;
        }
    }
    
    return res;
}
