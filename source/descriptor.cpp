/*
 * $Id: descriptor.cpp,v 1.11 2012-04-21 15:47:52 laci Exp $
 *
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

#include "descriptor.h"
#include "types.h"
#include "generalfuncs.h"

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
        normalize_seq(line);
        desc_ok &= parse_desc_properties(line);
        if(! desc_ok) desc_ok = parse_desc_order(line);
        if(! desc_ok) {
            err_str="Invalid descriptor line: "+line+
                "\nInvalid syntax of an secondary structure element or of a reorder command";
            return;
        }
    }

    // check whether all sses are defined and descriptor contains no inconsistencies
    if(! check_consistency() ) {
        err_str="The descriptor is not complete or contains a logical error";
        return;
    }

    // compute sizes of SSEs and distances between strands of helices
    compute_sizes();
    
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

// consistency check + wildcars exansion
bool Descriptor::check_consistency(){
    /*
    for(int i=1;i<sses.size();i++){
        cout<<i<<": \n";
        cout<<sses[i].pattern<<" "<<sses[i].complement<<endl;
        cout<<sses[i].num_mismatches<<" "<<sses[i].num_mispairings<<endl<<endl;
    }
    for(int i=0;i<srch_order.size();i++) cout<<srch_order[i]<<" "; cout<<endl;
    */

    for(int i=1;i<sses.size();i++){
        switch (sses[i].is_helix){
            case true: //a helix
                if( !expand_wildcards(sses[i].pattern) ) return false;
                if( !expand_wildcards(sses[i].complement) ) return false;
                if( sses[i].pattern.size() != sses[i].complement.size() ) return false;
                if( sses[i].pattern.size() == 0) return false;
                if( sses[i].num_mismatches > sses[i].pattern.size()) return false;
                if( sses[i].num_mispairings > sses[i].pattern.size()) return false;

                reverse(sses[i].complement.begin(), sses[i].complement.end()); //reverse the complement
                for(int j=0;j<sses[i].pattern.size();j++){ //complementarity check
                    if(!is_complemntary(sses[i].pattern[j], sses[i].complement[j], sses[i].transf_matrix)) return false;
                }
                break;
            case false: //a single strand
                if( !expand_wildcards(sses[i].pattern) ) return false;
                if( sses[i].pattern.size() == 0) return false;
                if( sses[i].num_mismatches > sses[i].pattern.size()) return false;
                break;
        }
    }

    if( !has_no_duplicates(predef_srch_order) ) return false;
    return true;
}

void Descriptor::compute_sizes(){
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

    /*
    for(int i=1;i<sses.size();i++){
        cout<<i<<": \n";
        cout<<"size_range "<<sses[i].size_range.first<<" "<<sses[i].size_range.second<<endl;
        cout<<"strand_dist "<<sses[i].strand_dist.first<<" "<<sses[i].strand_dist.second<<endl;
        cout<<"seq(comp) "<<sses[i].pattern<<" "<<sses[i].complement<<endl;
        cout<<"mismach(pair) "<<sses[i].num_mismatches<<" "<<sses[i].num_mispairings<<endl;
        cout<<"insert. "<<sses[i].num_insertions<<" "<<sses[i].allowed_insertion<<endl<<endl;
    }
    */
}

//compute information content of all elements
void Descriptor::compute_inf_contents(){
    string bases = "ACGT";
    string one="ACGTU";
    string two="MRWSYK";
    string three="VHBD";
    
    for(int i=1;i<sses.size();i++){
        double ic = 0;
        ///helix
        if (sses[i].is_helix){
            map< pair<char,char>, double> cacheIC;
            for(int j=0; j<sses[i].pattern.size(); j++){
                char c1 = sses[i].pattern[j];
                char c2 = sses[i].complement[j];
                pair<char,char> cc = make_pair(c1, c2);
                
                double pairIC = 0;
                
                //if we have calculated this pair before
                if (cacheIC.count(cc) > 0){
                    pairIC = cacheIC[cc];
                } else { //else calculate it and store it
                    //first, calculate the number of base pairs that match the motif pair
                    int num = 0;
                    for(int p1=0; p1<4; ++p1){
                        for(int p2=0; p2<4; ++p2){
                            if( fits(bases[p1], c1) &&
                                fits(bases[p2], c2) &&
                                is_complemntary(bases[p1], bases[p2], sses[i].transf_matrix)
                            ){
                                ++num;
                            }
                        }
                    }
                    //if the pair is ('*','*'), count in also the option of skipping that pair 
                    if(c1 == '*' && c1 == c2){ 
                        ++num;
                    }
                    
                    //information content is:
                    pairIC = 4 - log(num)/log(2);
                    cacheIC[cc] = pairIC;
                }
                
                //cout<<"PAIR: "<<c1<<","<<c2<<" pairIC= "<<pairIC<<endl;
                ic += pairIC;
            }
            ic -= 2*sses[i].num_mispairings;
        ///single strand
        } else {
            for(int j=0; j<sses[i].pattern.size(); j++){
                char c = sses[i].pattern[j];
                
                double strandIC = 0;
                if(c == 'N'){
                    strandIC += 0; // 2 - log(4)/log(2)
                } else if(c == '*'){
                    strandIC += 2 - log(5)/log(2);
                } else if(one.find(c) != string::npos){
                    strandIC += 2; // 2 - log(1)/log(2)
                } else if(two.find(c) != string::npos){
                    strandIC += 1; // 2 - log(2)/log(2)
                } else if(three.find(c) != string::npos){
                    strandIC += 2 - log(3)/log(2);
                }
                
                //cout<<"S: "<<c<<" strandIC= "<<strandIC<<endl;
                ic += strandIC;
            }
        }
        
        ic -= 2*sses[i].num_mismatches;
        
        if(sses[i].num_insertions > 0){
            char c = sses[i].allowed_insertion;
            
            double insEntropy = 0;
            if(c == 'N'){
                insEntropy += 2; // log(4)/log(2)
            } else if(c == '*'){
                insEntropy += log(5)/log(2);
            } else if(one.find(c) != string::npos){
                insEntropy += 0; // log(1)/log(2)
            } else if(two.find(c) != string::npos){
                insEntropy += 1; // log(2)/log(2)
            } else if(three.find(c) != string::npos){
                insEntropy += log(3)/log(2);
            }
            
            if(sses[i].is_helix){
                ic -= sses[i].num_insertions*( log(2*sses[i].pattern.size()-2)/log(2) + insEntropy );
            } else {
                ic -= sses[i].num_insertions*( log(sses[i].pattern.size()-1)/log(2) + insEntropy );
            }
        }
        
        sses[i].infContent = ic;
        
        //debug output
        /*
        for(map<string, int>::iterator it=transl.begin();it!=transl.end();++it){
            if( it->second == sses[i].id ){
                cout<<it->first<<"("<<sses[i].id<<") IC = "<<ic<<endl<<endl;
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
