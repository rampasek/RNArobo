/*
 * Project      : RNA motif searching in genomic sequences
 * Description  : the header file for generalfuncs.cpp
 *
 * Author       : Ladislav Rampasek <rampasek@gmail.com>
 * Institution  : Comenius University in Bratislava
 *
 */

#ifndef GENERALFUNCS_H
#define GENERALFUNCS_H

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <algorithm>
using namespace std;

namespace GF {

#define _A 0
#define _C 1
#define _G 2
#define _T 3
    
// end program with an error message
inline int die(string errmsg){
    fprintf(stderr, "ERROR: %s\n", errmsg.c_str());
    //cerr<<"ERROR: "<<errmsg<<"\n";
    //cerr<<"The program finished with an error.\n";
    exit(1);
}

// get next noncomment nonempty line from given ifstream
inline bool get_valuable_line(ifstream &fin, string &str){
    string line;
    string::iterator it;
    while(getline(fin, line)){
        it=line.begin();
        while(isspace(*it)) ++it;

        if(it == line.end()) continue;

        if(*it != '#') {
            str = line;
            return true;
        }
    }

    str = "";
    return false;
}

//normalize the given string = bring to uppercase and replace 'U's by 'T's
inline void normalize_seq(string::iterator begin, string::iterator end){
    for(string::iterator it=begin; it!=end; ++it){
        if('a' <= *it && *it <= 'z') *it += ('A'-'a');
        if(*it == 'U') *it = 'T';
    }
}


const static string iupacA = "WMRDHV";  //these contain A
const static string iupacC = "SMYBHV";  //these contain C
const static string iupacG = "SKRBDV";  //these contain G
const static string iupacT = "WKYBDH";  //these contain T

// evaluate if @ch fits @patt according to IUPAC Notation
inline bool fits(char ch, char &patt){
    //handle most often cases first
    if(patt=='N' || ch=='N' || patt=='*') return true;
    //if(ch=='U') ch='T';  
    if(ch==patt) return true;
    
    //the more complicated cases
    switch (patt){
        case 'A': //Adenine
            return iupacA.find(ch)!=string::npos;
        case 'C': //Cytosine
            return iupacC.find(ch)!=string::npos;
        case 'G': //Guanine
            return iupacG.find(ch)!=string::npos;
        //case 'U':
        case 'T': //Thymine in DNA; Uracil in RNA
            return iupacT.find(ch)!=string::npos;
        case 'M': //aMino
            return ch=='A' || ch=='C' || iupacA.find(ch)!=string::npos || iupacC.find(ch)!=string::npos;
        case 'R': //puRine
            return ch=='A' || ch=='G' || iupacA.find(ch)!=string::npos || iupacG.find(ch)!=string::npos;
        case 'W': //Weak
            return ch=='A' || ch=='T' || iupacA.find(ch)!=string::npos || iupacT.find(ch)!=string::npos;
        case 'S': //Strong
            return ch=='C' || ch=='G' || iupacC.find(ch)!=string::npos || iupacG.find(ch)!=string::npos;
        case 'Y': //pYrimidine
            return ch=='C' || ch=='T' || iupacC.find(ch)!=string::npos || iupacT.find(ch)!=string::npos;
        case 'K': //Keto
            return ch=='G' || ch=='T' || iupacG.find(ch)!=string::npos || iupacT.find(ch)!=string::npos;
        case 'V': //not T
            return ch!='T' && iupacT.find(ch)==string::npos;
        case 'H': //not G
            return ch!='G' && iupacG.find(ch)==string::npos;
        case 'D': //not C
            return ch!='C' && iupacC.find(ch)==string::npos;
        case 'B': //not A
            return ch!='A' && iupacA.find(ch)==string::npos;
    }
    return false;
}

// evaluate if @ch_s1 is complementary to @ch_s2 according to @transf_matrix and IUPAC Notation
inline bool is_complemntary(char &ch_s1, char &ch_s2, string &transf_matrix){
    //ch_s1=toupper(ch_s1);
    //ch_s2=toupper(ch_s2);

    if(ch_s1=='N' || ch_s2=='N') return true;
    if(ch_s1=='*' && ch_s2=='*') return true; //if both are wildcards - ok
    if(ch_s1=='*' || ch_s2=='*') return false; //if only one is wildcard - bad

    switch (ch_s1){
        case 'A': //Adenine
            return fits(ch_s2, transf_matrix[_A]);
        case 'C': //Cytosine
            return fits(ch_s2, transf_matrix[_C]);
        case 'G': //Guanine
            return fits(ch_s2, transf_matrix[_G]);
        //case 'U':
        case 'T': //Thymine in DNA; uracil in RNA
            return fits(ch_s2, transf_matrix[_T]);
        case 'M': //aMino
            return fits(ch_s2, transf_matrix[_A]) || fits(ch_s2, transf_matrix[_C]);
        case 'R': //puRine
            return fits(ch_s2, transf_matrix[_A]) || fits(ch_s2, transf_matrix[_G]);
        case 'W': //Weak
            return fits(ch_s2, transf_matrix[_A]) || fits(ch_s2, transf_matrix[_T]);
        case 'S': //Strong
            return fits(ch_s2, transf_matrix[_C]) || fits(ch_s2, transf_matrix[_G]);
        case 'Y': //pYrimidine
            return fits(ch_s2, transf_matrix[_C]) || fits(ch_s2, transf_matrix[_T]);
        case 'K': //Keto
            return fits(ch_s2, transf_matrix[_G]) || fits(ch_s2, transf_matrix[_T]);
        case 'V': //not T
            return fits(ch_s2, transf_matrix[_A]) || fits(ch_s2, transf_matrix[_C]) || fits(ch_s2, transf_matrix[_G]);
        case 'H': //not G
            return fits(ch_s2, transf_matrix[_A]) || fits(ch_s2, transf_matrix[_C]) || fits(ch_s2, transf_matrix[_T]);
        case 'D': //not C
            return fits(ch_s2, transf_matrix[_A]) || fits(ch_s2, transf_matrix[_G]) || fits(ch_s2, transf_matrix[_T]);
        case 'B': //not A
            return fits(ch_s2, transf_matrix[_C]) || fits(ch_s2, transf_matrix[_G]) || fits(ch_s2, transf_matrix[_T]);
    }
    return false;
}

// make complement of the given char
inline void in_place_complement(string::iterator i) {
    switch (*i){
        case 'A':
            *i = 'T';
            break;
        case 'C':
            *i = 'G';
            break;
        case 'G':
            *i = 'C';
            break;
        case 'T':
            *i = 'A';
            break;
        case 'N':
            //N doesn't change
            break;
        case 'M': //aMino
            *i = 'K';
            break;
        case 'R': //puRine
            *i = 'Y';
            break;
        case 'W': //Weak
            *i = 'S';
            break;
        case 'S': //Strong
            *i = 'W';
            break;
        case 'Y': //pYrimidine
            *i = 'R';
            break;
        case 'K': //Keto
            *i = 'M';
            break;
        case 'V': //not T
            *i = 'T';
            break;
        case 'H': //not G
            *i = 'G';
            break;
        case 'D': //not C
            *i = 'C';
            break;
        case 'B': //not A
            *i = 'A';
            break;
    }
}

// reverse and make complement of the given substring
inline void reverse_complement(string::iterator begin, string::iterator end) {
    while ((begin != end) && (begin != end-1)) {
        in_place_complement(begin);
        in_place_complement(--end);
        std::swap(*begin++, *end);
    }
    if(begin != end) in_place_complement(begin); //negate the middle nucleotide
}

// filter out all whitespaces from the given string @s
inline void filter_whitespaces(string &s){
    string tmp;
    unsigned int last = 0;
    unsigned int i = 0;
    unsigned int size = s.size();
    //tmp.reserve(size);
    
    while(i<size){
        while(i<size && !isspace(s[i])) ++i;
        if(i==size && last==0) return; //no whitespaces in s
        tmp += s.substr(last, i-last);
        
        while(i<size && isspace(s[i])) ++i;
        last = i;
    }
    s=tmp;
    
    /* slow version
    string tmp;
    tmp.reserve(s.size());
    for(int i=0;i<s.size();++i){
        if(!isspace(s[i])) tmp+=s[i];
    }
    s=tmp;
    */
}

}
#endif
