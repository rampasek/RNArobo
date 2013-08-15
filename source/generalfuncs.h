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

// evaluate if @ch fits @patt according to IUPAC Notation
inline bool fits(char &ch, char &patt){
    //ch=toupper(ch);
    //patt=toupper(patt);
    
    if(patt=='N' || ch=='N' || patt=='*') return true;
    
    //if(ch=='U') ch='T';  
    if(ch==patt) return true;
    
    switch (patt){
        /*case 'A': //Adenine
            return ch=='A';
        case 'C': //Cytosine
            return ch=='C';
        case 'G': //Guanine
            return ch=='G';
        case 'U':
        case 'T': //Thymine in DNA; uracil in RNA
            return ch=='T';*/
        case 'M': //aMino
            return ch=='A' || ch=='C';
        case 'R': //puRine
            return ch=='A' || ch=='G';
        case 'W': //Weak (2 H bonds)
            return ch=='A' || ch=='T';
        case 'S': //Strong (3 H bonds)
            return ch=='C' || ch=='G';
        case 'Y': //pYrimidine
            return ch=='C' || ch=='T';
        case 'K': //Keto
            return ch=='G' || ch=='T';
        case 'V': //not T
            return ch!='T';
        case 'H': //not G
            return ch!='G';
        case 'D': //not C
            return ch!='C';
        case 'B': //not A
            return ch!='A';
    }
    return false;
}

// evaluate if @ch_strand1 is complementary to @ch_strand2 according to @transf_matrix and IUPAC Notation
inline bool is_complemntary(char &ch_strand1, char &ch_strand2, string &transf_matrix){
    //ch_strand1=toupper(ch_strand1);
    //ch_strand2=toupper(ch_strand2);

    if(ch_strand1=='N' || ch_strand2=='N') return true;
    if(ch_strand1=='*' && ch_strand2=='*') return true; //if both are wildcards - ok
    if(ch_strand1=='*' || ch_strand2=='*') return false; //if only one is wildcard - bad

    switch (ch_strand1){
        case 'A': //Adenine
            return fits(ch_strand2, transf_matrix[0]);
        case 'C': //Cytosine
            return fits(ch_strand2, transf_matrix[1]);
        case 'G': //Guanine
            return fits(ch_strand2, transf_matrix[2]);
        //case 'U':
        case 'T': //Thymine in DNA; uracil in RNA
            return fits(ch_strand2, transf_matrix[3]);
    }
    return true;
}

// reverse and make complement of the given substring
inline void reverse_complement(string::iterator begin, string::iterator end) {
    reverse(begin, end);

    for(string::iterator i=begin; i!=end; ++i) {
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
            //case 'U':
            case 'T':
                *i = 'A';
                break;
            default :
                *i = 'N';
        }
    }
}

// filter out all whitespaces from the given string @s
inline void filter_whitespaces(string &s){
    string tmp;
    unsigned int last = 0;
    unsigned int i = 0;
    unsigned int size = s.size();
    tmp.reserve(size);
    
    while(i<size){
        while(i<size && !isspace(s[i])) ++i;
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
