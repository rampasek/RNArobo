/*
 * $Id: matrix.h,v 1.5 2012-03-14 14:44:06 laci Exp $
 *
 * Project      : RNA motif searching in genomic sequences
 * Description  : the header file for matrix.cpp
 *
 * Author       : Ladislav Rampasek <rampasek@gmail.com>
 * Institution  : Comenius University in Bratislava
 *
 */

#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
//#include <cstdarg>
//#include <cstdio>
#include <set>
#include <tr1/array>

using namespace std;

//data structure for dynamic programming
class Matrix{
    public:
        Matrix():access_counter(0ULL){ Matrix(1); };
        Matrix(unsigned int num_dimensions):num_dimensions(num_dimensions), access_counter(0ULL){};

        //bool get2(int num_args, ...);
        bool get(tr1::array<int, 7> &key);
        //void set2(int num_args, ...);
        void set(tr1::array<int, 7> &key);
        void clear(){ matrix.clear(); };
        bool empty(){ return matrix.empty(); };
        void set_dimensions(unsigned int new_num_dimensions){ num_dimensions=new_num_dimensions; clear(); };
        unsigned long long int ops(){ return access_counter; }
        void resetOpsCounter(){ access_counter=0; }
    private:
        unsigned int num_dimensions;
        unsigned long long int access_counter;
        std::set < unsigned long long int > matrix;
};

#endif
