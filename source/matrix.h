/*
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
#include <set>
#include <tr1/array>

using namespace std;

//data structure for dynamic programming
class Matrix{
    public:
        Matrix(){ Matrix(1); };
        Matrix(unsigned int num_dimensions) : num_dimensions(num_dimensions) {};

        bool get(tr1::array<unsigned int, 7> &key);
        void set(tr1::array<unsigned int, 7> &key);
        void clear(){ matrix.clear(); };
        bool empty(){ return matrix.empty(); };
        void set_dimensions(unsigned int new_num_dimensions){ num_dimensions=new_num_dimensions; clear(); };
    private:
        unsigned int num_dimensions;
        std::set<unsigned long long> matrix;
};

#endif
