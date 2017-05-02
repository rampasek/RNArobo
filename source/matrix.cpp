/*
 * Project      : RNA motif searching in genomic sequences
 * Description  : the data structure for dynamic programming
 *
 * Author       : Ladislav Rampasek <rampasek@gmail.com>
 * Institution  : Comenius University in Bratislava
 *
 */

#include <vector>
#include <set>
#include <array>
#include <cstring>

#include "matrix.h"

using namespace std;

//get value of a cell specified by given key of type vector<int>
bool Matrix::get(array<unsigned int, 7> &key){
    //if(key.size()!=num_dimensions) return 0;
    
    unsigned long long int x = 0;
    for(int i=0;i<3;i++){
        x <<= 15;
        x |= key[i];
    }
    for(int i=3;i<6;i++){
        x <<= 5;
        x |= key[i];
    }
    x <<= 1;
    x |= key[6];

    return (matrix.find(x) != matrix.end());
}

//set value of a cell specified by given key of type vector<int>
void Matrix::set(array<unsigned int, 7> &key){
    //if(key.size()!=num_dimensions) return;
    
    unsigned long long int x = 0;
    for(int i=0;i<3;i++){
        x <<= 15;
        x |= key[i];
    }
    for(int i=3;i<6;i++){
        x <<= 5;
        x |= key[i];
    }
    x <<= 1;
    x |= key[6];
    
    matrix.insert(x);
}

/*
//get value of a cell specified by arguments passed by "..."
bool Matrix::get2(int num_args, ...){
    //if(num_args!=num_dimensions) return 0;

    vector<int> key(num_dimensions);

    //get arguments
    va_list argp;
    va_start(argp,num_args);
        for(int i=0; i<num_args; ++i) {
            key[i]=va_arg(argp, int);
        }
    va_end(argp);

    return (matrix.find(key) != matrix.end());
}*/

/*
//set value of a cell specified by arguments passed by "..."
void Matrix::set2(int num_args, ...){
    //if(num_args!=num_dimensions) return;

    vector<int> key(num_dimensions);

    //get arguments
    va_list argp;
    va_start(argp,num_args);
        for(int i=0; i<num_args; ++i) {
            key[i]=va_arg(argp, int);
        }
    va_end(argp);

    matrix.insert(key);
}
*/
