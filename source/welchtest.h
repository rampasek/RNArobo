/*
 * Project      : RNA motif searching in genomic sequences
 * Description  : Welch Two Sample t-test
 *                (performs the upper, one-sided test)
 *
 * Author       : Ladislav Rampasek <rampasek@gmail.com>
 * Institution  : Comenius University in Bratislava
 *
 */

#ifndef WELCHTEST_H
#define WELCHTEST_H

#include <vector>

using namespace std;

namespace WT {

    
extern double criticalValuesAt80[102];
extern double criticalValuesAt90[102];
extern double criticalValuesAt95[102];
extern double criticalValuesAt975[102];
extern double criticalValuesAt99[102];

int welchTest(vector<double> &x, vector<double> &y, double* criticalValues);
bool test(void);

}
#endif
