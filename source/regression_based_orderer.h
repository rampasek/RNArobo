/*
 * Project      : RNA motif searching in genomic sequences
 * Description  : the header file for regression_based_orderer.cpp
 *
 * Author       : Ladislav Rampasek <rampasek@gmail.com>, Martin Kralik <majak47@gmail.com>
 * Institution  : Comenius University in Bratislava
 *
 */

#ifndef REGRESSION_BASED_ORDERER_H
#define REGRESSION_BASED_ORDERER_H

#include <vector>
#include <set>
#include <map>
#include <sstream>

#include "descriptor.h"

using namespace std;

class RegressionBasedOrderer{
    public:
        unsigned int K;
        vector< struct TupleStats > tupleStats;   // statistic data for tried k-tuples of elements
        set<int> activeTuples;   // tuples that have not been rejected yet (set of indices to tupleStats)
        int currentKTupleID;     // index of currently selected k-tuple in tupleStats
        vector<int> fixedOrder;     // the ordered list of SSEs that are already fixed, it's a prefix of srchOrder
        vector<int> searchOrder;     // the ordered list of SSEs to be searched (indices to the vector "sses")
        bool doIterativeTraining;
        bool samplingSearchDone;
        vector<pair<pair<double, double>, double> > samples;
        
        //stringstream fout;
        unsigned long long scannedWindows;
        unsigned long long scannedBases;
        
        RegressionBasedOrderer(){};
        RegressionBasedOrderer(Descriptor &dsc, unsigned int k, unsigned int trainSL, unsigned int alphaID, bool doIT);
        void setNewSearchOrder(int seqSize);
        
    private:
        Descriptor *desc;
        int currentSeqSize;
        unsigned int trainSetLimit;
        double *criticalValues;
        
        void prepareKTuples();
        int getTupleScore(vector<int> &tuple);
        int calculateDomainFlexibility(vector<bool> &fixed, int element);
        bool storeKTupleStats(int tupleID);
        void eliminateAgainstKTuple(int tupleID);
        int sampleKTuple();
        void completeSearchOrder();
};

#endif
