/*
 * Project      : RNA motif searching in genomic sequences
 * Description  : the header file for orderer.cpp
 *
 * Author       : Ladislav Rampasek <rampasek@gmail.com>
 * Institution  : Comenius University in Bratislava
 *
 */

#ifndef ORDERER_H
#define ORDERER_H

#include <vector>
#include <set>
#include <map>
#include <sstream>

#include "descriptor.h"

#define _IC_H  0
#define _IC_SS 1
#define _DF    2

using namespace std;

class Orderer{
    public:
        unsigned int K;
        vector< struct TupleStats > tupleStats;   // statistic data for tried k-tuples of elements
        set<int> activeTuples;   // tuples that have not been rejected yet (set of indices to tupleStats)
        int currentKTupleID;     // index of currently selected k-tuple in tupleStats
        vector<int> fixedOrder;     // the ordered list of SSEs that are already fixed, it's a prefix of srchOrder
        vector<int> searchOrder;     // the ordered list of SSEs to be searched (indices to the vector "sses")
        vector< vector<double> > heuristicParams; //coefficients for combining heuristics to score
        bool doIterativeTraining;
        
        //stringstream fout;
        unsigned long long scannedWindows;
        unsigned long long scannedBases;
        
        Orderer(){};
        Orderer(Descriptor &dsc, unsigned int k, unsigned int trainSL, int alphaID, bool doIT, vector< vector<double> > hParams);
        void setNewSearchOrder(int seqSize);
        bool isSearchDone(){return samplingSearchDone;};
        
    private:
        Descriptor *desc;
        int currentSeqSize;
        unsigned int trainSetLimit;
        double *criticalValues;
        bool samplingSearchDone;
        
        void prepareKTuples();
        int getTupleScore(vector<int> &tuple);
        int calculateDomainFlexibility(vector<bool> &fixed, int element);
        bool storeKTupleStats(int tupleID);
        void eliminateAgainstKTuple(int tupleID);
        int sampleKTuple();
        void completeSearchOrder();
};

#endif
