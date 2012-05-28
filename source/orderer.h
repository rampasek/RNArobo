/*
 * $Id: orderer.h,v 1.9 2012-04-23 08:08:03 laci Exp $
 *
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

using namespace std;

class Orderer{
    public:
        unsigned int K;
        vector< struct TupleStats > tupleStats;   // statistic data for tried k-tuples of elements
        set<int> activeTuples;   // tuples that have not been rejected yet (set of indices to tupleStats)
        int currentKTupleID;     // index of currently selected k-tuple in tupleStats
        vector<int> fixedOrder;     // the ordered list of SSEs that are already fixed, it's a prefix of srchOrder
        vector<int> searchOrder;     // the ordered list of SSEs to be searched (indices to the vector "sses")
				pair<vector<double>, vector<double> > heuristicParams;
        bool doIterativeTraining;
        bool samplingSearchDone;
        
        //stringstream fout;
        unsigned long long scannedWindows;
        unsigned long long scannedBases;
        
        Orderer(){};
        Orderer(Descriptor &dsc, unsigned int k, unsigned int trainSL, unsigned int alphaID, bool doIT, pair<vector<double>, vector<double> > hParams);
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
