/*
 * Project      : RNA motif searching in genomic sequences
 * Description  : implementation of Orderer class, which is to find a "good"
 *                element ordering for the backtrack in Simple_Search
 *
 * Author       : Ladislav Rampasek <rampasek@gmail.com>
 * Institution  : Comenius University in Bratislava
 *
 */

#include <vector>
#include <set>
#include <string>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include <iostream>
#include <iomanip>
#include "orderer.h"
#include "welchtest.h"

using namespace std;

Orderer::Orderer(Descriptor &dsc, unsigned int k, unsigned int trainSL, unsigned int alphaID, bool doIT, pair<vector<double>, vector<double> > hParams){
    srand ( time(NULL) );
    
    heuristicParams = hParams;
    desc = &dsc;
    K = min((int)k, (int)desc->sses.size() -1 -(int)desc->predef_srch_order.size());
    trainSetLimit = trainSL;
    doIterativeTraining = doIT;
    
    double *cVals[5] = { WT::criticalValuesAt80, WT::criticalValuesAt90, WT::criticalValuesAt95,
                         WT::criticalValuesAt975, WT::criticalValuesAt99
                       };
    if(alphaID >= 0 && alphaID < 5){
        criticalValues = cVals[alphaID];
    } else {
        criticalValues = cVals[2];
    }
    
    
    scannedWindows = 0;
    scannedBases = 0;
    currentKTupleID = -1;
    currentSeqSize = 0;
    samplingSearchDone = false;
    fixedOrder = desc->predef_srch_order;
    searchOrder = fixedOrder;
    completeSearchOrder();
    
    //fout<<"doIterativeTraining: "<<doIterativeTraining<<endl;
    //fout<<"signifLevel: 95"<<endl;
    
    prepareKTuples();
    //cerr<<"WT: "<<WT::test()<<endl<<endl;
}

///find and prepare good looking K-tuples for further empirical testing
void Orderer::prepareKTuples(){
    //find out which elements where not ordered by a user (those we have to order)
    vector<bool> has(desc->sses.size(), false);
    for(int i=0;i<fixedOrder.size();i++) has[fixedOrder[i]] = true;
  
    vector<int> elementsToOrder;
    elementsToOrder.reserve(desc->sses.size() - fixedOrder.size());
    for(int i=1;i<desc->sses.size();i++){
        if(!has[i]){
            elementsToOrder.push_back( i );
        }
    }

    //generate (by brute-force) and score all K-tuples
    vector< pair<int, vector<int> > > allTuples;
    allTuples.reserve( (unsigned int) pow((double)elementsToOrder.size(), (int)K) );
    vector< vector<int>::iterator > iters(K, elementsToOrder.begin());

    vector<int> tuple(K);
    int min_tscore = 1;
    while(true){
        //check whether the current configuration has no repetitions
        bool rep = false;
        set<int> lock;
        for(int i=0; i<iters.size(); ++i){
            if( lock.count(*iters[i]) > 0 ) {
                rep = true;
                break;
            }
            lock.insert(*iters[i]);
        }
            
        //if we have a new k-variation of elements, then score it and save it
        if(!rep){
            for(int i=0; i<iters.size(); ++i){
                tuple[i] = *iters[i];
                //cout<<tuple[i]<<" ";
            } //cout<<endl;

            //score the tuple and save it
            int tscore = getTupleScore(tuple);
            min_tscore = min(min_tscore, tscore);
            allTuples.push_back( make_pair(tscore, tuple) );
        }
        
        //move iterators to the next configuration
        int iiter = 0;
        while(iiter < iters.size()){
            ++iters[iiter];
            if(iters[iiter] == elementsToOrder.end()){
                iters[iiter] = elementsToOrder.begin();
                ++iiter;
            } else {
                break;
            }
        }
        //if there isn't next configuration, then end the generation
        if(iiter >= iters.size()){
            break;
        }
    }
    
    //to avoid problems with negative and zero scores, raise all scores above 0
    if(min_tscore <= 0){
        for(int i=0; i<allTuples.size(); ++i) allTuples[i].first += (-min_tscore)+1;
    }
    
    //sort the K-tuples according to their score and take only those "best looking"
    sort(allTuples.begin(), allTuples.end());

    int bestScore = allTuples.back().first;
    int scoreLimit = (int)(bestScore*.85); //TODO: change to 0.8
    
    struct TupleStats tmpStats;
    tmpStats.heuristicScore = 0;
    tmpStats.memOps.resize(K);
    tmpStats.basesScanned = 0;
    for(int i = allTuples.size()-1; i>-1 && allTuples[i].first >= scoreLimit && tupleStats.size() < trainSetLimit; i--){
        tmpStats.tuple = allTuples[i].second;
        tmpStats.heuristicScore = allTuples[i].first;
        
        //store to object variables
        activeTuples.insert( tupleStats.size() );
        tupleStats.push_back( tmpStats );
    }

    //debug log - print all tuples that made it to the list, and their scores
    /*for(int tt = 0; tt < tupleStats.size(); ++tt){
        cout<<"TupleID: "<<tt<<endl;
        cout<<"  the tuple: ";
        for(int i=0; i<tupleStats[tt].tuple.size(); ++i){
            cout<<tupleStats[tt].tuple[i]<<" ";
        }
        cout<<"  score = "<<tupleStats[tt].heuristicScore<<endl;
    }*/
    
}

///calculate heuristic score for given tuple of elements
///since part of the score is order dependant, fixedOrder is considered
///to be the beginning of the corresponding partial search order
int Orderer::getTupleScore(vector<int> &tuple){
    //cout<<"size: "<<tuple.size()<<" "<<fixedOrder.size()<<endl;
    //TODO: kombinacia HF do skore
    double tupleScore = 0;
        
    vector<bool> alreadyFixed(desc->sses.size(), false);
    for(int i=0; i<fixedOrder.size(); ++i){
        alreadyFixed[fixedOrder[i]] = true;
    }
    
    for(int i=0; i<tuple.size(); i++){
        //get the information content of the sse
        double ic = desc->sses[tuple[i]].infContent;
        
        //calculate approximate flexibility of the sse's search domain size
        int domainFlexibility = calculateDomainFlexibility(alreadyFixed, tuple[i]);
        if( desc->sses[tuple[i]].is_helix ){
            int secondStrandFlex = calculateDomainFlexibility(alreadyFixed, -tuple[i]);
            domainFlexibility = domainFlexibility * secondStrandFlex;
        }
        
        alreadyFixed[tuple[i]] = true;
        
        //calculate heuristic score for this sse and add it to the tuple's score using heuristicParams
        double elementScore = heuristicParams.first[i]*ic;
        elementScore += heuristicParams.second[i]*domainFlexibility;
        
        //weight scaling should be included in heuristicParams
				tupleScore += elementScore;
        
       //cout<<tuple[i]<<">  ic= "<<ic<<"\n    apxDF= "<<domainFlexibility<<endl;
    }
    //cout<<"TS: "<<tupleScore<<endl;
    return (int)round(tupleScore);
}

///calculate approximate flexibility of the sse's search domain size
///(if helix, then only of the first strand, "-element" means second strand)
int Orderer::calculateDomainFlexibility(vector<bool> &fixed, int element){
    int approxDFleft = -1;
    int approxDFright = -1;
    int tmpDFleft = 0;
    int tmpDFright = 0;
    bool wasSeen = false;
    for(int j=0; j<desc->motif.size(); ++j){
        int sseID = abs(desc->motif[j]);
        //if this motif element is fixed by the time sseID is search for (or it's the other strand)
        if(fixed[sseID] || (desc->sses[abs(element)].is_helix && desc->motif[j] == -element)){ 
            if(wasSeen){ //if @element has been seen already => we are at the end of the right side of its domain 
                approxDFright = tmpDFright;
                break;
            } else {
                tmpDFleft = 0;
            }
        } else { //else => this element is not fixed, so add in its flexibility
            if(wasSeen){ //right side of the domain
                tmpDFright += desc->sses[sseID].size_range.second - desc->sses[sseID].size_range.first;
            } else { //left side of the domain
                tmpDFleft += desc->sses[sseID].size_range.second - desc->sses[sseID].size_range.first;
            }
        }
        
        if(desc->motif[j] == element){
            wasSeen = true;
            approxDFleft = tmpDFleft;
            tmpDFright += desc->sses[sseID].size_range.second - desc->sses[sseID].size_range.first;
        }
    }
    if(tmpDFright == -1){
        approxDFright = tmpDFright;
    }
    
    return max(1, min(approxDFleft, approxDFright)+1);
}


///process the current k-tuple, then choose a new one + augment it to a complete ordering
void Orderer::setNewSearchOrder(int seqSize){
    //if we already have the best k-tuple then exit else continue to train the following k-tuples
    if(samplingSearchDone) return;
    
    ++scannedWindows;
    scannedBases += currentSeqSize;   
    
    //all but one k-tuple has been eliminated => we have the winner
    if(activeTuples.size() == 1) {
        ///FOUT begin
        /*fout<<"K: "<<K<<endl;
        fout<<"tuples sampled means (memOPs per window):\n";
        for(int i=0;i<tupleStats.size();++i){
            double sum=0;
            for(int j=0;j<tupleStats[i].sampledOpsPerWindow.size(); ++j){
                sum += tupleStats[i].sampledOpsPerWindow[j];
            }
            sum = sum/tupleStats[i].sampledOpsPerWindow.size();
            fout.setf(ios::floatfield); 
            fout<<setprecision(5)<<sum<<" ";
        } fout<<endl;
        
        fout<<"number of samples per tuple: ("<<tupleStats.size()<<")\n";
        for(int i=0;i<tupleStats.size();++i){
            fout<<tupleStats[i].sampledOpsPerWindow.size()<<" ";
        } fout<<endl;
        
        fout<<"WINNER: "<<*activeTuples.begin()<<endl;
        fout<<"COUNT: "<<scannedWindows<<"\nAVG.WIN: "<<scannedBases/(double)scannedWindows<<endl;
        */
        ///FOUT end
        
        //add the victorious k-tuple to fixedOrder
        if(currentKTupleID > -1){
            fixedOrder.insert(fixedOrder.end(),
                            tupleStats[currentKTupleID].tuple.begin(),
                            tupleStats[currentKTupleID].tuple.end()
                            );
        }
        
        //start training next K elements following in search ordering (if there are still some more to train)
        if(doIterativeTraining && fixedOrder.size() < searchOrder.size()-1){
            K = min((int)K, (int)desc->sses.size() -1 -(int)fixedOrder.size());
            
            currentKTupleID = -1;
            activeTuples.clear();
            tupleStats.clear();
            
            prepareKTuples();
        } else {
            samplingSearchDone = true;
            //cout<<"\nSAMPLING SEARCH DONE: "<<scannedWindows<<"\n\n";
            return;
        }
    }
    
    //store gathered info for the current k-tuple and run an elimination round
    bool didGatherSample = storeKTupleStats(currentKTupleID);
    if(didGatherSample){ eliminateAgainstKTuple(currentKTupleID); }

    //reset Ops counters for all elements
    for(int i=1; i<desc->sses.size(); i++){
        desc->sses[i].table.resetOpsCounter();
    }

    //test sampling
    /*ofstream ffout ("data.txt", ios::out | ios::app);
    for(int i=0;i<50000;++i){
        ffout<<sampleKTuple()<<", ";
    }
    ffout<<endl;
    ffout.close();
    */

    //sample another k-tuple and set a new search ordering
    searchOrder = fixedOrder;
    currentKTupleID = sampleKTuple();
    vector<int> nextTuple = tupleStats[currentKTupleID].tuple;
    searchOrder.insert(searchOrder.end(), nextTuple.begin(), nextTuple.end());
    
    completeSearchOrder();
    currentSeqSize = seqSize;
}

///store the gathered data for the current K-tuple (used memOps)
bool Orderer::storeKTupleStats(int tupleID){
    if (tupleID == -1 ) return false;
    
    unsigned long long allOps = 0;
    for(int i=0; i < tupleStats[tupleID].tuple.size(); ++i){
        //tupleStats[tupleID].memOps[i] += desc->sses[ tupleStats[tupleID].tuple[i] ].table.ops();
        allOps += desc->sses[ tupleStats[tupleID].tuple[i] ].table.ops();
    }
    
    //all ops for this tuple + the rest of the search order
    /*for(int i=fixedOrder.size(); i < searchOrder.size(); ++i){
        allOps += desc->sses[ searchOrder[i] ].table.ops();
    }*/

    tupleStats[tupleID].basesScanned += currentSeqSize;
    
    //don't add if the search didn't make it so far in the search order (the tuple wasn't searched)
    if(allOps <= 0){
        return false;
    }
    
    tupleStats[tupleID].sampledOpsPerWindow.push_back(allOps/double(currentSeqSize));
    
    return true;
}

///eliminate all k-tuples whose mean is significantly larger than that of the given k-tuple
///(or eliminate the given k-tuple if another k-tuple has a significantly smaller mean)
void Orderer::eliminateAgainstKTuple(int tupleID){
    if (tupleID == -1 ) return;
    //vector<int> tuple = tupleStats[tupleID].tuple;
    
    int minMeasurementsCount = 2;
    //if we have at least some measurements with this k-tuple
    if(tupleStats[tupleID].sampledOpsPerWindow.size() >= minMeasurementsCount){
        vector<int> eliminated;
        bool eliminateSelf = false;
        
        //loop through all active tuples and try to eliminate some
        for(set<int>::iterator it=activeTuples.begin(); it!=activeTuples.end(); ++it){
            if(*it == tupleID) continue;
            if(tupleStats[*it].sampledOpsPerWindow.size() >= minMeasurementsCount){
                //if the tupleID k-tuple is stat. sigificantly "faster" (smaller mean val.)
                // than *it k-tuple, then eliminate the *it k-tuple
                int testRes = WT::welchTest(tupleStats[*it].sampledOpsPerWindow,
                                 tupleStats[tupleID].sampledOpsPerWindow,
                                 criticalValues
                                );
                
                //decide about elimination
                if(testRes == 1){ //*it has higher mean
                    eliminated.push_back(*it);
                } else if(!eliminateSelf && testRes == -1){ //tupleID has higher mean
                    eliminated.push_back(tupleID);
                    eliminateSelf = true;
                } else if(!eliminateSelf &&
                         tupleStats[*it].sampledOpsPerWindow.size() > 100 &&
                         tupleStats[tupleID].sampledOpsPerWindow.size() > 100
                        ) //break tie between these tuples
                {
                    double sum1=0;
                    for(int j=0;j<tupleStats[tupleID].sampledOpsPerWindow.size(); ++j){
                        sum1 += tupleStats[tupleID].sampledOpsPerWindow[j];
                    }
                    sum1 = sum1/tupleStats[tupleID].sampledOpsPerWindow.size();
                    
                    double sum2=0;
                    for(int j=0;j<tupleStats[*it].sampledOpsPerWindow.size(); ++j){
                        sum2 += tupleStats[*it].sampledOpsPerWindow[j];
                    }
                    sum2 = sum2/tupleStats[*it].sampledOpsPerWindow.size();
                    
                    int el = (sum1 < sum2) ? *it : tupleID;
                    eliminated.push_back(el);
                    eliminateSelf = (el==tupleID) ? true : false;
                }
            }
        }
        
        //drop the eliminated k-tuples from the active set
        for(int i=0;i<eliminated.size();++i){
            activeTuples.erase(eliminated[i]);
            
            /*cout<<tupleID<<" ELIMINTATES: "<<eliminated[i]<<" after: "<<
                tupleStats[tupleID].sampledOpsPerWindow.size()<<" vs "<<
                tupleStats[eliminated[i]].sampledOpsPerWindow.size()<<endl;
            */
        }
    }

}

///sample another K-tuple from the list, return its index in tupleStats
int Orderer::sampleKTuple(){
    //sample from the K-tuples according to distribution of their score
    /*int minimum = 1;
    int scoreSum = 0;
    for(set<int>::iterator it=activeTuples.begin(); it!=activeTuples.end(); ++it){
        scoreSum += tupleStats[*it].score;
        if(tupleStats[*it].score < minimum){
            minimum = tupleStats[*it].score;
        }
    }
    */
    
    //prefer the K-tuple with the highest heuristic score
    //if(activeTuples.count(0)>0 && tupleStats[0].sampledOpsPerWindow.size() < 3) return 0;
    
    //take an i.i.d. sample
    int scoreSum = activeTuples.size()*10;
    int randNum = rand()%scoreSum +1;

    int sample = -1;
    for(set<int>::iterator it=activeTuples.begin(); it!=activeTuples.end(); ++it){
        randNum -= 10;
        if(randNum <= 0) {
            sample = *it;
            break;
        }
    }
    //cout<<"SAMPLE: "<<sample<<"   n"<<randNum<<"   s"<<scoreSum<<endl;
    return sample;
}

///augument the current search ordering to a complete one
void Orderer::completeSearchOrder(){
    vector<bool> has(desc->sses.size(),false);
    vector< pair<double, int> > auto_order;
    for(int i=0;i<searchOrder.size();i++) has[searchOrder[i]] = true;

    searchOrder.reserve(desc->sses.size()-1);

    //use information content heuristic to order elements according to their specificity
    auto_order.reserve(desc->sses.size() - searchOrder.size());
    for(int i=1;i<desc->sses.size();i++){
        if(!has[i]){
            double score = desc->sses[i].infContent;
            auto_order.push_back( make_pair(score, i) );
        }
    }
    
    sort(auto_order.begin(), auto_order.end());
    for(int i=auto_order.size()-1; i>-1; i--){
       searchOrder.push_back(auto_order[i].second);
    }
}
