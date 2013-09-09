/*
 * Project      : RNA motif searching in genomic sequences
 * Description  : the main program
 *
 * Author       : Ladislav Rampasek <rampasek@gmail.com>
 * Institution  : Comenius University in Bratislava
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <cstring>
#include <cmath>
#include <getopt.h>

#include "search.h"

using namespace std;
using namespace GF;

#ifndef RELEASE
    #define DEF_RELEASE     "2.1.0"
    #define DEF_RELEASEDATE "September 2013"
#endif


static string txt_banner = "RNARobo - RNA motif searching program";
static string txt_usage =
"\nrnarobo: version "+string(DEF_RELEASE)+", "+string(DEF_RELEASEDATE)+"\n\
Usage: rnarobo [OPTIONS] <descriptor-file> <sequence-file>\n\
\n\
  Available options: \n\
     -c              search both strands of database\n\
     -u              report only non-overlapping occurrences\n\
     -f              print output in plain FASTA format\n\
     -s              print output in FASTA format with element separators\n\
     --nratio FLOAT  set max allowed ratio of \"N\"s in occurrences to report\n\
     \n\
  To override default search order training parameters:\n\
  (defaults: --k 3 --limit 50 --alpha 0.01 --iterative TRUE )\n\
    --k INT            set length of tuples used in training\n\
    --limit INT        set max size of training set (max num of tuples)\n\
    --alpha FLOAT      set significance level for statistic tests,\n\
                         must be one of: 0.2, 0.1, 0.05, 0.025, 0.01 \n\
    --iterative BOOL   iteratively train whole the ordering TRUE / FALSE\n\
    --tonly            perform only the order training itself\n";
//     -d:    show dot-bracked represenation for every match (!!!CAN'T HANDLE INSERTIONS YET!!!)\n

inline void print_formated_output(bool fasta, unsigned int begin, unsigned int end, string name, string details, string seq){
    if (fasta == true) {
        printf(">%u-%u_%s %s\n%s\n", begin, end, name.c_str(), details.c_str(), seq.c_str());
    } else {
        printf("%6u %6u %12s %12s\n%s\n", begin, end, name.c_str(), details.c_str(), seq.c_str());
    }
}

int main(int argc, char* argv[]){
    clock_t tStart = clock(); //to measure runtime
    
    string file_db;       // name of sequence file to search
    string file_desc;     // name of file containing descriptor
    ifstream dbin, descin;

    string element_separator = "|"; // the string used to separate individual elements in the output

  /***********************************************
   * Parse command line arguments
   ***********************************************/
    if (argc < 3) die("Not enough parameters\n" + txt_usage);
    //for(int i=0;i<argc;++i) cout<<argv[i]<<endl;
    
    int opt_be_quiet = false;       /* OPTION: TRUE to silence verbosity*/    
    int opt_searchcomp = false;     /* OPTION: TRUE to search both strands of database*/
    int opt_uniq = false;           /* OPTION: TRUE report only non-overlapping motif occurrences*/
    int opt_fasta = false;          /* OPTION: TRUE to print output in FASTA format*/
    int opt_dotbracked = false;     /* OPTION: TRUE to print also dot-bracked representation*/
    float opt_max_nratio = 1.;      // max allowed ratio of "N"s in reported occurrences to its length
    
    int opt_tonly = false;   /* OPTION: TRUE to perform only the training of search ordering */
    int opt_k = -1;
    int opt_limit = -1;
    int opt_alpha = -1;
    int opt_iterative = -1;
    string opt_paramsfile = "";
    const char *alphas[] = {"0.2", "0.1", "0.05", "0.025", "0.01"};
    
    static struct option long_options[] = {
        {"k", 1, 0, 'k'},
        {"limit", 1, 0, 'l'},
        {"alpha", 1, 0, 'a'},
        {"iterative", 1, 0, 'i'},
        {"tonly", 0, 0, 't'},
        {"params", 1, 0, 'p'},
        {"complement", 0, 0, 'c'}, //long alias for -c
        {"uniq", 0, 0, 'u'}, //long alias for -u
        {"fasta", 0, 0, 'f'}, //long alias for -f
        {"sep", 0, 0, 's'}, //long alias for -s
        {"nratio", 1, 0, 'n'},
        {NULL, 0, NULL, 0}
    };
    extern int optind;
    int c, option_index = 0;
    while ((c = getopt_long(argc, argv, "cufs", long_options, &option_index)) != -1) {
        string par;
        
        switch(c) {
        case 0:
            //printf ("option %s %d", long_options[option_index].name, option_index);
            //if(optarg){ printf (" with arg %s", optarg); }
            //printf ("\n");
            break;
        case 'k':
            opt_k = atoi(optarg);
            if(opt_k <= 0) die("Illegal value of parameter 'k' : must be > 0.");
            break;
        case 'l':
            opt_limit = atoi(optarg);
            if(opt_limit <= 0) die("Illegal value of parameter 'limit': must be > 0.");
            break;
        case 'a':
            opt_alpha = (int)round(atof(optarg)*1000);
            switch(opt_alpha){
            case 200:
                opt_alpha = 0;
                break;
            case 100:
                opt_alpha = 1;
                break;
            case 50:
                opt_alpha = 2;
                break;
            case 25:
                opt_alpha = 3;
                break;
            case 10:
                opt_alpha = 4;
                break;
            default:
                die("Illegal 'alpha': "+string(optarg));
            }
            break;
        case 'i':
            par.assign(optarg);
            transform(par.begin(), par.end(), par.begin(), ::toupper);
            if(par == "TRUE"){
                opt_iterative = true;
            } else if(par == "FALSE"){
                opt_iterative = false;
            } else {
                opt_iterative = (bool)atoi(optarg);
            }
            break;
        case 't':
            opt_tonly = true;
            opt_searchcomp = true;
            break;
        case 'c':
            opt_searchcomp = true;
            break;
        case 'u':
            opt_uniq = true;
            break;
        case 'f':
            opt_fasta = true;
            opt_be_quiet = true;
            element_separator = "";
            break;
        case 's':
            opt_fasta = true;
            opt_be_quiet = true;
            element_separator = "|";
            break;
        case 'n':
            opt_max_nratio = max(0., min(1., atof(optarg)));
            break;
        case 'd':
            //opt_dotbracked = true;
            break;
        case 'p':
            opt_paramsfile = optarg;
            break;
        case '?':
            char opt[20];
            sprintf(opt, "Unknown arg \"-%c\"\n", optopt);
            die(string(opt) + txt_usage);
            break;
        }
    }
    
    //default parameters for Orderer
    unsigned int param_k = 3;
    unsigned int param_limit = opt_tonly?60:50;
    unsigned int param_alpha = opt_tonly?4:4;
    bool param_iterative = true;
    vector< vector<double> > params(3);
    double IC_H = 1;
    double IC_SS = 3;
    double DF = -0.2;
    
    //override default parameters by the user-defined
    if(opt_k != -1) param_k = opt_k;
    if(opt_limit != -1) param_limit = opt_limit;
    if(opt_alpha != -1) param_alpha = opt_alpha;
    if(opt_iterative != -1) param_iterative = opt_iterative;
    
    //set coeficients for mixing IC and DF heuristic into a score
    for(int i=param_k-1; i>=0; --i){
        params[_IC_H].push_back(pow(2.,i) * IC_H);
        params[_IC_SS].push_back(pow(2.,i) * IC_SS);
        params[_DF].push_back(pow(2.,i) * DF);
    }
    if(opt_paramsfile != "") {
        ifstream pfile(opt_paramsfile.c_str(), ifstream::in);
        for (int i=0; i<param_k; i++){
            pfile >> params[_IC_H][i];
            pfile >> params[_IC_SS][i];
            pfile >> params[_DF][i];
        }
    }
    /*cout<<optind<<endl;
    for ( ; optind < argc; optind++) {
            if (argv[optind]) {
                cout<<"mame "<<argv[optind]<<endl;
            }
    }*/
    if (argv[optind]){
        file_desc = argv[optind];
    }
    if (argv[optind+1]) {
        file_db = argv[optind+1];
    }
    //cout<<opt_searchcomp<<" "<<file_desc<<' '<<file_db<<endl;

  /***********************************************
   * Initialize sequence file, read descriptor
   ***********************************************/
    /// read descriptor
    descin.open( file_desc.c_str() );
    if(! descin.good()) die("Descriptor file "+file_desc+" could not be opened or it does not exist");
    Descriptor desc(descin);
    if(! desc.is_initialized()) die(desc.error_str());

    /// open sequence file
    dbin.open( file_db.c_str() );
    if(! dbin.good()) die("Sequence file "+file_db+" could not be opened or it does not exist");

  /***********************************************
   * Print header and do the search
   ***********************************************/
  
    if (! opt_be_quiet){
        printf("Starting rnarobo: version %s, %s\n", DEF_RELEASE, DEF_RELEASEDATE);
        printf("------------------------------------------------------------------\n");
        printf("Database file:                  %s\n", file_db.c_str());
        printf("Descriptor file:                %s\n", file_desc.c_str());
        printf("Complementary strand searched:  %s\n", (opt_searchcomp? "yes":"no"));
        printf("User-predefined search order:   %s\n",
               (desc.predef_srch_order.size()==0
               ? "none" 
               : desc.search_order_to_str(desc.predef_srch_order).c_str()) );
        printf("Order training params:          k=%u limit=%u alpha=%s iter=%s\n",
               param_k, param_limit, alphas[param_alpha], (param_iterative? "true":"false") );
        printf("------------------------------------------------------------------\n");
        
        if(opt_tonly){
            printf("\nTRAINING SEARCH ORDER ONLY\n");
        } else {
            printf("\n seq-f  seq-t     name      description\n");
            printf(  "------ ------ ------------ ------------\n");
        }
    }

    Orderer orderer(desc, param_k, param_limit, param_alpha, param_iterative, params);
    Simple_Search ssearch(desc, orderer);
    string line;
    long long total_bases_scanned = 0;
    long long total_matches = 0;
    long long reported_matches = 0;
    
    //long long evalBases = 0;
    //unsigned int evalWindows = 0;
    //double opsPerBase = -1;
    bool doSearch = true;
    
    while( get_valuable_line(dbin,line) && doSearch){
        if(line[0] != '>') die("Incorrect sequence format");
        line[0]=' ';

        string sq_name, sq_details, sq;
        istringstream sin(line);
        if( !(sin>>sq_name)) die("Incorrect sequence format");
        while( sin>>line ) sq_details+=" "+line;

        while('A'<=dbin.peek() && dbin.peek()<='z') {
            getline(dbin,line);
            normalize_seq(line.begin(), line.end());
            sq+=line;
        }
        filter_whitespaces(sq);
        //cout<<sq_name<<":"<<sq<<endl<<endl;
            
        int max_motif_length = ssearch.desc->get_max_motif_length();
        //to how long pieces we will chop up the sequence
        int max_seq_length = max(3000, 20*max_motif_length);
        if(ssearch.orderer->isSearchDone()) max_seq_length = max(10000, 50*max_motif_length);
        
        //to store beginnings of found matches (in both strands) - for filtering repeating matches
        set <pair <unsigned int, unsigned int> > found_matches, found_op_matches;
        found_matches.clear();
        found_op_matches.clear();
        //position of previous matches in the + and - strands
        pair <unsigned int, unsigned int> prev_match, prev_op_match; 

        //cut out "N"regions longer than 10 -> the sequence is divided into blocks
        vector< pair<string, unsigned int> > seq_blocks; //pair<block's sequence, block's beginning position in original seqence>
        seq_blocks.reserve(4);
        
        int prev_found = 0;
        char *p = NULL;
        p = strstr(const_cast<char*>(sq.c_str()), "NNNNNNNNNN");
        while(p){
            int found_pos = p - sq.c_str();
            seq_blocks.push_back(make_pair(sq.substr(prev_found, found_pos-prev_found+1), prev_found));
            while(found_pos+10 < sq.size() && sq[found_pos+10] == 'N'){
                ++found_pos;
                ++p;
            }
            prev_found = found_pos;
            p = strstr(p+1, "NNNNNNNNNN");
        }
        
        seq_blocks.push_back(make_pair(sq.substr(prev_found, sq.size()), prev_found));
        
        for(int k=0; k<seq_blocks.size() && doSearch; k++){
            //count number of bases scanned
            total_bases_scanned += seq_blocks[k].first.size();
            if(opt_searchcomp){ total_bases_scanned += seq_blocks[k].first.size(); }
            
            //chop up the given block to partitions of length at most max_seq_length
            // ! two partitions must have overlap of max_motif_length !
            int pointer = 0;
            vector< pair<string, unsigned int> > seq_partitions;
            seq_partitions.reserve(sq.size() / max_seq_length +1);
            
            while(pointer < seq_blocks[k].first.size()){
                int start = max(0, pointer-max_motif_length+1);

                string s = seq_blocks[k].first.substr(start, max_seq_length);
                seq_partitions.push_back(make_pair(s, start));
                pointer = start + max_seq_length;
            }

            //search for the motif in partitions
            for(int j=0; j<seq_partitions.size() && doSearch; j++){
                ssearch.search(seq_partitions[j].first);

                unsigned int offset = seq_partitions[j].second + seq_blocks[k].second;

                //print found hits
                for(int i=0;i<ssearch.solutions.size() && !opt_tonly;i++){
                    unsigned int match_begin = ssearch.solutions[i][0].first +1 +offset;
                    unsigned int match_end = ssearch.solutions[i].back().second +offset;
                    pair<unsigned int, unsigned int> match_pos = make_pair(match_begin, match_end);
                    
                    //filter duplicates (serves also as counter for number of found occurrences)
                    if(found_matches.count(match_pos) != 0){
                        continue;
                    }
                    found_matches.insert(match_pos);
                    
                    //filter out overlapping matches if the option is set
                    if(opt_uniq){
                        //this match must not have any overlap with the previous one
                        if( (prev_match.first <= match_begin && match_begin <= prev_match.second) ||
                            (prev_match.first <= match_end && match_end <= prev_match.second) ||
                            (match_begin <= prev_match.first && prev_match.second <= match_end)
                        ){
                            continue;
                        }
                        //cout<<prev_match.first<<" "<<prev_match.second<<" vs "<<match_begin<<" "<<match_end<<endl;
                    }
                    
                    //filter by the number of Ns ratio
                    int num_ns = 0;
                    for(int x=ssearch.solutions[i][0].first; x<ssearch.solutions[i].back().second; ++x){
                        if(seq_partitions[j].first[x] == 'N') ++num_ns;
                    }
                    if(float(num_ns)/(match_end-match_begin+1) > opt_max_nratio) continue;
                    
                    //report the found motif occurrence
                    prev_match = match_pos;
                    ++reported_matches;
                    print_formated_output(opt_fasta,
                                        match_begin,
                                        match_end,
                                        sq_name,
                                        sq_details,
                                        ssearch.solution_to_str(i, seq_partitions[j].first, element_separator)
                                    );
                    if (opt_dotbracked == true) {
                        printf("%s\n", ssearch.solution_to_dotbracket(i, element_separator).c_str());
                    }
                    
                }

                //search the partition also in the opposite direction (if opt_searchcomp == true)
                if(opt_searchcomp){
                    reverse_complement(seq_partitions[j].first.begin(), seq_partitions[j].first.end());
                    ssearch.search(seq_partitions[j].first);

                    //print found hits
                    for(int i=0;i<ssearch.solutions.size() && !opt_tonly;i++){
                        unsigned int match_begin = seq_partitions[j].first.size() -ssearch.solutions[i][0].first + offset;
                        unsigned int match_end = seq_partitions[j].first.size() -ssearch.solutions[i].back().second +1 +offset;
                        pair<unsigned int, unsigned int> match_pos = make_pair(match_begin, match_end);
                        
                        //filter duplicates (serves also as counter for number of found occurrences)
                        if(found_op_matches.count(match_pos) != 0){
                            continue;
                        }
                        found_op_matches.insert(match_pos);
                        
                        //filter out overlapping matches if the option is set
                        if(opt_uniq){
                            //this match must not have any overlap with the previous one
                            if( (prev_op_match.second <= match_begin && match_begin <= prev_op_match.first) ||
                                (prev_op_match.second <= match_end && match_end <= prev_op_match.first) ||
                                (match_begin >= prev_op_match.first && prev_op_match.second >= match_end)
                            ){
                                continue;
                            }
                            //cout<<prev_op_match.first<<" "<<prev_op_match.second<<" 'vs' "<<match_begin<<" "<<match_end<<endl;
                        }
                        
                        //filter by the number of Ns ratio
                        int num_ns = 0;
                        for(int x=ssearch.solutions[i][0].first; x<ssearch.solutions[i].back().second; ++x){
                            if(seq_partitions[j].first[x] == 'N') ++num_ns;
                        }
                        if(float(num_ns)/(match_begin-match_end+1) > opt_max_nratio) continue;
                        
                        //report the found motif occurrence
                        prev_op_match = match_pos;
                        ++reported_matches;
                        print_formated_output(opt_fasta,
                                        match_begin,
                                        match_end,
                                        sq_name,
                                        sq_details,
                                        ssearch.solution_to_str(i, seq_partitions[j].first, element_separator)
                                    );
                        if (opt_dotbracked == true) {
                            printf("%s\n", ssearch.solution_to_dotbracket(i, element_separator).c_str());
                        }
                        
                    }
                }
                
                //if order training is done, then evaluate the order on up to 100 windows
                /*if(ssearch.orderer->isSearchDone() && evalWindows < 101){
                    evalBases += seq_partitions[j].first.size();
                    if(opt_searchcomp){ evalBases += seq_partitions[j].first.size(); }
                    
                    if(evalWindows == 0){ //reset ops counters for evaluation
                        for(int i=1; i<ssearch.desc->sses.size(); i++){
                            ssearch.desc->sses[i].table.resetOpsCounter();
                        }
                        evalBases = 0;
                    }

                    //update statistics
                    unsigned long long ops = 0;
                    for(int i=1; i<ssearch.desc->sses.size(); i++){
                        ops += ssearch.desc->sses[i].table.ops();
                    }
                    opsPerBase = ops/(double)evalBases;

                    //if training only, then finish once the evaluation is done
                    if(evalWindows == 100 && opt_tonly){
                        doSearch = false;
                    }
                    
                    ++evalWindows;
                }*/
                
            }

        }
        total_matches += found_matches.size();
        total_matches += found_op_matches.size();

    }
    
    //print final messages
    if(opt_tonly && !ssearch.orderer->isSearchDone()){
        printf("\nWARNING! The database was not long enough to complete the training.\n\n");
    } else {
        printf("\n----- %s DONE -----\n", (opt_tonly?"TRAINING":"SEARCH"));
    }
    printf("Total bases scanned: %lld\n", total_bases_scanned);
    if(!opt_tonly){
        printf("Reported matches:    %lld (%lld total)\n", reported_matches, total_matches);
    }
    double elapsed = (double)(clock() - tStart)/CLOCKS_PER_SEC;
    printf("Time since start:    %02.0fh %02.0fm %02.0fs (%.2fs)\n",
           floor(elapsed/3600.0),
           floor(fmod(elapsed,3600.0)/60.0),
           fmod(elapsed,60.0),
           elapsed
          );
    printf("Final search order%s %s\n",
           (ssearch.orderer->isSearchDone() ? ": " : "*:"),
           ssearch.desc->search_order_to_str(ssearch.orderer->searchOrder).c_str()
          );
    /*printf("Est. difficulty:     ");
        if(opsPerBase == -1){
            printf("-\n");
        } else {
            printf("%.3f ops/base\n", opsPerBase);
        }*/
    
    //DEBUGGING AND DDEO PERFORMANCE MEASURE LOG 
    /*cout<<"\n";
    cout<<ssearch.orderer->fout.str()<<endl;
    if(ssearch.orderer->activeTuples.size() > 1) {
        cout<<"K: "<<ssearch.orderer->K<<endl;
        cout<<"tuples sampled means (ops per window):\n";
        for(int i=0;i<ssearch.orderer->tupleStats.size();++i){
            double sum=0;
            for(int j=0;j<ssearch.orderer->tupleStats[i].sampledOpsPerWindow.size(); ++j){
                sum += ssearch.orderer->tupleStats[i].sampledOpsPerWindow[j];
            }
            sum = sum/ssearch.orderer->tupleStats[i].sampledOpsPerWindow.size();
            cout.setf(ios::floatfield); 
            cout<<setprecision(5)<<sum<<" ";
        } cout<<endl;
        
        cout<<"number of samples per tuple: ("<<ssearch.orderer->tupleStats.size()<<")\n";
        for(int i=0;i<ssearch.orderer->tupleStats.size();++i){
            cout<<ssearch.orderer->tupleStats[i].sampledOpsPerWindow.size()<<" ";
        } cout<<endl;
        
        cout<<"FINALISTS: ";
        for(set<int>::iterator it=ssearch.orderer->activeTuples.begin(); it!=ssearch.orderer->activeTuples.end() ;++it){
            cout<<*it<<" ";
        }cout<<endl;
        cout<<"COUNT: "<<ssearch.orderer->scannedWindows<<"\nAVG.WIN: "
            <<ssearch.orderer->scannedBases/(double)ssearch.orderer->scannedWindows<<endl;
    }*/
    
    /*cout<<"ops stats from table: \n";
    unsigned long long ops = 0ULL;
    for(int i=0; i<ssearch.orderer->searchOrder.size(); i++){
        int id = ssearch.orderer->searchOrder[i];
        string name;
        for(map<string, int>::iterator it=ssearch.desc->transl.begin(); it!=ssearch.desc->transl.end(); ++it){
            if(it->second == id){
                name = it->first;
                break;
            }
        }
        cout<<name<<": "<<ssearch.desc->sses[id].table.ops()<<endl;
        ops += ssearch.desc->sses[id].table.ops();
    }
    cout<<"AllOps: "<<ops<<endl;
    cout<<"ops/second: "<<ops/elapsed<<endl;

    cout<<"\n\nops stats from tuples: \n";
    unsigned long long ops2 = 0ULL;
    for(int i=0; i<ssearch.orderer->tupleStats.size(); i++){
        struct TupleStats stats = ssearch.orderer->tupleStats[i];
        unsigned long long tupleOps = 0ULL;
        cout<<i<<": ";
        for(int j=0; j<stats.tuple.size(); ++j){
            cout<<stats.tuple[j]<<"$";
            string name;
                for(map<string, int>::iterator it=ssearch.desc->transl.begin(); it!=ssearch.desc->transl.end(); ++it){
                    if(it->second == stats.tuple[j]){
                        name = it->first;
                        break;
                    }
                }
            cout<<name<<"@"<<stats.memOps[j]<<" ";
            tupleOps += stats.memOps[j];
        }
        cout<<"\n   OpsForThisTuple: "<<tupleOps<<endl;
        ops2 += tupleOps;
        
        cout<<"   SampledOpsPerWindow: ";
        for(int j=0; j<stats.sampledOpsPerWindow.size(); ++j){
            cout<<stats.sampledOpsPerWindow[j]<<" ";
        }cout<<endl;
    }
    cout<<"AllOps: "<<ops2<<endl;
    cout<<"ops/second: "<<ops2/elapsed<<endl;

    cout<<"\nAllTOGETHER: "<<ops+ops2<<endl;
    cout<<"ops/second: "<<(ops+ops2)/elapsed<<endl;
    */
    return 0;
}
