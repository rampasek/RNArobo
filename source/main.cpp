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
#include <cmath>
#include <getopt.h>

#include "search.h"

using namespace std;
using namespace GF;

#ifndef RELEASE
    #define DEF_RELEASE     "2.0.1"
    #define DEF_RELEASEDATE "August 2013"
#endif


static string txt_banner = "RNARobo - RNA motif searching program";
static string txt_usage =
"\nrnarobo: version "+string(DEF_RELEASE)+", "+string(DEF_RELEASEDATE)+"\n\
Usage: rnarobo [OPTIONS] <descriptor-file> <sequence-file>\n\
\n\
  Available options: \n\
     -c             search both strands of database\n\
     -f             print output in plain FASTA format\n\
     -s             print output in FASTA format with element separators\n\
     \n\
  To override default search order training parameters:\n\
  (defaults: --k 3 --limit 40 --alpha 0.025 --iterative TRUE )\n\
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
    int opt_fasta = false;          /* OPTION: TRUE to print output in FASTA format*/
    int opt_dotbracked = false;     /* OPTION: TRUE to print also dot-bracked representation*/
    
    int opt_tonly = false;   /* OPTION: TRUE to perform only the training of search ordering */
    int opt_k = -1;
    int opt_limit = -1;
    int opt_alpha = -1;
    int opt_iterative = -1;
    string opt_paramsfile = "";
    const char *alphas[] = {"0.2", "0.1", "0.05", "0.025", "0.01"};
    pair<vector<double>, vector<double> > params;
    double IC = 1;
    double DF = -0.3;
    params.first.push_back(4*IC);
    params.first.push_back(2*IC);
    params.first.push_back(1*IC);
    params.second.push_back(4*DF);
    params.second.push_back(2*DF);
    params.second.push_back(1*DF);

    static struct option long_options[] = {
        {"k", 1, 0, 'k'},
        {"limit", 1, 0, 'l'},
        {"alpha", 1, 0, 'a'},
        {"iterative", 1, 0, 'i'},
        {"tonly", 0, 0, 't'},
        {"params", 1, 0, 'p'},
        {"complement", 0, 0, 'c'}, //long alias for -c
        {"fasta", 0, 0, 'f'}, //long alias for -f
        {"sep", 0, 0, 's'}, //long alias for -s
        {NULL, 0, NULL, 0}
    };
    extern int optind;
    int c, option_index = 0;
    while ((c = getopt_long(argc, argv, "cfs", long_options, &option_index)) != -1) {
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
    unsigned int param_limit = opt_tonly?60:40;
    unsigned int param_alpha = opt_tonly?3:3;
    bool param_iterative = true;
    //override default parameters by the user-defined
    if(opt_k != -1) param_k = opt_k;
    if(opt_limit != -1) param_limit = opt_limit;
    if(opt_alpha != -1) param_alpha = opt_alpha;
    if(opt_iterative != -1) param_iterative = opt_iterative;
    if(opt_paramsfile != "") {
        ifstream pfile(opt_paramsfile.c_str(), ifstream::in);
        for (int i=0; i<param_k; i++) pfile >> params.first[i];
        for (int i=0; i<param_k; i++) pfile >> params.second[i];
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
        //printf("IC parameters for heuristic:"); for (int i=0; i<param_k; i++) printf(" %f", params.first[i]); printf("\n");
        //printf("DF parameters for heuristic:"); for (int i=0; i<param_k; i++) printf(" %f", params.second[i]); printf("\n");
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
    
    long long evalBases = 0;
    int evalWindows = 0;
    double opsPerBase = -1;
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
            normalize_seq(line);
            sq+=line;
        }
        filter_whitespaces(sq);
        //cout<<sq_name<<":"<<sq<<endl<<endl;
        
        //count number of bases scanned
        total_bases_scanned += sq.size();
        if(opt_searchcomp){ total_bases_scanned += sq.size(); }
            
        int max_motif_length = ssearch.desc->get_max_motif_length();
        //to how long pieces we will chop up the sequence
        int max_seq_length = max(5000, 20*max_motif_length);
        
        //to store beginnings of found matches (in both strands) - for filtering repeating matches
        set <pair <unsigned int, unsigned int> > found_matches, found_op_matches;
        found_matches.clear();
        found_op_matches.clear();

        //cut out "N"regions longer than 10 -> the sequence is divided into blocks
        vector< pair<string, unsigned int> > seq_blocks; //first is block's sequence, second is block's beginning position in original seqence
        int prev_found = 0;
        int found = sq.find("NNNNNNNNNN");
        while(found!=string::npos){
            seq_blocks.push_back(make_pair(sq.substr(prev_found, found-prev_found+1), prev_found));
            while(sq[found+10]=='N') ++found;
            prev_found = found;
            found=sq.find("NNNNNNNNNN",found+1);
        }
        seq_blocks.push_back(make_pair(sq.substr(prev_found, sq.size()), prev_found));

        for(int k=0; k<seq_blocks.size() && doSearch; k++){
            //chop up the given block to partitions of length at most max_seq_length
            // ! two partitions must have overlap of max_motif_length !
            int pointer = 0;
            vector< pair<string, unsigned int> > seq_partitions;
            while(pointer < seq_blocks[k].first.size()){
                int start = max(0, pointer-max_motif_length+1);

                string s = seq_blocks[k].first.substr(start, max_seq_length);
                seq_partitions.push_back(make_pair(s, start));
                pointer = start + max_seq_length;
            }

            //search for the motif in partitions
            for(int j=0; j<seq_partitions.size() && doSearch; j++){
                ssearch.search(seq_partitions[j].first);

                unsigned int offset = seq_partitions[j].second+ seq_blocks[k].second;

                //print found hits
                for(int i=0;i<ssearch.solutions.size() && !opt_tonly;i++){
                    unsigned int match_begin = ssearch.solutions[i][0].first +1 +offset;
                    unsigned int match_end = ssearch.solutions[i].back().second +offset;
                    pair<unsigned int, unsigned int> match_pos = make_pair(match_begin, match_end);
                    
                    if(found_matches.find(match_pos) != found_matches.end()){
                        continue;
                    } else {
                        found_matches.insert(match_pos);

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

                //search the partition also in the opposite direction (if opt_searchcomp == true)
                if(opt_searchcomp){
                    reverse_complement(seq_partitions[j].first.begin(), seq_partitions[j].first.end());
                    ssearch.search(seq_partitions[j].first);

                    //print found hits
                    for(int i=0;i<ssearch.solutions.size() && !opt_tonly;i++){
                        unsigned int match_begin = seq_partitions[j].first.size() -ssearch.solutions[i][0].first + offset;
                        unsigned int match_end = seq_partitions[j].first.size() -ssearch.solutions[i].back().second +1 +offset;
                        pair<unsigned int, unsigned int> match_pos = make_pair(match_begin, match_end);
                        
                        if(found_op_matches.find(match_pos) != found_op_matches.end()){
                            continue;
                        } else {
                            found_op_matches.insert(match_pos);
                        
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
                }
                
                //if order training is done, then evaluate the order
                if(ssearch.orderer->samplingSearchDone && evalWindows < 302){
                    evalBases += seq_partitions[j].first.size();
                    if(opt_searchcomp){ evalBases += seq_partitions[j].first.size(); }
                    
                    if(evalWindows == 0){
                        for(int i=1; i<ssearch.desc->sses.size(); i++){
                            ssearch.desc->sses[i].table.resetOpsCounter();
                        }
                        
                        evalBases = 0;
                    }
                    
                    if(evalWindows == 300){
                        unsigned long long ops = 0;
                        for(int i=1; i<ssearch.desc->sses.size(); i++){
                            ops += ssearch.desc->sses[i].table.ops();
                        }
                        opsPerBase = ops/(double)evalBases;
                        
                        //if training only, then finish
                        if(opt_tonly) doSearch = false;
                    }
                    
                    ++evalWindows;
                }
                
            }

        }
        total_matches += found_matches.size();
        total_matches += found_op_matches.size();

    }
    
    //print final messages
    if(opt_tonly && !ssearch.orderer->samplingSearchDone){
        printf("\nWARNING! The database was not long enough to complete the training.\n\n");
    } else {
        printf("\n----- %s DONE -----\n", (opt_tonly?"TRAINING":"SEARCH"));
    }
    printf("Total bases scanned: %lld\n", total_bases_scanned);
    if(!opt_tonly){
        printf("Found matches:       %lld\n", total_matches);
    }
    double elapsed = (double)(clock() - tStart)/CLOCKS_PER_SEC;
    printf("Time since start:    %02.0fh %02.0fm %02.0fs (%.2fs)\n",
           floor(elapsed/3600.0),
           floor(fmod(elapsed,3600.0)/60.0),
           fmod(elapsed,60.0),
           elapsed
          );
    printf("Final search order%s %s\n",
           (opsPerBase==-1 ? "*:" : ": "),
           ssearch.desc->search_order_to_str(ssearch.orderer->searchOrder).c_str()
          );
    printf("Est. difficulty:     ");
        if(opsPerBase == -1){
            printf("-\n");
        } else {
            printf("%.3f ops/base\n", opsPerBase);
        }
    
    //DEBUGGING AND PERFORMANCE MEASURE LOG 
    /*cout<<"\n";
    cout<<ssearch.orderer->fout.str()<<endl;
    if(ssearch.orderer->activeTuples.size() > 1) {
        cout<<"K: "<<ssearch.orderer->K<<endl;
        cout<<"tuples sampled means (memOPs per window):\n";
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
    
    
/*
    cout<<"memOps stats from table: \n";
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
    cout<<"AllMemOps: "<<ops<<endl;
    cout<<"MemOps/second: "<<ops/elapsed<<endl;

    cout<<"\n\nmemOps stats from tuples: \n";
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
    cout<<"AllMemOps: "<<ops2<<endl;
    cout<<"MemOps/second: "<<ops2/elapsed<<endl;

    cout<<"\nAllTOGETHER: "<<ops+ops2<<endl;
    cout<<"MemOps/second: "<<(ops+ops2)/elapsed<<endl;
*/
    return 0;
}
