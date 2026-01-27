#include <iostream>
#include <istream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <vector>
#include <istream>
#include <sstream>
#include <set>
#include <map>
#include <queue>
#include <unistd.h>

#include "telseq.h"
#include "Timer.h"
#include "math.h"
#include "Util.h"
#include "prettyprint.h"
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <mutex>
#include <thread>
#include <htslib/thread_pool.h>

//
// Getopt
//

#define PROGRAM_BIN "telseq"
#define AUTHOR "Tetsushi Nakao"

static const char *TELSEQ_VERSION_MESSAGE =
"TelSeq Version " PACKAGE_VERSION "\n"
"Original version from github.com/zd1/telseq\n"
"Modified by " AUTHOR ".\n";

static const char *TELSEQ_USAGE_MESSAGE =
"Program: " PACKAGE_NAME "\n"
"Version: " PACKAGE_VERSION "\n"
"Contact: " AUTHOR " [" PACKAGE_BUGREPORT "]\n\n"
"This package is modified TelSeq to accept CRAM file and multi-threading."
"Usage: " PROGRAM_BIN " [OPTION] <in.1.bam> <in.2.bam> <...> \n"
"Scan BAM and estimate telomere length. \n"
"   <in.bam>                 one or more BAM files to be analysed. File names can also be passed from a pipe, \n "
"                            with each row containing one BAM path.\n"
"   -f, --bamlist=STR        a file that contains a list of file paths of BAMs. It should has only one column, \n"
"                            with each row a BAM file path. -f has higher priority than <in.bam>. When specified, \n"
"                            <in.bam> are ignored.\n"
"   -o, --output_dir=STR     output file for results. Ignored when input is from stdin, in which case output will be stdout. \n"
"   -H                       remove header line, which is printed by default.\n"
"   -h                       print the header line only. The text can be used to attach to result files, useful\n"
"                            when the headers of the result files are suppressed. \n"
"   -k                       threshold of the amount of TTAGGG/CCCTAA repeats in read for a read to be considered telomeric. default = 7.\n"
"\nTesting functions\n------------\n"
"   -r                       read length. default = 100\n"
"   -z                       use user specified pattern for searching [ATGC]*.\n"
"   -e, --exomebed=STR       specifiy exome regions in BED format. These regions will be excluded \n"
"   -w,                      consider BAMs in the speicfied bamlist as one single BAM. This is useful when \n"
"                            the initial alignemt is separated for some reason, such as one for mapped and one for ummapped reads. \n"
"   -g, --reference_fasta    reference file is required for CRAM processing. .fai is also required.\n"
"   -p, --NUM_THREADS        number of threads. default = 4\n"
"   --help                   display this help and exit\n"

"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static StringVector bamlist;
    static std::string outputfile = "";
    static std::string exomebedfile = "";
    static std::string reference = "";
    static std::map< std::string, std::vector<range> > exomebed;
    static bool writerheader = true;
    static bool mergerg = false;
    static bool ignorerg = false;
    static bool onebam = false; // whether to consider all bams as one bam
    static int tel_k= ScanParameters::TEL_MOTIF_CUTOFF;
    static std::string unknown = "UNKNOWN";
    static std::string PATTERN;
    static std::string PATTERN_REV;
    static int NUM_THREADS = 4;

}

static const char* shortopts = "f:o:k:z:e:r:p:g:s:Hhvmuw";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "bamlist",		optional_argument, NULL, 'f' },
    { "output-dir",		optional_argument, NULL, 'o' },
    { "exomebed",		optional_argument, NULL, 'e' },
    { "reference",		optional_argument, NULL, 'g' },
    { "NUM_THREADS",         optional_argument, NULL, 'p' },
    { "help",               no_argument,       NULL, OPT_HELP },
    { "version",            no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

// combine counts in two ScanResults objects
void add_results(ScanResults& x, ScanResults& y){

    x.numTotal += y.numTotal;
    x.numMapped += y.numMapped;
    x.numDuplicates += y.numDuplicates;
    x.n_exreadsExcluded += y.n_exreadsExcluded;
    x.n_exreadsChrUnmatched += y.n_exreadsChrUnmatched;
    x.n_totalunfiltered += y.n_totalunfiltered;

    for (std::size_t j = 0, max = x.telcounts.size(); j != max; ++j){
        x.telcounts[j] +=y.telcounts[j];
    }
    for (std::size_t k = 0, max = x.gccounts.size(); k != max; ++k){
        x.gccounts[k] += y.gccounts[k];
    }

}

// merge results in result list into one
void merge_results_by_readgroup(

    std::vector< std::map<std::string, ScanResults> >& resultlist){
    std::vector< std::map<std::string, ScanResults> > mergedresultslist;
    std::map<std::string, ScanResults> mergedresults;

    for(size_t i =0; i< resultlist.size(); i++){
        auto rmap = resultlist[i];
        for(std::map<std::string, ScanResults>::iterator it= rmap.begin();
                it != rmap.end(); ++it){

            std::string rg = it ->first;
            ScanResults result = it -> second;

            if(mergedresults.find(rg) == mergedresults.end()){
                mergedresults[rg] = result;
            }else{
                add_results(mergedresults[rg], result);
            }
        }
    }
    mergedresultslist.push_back(mergedresults);
    resultlist = mergedresultslist;

}

// Prep for multi-threading
struct ThreadData {
    int thread_id;
    std::queue<bam1_t*> *queue;
    pthread_mutex_t *queue_mutex;
    pthread_mutex_t *count_mutex;
    pthread_cond_t *queue_not_empty;
    pthread_cond_t *queue_not_full;
    std::map<std::string, ScanResults> *resultmap;
};

void *processReads(void *threadarg) {
    struct ThreadData *my_data;
    my_data = (struct ThreadData *) threadarg;
    std::cerr << "Starting thread with ID " << my_data->thread_id << "\n";
    while (true) {
        pthread_mutex_lock(my_data->queue_mutex);
        while (my_data->queue->empty()) {
            pthread_cond_wait(my_data->queue_not_empty, my_data->queue_mutex);
        }
        bam1_t *b = my_data->queue->front();
        my_data->queue->pop();
        pthread_cond_signal(my_data->queue_not_full);
        pthread_mutex_unlock(my_data->queue_mutex);
        if (b == nullptr) {
            break;
        }
        // Process the read
        std::string tag = "unknown";

        pthread_mutex_lock(my_data->count_mutex);
          my_data->resultmap->operator[](tag).numTotal += 1;
        pthread_mutex_unlock(my_data->count_mutex);

        uint8_t *sequence = bam_get_seq(b); // Get a pointer to the sequence data
        std::string sequenceString;
        sequenceString.reserve(b->core.l_qseq);  // Pre-allocate memory for the sequence string
        for (int i = 0; i < b->core.l_qseq; ++i) {
	      sequenceString += seq_nt16_str[bam_seqi(sequence, i)]; // Convert integer to base character
	    }

        double gc = calcGC(sequenceString, b->core.l_qseq);
        int ptn_count = countMotif(sequenceString, opt::PATTERN, opt::PATTERN_REV);
        // when the read length exceeds 100bp, number of patterns might exceed the boundary
        if (ptn_count > ScanParameters::TEL_MOTIF_N - 1) {
          bam_destroy1(b);
          continue;
        }

        pthread_mutex_lock(my_data->count_mutex);
          (*my_data->resultmap)[tag].telcounts[ptn_count] += 1;
        pthread_mutex_unlock(my_data->count_mutex);

        if (gc >= ScanParameters::GC_LOWERBOUND && gc <= ScanParameters::GC_UPPERBOUND) {
          // get index for GC bin.
          int idx = floor((gc - ScanParameters::GC_LOWERBOUND) / ScanParameters::GC_BINSIZE);
          assert(idx >= 0 && idx <= ScanParameters::GC_BIN_N - 1);
          if (idx > ScanParameters::GC_BIN_N - 1) {
            std::cerr << " GC:{" << gc << "} telcounts:{" << ptn_count << "} GC bin index out of bound:" << idx << "\n";
            exit(EXIT_FAILURE);
          }

        pthread_mutex_lock(my_data->count_mutex);
          (*my_data->resultmap)[tag].gccounts[idx] += 1;
        pthread_mutex_unlock(my_data->count_mutex);

        }

        bam_destroy1(b);
    }
    pthread_exit(NULL);
}


// Does accept one argument only.
int scanCram() {
    std::vector<std::map<std::string, ScanResults>> resultlist;

    std::cout << opt::bamlist << "\n";
    std::cout << opt::bamlist.size() << " CRAMs" << std::endl;

    // Assuming opt::bamlist is the variable you're using from the opt namespace
    auto& bamlist = opt::bamlist;
    
    // Declare the mutex and errorMessages before the lambda
    std::mutex errorMutex;
    std::vector<std::string> errorMessages;

    // Initialize thread pool
    htsThreadPool p = {NULL, 0};
    p.pool = hts_tpool_init(opt::NUM_THREADS); // assuming opt::NUM_THREADS is the number of threads you want

    std::vector<std::thread> threads;
    std::mutex resultMutex;

    for (std::size_t i = 0; i < opt::bamlist.size(); i++) {
        std::thread t([&resultMutex, &resultlist, &bamlist, &errorMutex, &errorMessages, &p](std::size_t index) {
            // lambda function to capture the current loop variable
            // The code that processes each CRAM file
            // storing results for each read group (RG tag). use
            // read group ID as key.
            std::map<std::string, ScanResults> resultmap;
            // store where the overlap was last found in the case of exome seq
            std::map<std::string, std::vector<range>::iterator> lastfound;
            std::vector<range>::iterator searchhint;

            std::string tag = "unknown";

            // Extract sample name
            size_t lastSlash = opt::bamlist[index].find_last_of("/"); // Find last '/'
            size_t startPos = (lastSlash == std::string::npos) ? 0 : lastSlash + 1; // If there's no '/', start searching from the beginning of the string
            size_t dotCram = opt::bamlist[index].find(".cram", startPos); // Find '.cram' after the last slash or from the beginning
        
            if (dotCram != std::string::npos) {
                resultmap[tag].sample = opt::bamlist[index].substr(startPos, dotCram - startPos);
            } else {
                std::cout << "Could not find sample name in the file path provided." << std::endl;
            }

            // Open the cram files for reading/writing
            samFile *cramFile = sam_open(opt::bamlist[index].c_str(), "r");
            if (!cramFile) {
                std::lock_guard<std::mutex> lock(errorMutex);
                errorMessages.push_back("Failed to open CRAM file: " + opt::bamlist[index]);
                return;
            }

            // Get cram headers
            bam_hdr_t *header = sam_hdr_read(cramFile);
            bam1_t *read = bam_init1();

            // Load index file
            std::string crai = opt::bamlist[index] + ".crai";

            
            // Load reference
            if (hts_set_opt(cramFile, CRAM_OPT_REFERENCE, opt::reference.c_str()) < 0) {
                std::cerr << "Failed to set reference for " << opt::bamlist[index] << std::endl;
                
                // Cleanup resources on failure before returning or exiting
                bam_destroy1(read);
                bam_hdr_destroy(header);
                sam_close(cramFile);
                
                // Handle the error (either propagate it or just log it)
                return;
            }


            // Set the number of threads for reading
            hts_set_thread_pool(cramFile, &p);

            while ((sam_read1(cramFile, header, read)) >= 0) {

              resultmap[tag].numTotal += 1;

              uint8_t *sequence = bam_get_seq(read); // Get a pointer to the sequence data

              std::string sequenceString;
              int seq_len = read->core.l_qseq;
              sequenceString.reserve(seq_len);
              for (int i = 0; i < seq_len; ++i) {
                sequenceString += seq_nt16_str[bam_seqi(sequence, i)]; // Convert integer to base character
              }

              // std::cerr << "Sequence is " << sequenceString << "\n";

              double gc = calcGC(sequenceString, seq_len);
              int ptn_count = countMotif(sequenceString, opt::PATTERN, opt::PATTERN_REV);
              // when the read length exceeds 100bp, number of patterns might exceed the boundary
              if (ptn_count > ScanParameters::TEL_MOTIF_N - 1) {
                continue;
              }
              resultmap[tag].telcounts[ptn_count] += 1;

              if (gc >= ScanParameters::GC_LOWERBOUND && gc <= ScanParameters::GC_UPPERBOUND) {
                // get index for GC bin.
                int idx = floor((gc - ScanParameters::GC_LOWERBOUND) / ScanParameters::GC_BINSIZE);
                assert(idx >= 0 && idx <= ScanParameters::GC_BIN_N - 1);
                resultmap[tag].gccounts[idx] += 1;
              }

            }

            bam_destroy1(read);
            bam_hdr_destroy(header);
            sam_close(cramFile);

            std::lock_guard<std::mutex> lock(resultMutex);
            resultlist.push_back(resultmap);

        }, i);
        threads.push_back(std::move(t));
    }

    // Join threads after processing
    for (auto& thread : threads) {
        if (thread.joinable()) {
            thread.join();
        }
    }

    hts_tpool_destroy(p.pool); // Cleanup the thread pool

    outputresults(resultlist);
    std::cerr << "Completed writing results\n";

    return 0;
}

void printlog(std::vector< std::map<std::string, ScanResults> > resultlist){

	for(size_t i =0; i< resultlist.size(); i++){
		auto rmap = resultlist[i];
		for(std::map<std::string, ScanResults>::iterator it= rmap.begin();
				it != rmap.end(); ++it){

			std::string rg = it ->first;
			ScanResults result = it -> second;
			std::cout << "BAM:" << rg << std::endl;
			std::cout << "	chr ID unmatched reads: " << result.n_exreadsChrUnmatched << std::endl;
			std::cout << "	exome reads excluded: " << result.n_exreadsExcluded << std::endl;
		}
	}

}

void printout(std::string rg, ScanResults result, std::ostream* pWriter){

	*pWriter << rg << ScanParameters::FIELD_SEP;
	*pWriter << result.lib << ScanParameters::FIELD_SEP;
	*pWriter << result.sample << ScanParameters::FIELD_SEP;
	*pWriter << result.numTotal << ScanParameters::FIELD_SEP;
	*pWriter << result.numMapped << ScanParameters::FIELD_SEP;
	*pWriter << result.numDuplicates << ScanParameters::FIELD_SEP;

	result.telLenEstimate = calcTelLength(result);
	if(result.telLenEstimate==-1){
		std::cerr << "Telomere length estimate unknown. No read was found with telomeric GC composition.\n";
		*pWriter << opt::unknown << ScanParameters::FIELD_SEP;
	}else if(result.telLenEstimate>1000000){
		std::cerr << "Telomere length estimate unknown. Possibly due to not enough representation of genome.\n";
		*pWriter << opt::unknown << ScanParameters::FIELD_SEP;
	}else if(result.telLenEstimate==0){
		std::cerr << "Telomere length estimate unknown. No read contains more than " << opt::tel_k << " telomere repeats.\n";
		*pWriter << opt::unknown << ScanParameters::FIELD_SEP;
	}
	else{
		*pWriter << result.telLenEstimate << ScanParameters::FIELD_SEP;
	}

	for (std::size_t j = 0, max = result.telcounts.size(); j != max; ++j){
		*pWriter << result.telcounts[j] << ScanParameters::FIELD_SEP;
	}
	for (std::size_t k = 0, max = result.gccounts.size(); k != max; ++k){
		*pWriter << result.gccounts[k] << ScanParameters::FIELD_SEP;
	}
	*pWriter << "\n";
}

int outputresults(std::vector< std::map<std::string, ScanResults> > resultlist){

	std::ostream* pWriter;
	bool tostdout = opt::outputfile.empty() ? true:false;
	if(tostdout){
		pWriter = &std::cout;
	}else{
		pWriter = createWriter(opt::outputfile);
	}

	if(opt::writerheader){
		Headers hd;
		for(size_t h=0; h<hd.headers.size();h++){
			*pWriter << hd.headers[h] << ScanParameters::FIELD_SEP;
		}
		*pWriter << "\n";
	}

	for(size_t i=0; i < resultlist.size();++i){
		ScanResults mergedrs;
		std::string grpnames = "";

		std::map<std::string, ScanResults> resultmap = resultlist[i];

		// if merge read groups, take weighted average for all measures
		bool domg = opt::mergerg && resultmap.size() > 1? true:false;

		for(std::map<std::string, ScanResults>::iterator it= resultmap.begin();
				it != resultmap.end(); ++it){

			std::string rg = it ->first;
			ScanResults result = it -> second;

			if(domg){
				if(grpnames.size()==0){
					grpnames += rg;
				}else{
					grpnames += "|"+rg;
				}

				mergedrs.sample = result.sample;
				mergedrs.numTotal += result.numTotal;
				mergedrs.numMapped += result.numMapped * result.numTotal;
				mergedrs.numDuplicates += result.numDuplicates * result.numTotal;
				mergedrs.telLenEstimate += calcTelLength(result)* result.numTotal;

				for (std::size_t j = 0, max = result.telcounts.size(); j != max; ++j){
					mergedrs.telcounts[j] += result.telcounts[j]* result.numTotal;
				}
				for (std::size_t k = 0, max = result.gccounts.size(); k != max; ++k){
					mergedrs.gccounts[k] += result.gccounts[k]* result.numTotal;
				}
				continue;
			}else{
				printout(rg, result, pWriter);
			}
		}

		//in this case calculate weighted average
		if(domg){

			mergedrs.numMapped /= mergedrs.numTotal;
			mergedrs.numDuplicates /= mergedrs.numTotal;
			mergedrs.telLenEstimate /= mergedrs.numTotal;

			for (std::size_t j = 0, max = mergedrs.telcounts.size(); j != max; ++j){
				mergedrs.telcounts[j] /= mergedrs.numTotal;
			}
			for (std::size_t k = 0, max = mergedrs.gccounts.size(); k != max; ++k){
				mergedrs.gccounts[k] /= mergedrs.numTotal;
			}

			mergedrs.numTotal  /= resultmap.size();

			printout(grpnames, mergedrs, pWriter);
		};
	}

	if(!tostdout){
		delete pWriter;
	}
	return 0;
}

double calcTelLength(ScanResults results){

	float acc = 0;
	for(std::size_t i=opt::tel_k, max = results.telcounts.size(); i !=max; ++i){
		acc += results.telcounts[i];
	}

	float gc_tel = 0;
	for(std::size_t i=0, max = results.gccounts.size(); i !=max; ++i){
		float gc1=ScanParameters::GC_LOWERBOUND + ScanParameters::GC_BINSIZE*i;
		float gc2=ScanParameters::GC_LOWERBOUND + ScanParameters::GC_BINSIZE*(i+1);
		if(gc1 >= ScanParameters::GC_TELOMERIC_LOWERBOUND && gc2 <= ScanParameters::GC_TELOMERIC_UPPERBOUND ){
			gc_tel += results.gccounts[i];
		}
	}

	if(gc_tel == 0){
		return -1;
	}

	return (acc/gc_tel)*float(ScanParameters::GENOME_LENGTH_AT_TEL_GC)/ScanParameters::LENGTH_UNIT/ScanParameters::TELOMERE_ENDS;
}

int countMotif(const std::string &read, const std::string &pattern, const std::string &pattern_rev){
    int motifcount1 = 0;
    int motifcount2 = 0;

    size_t p1 = 0, p2 = 0;

    while (p1 < read.size() || p2 < read.size()) {
        p1 = read.find(pattern, p1);
        if (p1 != std::string::npos) {
            motifcount1++;
            p1 += pattern.size();
        }

        p2 = read.find(pattern_rev, p2);
        if (p2 != std::string::npos) {
            motifcount2++;
            p2 += pattern_rev.size();
        }
    }

    return std::max(motifcount1, motifcount2);
}

double calcGC(const std::string& seq, int length)
{
    int num_gc = 0;
    for(size_t i = 0; i < length; ++i)
    {
        char c = seq[i];
        if(c == 'C' || c == 'G')
            ++num_gc;
    }
    return static_cast<double>(num_gc) / length;
}

void update_pattern(){

	opt::PATTERN = ScanParameters::PATTERN;
	opt::PATTERN_REV = ScanParameters::PATTERN_REVCOMP;

    // update total motif counts when read length and/or pattern have been specified by user
    ScanParameters::TEL_MOTIF_N = ScanParameters::READ_LENGTH/ScanParameters::PATTERN.size() +1;
	if(opt::tel_k < 1 or opt::tel_k > ScanParameters::TEL_MOTIF_N-1){
		std::cerr << "k out of bound. k must be an integer from 1 to " <<  ScanParameters::TEL_MOTIF_N-1 << "\n";
		exit(EXIT_FAILURE);
	}

}


//
// Handle command line arguments
//
void parseScanOptions(int argc, char** argv)
{

	std::string bamlistfile =  "";
	std::string rev = "";

	Headers hd;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case 'f':
            	arg >> bamlistfile; break;
            case 'o':
            	arg >> opt::outputfile; break;
            case 'g':
            	arg >> opt::reference; break;
            case 'H':
            	opt::writerheader=false; break;
            case 'm':
                opt::mergerg = true; break;
            case 'u':
                opt::ignorerg = true; break;
            case 'w':
                opt::onebam = true; break;
            case 'h':
        		for(size_t h=0; h<hd.headers.size();h++){
        			std::cout << hd.headers[h] << ScanParameters::FIELD_SEP;
        		}
        		std::cout << "\n";
        		exit(EXIT_SUCCESS);
            case 'k':
            	arg >> opt::tel_k;
            	break;
            case 'r':
				arg >> ScanParameters::READ_LENGTH;
				if(ScanParameters::READ_LENGTH <= 0 || ScanParameters::READ_LENGTH > 100000){
					std::cerr << "please specify valid read length that is greater than 0 and length than 100kb" << "\n";
					exit(EXIT_FAILURE);
				}
				break;
            case 'p':
                arg >> opt::NUM_THREADS;
				break;
            case 'z':
            	arg >> ScanParameters::PATTERN;
				ScanParameters::PATTERN_REVCOMP = reverseComplement(ScanParameters::PATTERN);
				std::cerr << "use user specified pattern " <<  ScanParameters::PATTERN << "\n";
				std::cerr << "reverse complement " <<  ScanParameters::PATTERN_REVCOMP << "\n";
            	break;
            case 'e':
            	arg >> opt::exomebedfile;
            	opt::exomebed = readBedAsVector(opt::exomebedfile);
            	std::cout << "loaded "<< opt::exomebed.size() << " exome regions \n"<< std::endl;
//            	std::cout << opt::exomebed << "\n";
            	break;

            case OPT_HELP:
                std::cout << TELSEQ_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << TELSEQ_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    update_pattern();

    // deal with cases of API usage:
    // telseq a.bam b.bam c.bam ...
    // | telseq
    // telseq -f

    if (argc - optind < 1) // no argument specified
    {
    	// check if it is from pipe
		if(!isatty(fileno(stdin))){
			std::string line;
			while (std::getline(std::cin, line))
			{
				if(line.empty()){
					continue;
				}
				opt::bamlist.push_back(line);
			}
		}else if(bamlistfile.empty() ){ // check if not from a pipe, -f must be spceified
//    		std::cerr << SUBPROGRAM ": No BAM specified. Please specify BAM either directly, by using -f option or piping BAM file path.\n";
    		std::cout << TELSEQ_USAGE_MESSAGE;
    		exit(EXIT_FAILURE);
    	}
    }

    else if (argc - optind >= 1) // if arguments are specified
    {
    	// -f has higher priority, when specified, ignore arguments.
    	if(bamlistfile.empty()){
    		for(int i = optind; i < argc; ++i ){
    		    opt::bamlist.push_back(argv[i]);
    		}
    	}
    }

    // read in bamlist
    if(!bamlistfile.empty()){
        size_t filesize = getFilesize(bamlistfile);
        if(filesize == 0)
        {
            std::cerr << PROGRAM_BIN ": BAMLIST file specified by -f is empty\n";
            exit(EXIT_FAILURE);
        }

        std::istream* pReader = createReader(bamlistfile);
        std::string line;

        while(getline(*pReader, line))
        {
        	if(line.empty()){
        		continue;
        	}
            opt::bamlist.push_back(line);
        }

        size_t bamsize = opt::bamlist.size();
        if(bamsize == 0 ){
            std::cerr << PROGRAM_BIN ": Could not find any sample in BAMLIST file specified.\n";
            exit(EXIT_FAILURE);
        }
        delete pReader;
    }
}



int main(int argc, char** argv)
{
	Timer* pTimer = new Timer("scan BAM");
	parseScanOptions(argc, argv);
	scanCram();
    delete pTimer;
    return 0;
}