/*****************************************************************************
 *   TakeABreak: Breakpoint detection from raw unassembled NGS reads
 *   A tool from the GATB (Genome Assembly Tool Box)
 *   Copyright (C) 2014  INRIA
 *   Authors: C.Lemaitre, P.Peterlongo, E.Drezen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#include <gatb/gatb_core.hpp>
#include <DBGWalker.hpp>
#include <TakeABreak.hpp>
#include <iostream>
#include <algorithm>
#include <LCS.hpp>
#include <Solution.hpp>
#include <time.h>
#include <limits> //for std::numeric_limits



//#define check_memory // http://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
#ifdef check_memory
#include<mach/mach.h>
struct task_basic_info t_info;
mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
#endif


// We use the required packages

using namespace std;

//#define debug
//#define only_canonical //no longer usefull
//#define COMBINATORIAL_FACTOR 1000
#define DEBUG(a)  printf a

#include <gatb/system/api/config.hpp>
//to get the version number of gatb-core
// WARNING when doing cmake -DMAJOR=3 -DMINOR=4 -DPATCH=11 .., this will override the gatb-core version (but only done when making delivery, the user should classically build the project with cmake .. (except for binary versions...)
// variable name is STR_LIBRARY_VERSION

char * getVersion(){
	return (char *)"1.1.2"; // TO UPDATE AT EACH RELEASE
    //return (char *)STR_LIBRARY_VERSION;
}




/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
TakeABreak::TakeABreak () : Tool("TakeABreak"), _kmerSize(27), shannon_limit(1.7), synchro(0)
{

    /* Commented this until there is a function hide() for parser
     because we want some of the options of Graph to be masked to the user
     plus it enables to master the order of appearance of options in the Help message
    //Getting the options of graph (ie. dbgh5)
    //OptionsParser parser = Graph::getOptionsParser(false);
    //getParser()->add (parser);
     */
	//getParser()->push_back (Graph::getOptionsParser(false),0);

//	IOptionsParser* generalParser = getParser();
//	generalParser->setName("general");
	IOptionsParser* generalParser = new OptionsParser("General");
	generalParser->push_back (new OptionNoParam (STR_HELP, "help", false));
	generalParser->push_back (new OptionNoParam (STR_VERSION, "version", false));
	generalParser->push_back (new OptionOneParam (STR_NB_CORES,    "number of cores",      false, "0"  ));
	generalParser->push_back (new OptionOneParam (STR_VERBOSE,     "verbosity level",      false, "1"  ));
//
	setParser (new OptionsParser ("TakeABreak"));

	IOptionsParser* inputParser = new OptionsParser("Input / output");
	inputParser->push_front (new OptionOneParam (STR_URI_OUTPUT, "prefix for output files", false, ""));
	inputParser->push_front (new OptionOneParam (STR_URI_GRAPH, "input graph file (likely a hdf5 file)",  false, ""));
	inputParser->push_front (new OptionOneParam (STR_URI_INPUT, "input read file(s)",  false, ""));

	IOptionsParser* inversionParser = new OptionsParser("Inversion");
	inversionParser->push_front (new OptionOneParam (STR_LCT, "local complexity threshold (LCT)", false, "100"));
	inversionParser->push_front (new OptionOneParam (STR_TOLERANCE_RC, "maximal repeat size at the breakpoint (longest common suffix between a and b')", false, "8"));
	inversionParser->push_front (new OptionOneParam (STR_MAX_SIM, "max similarity percentage between a and b' and between u and v'", false, "80"));

	IOptionsParser* graphParser = new OptionsParser("Graph");
	graphParser->push_front (new OptionOneParam (STR_MAX_MEMORY, "max memory (in MBytes)", false, "2000"));
	graphParser->push_front (new OptionOneParam (STR_MAX_DISK, "max disk   (in MBytes)", false, "0"));
	graphParser->push_front (new OptionOneParam (STR_SOLIDITY_KIND, "way to consider a solid kmer with several datasets (sum, one, or all)", false, "one"));
	string abundanceMax = Stringify::format("%ld", std::numeric_limits<CountNumber>::max()); //to be sure in case CountNumber definition changes
	graphParser->push_front (new OptionOneParam (STR_KMER_ABUNDANCE_MAX, "maximal abundance threshold for solid kmers", false, abundanceMax));
	graphParser->push_front (new OptionOneParam (STR_KMER_ABUNDANCE_MIN, "minimal abundance threshold for solid kmers", false, "auto"));
	graphParser->push_front (new OptionOneParam (STR_KMER_SIZE, "size of a kmer", false, "31"));

	getParser()->push_front(generalParser);
	getParser()->push_front(inversionParser);
	getParser()->push_front(graphParser);
	getParser()->push_front(inputParser);

	//OptionFailure(getParser(), "options -graph and -in are incompatible, but at least one of these is mandatory");

//	IOptionsParser* graphParser = Graph::getOptionsParser(false);
//
//	if (IOptionsParser* generalParser = graphParser->getParser ("general"))  {  generalParser->setVisible(false);  }
//	if (IOptionsParser* branchingParser = graphParser->getParser ("general"))  {  branchingParser->setVisible(false);  }
//	if (IOptionsParser* inParser = graphParser->getParser (STR_URI_INPUT))  {  inParser->setVisible(false);  }
//	if (IOptionsParser* inParser = graphParser->getParser (STR_URI_INPUT))  {  inParser->setVisible(false);  }

	// Here getParser() has already the options inherited from Tool : nb-cores, verbose et help
    // push_front options, in the order of more and more importance
//    getParser()->push_back (new OptionNoParam (STR_VERSION, "version", false));
//
//    // from Graph.cpp
//    getParser()->push_front (new OptionOneParam (STR_MAX_MEMORY, "max memory (in MBytes)", false, "2000"));
//    getParser()->push_front (new OptionOneParam (STR_MAX_DISK, "max disk   (in MBytes)", false, "0"));
//
//    getParser()->push_front (new OptionOneParam (STR_LCT, "local complexity threshold (LCT)", false, "100"));
//    getParser()->push_front (new OptionOneParam (STR_TOLERANCE_RC, "maximal repeat size at the breakpoint (longest common suffix between a and b')", false, "8"));
//    getParser()->push_front (new OptionOneParam (STR_MAX_SIM, "max similarity percentage between a and b' and between u and v'", false, "80"));
//
//    getParser()->push_front (new OptionOneParam (STR_SOLIDITY_KIND, "way to consider a solid kmer with several datasets (sum, min or max)", false, "sum"));
//    getParser()->push_front (new OptionOneParam (STR_KMER_ABUNDANCE_MAX, "maximal abundance threshold for solid kmers", false, "4294967295"));
//    getParser()->push_front (new OptionOneParam (STR_KMER_ABUNDANCE_MIN, "minimal abundance threshold for solid kmers", false, "3"));
//    getParser()->push_front (new OptionOneParam (STR_KMER_SIZE, "size of a kmer", false, "31"));
//    getParser()->push_front (new OptionOneParam (STR_URI_OUTPUT, "prefix for output files", false, ""));
//    getParser()->push_front (new OptionOneParam (STR_URI_GRAPH, "input graph file (likely a hdf5 file)",  false, ""));
//    getParser()->push_front (new OptionOneParam (STR_URI_INPUT, "input read file(s)",  false, ""));
//
    //remove unused options -- should be replaced by future function hide()
//    getParser()->remove(STR_BLOOM_TYPE);
//    getParser()->remove(STR_DEBLOOM_TYPE);
//    getParser()->remove(STR_MPHF_TYPE);
//    getParser()->remove(STR_BANK_CONVERT_TYPE);
//    getParser()->remove(STR_URI_SOLID_KMERS);
//    getParser()->remove(STR_BRANCHING_TYPE);
//    getParser()->remove(STR_URI_OUTPUT_DIR);
    
    // remove redundant options -- should not be useful anymore (if dealt in gatb-core directly)
//    getParser()->remove(STR_NB_CORES);
//    getParser()->remove(STR_VERBOSE);
//    getParser()->remove(STR_HELP);
    
    
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
TakeABreak::~TakeABreak ()
{
    if (synchro)  { delete synchro; }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void TakeABreak::displayVersion(std::ostream& os){

	os << "* * * * * * * * * * * * * * * * * * * * * *" << endl;
	os << "* TakeABreak version "<< getVersion() << "                *" << endl; //<< " AGPL licence" <<endl;
	os << "* Using gatb-core version "<< STR_LIBRARY_VERSION <<  "           *" << endl;
	os << "* Supported kmer sizes <" << KSIZE_STRING << "   *" << endl;
	os << "* * * * * * * * * * * * * * * * * * * * * *" << endl;
	//os << "*********************************************************" << endl;


	//os << "built " << gatb::core::system::impl::System::info().getBuildDate().c_str() << endl; // build date of gatb-core, not TakeABreak
	//        cout<< Stringify::format ("* version %s (%s)\n* built on %s with compiler '%s'\n* supported kmer sizes %d %d %d %d",
	//                    gatb::core::system::impl::System::info().getVersion().c_str(),
	//                    gatb::core::system::impl::System::info().getBuildDate().c_str(),
	//                    gatb::core::system::impl::System::info().getBuildSystem().c_str(),
	//                    gatb::core::system::impl::System::info().getBuildCompiler().c_str(),
	//                    KSIZE_1, KSIZE_2, KSIZE_3, KSIZE_4
	//                    )
	//                    << endl;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void TakeABreak::execute ()
{

	//option VERSION is managed by Tool at the beginning of method run(), before calling execute()
	// To customize the version message : override the method displayVersion()
//    if (getInput()->get(STR_VERSION) != 0){
//        cout << "TakeABreak version "<< getVersion() << " AGPL licence" <<endl;
//        cout << "Using gatb-core version "<< STR_LIBRARY_VERSION << endl;
//        return;
//    }

	if (getInput()->get(STR_HELP) != 0){
		throw OptionFailure(getParser(), "");
	}

    if ((getInput()->get(STR_URI_GRAPH) != 0 && getInput()->get(STR_URI_INPUT) != 0) || (getInput()->get(STR_URI_GRAPH) == 0 && getInput()->get(STR_URI_INPUT) == 0))
    {
        throw OptionFailure(getParser(), "Error: options -graph and -in are incompatible, but at least one of these is mandatory");
    }

    // If outputPrefix is not provided we create one using the current date-time
    if (getInput()->get(STR_URI_OUTPUT) == 0)
    {
        time_t     now = time(0);
        struct tm  tstruct;
        char       buf[80];
        tstruct = *localtime(&now);
        strftime(buf, sizeof(buf), "%Y-%m-%d.%I:%M", &tstruct);
        string outputPrefix="TakeABreak_Expe-"+string(buf);

        getInput()->add (0, STR_URI_OUTPUT, outputPrefix);

    }
    
    
//    _log_file=getInput()->getStr(STR_URI_OUTPUT)+".log";
//    FILE * log = fopen(_log_file.c_str(), "w");
//    if(log != NULL){
//        cout << "Log info are dumped in file " << _log_file << endl;
//    }
    
    
    // Getting the graph
    
    // Case 1 : -in option, we create the graph from read files
    if (getInput()->get(STR_URI_INPUT) != 0)
    {
        //fprintf(log,"Creating the graph from file(s) %s\n",getInput()->getStr(STR_URI_INPUT).c_str());
        
        // We need to add the options of dbgh5/Graph that were masked to the user (or we could create a new Properties object)
        getInput()->add(0,STR_BANK_CONVERT_TYPE,"tmp");
        getInput()->add(0,STR_URI_OUTPUT_DIR, ".");
        getInput()->add(0,STR_URI_OUTPUT_TMP, ".");
        getInput()->add(0,STR_BLOOM_TYPE, "neighbor"); //neighbor
        getInput()->add(0,STR_DEBLOOM_TYPE, "cascading");
        getInput()->add(0,STR_DEBLOOM_IMPL, "minimizer"); //minimizer => STR_BLOOM_TYPE = neighbor
        getInput()->add(0,STR_BRANCHING_TYPE, "stored");
        getInput()->add(0,STR_INTEGER_PRECISION, "0");
        getInput()->add(0,STR_MPHF_TYPE, "none");
        getInput()->add(0,STR_MINIMIZER_SIZE, "8");
        getInput()->add(0,STR_REPARTITION_TYPE, "0");
        getInput()->add(0,STR_MINIMIZER_TYPE, "0");
        getInput()->add(0,STR_HISTOGRAM_MAX, "10000");
        getInput()->add(0,STR_KMER_ABUNDANCE_MIN_THRESHOLD,"3");
        //getInput()->add(0,STR_URI_SOLID_KMERS, ""); //surtout ne pas decommenter cette ligne, sinon les kmers solids sont stockes dans le fichier ./.h5 et les infos ne sont plus dans le output.h5
        
        //Warning if kmer size >128 cascading debloom does not work
        if(getInput()->getInt(STR_KMER_SIZE)>128){
            getInput()->get(STR_DEBLOOM_TYPE)->value="original";
        }
        
        _graph = Graph::create (getInput());
        //cout << _graph.getInfo() << endl; // to debug can be useful
        _kmerSize = getInput()->getInt(STR_KMER_SIZE);
        
    }
    
    // Case 2 : -graph option, we load the graph from a .h5 file
    if (getInput()->get(STR_URI_GRAPH) != 0)
    {
        //fprintf(log,"Loading the graph from file %s\n",getInput()->getStr(STR_URI_GRAPH).c_str());
        
        _graph = Graph::load (getInput()->getStr(STR_URI_GRAPH));
        //cout << _graph.getInfo() << endl; // to debug can be useful
        _kmerSize = _graph.getKmerSize();
    }
    
//    stringstream infoGraph;
//    infoGraph << _graph.getInfo();
//    fprintf(log,"%s",infoGraph.str().c_str());
    
    
    _output_file=getInput()->getStr(STR_URI_OUTPUT)+".fasta";
    FILE * out = fopen(_output_file.c_str(), "w");
    if(out == NULL){
        //cerr <<" Cannot open file "<< _output_file <<" for writting" << endl;
        string message = "Cannot open file "+ _output_file + " for writting";
        throw Exception(message.c_str());

    }

    
    //Getting Inversion parameters
    _tolerance_rc = getInput()->getInt(STR_TOLERANCE_RC);
    _nbCores = getInput()->getInt(STR_NB_CORES);
    _max_sim = getInput()->getInt(STR_MAX_SIM);
    _LCT = getInput()->getInt(STR_LCT);
    
    if(_tolerance_rc>_kmerSize-1){
        _tolerance_rc=_kmerSize-1;
        fprintf(stderr," Warning : reverse tolerance can't be bigger than k-2, setting reverse tolerance to k-1=%d \n", _tolerance_rc);
    }

   
    // Inversion detection
    LCS lcs_instance(_kmerSize, _max_sim);
    synchro = System::thread().newSynchronizer();

    time_t start = time(0);
    find_ALL_occurrences_of_inversion_pattern(lcs_instance);
    size_t number_inv_found=writeResults(out);
    time_t end = time(0);
    fclose(out);

    
#ifdef PRINTALL
    string output_file2=getInput()->getStr(STR_URI_OUTPUT)+"_all.fasta";
    FILE * out2 = fopen(output_file2.c_str(), "w");
    size_t number_untrunc=writeUntruncResults(out2);
    cout << "Nb occurrences found = " << _nbOccurrences << endl;
    cout << "Nb untruncated solutions found = " << number_untrunc << endl;
    cout << "Nb canonical solutions found = " << number_inv_found << endl;
#endif

    double seconds=difftime(end,start);
    //cout << seconds << endl;
    
    // Printing info on the run
    //getInfo()->get(getName())->value="done"; // par defaut une cle = getName() sans valeur
    getInfo()->add(1,"version",getVersion());
    getInfo()->add(1,"gatb-core-library",STR_LIBRARY_VERSION);
    getInfo()->add(1,"supported_kmer_sizes","%s", KSIZE_STRING);
    resumeParameters();
    resumeResults(number_inv_found,seconds);
    
//    stringstream infoResults;
//    infoResults << *getInfo();
//    fprintf(log,"%s",infoResults.str().c_str());
//    fclose(log);


}


void TakeABreak::resumeParameters(){
    
    //Properties resumeParams;
	//getInfo()->add(0,"test",getInput()->getStr(STR_MINIMIZER_SIZE).c_str());
    getInfo()->add(0,"Parameters");
    getInfo()->add(1,"Input data");
    if (getInput()->get(STR_URI_INPUT) != 0){
        getInfo()->add(2,"Reads",getInput()->getStr(STR_URI_INPUT).c_str());
    }
    if (getInput()->get(STR_URI_GRAPH) != 0){
        getInfo()->add(2,"Graph",getInput()->getStr(STR_URI_GRAPH).c_str());
    }
    getInfo()->add(1,"Graph");
    getInfo()->add(2,"kmer-size","%i", _kmerSize);

    try { // version actuelle info manquante si -graph
    	getInfo()->add(2,"nb_banks",_graph.getInfo().getStr("nb_banks").c_str());
    } catch (Exception e) {
    	// doing nothing
    }
    getInfo()->add(2,"solidity_kind",_graph.getInfo().getStr("solidity_kind").c_str());
    //getInfo()->add(2,"abundance_min",_graph.getInfo().getStr("abundance_min").c_str()); // to change en thresholds
    try { // entour try/catch ici au cas ou le nom de la cle change dans gatb-core
    	getInfo()->add(2,"abundance_min (auto inferred)",_graph.getInfo().getStr("cutoffs_auto.values").c_str());
    } catch (Exception e) {
    	// doing nothing
    }
    string min_abundance;
    if(_graph.getInfo().getStr("solidity_kind")=="sum"){
    	int thre = _graph.getInfo().getInt("thresholds"); //with getInt obtains the first number (if sum : threshold = x x x ...)
    	stringstream ss;
    	ss << thre;
    	min_abundance = ss.str();
    }
    else{
    	min_abundance = _graph.getInfo().getStr("thresholds").c_str();
    }
    getInfo()->add(2,"abundance_min (used)",min_abundance);
    try { // version actuelle info manquante si -graph
    	getInfo()->add(2,"abundance_max",_graph.getInfo().getStr("abundance_max").c_str());
    } catch (Exception e) {
    	// doing nothing
    }
    getInfo()->add(2,"nb_solid_kmers",_graph.getInfo().getStr("kmers_nb_solid").c_str());
    getInfo()->add(2,"nb_branching_nodes",_graph.getInfo().getStr("nb_branching").c_str());

    getInfo()->add(1,"Inversions");
    getInfo()->add(2,"repeat","%i", _tolerance_rc);
    getInfo()->add(2,"max_sim","%i", _max_sim);
    getInfo()->add(2,"LCT","%i", _LCT);
    
    //return resumeParams;
}

void TakeABreak::resumeResults(size_t number_inv_found, double seconds){
    
    getInfo()->add(0,"Results");
    getInfo()->add(1,"nb_occurrences", "%i",_nbOccurrences);
    getInfo()->add(1,"nb_distinct_solutions", "%i",number_inv_found);
    getInfo()->add(1,"time", "%.1f s",seconds);
    getInfo()->add(1,"output_files");
    if(getInput()->get(STR_URI_INPUT) != 0){
        getInfo()->add(2,"graph_file", "%s.h5",getInput()->getStr(STR_URI_OUTPUT).c_str());
    }
    getInfo()->add(2,"solution_file", "%s",_output_file.c_str());
    //getInfo()->add(2,"log_file", "%s",_log_file.c_str());

}

// dirty but efficient: declare
char int2char['z'];

// NO LONGER USED : for Shannon index
void init_int2char(){
    
    int i;
    for(i=0;i<'U';i++) int2char[i]=-1;
    
    int2char[(int)'a']=0;
    int2char[(int)'A']=0;
    int2char[(int)'c']=1;
    int2char[(int)'C']=1;
    int2char[(int)'g']=2;
    int2char[(int)'G']=2;
    int2char[(int)'t']=3;
    int2char[(int)'T']=3;
    int2char[(int)'n']=4;
    int2char[(int)'N']=4;
    
}

// NO LONGER USED : for Shannon index
inline int get_value_char(char val)
{
	return (int)int2char[(int)val];
}

// NO LONGER USED
// Compute the Shannon index of a sequence //Fixme: optimizer en lisant plus d'un caractere a la fois
float shannon_index(const string& read)
{
	float index = 0;
	float freq[5];
	int i;
    const int size_read = read.length();
	for(i=0; i<5; i++) freq[i]=0;
    
	// Frequency of each letter (A, C, G, T or N)
	for(i=0; i<read.length(); i++)
	{
		freq[get_value_char(read[i])]++;
	}
	// Shannon index
	for(i=0; i<5; i++)
	{
		freq[i] = freq[i]/size_read;
		if(freq[i] != 0) index += freq[i]*(log(freq[i])/log(2));
	}
    
	return fabs(index);
}

// NO LONGER USED Returns True if the Shannon entropy of the node is above the limit (used as a filter)
bool conserve_node(const string& read, const float shannon_limit){
    // if the first or last X characters are the same, discard the node
    int i;
//    int X = read.length()/2;
    int X = 8; // TODO: parameter
    
   
    for(i=1;i<X;i++) if(read[i-1]!=read[i]) break;
    if(i==X) return false;
    
    for(i=read.length()-X+1;i<read.length();i++) if(read[i-1]!=read[i]) break;
    if(i==read.length()) return false;
    
    return shannon_index(read)>shannon_limit;
}

// NO LONGER USED Counts the number of substitutions between 2 kmers (when aligned naively)
inline int number_substitutions(const string& a, const string& b){
    int res=0;
    for(int i=0;i<a.length();i++) if(a[i]!=b[i]) res++;
//    cout<<a<<endl<<b<<" = "<<res<<" subst"<<endl;
    return res;
}


// USED ONLY AS A TEST : Returns True if prefixes of size _tolerance_rc of both strings are different
bool TakeABreak::check_tolerance(const Node& s1, const Node& s2, const int _tolerance_rc){
    if(_graph.toString(s1).substr(0,_tolerance_rc+1).compare(_graph.toString(s2).substr(0,_tolerance_rc+1))==0) return false;
    return true;
}

// NO LONGER USED : Returns True if both pairs a-bbar and u-vbar have more than 10 substutions (when aligned naively)
bool TakeABreak::conserve_inversion(const Node& a, const Node& u, const Node& v, const Node& b){
    const int limit_subsitutions=10; // TODO parameter?
    // checks that a is different enough from rc(b)
    if(number_substitutions(_graph.toString(a),_graph.toString(_graph.reverse(b)))<=limit_subsitutions) return false;
    // checks that u is different enough from rc(v)
    if(number_substitutions(_graph.toString(u),_graph.toString(_graph.reverse(v)))<=limit_subsitutions) return false;
    
    return true;
}

// Returns the position of the smallest string
inline int which_min_8(const string strings[]){
    int min=0;
    for(int i=1;i<8;i++) if (strings[i].compare(strings[min])<0) min=i;
    return min;
    
}

// Computes longest common prefix between two sequences
inline int lcp(const string& a, const string& b){
    int lcp=0;
    while(true){
        if(a.at(lcp)!=b.at(lcp)) return lcp;
        lcp++;
    }
    
}



// Returns the canonical string of an occurence (4 sequences of size at most 2k-2 each)
// Algo

// 1: case without small repeats :
// writing the canonical version:
// min(au, av', u'a', u'b, b'v', b'u, vb, va') --> indicates which solution we choose a, u', b', or v (liviu = 0, 1, 2, 3)
// 0:0 =    au  ; vb  ; av' ; u'b  (if au==min)
// 1:1 =    u'a'; b'v'; u'b ; av'  (if u'a'==min)
// 2:2 =    vb  ; au  ; va' ; b'u  (if vb==min)
// 3:3 =    b'v'; u'a'; b'u ; va'  (if b'v'==min)
// 4:0bis = av' ; u'b ; au  ; vb   (if av'==min)
// 5:1bis = u'b ; av' ; u'a'; b'v' (if u'b==min)
// 6:2bis = va' ; b'u ; vb  ; au   (if va'==min)
// 7:3bis = b'u ; va' ; b'v'; u'a' (if b'u==min)
// code: create the 8 sequences, 2k-2 sequences : removes first and last character, and writes the smallest

// 2: general case, with repeats. We need to keep only the subsequences common to all 4 possibilities of detection (starting node)
// x=lcp(a',b)
// u=u[0,|u|-x] // (end of a = beginning of rev comp of b)
// v=v[x,|v|]
Solution TakeABreak::get_canonicalSolution(const Node& a, const Node& u, const Node& v, const Node& b){
    string strings[8];
    
    //cout << "a=" << _graph.toString(a) << " u=" << _graph.toString(u) << " v=" << _graph.toString(v) << " b=" << _graph.toString(b) << endl;
    
    Solution canon;
    
    Node abar=_graph.reverse(a);
    Node bbar=_graph.reverse(b);
    Node ubar=_graph.reverse(u);
    Node vbar=_graph.reverse(v);
    
    
    int x=lcp(_graph.toString(abar),_graph.toString(b));
    
    string new_u_str = _graph.toString(u).substr(0, _kmerSize-x);
    string new_ubar_str = _graph.toString(ubar).substr(x, _kmerSize-x);
    string new_v_str = _graph.toString(v).substr(x, _kmerSize-x);
    string new_vbar_str = _graph.toString(vbar).substr(0, _kmerSize-x);
    
    strings[0]=_graph.toString(a)+new_u_str; // 0
    strings[0]=strings[0].substr(1,strings[0].length()-2);
    
    strings[1]=new_ubar_str+_graph.toString(abar); // 1
    strings[1]=strings[1].substr(1,strings[1].length()-2);
    
    strings[2]=new_v_str+_graph.toString(b); // 2
    strings[2]=strings[2].substr(1,strings[2].length()-2);
    
    strings[3]=_graph.toString(bbar)+new_vbar_str; // 3
    strings[3]=strings[3].substr(1,strings[3].length()-2);
    
    strings[4]=_graph.toString(a)+new_vbar_str; // 0bis
    strings[4]=strings[4].substr(1,strings[4].length()-2);
    
    strings[5]=new_ubar_str+_graph.toString(b); // 1bis
    strings[5]=strings[5].substr(1,strings[5].length()-2);
    
    strings[6]=new_v_str+_graph.toString(abar); // 2bis
    strings[6]=strings[6].substr(1,strings[6].length()-2);
    
    strings[7]=_graph.toString(bbar)+new_u_str; // 3bis
    strings[7]=strings[7].substr(1,strings[7].length()-2);
    
    int min = which_min_8(strings);
    
    switch (min) {
        case 0:
            canon.setSequences(strings[0],strings[2]);
            break;
            
        case 4:
            canon.setSequences(strings[4],strings[5]);
            break;
#ifndef only_canonical
        case 1:
            canon.setSequences(strings[1],strings[3]);
            break;
            
        case 2:
            canon.setSequences(strings[2],strings[0]);
            break;
            
        case 3:
            canon.setSequences(strings[3],strings[1]);
            break;
            
        case 5:
            canon.setSequences(strings[5],strings[4]);
            break;
            
        case 6:
            canon.setSequences(strings[6],strings[7]);
            break;
            
        case 7:
            canon.setSequences(strings[7],strings[6]);
            break;
#endif
    }
    //cout << "canon=" << canon._auvb << endl;
    return canon;
}


size_t TakeABreak::writeResults(FILE * out)
{
    size_t count=0;
    while (!_solutions.empty()) {
        Solution sol=*_solutions.begin();
        size_t res=sol.writeFastaBreakpoints(out,count,_kmerSize);
        _solutions.erase(_solutions.begin());
        count=count+res;
    }
    return count;
}

// to improve : duplication of code
// Function used only if PRINTALL=1
size_t TakeABreak::writeUntruncResults(FILE * out)
{
    size_t count=0;
    while (!_untruncSolutions.empty()) {
        Solution sol=*_untruncSolutions.begin();
        size_t res=sol.writeFastaBreakpoints(out,count,_kmerSize);
        //size_t res=sol.writeFastaNodes(out,count); // here each fasta entry is the sequence of a node, not the concatenation as a breakpoint sequence
        _untruncSolutions.erase(_untruncSolutions.begin());
        count=count+res;
    }
    return count;
}

// Same as get_canonicalSolution except that nodes are not truncated if there is a repeated sequence
// Enables to see all breakpoint solutions (when there are variations at the extremities of the breakpoints)
// BUT : solutions with a repeat are very likely to be output twice, depending if the repeat is in a/b' or in u/v'
// Only used when PRINTALL=1
Solution TakeABreak::get_untruncatedCanonicalSolution(const Node& a, const Node& u, const Node& v, const Node& b){
    string strings[8];
    
    //cout << "a=" << _graph.toString(a) << " u=" << _graph.toString(u) << " v=" << _graph.toString(v) << " b=" << _graph.toString(b) << endl;
    
    Solution canon;
    
    Node abar=_graph.reverse(a);
    Node bbar=_graph.reverse(b);
    Node ubar=_graph.reverse(u);
    Node vbar=_graph.reverse(v);
    
    
//    int x=lcp(_graph.toString(abar),_graph.toString(b));
//    string new_u_str = _graph.toString(u).substr(0, _kmerSize-x);
//    string new_ubar_str = _graph.toString(ubar).substr(x, _kmerSize-x);
//    string new_v_str = _graph.toString(v).substr(x, _kmerSize-x);
//    string new_vbar_str = _graph.toString(vbar).substr(0, _kmerSize-x);
    
    strings[0]=_graph.toString(a)+_graph.toString(u); // 0
//    strings[0]=strings[0].substr(1,strings[0].length()-2);
    
    strings[1]=_graph.toString(ubar)+_graph.toString(abar); // 1
//    strings[1]=strings[1].substr(1,strings[1].length()-2);
    
    strings[2]=_graph.toString(v)+_graph.toString(b); // 2
//    strings[2]=strings[2].substr(1,strings[2].length()-2);
    
    strings[3]=_graph.toString(bbar)+_graph.toString(vbar); // 3
//    strings[3]=strings[3].substr(1,strings[3].length()-2);
    
    strings[4]=_graph.toString(a)+_graph.toString(vbar); // 0bis
//    strings[4]=strings[4].substr(1,strings[4].length()-2);
    
    strings[5]=_graph.toString(ubar)+_graph.toString(b); // 1bis
//    strings[5]=strings[5].substr(1,strings[5].length()-2);
    
    strings[6]=_graph.toString(v)+_graph.toString(abar); // 2bis
//    strings[6]=strings[6].substr(1,strings[6].length()-2);
    
    strings[7]=_graph.toString(bbar)+_graph.toString(u); // 3bis
//    strings[7]=strings[7].substr(1,strings[7].length()-2);
    
    int min = which_min_8(strings);
    
    switch (min) {
        case 0:
            canon.setSequences(strings[0],strings[2]);
            break;
            
        case 4:
            canon.setSequences(strings[4],strings[5]);
            break;
        case 1:
            canon.setSequences(strings[1],strings[3]);
            break;
            
        case 2:
            canon.setSequences(strings[2],strings[0]);
            break;
            
        case 3:
            canon.setSequences(strings[3],strings[1]);
            break;
            
        case 5:
            canon.setSequences(strings[5],strings[4]);
            break;
            
        case 6:
            canon.setSequences(strings[6],strings[7]);
            break;
            
        case 7:
            canon.setSequences(strings[7],strings[6]);
            break;
    }
    //cout << "canon=" << canon._auvb << endl;
    return canon;
}


/********************************************************************************/
//For each branching kmer, on each strand -> a:
//  If a->out_neighbor <2: continue (avoids dead end "branching kmers")
//  For each (a+1) (neigbor de a)
//    Ualpha = (k-1)-extension of (a+1), without any constraint
//    For each master_group from 0 to a->out_neighbor-2
//      For each u in Umaster_group:
//          B = give_me_B(a,u,_tolerance_rc)
//          For each b of B:
//             For each alpha > master_group // groups for which we have not yet computed the b's
//                For each vbar in Ualpha:
//                   If check_path(v,b):
//                      output (auvb)
void TakeABreak::MainLoopFunctor::operator() (Node& nodeA)
{
    LCS& lcs_instance = lcs();

//cout << "node " << ref._graph.toString(nodeA) << "  outdeg=" << ref._graph.outdegree (nodeA) << "  thread=" << System::thread().getThreadSelf() << endl;
//cout.flush();

#ifdef check_memory
    if (KERN_SUCCESS != task_info(mach_task_self(),
            TASK_BASIC_INFO, (task_info_t)&t_info,
            &t_info_count))
    {
        exit(1);
    }
    if(t_info.resident_size>max_mem) max_mem=t_info.resident_size;
#endif

    //if(!conserve_node(_graph.toString(nodeA).c_str(),shannon_limit)) continue;

    /** We iterate the strands */
    for(int strand=0;strand<2;strand++,  nodeA=ref._graph.reverse(nodeA)){ // each strand
#ifdef debug
        cout<<"strand = "<<strand<<" Node A = "<<ref._graph.toString(nodeA)<<" -- out degree ="<<ref._graph.outdegree (nodeA)<<endl;
#endif
        // test in nodeA is branching
        if(ref._graph.outdegree (nodeA)<2) {
            continue;
        }
        // get all immediate out-neighbors of node A
        Graph::Vector<Node> neighborsA = ref._graph.neighbors<Node> (nodeA, DIR_OUTCOMING);

#ifdef debug
        cout<<neighborsA.size()<<" immediate neighbor"<<endl;
#endif

        // will store a DBGwalker containing all (k-1) neighbors for each immediate neighbor (alpha) (U sets)
        vector<DBGWalker> Ualpha (neighborsA.size());

        // retreiving the U's for each alpha
        size_t maxUSize = 0; // remind the max U size encountered so far, to limit the recursive search with the LCT parameter
        for(int alpha=0;alpha<neighborsA.size();alpha++)
        {
            Ualpha[alpha].setInfo (ref._graph, ref._LCT); // initialize the DBGwalker object

            Ualpha[alpha].find_all_at_depth(neighborsA[alpha], ref._graph.getKmerSize()-1, maxUSize);

            if (maxUSize < Ualpha[alpha].size())  { maxUSize = Ualpha[alpha].size(); }

        }// end retreiving the U's for each alpha

        // sort the U's in increasing size order : heuristic to reduce the loop sizes
        std::sort (Ualpha.begin(), Ualpha.end());


        // checks that at least 2 Ualpha are not empty (the two largest ones)
        if(Ualpha[neighborsA.size()-1].size()==0 || Ualpha[neighborsA.size()-2].size()==0)
        {
            continue;
        }
        // filtering on the cardinality of number of kmers reachable from node A
        if(Ualpha[neighborsA.size()-1].size() * Ualpha[neighborsA.size()-2].size() > ref._LCT)
        {
            continue;
        }

        for(int master_group=0;master_group<neighborsA.size()-1;master_group++){ // choosing the master group
            //cout<<"master group "<<master_group<<endl;
            vector<Node>& current_master_u = Ualpha[master_group].reachable_neighbor; // get the U set of master group
            for(int id_u=0;id_u<Ualpha[master_group].size();id_u++){ // for each u of the Umaster

                DBGWalker B;
                B.setInfo(ref._graph, ref._LCT);

                //if(!conserve_node(_graph.toString(current_master_u[id_u]),shannon_limit)) continue;
                B.find_B(nodeA,current_master_u[id_u],ref._tolerance_rc, lcs_instance, maxUSize); // get the set of B

                // if number of comparisons of size of set B with the largest U is too big, we stop.
                if(maxUSize*B.size()>ref._LCT) continue;

                for(int id_b=0;id_b<B.size();id_b++){ // for each b in B
                    //                        if(!conserve_node(_graph.toString(set_B[id_b]),shannon_limit)) continue;
                    for(int alpha_other_group=master_group+1;alpha_other_group<neighborsA.size();alpha_other_group++){ // for each other group
                        for(int id_v_other_group=0;id_v_other_group<Ualpha[alpha_other_group].size();id_v_other_group++){ // each u (== vbar) in the other groups
                            Node v=ref._graph.reverse(Ualpha[alpha_other_group].reachable_neighbor[id_v_other_group]);
                            //                                if(!conserve_node(_graph.toString(v),shannon_limit)) continue;
#ifdef debug
                            cout<<"checking:"<<endl;
                            cout<<"  "<<_graph.toString(nodeA)<<" shannon "<<shannon_index(_graph.toString(nodeA))<<endl;
                            cout<<"  "<<_graph.toString(current_master_u[id_u])<<" shannon "<<shannon_index(_graph.toString(current_master_u[id_u]))<<endl;
                            cout<<"  "<<_graph.toString(v)<<" shannon "<<shannon_index(_graph.toString(v))<<endl;
                            cout<<"  "<<_graph.toString(set_B[id_b])<<" shannon "<<shannon_index(_graph.toString(set_B[id_b]))<<endl;;
#endif
                            if(lcs_instance.size_lcs(ref._graph.toString(current_master_u[id_u]), ref._graph.toString(ref._graph.reverse(v)))<lcs_instance.lcs_threshold && ref.checkPath(v,B.reachable_neighbor[id_b]))
                                //                                   &&
                                //                                   conserve_inversion(nodeA,current_master_u[id_u], v, set_B[id_b])
                                //                                   &&
                                //                                   check_tolerance(_graph.reverse(nodeA), set_B[id_b], _tolerance_rc) // tolerance of _tolerance_rc characters similar between a and b'
                                //                                   &&
                                //                                   check_tolerance(v, _graph.reverse(current_master_u[id_u]), _tolerance_rc) // tolerance of _tolerance_rc characters similar between v and u'
                            { // found one path
                                //                                    cout<<"I found inversion (a,u,v,b) : ("<<number_substitutions(_graph.toString(nodeA),_graph.toString(_graph.reverse(set_B[id_b])))<<" subst a, rcb) - combinatorial_factor="<<combinatorial_factor<<endl;
                                //                                    cout<<"  "<<_graph.toString(nodeA)<<" shannon "<<shannon_index(_graph.toString(nodeA))<<endl;
                                //                                    cout<<"  "<<_graph.toString(current_master_u[id_u])<<" shannon "<<shannon_index(_graph.toString(current_master_u[id_u]))<<endl;
                                //                                    cout<<"  "<<_graph.toString(v)<<" shannon "<<shannon_index(_graph.toString(v))<<endl;
                                //                                    cout<<"  "<<_graph.toString(set_B[id_b])<<" shannon "<<shannon_index(_graph.toString(set_B[id_b]))<<endl;
                                
                                //ref.print_canonical(nodeA, current_master_u[id_u], v, B.reachable_neighbor[id_b], number_inv_found, out);
                                Solution canon=ref.get_canonicalSolution(nodeA, current_master_u[id_u], v, B.reachable_neighbor[id_b]);
                                if (canon._auvb.size()>0){ // en mode onlycanon on n'insert pas les Solutions vides
                                    //cout << "lcs u-v'" << lcs_instance.size_lcs(ref._graph.toString(current_master_u[id_u]), ref._graph.toString(ref._graph.reverse(v))) << endl;
                                    //cout << "lcs a-b'" << lcs_instance.size_lcs(ref._graph.toString(nodeA), ref._graph.toString(ref._graph.reverse(B.reachable_neighbor[id_b]))) << endl;
                                    ref.getSynchro()->lock();
                                    ref._nbOccurrences=ref._nbOccurrences+1;
                                    ref._solutions.insert(canon);
#ifdef PRINTALL
                                    Solution occ=ref.get_untruncatedCanonicalSolution(nodeA, current_master_u[id_u], v, B.reachable_neighbor[id_b]);
                                    ref._untruncSolutions.insert(occ);
#endif
                                    ref.getSynchro()->unlock();

                                }


                            } // end found one path
                        } // end each u (== vbar) in the other groups
                    } // end each other group
                } // end each b in B
            } // end each u of the Umaster
        } // end choosing the master group
    } // end each strand

}

/********************************************************************************/
// Algorithm is in MainLoopFunctor()
void TakeABreak::find_ALL_occurrences_of_inversion_pattern (LCS& lcsParam)
{
    // dirty:
    init_int2char();
    
    
#ifdef check_memory
    int64_t max_mem=0;
#endif
    
    //int number_inv_found=0;

    /** We use a dispatcher object for parallelization, with a given number of cores to be used. */
    Dispatcher dispatcher (_nbCores);

    /** We define an iterator over the branching nodes of the graph. We use also progress information. */
    ProgressGraphIterator<BranchingNode, ProgressTimerAndSystem> branchingNodes (_graph.iterator<BranchingNode>(), "looping nodes");
    
    /** We use one LCS object per thread => this is done thanks to the ThreadObject class. */
    ThreadObject<LCS> lcs (lcsParam);

    /** Functor called for each branching node. */
    MainLoopFunctor functor (*this, lcs);

    _nbOccurrences =0;
    /** We iterate all the branching nodes of the graph. */
    IDispatcher::Status status = dispatcher.iterate (branchingNodes, functor);
    
    
#ifdef check_memory
    printf("maximal memory used = %llu\n", max_mem);
#endif
    
//    cout<<"\r"<<100*nb_treated/branchingNodes.size()<<" % done ("<<nb_treated<<" branching nodes treated)"<<endl;
}

// Function to check if there is a path of size k in the DBG to get from nodeV to nodeB
// We test if each edge exists
bool TakeABreak::checkPath (Node nodeV, Node nodeB)
{

    bool check=true;
    
    // First getting the string-kmer of target node
    string strB = _graph.toString(nodeB);
    // cout << "node B to string =" << strB << endl;
    
    // Then we will walk the graph from nodeV towards nodeB, constraining the edges to the correct nucleotides
    Node currentNode=nodeV;
    int i=0;
    while (check && i<_kmerSize){
    //for(int i=0;i<_kmerSize && check;i++){
        check=false;
        // We retrieve all outcoming edges from current nodes
        Graph::Vector<Edge> edges = _graph.neighbors<Edge> (currentNode, DIR_OUTCOMING);
        int j=0;
        // cout << "step i=" << i << endl;
        
        // We check among the outcoming nodes if one has the correct nucleotide (ith nt of strB)
        while (!check && j<edges.size()){
        //for(int j=0;j<edges.size() && !check;j++){
            // cout << _graph.debugString(edges[j]) << endl;
            if(ascii(edges[j].nt) == strB[i]){
                // correct nucleotide, we change the current node and go to the next position i
                //cout << "nt ok=" << ascii(edges[j].nt) << endl;
                currentNode = edges[j].to;
                check=true;
            }
            j++;
        }
        i++;
    }
    
    return check;
}


/********************************************************************************/
// Old functions (no longer used)


void TakeABreak::printParameters(FILE * log){
    stringstream allParams;
    allParams << "*******Parameter values*******" <<endl;
    //allParams << "Output prefix: " << getInput()->getStr(STR_URI_OUTPUT) << endl;
    allParams << "de Bruijn Graph: " << endl;
    allParams << "\tk=" << _kmerSize << endl;
    allParams << "\tmin abundance=" << _graph.getInfo().getStr(ATTR_KMER_ABUNDANCE) << endl;
    
    allParams <<  "Inversion detection:" << endl;
	allParams <<  "\trepeat=" << _tolerance_rc << endl;
	allParams <<  "\tLCT=" << _LCT << endl;
	allParams <<  "\tmax_sim="<< _max_sim << endl;
    allParams <<  "*****************************" << endl;
    
    fprintf(log,"%s",allParams.str().c_str());
    
}


// prints directly in a file 4 sequences of size at most 2k-2 each
// no longer used, see Solution TakeABreak::get_canonical(...)

// 1: case without small repeats :
// writing the canonical version:
// min(au, av', u'a', u'b, b'v', b'u, vb, va') --> indicates which solution we choose a, u', b', or v (liviu = 0, 1, 2, 3)
// 0:0 =    au  ; vb  ; av' ; u'b  (if au==min)
// 1:1 =    u'a'; b'v'; u'b ; av'  (if u'a'==min)
// 2:2 =    vb  ; au  ; va' ; b'u  (if vb==min)
// 3:3 =    b'v'; u'a'; b'u ; va'  (if b'v'==min)
// 4:0bis = av' ; u'b ; au  ; vb   (if av'==min)
// 5:1bis = u'b ; av' ; u'a'; b'v' (if u'b==min)
// 6:2bis = va' ; b'u ; vb  ; au   (if va'==min)
// 7:3bis = b'u ; va' ; b'v'; u'a' (if b'u==min)
// code: create the 8 sequences, 2k-2 sequences : removes first and last character, and writes the smallest

// 2: general case, with repeats. We need to keep only the subsequences common to all 4 possibilities of detection (starting node)
// x=lcp(a',b)
// u=u[0,|u|-x] // (end of a = beginning of rev comp of b)
// v=v[x,|v|]
bool TakeABreak::print_canonical(const Node& a, const Node& u, const Node& v, const Node& b, int& number_inv_found, FILE * out){
    string strings[8];
    bool output = false;
    
    LocalSynchronizer local(synchro);
    
    Node abar=_graph.reverse(a);
    Node bbar=_graph.reverse(b);
    Node ubar=_graph.reverse(u);
    Node vbar=_graph.reverse(v);
    
    
    int x=lcp(_graph.toString(abar),_graph.toString(b));
    
    
    string new_u_str = _graph.toString(u).substr(0, _kmerSize-x);
    string new_ubar_str = _graph.toString(ubar).substr(x, _kmerSize-x);
    string new_v_str = _graph.toString(v).substr(x, _kmerSize-x);
    string new_vbar_str = _graph.toString(vbar).substr(0, _kmerSize-x);
    
    strings[0]=_graph.toString(a)+new_u_str; // 0
    strings[0]=strings[0].substr(1,strings[0].length()-2);
    
    strings[1]=new_ubar_str+_graph.toString(abar); // 1
    strings[1]=strings[1].substr(1,strings[1].length()-2);
    
    strings[2]=new_v_str+_graph.toString(b); // 2
    strings[2]=strings[2].substr(1,strings[2].length()-2);
    
    strings[3]=_graph.toString(bbar)+new_vbar_str; // 3
    strings[3]=strings[3].substr(1,strings[3].length()-2);
    
    strings[4]=_graph.toString(a)+new_vbar_str; // 0bis
    strings[4]=strings[4].substr(1,strings[4].length()-2);
    
    strings[5]=new_ubar_str+_graph.toString(b); // 1bis
    strings[5]=strings[5].substr(1,strings[5].length()-2);
    
    strings[6]=new_v_str+_graph.toString(abar); // 2bis
    strings[6]=strings[6].substr(1,strings[6].length()-2);
    
    strings[7]=_graph.toString(bbar)+new_u_str; // 3bis
    strings[7]=strings[7].substr(1,strings[7].length()-2);
    
    int min = which_min_8(strings);
    
    //    cout<<" min = "<<min<<" (x= "<<x<<")"<<endl;
    //    cout<<strings[0]<<endl;
    //        cout<<strings[1]<<endl;
    //        cout<<strings[2]<<endl;
    //        cout<<strings[3]<<endl;
    //        cout<<strings[4]<<endl;
    //        cout<<strings[5]<<endl;
    //        cout<<strings[6]<<endl;
    //        cout<<strings[7]<<endl;
    switch (min) {
        case 0:
            // fasta mode (not removing first and last letters):
            //            fprintf(out, ">inv_%d au\n%s\n>vb\n%s\n>av'\n%s\n>u'b\n%s\n", number_inv_found, strings[0].c_str(), strings[2].c_str(), strings[4].c_str(), strings[5].c_str());
            
            // r2sv mode: (allows the sort -u)
            // Prints au vb av' u'b
            // 1/ removing the first character of au and av' and last character of vb and u'b
            // 2/ only if x==0, removing the last character of au and av' and the first character of vb and u'b (else, if x>0, it was already removed)
            //            fprintf(out, "au;%s;>vb;%s;>av';%s;>u'b;%s\n",
            //                    x>0?strings[0].substr(1).c_str():strings[0].substr(1,strings[0].length()-1).c_str(), // au, removing first character and maybe last
            //                    x>0?strings[2].substr(0,strings[2].length()-1).c_str():strings[2].substr(1,strings[2].length()-1).c_str(), // vb, maybe removeing first and removeing last
            //                    x>0?strings[4].substr(1).c_str():strings[4].substr(1,strings[4].length()-1).c_str(), // av' removing first character and maybe last
            //                    x>0?strings[5].substr(0,strings[5].length()-1).c_str():strings[5].substr(1,strings[5].length()-1).c_str()); // vb, maybe removeing first and removeing last
            
            
            // r2sv mode: (allows the sort -u)
            // Prints au vb av' u'b
            // removing the first and last character of au and av' and first and last character of vb and u'b
            fprintf(out, "au;%s;>vb;%s;>av';%s;>u'b;%s\n",
                    strings[0].c_str(),  // au removing first and  last
                    strings[2].c_str(),  // vb removing first and  last
                    strings[4].c_str(),  // av' removing first and  last
                    strings[5].c_str()); // v'b, removing first and  last
            
            output=true;
            break;
            
        case 4:
            
            fprintf(out, "au;%s;>vb;%s;>av';%s;>u'b;%s\n",
                    strings[4].c_str(),  // av' removing first and  last
                    strings[5].c_str(),  // vb  removing first and  last
                    strings[0].c_str(),  // au  removing first and  last
                    strings[2].c_str()); // vb  removing first and  last
            
            output=true;
            break;
            
            // Commented all this code, that is useless.
#ifndef only_canonical
        case 1:
            fprintf(out, "au;%s;>vb;%s;>av';%s;>u'b;%s\n",
                    strings[1].c_str(),  // av' removing first and  last
                    strings[3].c_str(),  // vb  removing first and  last
                    strings[5].c_str(),  // au  removing first and  last
                    strings[4].c_str()); // vb  removing first and  last
            
            output=true;
            break;
            
        case 2:
            fprintf(out, "au;%s;>vb;%s;>av';%s;>u'b;%s\n",
                    strings[2].c_str(),  // av' removing first and  last
                    strings[0].c_str(),  // vb  removing first and  last
                    strings[6].c_str(),  // au  removing first and  last
                    strings[7].c_str()); // vb  removing first and  last
            
            output=true;
            break;
            
            
        case 3:
            fprintf(out, "au;%s;>vb;%s;>av';%s;>u'b;%s\n",
                    strings[3].c_str(),  // av' removing first and  last
                    strings[1].c_str(),  // vb  removing first and  last
                    strings[7].c_str(),  // au  removing first and  last
                    strings[6].c_str()); // vb  removing first and  last
            
            output=true;
            break;
            
        case 5:
            fprintf(out, "au;%s;>vb;%s;>av';%s;>u'b;%s\n",
                    strings[5].c_str(),  // av' removing first and  last
                    strings[4].c_str(),  // vb  removing first and  last
                    strings[1].c_str(),  // au  removing first and  last
                    strings[3].c_str()); // vb  removing first and  last
            
            output=true;
            break;
            
        case 6:
            fprintf(out, "au;%s;>vb;%s;>av';%s;>u'b;%s\n",
                    strings[6].c_str(),  // av' removing first and  last
                    strings[7].c_str(),  // vb  removing first and  last
                    strings[2].c_str(),  // au  removing first and  last
                    strings[0].c_str()); // vb  removing first and  last
            
            output=true;
            break;
            
        case 7:
            fprintf(out, "au;%s;>vb;%s;>av';%s;>u'b;%s\n",
                    strings[7].c_str(),  // av' removing first and  last
                    strings[6].c_str(),  // vb  removing first and  last
                    strings[3].c_str(),  // au  removing first and  last
                    strings[1].c_str()); // vb  removing first and  last
            
            output=true;
            break;
            
#endif //not only_canonical
    }
    
    //        case 1:
    //            fprintf(out, ">inv_%d au\n%s\n>vb\n%s\n>av'\n%s\n>u'b\n%s\n", number_inv_found, strings[1].c_str(), strings[3].c_str(), strings[5].c_str(), strings[4].c_str());
    //
    ////            cout<<strings[1]<<" "<<strings[3]<<" "<<strings[5]<<" "<<strings[4]<<endl;
    //            output=true;
    //            break;
    //        case 2:
    //            fprintf(out, ">inv_%d au\n%s\n>vb\n%s\n>av'\n%s\n>u'b\n%s\n", number_inv_found, strings[2].c_str(), strings[0].c_str(), strings[6].c_str(), strings[7].c_str());
    //
    ////            cout<<strings[2]<<" "<<strings[0]<<" "<<strings[6]<<" "<<strings[7]<<endl;
    //            output=true;
    //            break;
    //        case 3:
    //            fprintf(out, ">inv_%d au\n%s\n>vb\n%s\n>av'\n%s\n>u'b\n%s\n", number_inv_found, strings[3].c_str(), strings[1].c_str(), strings[7].c_str(), strings[6].c_str());
    //
    //
    //
    ////            cout<<strings[3]<<" "<<strings[1]<<" "<<strings[7]<<" "<<strings[6]<<endl;
    //            output=true;
    //            break;
    //
    //        case 5:
    //            fprintf(out, ">inv_%d au\n%s\n>vb\n%s\n>av'\n%s\n>u'b\n%s\n", number_inv_found, strings[5].c_str(), strings[4].c_str(), strings[1].c_str(), strings[3].c_str());
    //
    ////            cout<<strings[5]<<" "<<strings[4]<<" "<<strings[1]<<" "<<strings[3]<<endl;
    //            output=true;
    //            break;
    //        case 6:
    //            fprintf(out, ">inv_%d au\n%s\n>vb\n%s\n>av'\n%s\n>u'b\n%s\n", number_inv_found, strings[6].c_str(), strings[7].c_str(), strings[2].c_str(), strings[0].c_str());
    //
    ////            cout<<strings[6], strings[7]<<" "<<strings[2]<<" "<<strings[0]<<endl;
    //            output=true;
    //            break;
    //        case 7:
    //            fprintf(out, ">inv_%d au\n%s\n>vb\n%s\n>av'\n%s\n>u'b\n%s\n", number_inv_found, strings[7].c_str(), strings[6].c_str(), strings[3].c_str(), strings[1].c_str());
    //
    ////            cout<<strings[7], strings[6]<<" "<<strings[3]<<" "<<strings[1]<<endl;
    //            output=true;
    //            break;        default:
    //            break;
    
    if(output){
        number_inv_found++;
        //        delete[] strings;
        return true;
    }
    return false;
}

// Returns the canonical string of an occurence (4 sequences of size at most 2k-2 each)
// no longer used, see Solution TakeABreak::get_canonical(...)
string TakeABreak::get_canonical(const Node& a, const Node& u, const Node& v, const Node& b){
    string strings[8];
    char canon[512];
    //LocalSynchronizer local(synchro);
    
    Node abar=_graph.reverse(a);
    Node bbar=_graph.reverse(b);
    Node ubar=_graph.reverse(u);
    Node vbar=_graph.reverse(v);
    
    
    int x=lcp(_graph.toString(abar),_graph.toString(b));
    
    
    string new_u_str = _graph.toString(u).substr(0, _kmerSize-x);
    string new_ubar_str = _graph.toString(ubar).substr(x, _kmerSize-x);
    string new_v_str = _graph.toString(v).substr(x, _kmerSize-x);
    string new_vbar_str = _graph.toString(vbar).substr(0, _kmerSize-x);
    
    strings[0]=_graph.toString(a)+new_u_str; // 0
    strings[0]=strings[0].substr(1,strings[0].length()-2);
    
    strings[1]=new_ubar_str+_graph.toString(abar); // 1
    strings[1]=strings[1].substr(1,strings[1].length()-2);
    
    strings[2]=new_v_str+_graph.toString(b); // 2
    strings[2]=strings[2].substr(1,strings[2].length()-2);
    
    strings[3]=_graph.toString(bbar)+new_vbar_str; // 3
    strings[3]=strings[3].substr(1,strings[3].length()-2);
    
    strings[4]=_graph.toString(a)+new_vbar_str; // 0bis
    strings[4]=strings[4].substr(1,strings[4].length()-2);
    
    strings[5]=new_ubar_str+_graph.toString(b); // 1bis
    strings[5]=strings[5].substr(1,strings[5].length()-2);
    
    strings[6]=new_v_str+_graph.toString(abar); // 2bis
    strings[6]=strings[6].substr(1,strings[6].length()-2);
    
    strings[7]=_graph.toString(bbar)+new_u_str; // 3bis
    strings[7]=strings[7].substr(1,strings[7].length()-2);
    
    int min = which_min_8(strings);
    
    switch (min) {
        case 0:
            sprintf(canon,"au;%s;>vb;%s;>av';%s;>u'b;%s\n",
                    strings[0].c_str(),  // au removing first and  last
                    strings[2].c_str(),  // vb removing first and  last
                    strings[4].c_str(),  // av' removing first and  last
                    strings[5].c_str()); // v'b, removing first and  last
            break;
            
        case 4:
            
            sprintf(canon,"au;%s;>vb;%s;>av';%s;>u'b;%s\n",
                    strings[4].c_str(),  // av' removing first and  last
                    strings[5].c_str(),  // vb  removing first and  last
                    strings[0].c_str(),  // au  removing first and  last
                    strings[2].c_str()); // vb  removing first and  last
            break;
            
        case 1:
            sprintf(canon,"au;%s;>vb;%s;>av';%s;>u'b;%s\n",
                    strings[1].c_str(),  // av' removing first and  last
                    strings[3].c_str(),  // vb  removing first and  last
                    strings[5].c_str(),  // au  removing first and  last
                    strings[4].c_str()); // vb  removing first and  last
            break;
            
        case 2:
            sprintf(canon,"au;%s;>vb;%s;>av';%s;>u'b;%s\n",
                    strings[2].c_str(),  // av' removing first and  last
                    strings[0].c_str(),  // vb  removing first and  last
                    strings[6].c_str(),  // au  removing first and  last
                    strings[7].c_str()); // vb  removing first and  last
            break;
            
            
        case 3:
            sprintf(canon,"au;%s;>vb;%s;>av';%s;>u'b;%s\n",
                    strings[3].c_str(),  // av' removing first and  last
                    strings[1].c_str(),  // vb  removing first and  last
                    strings[7].c_str(),  // au  removing first and  last
                    strings[6].c_str()); // vb  removing first and  last
            break;
            
        case 5:
            sprintf(canon,"au;%s;>vb;%s;>av';%s;>u'b;%s\n",
                    strings[5].c_str(),  // av' removing first and  last
                    strings[4].c_str(),  // vb  removing first and  last
                    strings[1].c_str(),  // au  removing first and  last
                    strings[3].c_str()); // vb  removing first and  last
            break;
            
        case 6:
            sprintf(canon,"au;%s;>vb;%s;>av';%s;>u'b;%s\n",
                    strings[6].c_str(),  // av' removing first and  last
                    strings[7].c_str(),  // vb  removing first and  last
                    strings[2].c_str(),  // au  removing first and  last
                    strings[0].c_str()); // vb  removing first and  last
            break;
            
        case 7:
            sprintf(canon,"au;%s;>vb;%s;>av';%s;>u'b;%s\n",
                    strings[7].c_str(),  // av' removing first and  last
                    strings[6].c_str(),  // vb  removing first and  last
                    strings[3].c_str(),  // au  removing first and  last
                    strings[1].c_str()); // vb  removing first and  last
            break;
            
    }
    
    return canon;
}

