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



#ifndef _TAKEABREAK_HPP_
#define _TAKEABREAK_HPP_

/********************************************************************************/

#include <gatb/gatb_core.hpp>

#include <stdio.h>
#include <string>
#include <sstream>
#include <LCS.hpp>
#include <Solution.hpp>
#include <set>
using namespace std;

static const char* STR_TOLERANCE_RC = "-repeat";
static const char* STR_LCT = "-lct";
static const char* STR_MAX_SIM = "-max-sim";


/********************************************************************************/
/* Class TakeABreak*/
/********************************************************************************/

class TakeABreak : public Tool
{
//private:
public:

    
    size_t _kmerSize;
    Graph _graph;
    int _tolerance_rc;
    int _max_sim;
    int _LCT;
    float shannon_limit;
    size_t _nbCores;
    string _output_file;
    //string _log_file;
    set<Solution> _solutions;
    set<Solution> _untruncSolutions; // only when PRINTALL is activated
    size_t _nbOccurrences;
	
public:

    /** */
    TakeABreak ();
    virtual ~TakeABreak();
    
    void printParameters(FILE * log); // no longer used
    
    /** fills getInfo() with parameters informations
     */
    void resumeParameters();
    /** fills getInfo() with results informations
     */
    void resumeResults(size_t number_inv_found, double seconds);
    
    
    /** Search for all occurrences of inversion by testing all starting branching nodes
     */
    void find_ALL_occurrences_of_inversion_pattern (LCS& lcs_instance);
    
    /** Returns true if there exists a path in the graph from node v to node b
     */
    bool checkPath (Node nodeV, Node nodeB);
    
    bool conserve_inversion(Node& a, Node& u, Node& v, Node& b);
    // the two sequences s1 and s2 are distinct enough (return true) if (s1[0:size_tolerance_rc] !=  s2[0:size_tolerance_rc])
    bool check_tolerance(Node& s1, Node& s2, const int size_tolerance_rc);
	
	/** Returns the canonical solution of the occurrence composed of the four nodes a, u, v and b
     */
    Solution get_canonicalSolution(Node& a, Node& u, Node& v, Node& b);
    /** Writes in the file the set of canonical solutions
     */
    size_t writeResults(FILE * out);
    
    /** Returns the canonical -untruncated- solution of the occurrence composed of the four nodes a, u, v and b
     */
    Solution get_untruncatedCanonicalSolution(Node& a, Node& u, Node& v, Node& b);
    /** Writes in the file the set of untruncated canonical solutions
     */
    size_t writeUntruncResults(FILE * out);
    
    //no longer used
    bool print_canonical(Node& a, Node& u, Node& v, Node& b, int& number_inv_found, FILE * out);
    string get_canonical(Node& a, Node& u, Node& v, Node& b);

    
    ISynchronizer* getSynchro(){return synchro;}
    
private:
    
    /** */
    void execute ();
    void displayVersion(std::ostream& os);
    ISynchronizer* synchro;
    
    /** Functor called for each node.
     * NOTE: Previous version used lambda expressions (which made the code clearer) BUT old compilers may not support it,
     * so we go back to the old way with a functor coded as a struct. */
    struct MainLoopFunctor
    {
        MainLoopFunctor (TakeABreak& ref, ThreadObject<LCS>& lcs)
            : ref(ref), lcs(lcs)  {}

        void operator() (Node& node);

        TakeABreak& ref;
        ThreadObject<LCS>& lcs;
    };
    
    friend struct MainLoopFunctor;
};

/********************************************************************************/

#endif /* _TAKEABREAK_HPP_ */

