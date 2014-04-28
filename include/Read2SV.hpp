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



#ifndef _READ2SV_HPP_
#define _READ2SV_HPP_

/********************************************************************************/

#include <gatb/gatb_core.hpp>

#include <stdio.h>
#include <string>
#include <sstream>
#include <LCS.hpp>

/********************************************************************************/
/* Class Read2SV*/
/********************************************************************************/

class Read2SV : public Tool
{
//private:
public:


//    static const char* STRING_KMER_SIZE;
//    static const char* STRING_URI_SOLID_KMERS;
//    static const char* STRING_KMER_CFP;
    size_t          _kmerSize;
    Graph _graph;
    int size_tolerance_rc;
    float shannon_limit;
    size_t nbCores;
	
public:

    /** */
    Read2SV (char* graphFile, const int size_tolerance_rc);
    virtual ~Read2SV();
    void find_ALL_occurrences_of_inversion_pattern (LCS& lcs_instance, FILE * out, const int& local_complexity_threshold);
    bool checkPath (Node nodeV, Node nodeB);
    bool conserve_inversion(const Node& a, const Node& u, const Node& v, const Node& b);
    // the two sequences s1 and s2 are distinct enought (return true) if (s1[0:size_tolerance_rc] !=  s2[0:size_tolerance_rc])
    bool check_tolerance(const Node& s1, const Node& s2, const int size_tolerance_rc);
	
	
    bool print_canonical(const Node& a, const Node& u, const Node& v, const Node& b, int& number_inv_found, FILE * out);
private:
    
    /** */
    void execute ();
    ISynchronizer* synchro;
    
    
};

/********************************************************************************/

#endif /* _READ2SV_HPP_ */

