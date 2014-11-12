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



#ifndef TakeABreak_Solution_hpp
#define TakeABreak_Solution_hpp

#include <gatb/gatb_core.hpp>

#include <string>
#include <sstream>


/********************************************************************************/
/* Class Solution*/
/********************************************************************************/
class Solution{
public:
    string _auvb;
    Solution();
    ~Solution();
    
    void setSequences(string au,string vb);
    
    // Prints the solution in the output file as a multi-fasta
    size_t writeFastaOutput(FILE * out, size_t id);
    
    // necessary to put several Solution objects in a set
    bool operator< (const Solution& other) const;
    
private:
    // Computes the reverse complement of a string
    string reverseComplement(const string& dna);
};
#endif
