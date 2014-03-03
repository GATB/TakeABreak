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



#ifndef read2SV_LCS_hpp
#define read2SV_LCS_hpp

#include <gatb/gatb_core.hpp>

#include <string>
#include <sstream>


/********************************************************************************/
/* Class LCS*/
/********************************************************************************/
class LCS{
public:
    size_t  _kmerSize;
    int lcs_threshold;
    int** opt;
    LCS(const size_t& k, const int& max_percentage);
    ~LCS();
    
    // Computes the lcs of two objects of size k
    const int size_lcs(const std::string& a, const std::string& b);
private:
    void reinit();
};
#endif
