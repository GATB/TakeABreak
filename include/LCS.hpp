//
//  LCS.hpp
//  read2SV
//
//  Created by Pierre Peterlongo on 14/01/2014.
//  Copyright (c) 2014 Pierre Peterlongo. All rights reserved.
//

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
