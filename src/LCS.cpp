//
//  LCS.cpp
//  read2SV
//
//  Created by Pierre Peterlongo on 14/01/2014.
//  Copyright (c) 2014 Pierre Peterlongo. All rights reserved.
//

#include "LCS.hpp"


// We use the required packages

using namespace std;

LCS::LCS(const size_t& k, const int& max_percentage){
    _kmerSize = k;
//    lcs_threshold=(2*_kmerSize/3);
    lcs_threshold=_kmerSize*max_percentage/100;
    opt = (int **) malloc(sizeof(int *)*(k+1));
    if(opt == NULL) {fprintf(stderr,"cannot allocate memory for LCS computation\n"); exit(1);}
    for(int i=0;i<k+1;i++) {
        opt[i]=(int *)malloc(sizeof(int)*(k+1));
        if(opt[i] == NULL) {fprintf(stderr,"cannot allocate memory for LCS computation (line %d)\n",i); exit(1);}
    }
}

LCS::~LCS(){
    for(int i=0;i<_kmerSize+1;i++) free(opt[i]);
    free(opt);
}

void LCS::reinit() {
//    memset(opt, 0, sizeof(int) * (_kmerSize) * (_kmerSize));
      for (int i = 0; i <=_kmerSize; i++){
        opt[_kmerSize][i]=0;
        opt[i][_kmerSize]=0;
    }
}

const int LCS::size_lcs(const string& a, const string& b){
    reinit();
    // compute length of LCS and all subproblems via dynamic programming
    for (int i = _kmerSize-1; i > -1; i--) {
        for (int j = _kmerSize-1; j > -1; j--) {
            if (a[i] == b[j])
                opt[i][j] = opt[i+1][j+1] + 1;
            else
                opt[i][j] = max(opt[i+1][j], opt[i][j+1]);
        }
    }
    return opt[0][0];

}