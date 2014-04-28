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

LCS::LCS(const LCS& lcs){
    _kmerSize = lcs._kmerSize;
    lcs_threshold = lcs.lcs_threshold;
    opt = (int **) malloc(sizeof(int *)*(_kmerSize+1));
    if(opt == NULL) {fprintf(stderr,"cannot allocate memory for LCS computation\n"); exit(1);}
    for(int i=0;i<_kmerSize+1;i++) {
        opt[i]=(int *)malloc(sizeof(int)*(_kmerSize+1));
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
