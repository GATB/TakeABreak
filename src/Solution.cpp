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


#include <iostream>
#include "Solution.hpp"


// We use the required packages

using namespace std;


Solution::Solution(){
    
}

void Solution::setSequences(string au,string vb){
    _auvb=au+vb;
}


Solution::~Solution(){
}

bool Solution::operator< (const Solution& other) const  {
    
    return _auvb<other._auvb;
}

size_t Solution::writeFastaBreakpoints(FILE * out, size_t id) {
    if(_auvb.size()>0){
        size_t l=_auvb.size(); // should be even
        size_t l2=l/2;
        size_t lmin=l2/2;
        size_t lmax=l2-lmin;
        // the size of breakpoint sequences can be odd
        //keeping a and b of length lmax, and u and v of length lmin
        string a=_auvb.substr(0,lmax);
        string u=_auvb.substr(lmax,lmin);
        string v=_auvb.substr(l2,lmin);
        string b=_auvb.substr(l2+lmin,lmax);
        stringstream res;
        res <<">inv_"<<id<<" au\n"<<a<<u<<"\n";
        res <<">vb\n"<<v<<b<<"\n";
        res <<">av'\n"<<a<<reverseComplement(v)<<"\n";
        res <<">u'b\n"<<reverseComplement(u)<<b<<"\n";
        
        fprintf(out,"%s",res.str().c_str());
        return 1;
    }
    else{
        return 0;
    }
}

size_t Solution::writeFastaNodes(FILE * out, size_t id) {
    if(_auvb.size()>0){
        size_t l=_auvb.size(); // should be even
        size_t l2=l/2;
        size_t lmin=l2/2;
        size_t lmax=l2-lmin;
        // the size of breakpoint sequences can be odd
        //keeping a and b of length lmax, and u and v of length lmin
        string a=_auvb.substr(0,lmax);
        string u=_auvb.substr(lmax,lmin);
        string v=_auvb.substr(l2,lmin);
        string b=_auvb.substr(l2+lmin,lmax);
        stringstream res;
        res <<">inv_"<<id<<" a\n"<<a<<"\n";
        res <<">u\n"<<u<<"\n";
        res <<">v\n"<<v<<"\n";
        res <<">b\n"<<b<<"\n";
        
        fprintf(out,"%s",res.str().c_str());
        return 1;
    }
    else{
        return 0;
    }
}


string Solution::reverseComplement(const string& dna) {
    
    //cout<<"dna="<<dna<<endl;
    string revComp= "";
    
    for (string::const_reverse_iterator it = dna.rbegin(); it != dna.rend(); it++) {
        switch(*it) { //each element in the temp vector to it's complement element.
            case 'a' :
                revComp += "t";
                break;
            case 't' :
                revComp += "a";
                break;
            case 'c' :
                revComp += "g";
                break;
            case 'g' :
                revComp += "c";
                break;
            case 'A' :
                revComp += "T";
                break;
            case 'T' :
                revComp += "A";
                break;
            case 'C' :
                revComp += "G";
                break;
            case 'G' :
                revComp += "C";
                break;
        }
    }
    
    //cout<<"revComp="<<revComp<<endl;
    return revComp;
}
