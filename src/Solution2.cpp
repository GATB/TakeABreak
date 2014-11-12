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

Solution::Solution(string a,string u,string v,string b){
    _a=a;
    _b=b;
    _u=u;
    _v=v;
}

Solution::Solution(){
    
}

void Solution::setSequences(string au,string vb){
    size_t l=au.size();
    size_t lmin=l/2;
    size_t lmax=l-lmin;
    // the size of breakpoint sequences can be odd
    //cerr<<l<<"  min="<<lmin<<"  max="<<lmax<<endl;
    //keeping a and b of length lmax, and u and v of length lmin
    _a=au.substr(0,lmax);
    _u=au.substr(lmax,lmin);
    _v=vb.substr(0,lmin);
    _b=vb.substr(lmin,lmax);
}


Solution::~Solution(){
}

bool Solution::operator< (const Solution& other) const  {
    
    if(_a==other._a){
        if(_u==other._u){
            if(_v==other._v){
                if(_b==other._b){
                    return false;
                }
                else{
                    return _b<other._b;
                }
            }
            else{
                return _v<other._v;
            }
        }
        else{
            return _u<other._u;
        }
    }
    else{
        return _a<other._a;
    }
}

void Solution::writeFastaOutput(FILE * out, size_t id) {
    //if(_a.size()>0){
        stringstream res;
        res <<">inv_"<<id<<" au\n"<<_a<<_u<<"\n";
        res <<">vb\n"<<_v<<_b<<"\n";
        res <<">av'\n"<<_a<<reverseComplement(_v)<<"\n";
        res <<">u'b\n"<<reverseComplement(_u)<<_b<<"\n";
        
        fprintf(out,"%s",res.str().c_str());
    //}
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
