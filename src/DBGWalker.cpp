/*****************************************************************************
 *   TakeABreak: Breakpoint detection from raw unassembled NGS reads
 *   A tool from the GATB (Genome Assembly Tool Box)
 *   Copyright (C) 2014  C.Lemaitre, P.Peterlongo, E.Drezen
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



//
//  BBGWalker.cpp
//  read2SV
//
//  Created by Pierre Peterlongo on 26/11/2013.
//  Copyright (c) 2013 Pierre Peterlongo. All rights reserved.
//



#include <DBGWalker.hpp>
#include <iostream>
using namespace std;



DBGWalker::~DBGWalker(){
}
    
/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
DBGWalker::DBGWalker (char* graphFile)
{
    _graph = Graph::load (graphFile);
}
DBGWalker::DBGWalker (Graph& _graph){
    this->_graph=_graph;
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
void DBGWalker::recursive_find_all_at_depth(const Node& nodeA, const int depth){
    
    // deal with possible eroneous call of the function
    if(depth<0) return;
    
    // finished, we store th reached node in the "reachable_neighbor“ vector
    if(depth==0) {
        reachable_neighbor.push_back(nodeA);
        return;
    }
//    cout<<"depth = "<<depth<<endl;
    // go on on all soons.
    Graph::Vector<Node> neighbor = _graph.neighbors<Node> (nodeA, DIR_OUTCOMING);
//    cout<<" n size =  "<<neighbor.size()<<endl;
    for(int i=0;i<neighbor.size();i++) recursive_find_all_at_depth(neighbor[i],depth-1);
//    cout<<"finish depth = "<<depth<<endl;
}


/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
void DBGWalker::find_all_at_depth(const Node& nodeA, const int depth){
    reachable_neighbor.clear();
    recursive_find_all_at_depth(nodeA, depth);
}




/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
void DBGWalker::recursive_find_B(const Node& cur, const int size_tolerance_rc, const int depth, const int distance_from_start, const string& forbiden_kmer,
                                 LCS& lcs_instance, const string& rca){
    // deal with possible eroneous call of the function
    if(depth<0) return;
//#define debug
#ifdef debug
    cout<<"DEBUG rec_find_B "<<depth<<" "<<_graph.toString(cur).c_str()<<" "<<forbiden_kmer<<" "<<(_graph.toString(cur).compare(forbiden_kmer)==0)<<endl;
#endif //debug
    // finished, we store the reached node in the "reachable_neighbor“ vector (this is thus bbar)
    if(depth==0) {
#ifdef debug
        cout<<"lcs a'b "<<lcs_instance.size_lcs(rca, _graph.toString(cur))<<endl;
#endif //debug
        if(lcs_instance.size_lcs(rca, _graph.toString(cur))<=lcs_instance.lcs_threshold)
           reachable_neighbor.push_back(cur);
        return;
    }
#ifdef debug
        if(size_tolerance_rc && distance_from_start==(size_tolerance_rc+1) )
            cout<<_graph.toString(cur).c_str()<<" "<<forbiden_kmer<<" "<<(_graph.toString(cur).compare(forbiden_kmer)==0)<<endl;
#endif
    // checks that we are not comming back "a" to after "size_tolerance" steps.
    if(distance_from_start==(size_tolerance_rc+1) && _graph.toString(cur).compare(forbiden_kmer)==0)
        return;
    
    // go on on all soons.
    Graph::Vector<Node> neighbor = _graph.neighbors<Node> (cur, DIR_OUTCOMING);
    for(int i=0;i<neighbor.size();i++) recursive_find_B(neighbor[i],
                                                        size_tolerance_rc,
                                                        depth-1,
                                                        distance_from_start+1,
                                                        forbiden_kmer,
                                                        lcs_instance,
                                                        rca);
    
}



/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
void DBGWalker::find_B(const Node& a,const Node& u, const int size_tolerance_rc, LCS& lcs_instance){
//#define debug
#ifdef debug
    cout<<"find B"<<endl;
    cout<<"a="<<_graph.toString(a)<<endl;
    cout<<"u="<<_graph.toString(u)<<endl;
#endif // debug
    
    Node rcu = _graph.reverse(u);
    Node rca = _graph.reverse(a);
    string rc_au = _graph.toString(rcu)+_graph.toString(rca);
    string forbiden_kmer = rc_au.substr(size_tolerance_rc+1, _graph.getKmerSize());
    
    
    reachable_neighbor.clear();
    recursive_find_B(rcu,size_tolerance_rc,_graph.getKmerSize(),0,forbiden_kmer, lcs_instance, _graph.toString(rca));
    
    // checks if the LCS is sufficient.
    
    
}



















