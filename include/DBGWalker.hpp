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

#ifndef read2SV_DBGWalker_hpp
#define read2SV_DBGWalker_hpp

#include <gatb/gatb_core.hpp>
#include <LCS.hpp>

#include <vector>
#include <string>
#include <sstream>



class DBGWalker {
public:
    Graph _graph; // needed to have the string of a node
    std::vector<Node> reachable_neighbor;
    size_t idxRecursion;
    size_t _lct; // local complexity threshold, size limit for reachable_neighbor (limit also for the product of the current size x a given one)


    DBGWalker()  { }
    DBGWalker(char* graphFile, size_t lct); // no longer used
    DBGWalker(Graph& graph, size_t lct); // no longer used
    ~DBGWalker();

    /**
     * Initialization of the object. Note that it is mandatory...
     */
    void setInfo (const Graph& graph, size_t lct)  { _graph = graph;  _lct=lct;  reachable_neighbor.resize (_lct+1); }
    
    /**
     * Gives the size of the filled reachable_neighbor 
     * Note that reachable_neighbor.size() is to forbid, since it will always output _lct + 1, which is not the actual number of neighbors)
     */
    size_t size () const { return idxRecursion; }
    /**
     * enables to sort a set of DBGwalker objects
     */
    bool operator< (const DBGWalker& other) const  {  return size() < other.size(); }

    /**
     * Finds all kmers reachable from nodeA, after k steps.
     * fills the reachable_neighbor vector
     * The last parameter : maxOtherSize enables to stop the recursion search earlier according to the _lct limit
     */
    void find_all_at_depth(const Node& nodeA, const int depth, size_t maxOtherSize);
    
    /**
     * Provides all kmers that can be reached from node rcu, after k steps
     * constraint: the kmer at step size_tolerance_rc(+1?) (=forbiden kmer) should be distinct from the 'equivalent' in au (= rc(au)[size_tolerance_rc, size_tolerance_rc+k]
     * fills the reachable_neighbor vector
     */
    void find_B(const Node& a, const Node& u, const int size_tolerance_rc, LCS& lcs_instance, size_t maxOtherSize);
    
    
private:
    void recursive_find_all_at_depth(const Node& nodeA, const int depth, size_t maxOtherSize);
    void recursive_find_B(const Node& cur, const int size_tolerance_rc, const int depth, const int distance_from_start, const std::string& forbiden_kmer,
                          LCS& lcs_instance, size_t maxOtherSize, const std::string& rca);
};


#endif
