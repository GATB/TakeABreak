/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  R.Chikhi, G.Rizk, E.Drezen
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
    
    DBGWalker(char* graphFile);
    DBGWalker(Graph& _graph);
    ~DBGWalker();
    
    /**
     * Finds all kmers reachable from nodeA, after k steps.
     * fills the reachable_neighbor vector
     */
    void find_all_at_depth(const Node& nodeA, const int depth);
    
    /**
     * Provides all kmers that can be reached from node rcu, after k steps
     * constraint: the kmer at step size_tolerance_rc(+1?) (=forbiden kmer) should be distinct from the 'equivalent' in au (= rc(au)[size_tolerance_rc, size_tolerance_rc+k]
     * fills the reachable_neighbor vector
     */
    void find_B(const Node& a, const Node& u, const int size_tolerance_rc, LCS& lcs_instance);
    
    
private:
    void recursive_find_all_at_depth(const Node& nodeA, const int depth);
    void recursive_find_B(const Node& cur, const int size_tolerance_rc, const int depth, const int distance_from_start, const std::string& forbiden_kmer,
                          LCS& lcs_instance, const std::string& rca);
};


#endif
