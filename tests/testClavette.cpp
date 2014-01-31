


#include <gatb/gatb_core.hpp>


#include <Read2SV.hpp>
#include <iostream>



// We use the required packages

using namespace std;




int main (int argc, char* argv[])
{
//    if (argc < 2)
//    {
//        cerr << "you must provide:" << endl;
//        cerr << "   1) graph file"  << endl;
//        return EXIT_FAILURE;
//    }
//    
//    char* graphFile = argv[1];
//     
//    Read2SV read2sv(graphFile, 3);
//    read2sv.find_ALL_occurrences_of_inversion_pattern();
//    
//    // Test pour la méthode checkPath (a bouger dans ../test/testUnit.cpp)
//    // idée : je prends un mot de taille > 2k de mon jeu de données (ici reads1.fa : ../../../thirdparty/gatb-core/gatb-core/test/db/reads1.fa)
//    //        je récupère le premier kmer (a), puis celui à k plus loin (b) et enfin un autre à k+2 plys loin que a (c)
//    //        checkPath(a,b) -> true, checkPath(a,c) -> false
//    // Pour l'instant ca marche pas
//    char* seq = (char*) "CATCATTGTTTATCAATGATAAAATATAATA";
//    int k = read2sv._kmerSize;
//    
//    Node nodeA = read2sv._graph.buildNode(Data(seq));
//    Node nodeB = read2sv._graph.buildNode(Data(seq),k);
//    Node nodeC = read2sv._graph.buildNode(Data(seq),k+3);
//    
//
//    // en attendant avec noeuds pris au hasard dans le graphe
//    //    Graph::Iterator<Node> itNodes =read2sv._graph.iterator<Node>();
//    //    itNodes.first();
//    //    Node nodeA=itNodes.item();
//    //    itNodes.next();
//    //    Node nodeB=itNodes.item();
//    //    itNodes.next();
//    //    Node nodeC=itNodes.item();
//    
//    cout << "node A = " << read2sv._graph.toString(nodeA) << endl;
//    cout << "node B = " << read2sv._graph.toString(nodeB) << endl;
//    cout << "node C = " << read2sv._graph.toString(nodeC) << endl;
//    
//    cout << "checkPath(A,B) = " << boolalpha << read2sv.checkPath(nodeA,nodeB) << endl;
//    cout << "checkPath(A,C) = " << boolalpha << read2sv.checkPath(nodeA,nodeC) << endl;
//    
//    
    return EXIT_SUCCESS;
}
