


#include <gatb/gatb_core.hpp>
#include <DBGWalker.hpp>
#include <LCS.hpp>
#include <Read2SV.hpp>
#include <iostream>

#include<mach/mach.h>



// We use the required packages

using namespace std;


char * getVersion(){
	return (char *)"1.0.0 TODO Machin Licence";
}
//#define VERBOSE


void print_usage_and_exit(char * name){
	fprintf (stderr, "NAME\n%s, version %s\n", name, getVersion());
	fprintf (stderr, "\nUSAGE\n%s -i input_graph [-t tolerence walk back] [-o name] [-h] \n", name);
	fprintf (stderr, "\nDESCRIPTION\n");
    
    
	fprintf (stderr, "\nMANDATORY\n");
	fprintf (stderr, "\t -i STRING: File name of the input graph (.h5), omitting the extension name\n");
    
	fprintf (stderr, "\nOPTIONS\n");
	fprintf (stderr, "\t -o STRING file_name for writing results. Default: standard output \n");
	fprintf (stderr, "\t -m INT: max similarity percentage: Inversions with a and b' (or u and v') whose longuest common subsequence size is bigger than k*(this value)/100 are discarded. Defaults: 5 \n");
	fprintf (stderr, "\t -r INT: max repeated size suffix of u and v': Defaults: 5 \n");
	fprintf (stderr, "\t -c INT: local complexity threshold: Defaults: 100 \n");
	fprintf (stderr, "\t -h prints this message and exit\n");
    
    
	exit(0);
}


int main (int argc, char* argv[])
{
    std::cout.setf(std::ios::unitbuf); // avoids the buffer on the cout.
    
   
    
    char* graphFile = NULL;
    int tolerance_rc= 5;
    int max_percentage = 0.8;
    int local_complexity_threshold = 100;
    
    char * output_file= NULL;
    
    // dealing with options
    while (1)
	{
        int witness = getopt (argc, argv, "hr:i:o:m:c:");
		if (witness == -1){
			break;
		}
		switch (witness)
		{
            case 'i':
                graphFile=strdup(optarg);
                break;
            case 'o':
                output_file=strdup(optarg);
                printf("will output results in %s\n", output_file);
                break;
            case 'm':
                max_percentage=atoi(optarg);
                break;
            case 'h':
                print_usage_and_exit(argv[0]);
                break;
            case 'r':
                tolerance_rc=atoi(optarg);
                break;
            case 'c':
                local_complexity_threshold=atoi(optarg);
                break;
            default:
                printf ("Unknown option %c\n", witness);
                print_usage_and_exit(argv[0]);
		}
	}
    
    if(graphFile == NULL){
        fprintf(stderr," Detected error: you must provide an input graph file \n");
        print_usage_and_exit(argv[0]);
    }
    
    
    
    
    Read2SV read2sv(graphFile, tolerance_rc);
    
    
    if(tolerance_rc>read2sv._kmerSize-1){
        tolerance_rc=read2sv._kmerSize-1;
        fprintf(stderr," Warning : tolerence can't be bigger than k-2, set tolerence to k-1=%d \n", tolerance_rc);
    }
    LCS lcs_instance(read2sv._kmerSize, max_percentage);
    FILE * out;
    if(output_file) out = fopen(output_file, "w");
    else out=stdout;
    if(out == NULL){
        fprintf(stderr," Cannot open file %s for writting \n", output_file);
        exit(1);
    }
    
    read2sv.find_ALL_occurrences_of_inversion_pattern(lcs_instance, out, local_complexity_threshold);
    fclose(out);
    
   
    
    return EXIT_SUCCESS;
    
//    DBGWalker dbgw(graphFile);
//    // en attendant avec noeuds pris au hasard dans le graphe
//    Graph::Iterator<Node> itNodes =dbgw._graph.iterator<Node>();
//    
//    for (itNodes.first(); !itNodes.isDone(); itNodes.next()) {
//        Node nodeA=itNodes.item();
//        
//        int depth=dbgw._graph.getKmerSize();
//        dbgw.find_all_at_depth(nodeA,depth);
//        
//        cout<<"from node "<<dbgw._graph.toString(nodeA)<<" I can reach these nodes at depth "<<depth<<endl;
//        for (int i=0;i<dbgw.reachable_neighbor.size();i++)
//            cout<<"     node "<<dbgw._graph.toString(dbgw.reachable_neighbor[i])<<endl;
//        
//        
//        if(dbgw.reachable_neighbor.empty()) continue;
//        Node nodeU = dbgw.reachable_neighbor[0];
//        dbgw.find_B(nodeA,nodeU,3);
//        
//        cout<<"from nodes au= "<<dbgw._graph.toString(nodeA)<<","<<dbgw._graph.toString(nodeU)<<" I can reach the following set B :"<<endl;
//        for (int i=0;i<dbgw.reachable_neighbor.size();i++)
//            cout<<"     node "<<dbgw._graph.toString(dbgw.reachable_neighbor[i])<<endl;
//        if(dbgw.reachable_neighbor.size()>0) return EXIT_SUCCESS;
//        
//    }
//    return EXIT_SUCCESS;
}

